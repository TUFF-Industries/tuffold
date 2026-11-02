//! What someone can send when Tuffold breaks: a crash log written on panic,
//! and a report of the machine, the models, a project's counts, and the last
//! crashes. Nothing leaves the machine by itself; the report is a text file
//! to read before sending. It holds counts, never sample, run, or protein
//! names, and shows the home folder as `~`. A crash message can still name a
//! file, so the report says to check. Panics caught on purpose (a guarded
//! spectra reader, one failed file in a search) are logged too, so a log means
//! a panic, not always a crash.

use std::fmt::Write as _;
use std::fs::{self, OpenOptions};
use std::io::Write as _;
use std::path::{Path, PathBuf};
use std::sync::atomic::{AtomicUsize, Ordering};
use std::time::{SystemTime, UNIX_EPOCH};

use tuffold_store::Project;

use crate::deep::models_status;

const CRASH_PREFIX: &str = "crash-";
/// Panic logs one process writes: a panic repeated per spectrum mustn't fill
/// the disk, and the first ones are the ones that explain it.
const CRASHES_PER_RUN: usize = 20;
/// Older crash logs are removed when a program starts.
const KEPT_CRASHES: usize = 20;
/// Crash logs quoted in a report.
const REPORTED_CRASHES: usize = 3;

/// Where crash logs go: `~/Library/Logs/Tuffold` on macOS, else
/// `$XDG_STATE_HOME/tuffold/logs` (`~/.local/state`).
pub fn log_dir() -> Option<PathBuf> {
    let home = std::env::var_os("HOME").map(PathBuf::from);
    if cfg!(target_os = "macos") {
        return home.map(|home| home.join("Library/Logs/Tuffold"));
    }
    std::env::var_os("XDG_STATE_HOME")
        .map(PathBuf::from)
        .or_else(|| home.map(|home| home.join(".local/state")))
        .map(|state| state.join("tuffold/logs"))
}

/// Keep writing panics to a crash log, on top of what the default hook does.
pub fn install_panic_log(program: &'static str) {
    let previous = std::panic::take_hook();
    static WRITTEN: AtomicUsize = AtomicUsize::new(0);
    std::panic::set_hook(Box::new(move |info| {
        let write = WRITTEN.fetch_add(1, Ordering::Relaxed) < CRASHES_PER_RUN;
        if let Some(dir) = log_dir().filter(|_| write) {
            let backtrace = std::backtrace::Backtrace::force_capture();
            let thread = std::thread::current();
            let text = format!(
                "{}\nthread: {}\n{info}\n\n{backtrace}\n",
                machine(program),
                thread.name().unwrap_or("unnamed"),
            );
            // A failing log mustn't hide the panic itself.
            let _ = write_crash(&dir, &text);
        }
        previous(info);
    }));
    if let Some(dir) = log_dir() {
        prune_crashes(&dir, KEPT_CRASHES);
    }
}

fn machine(program: &str) -> String {
    format!(
        "{program} {} on {} {}",
        env!("CARGO_PKG_VERSION"),
        std::env::consts::OS,
        std::env::consts::ARCH
    )
}

fn write_crash(dir: &Path, text: &str) -> std::io::Result<PathBuf> {
    fs::create_dir_all(dir)?;
    let seconds = SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .map_or(0, |t| t.as_secs());
    // Panics cascade within a second, and the first is the one that matters.
    for attempt in 0..CRASHES_PER_RUN {
        let path = dir.join(format!(
            "{CRASH_PREFIX}{seconds:012}-{}-{attempt:02}.log",
            std::process::id()
        ));
        match OpenOptions::new().write(true).create_new(true).open(&path) {
            Ok(mut file) => {
                file.write_all(text.as_bytes())?;
                return Ok(path);
            }
            Err(error) if error.kind() == std::io::ErrorKind::AlreadyExists => {}
            Err(error) => return Err(error),
        }
    }
    Err(std::io::ErrorKind::AlreadyExists.into())
}

/// Crash logs, newest first.
fn crash_logs(dir: &Path) -> Vec<PathBuf> {
    let mut logs: Vec<PathBuf> = fs::read_dir(dir)
        .into_iter()
        .flatten()
        .flatten()
        .map(|entry| entry.path())
        .filter(|path| {
            path.file_name()
                .and_then(|name| name.to_str())
                .is_some_and(|name| name.starts_with(CRASH_PREFIX))
        })
        .collect();
    logs.sort_unstable_by(|a, b| b.cmp(a));
    logs
}

fn prune_crashes(dir: &Path, keep: usize) {
    for old in crash_logs(dir).into_iter().skip(keep) {
        let _ = fs::remove_file(old);
    }
}

/// The report, with this machine's crash logs.
pub fn diagnostics_report(project: Option<&Path>) -> String {
    let home = std::env::var("HOME").ok();
    report(
        log_dir().as_deref(),
        project,
        home.as_deref(),
        &models_line(),
    )
}

fn models_line() -> String {
    let status = models_status();
    match status.problem {
        Some(problem) => format!("unusable ({problem})"),
        None if status.installed => "installed".to_string(),
        None => "not installed".to_string(),
    }
}

fn report(logs: Option<&Path>, project: Option<&Path>, home: Option<&str>, models: &str) -> String {
    let mut text = format!(
        "Tuffold diagnostics. Read this before sending it: crash messages can name a file.\n\n\
         {}\ndeep-learning models: {models}\n",
        machine("tuffold")
    );
    if let Some(path) = project {
        text.push_str(&project_section(path));
    }
    let crashes = logs.map(crash_logs).unwrap_or_default();
    let _ = writeln!(text, "\npanic logs: {}", crashes.len());
    for crash in crashes.iter().take(REPORTED_CRASHES) {
        let name = crash.file_name().unwrap_or_default().to_string_lossy();
        let body = fs::read_to_string(crash).unwrap_or_else(|e| format!("unreadable: {e}\n"));
        let _ = write!(text, "\n--- {name}\n{body}");
    }
    // As a whole path component: /srv/al isn't the start of /srv/alice.
    match home
        .map(|home| home.trim_end_matches('/'))
        .filter(|home| !home.is_empty())
    {
        Some(home) => text.replace(&format!("{home}/"), "~/"),
        None => text,
    }
}

fn project_section(path: &Path) -> String {
    let opened = Project::open_read_only(path).and_then(|project| project.summary());
    match opened {
        // The error may quote the path, whose file name may be a study's.
        Err(error) => format!(
            "\nproject: can't be read: {}\n",
            error
                .to_string()
                .replace(&path.display().to_string(), "<project>")
        ),
        Ok(summary) => format!(
            "\nproject: {} imports, {} biosamples, {} runs, {} spectra, {} PSMs ({} decoy), \
             {} peptides, {} proteins ({} decoy, {} with sequences)\n",
            summary.imports,
            summary.biosamples,
            summary.ms_runs,
            summary.spectra,
            summary.psms,
            summary.decoy_psms,
            summary.peptides,
            summary.proteins,
            summary.decoy_proteins,
            summary.proteins_with_sequence,
        ),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn crash_logs_are_kept_newest_first_and_pruned() {
        let dir = tempfile::tempdir().unwrap();
        for seconds in [5, 9, 7] {
            let name = format!("{CRASH_PREFIX}{seconds:012}-1.log");
            fs::write(dir.path().join(name), seconds.to_string()).unwrap();
        }
        fs::write(dir.path().join("other.txt"), "not a crash").unwrap();
        let first = crash_logs(dir.path());
        assert_eq!(first.len(), 3);
        assert_eq!(fs::read_to_string(&first[0]).unwrap(), "9");
        prune_crashes(dir.path(), 1);
        let left = crash_logs(dir.path());
        assert_eq!(left.len(), 1);
        assert_eq!(fs::read_to_string(&left[0]).unwrap(), "9");
        assert!(dir.path().join("other.txt").exists());
    }

    #[test]
    fn a_report_has_counts_not_names_and_hides_the_home_folder() {
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("demo.tuff");
        crate::demo::create_demo_project(&path).unwrap();
        let logs = dir.path().join("logs");
        write_crash(
            &logs,
            &format!("panicked at {}/thing", dir.path().display()),
        )
        .unwrap();
        let home = dir.path().to_string_lossy().into_owned();
        let text = report(Some(&logs), Some(&path), Some(&home), "not installed");
        assert!(
            text.contains("deep-learning models: not installed"),
            "{text}"
        );
        assert!(text.contains("4 biosamples"), "{text}");
        assert!(text.contains("panic logs: 1"), "{text}");
        assert!(text.contains("panicked at ~/thing"), "{text}");
        assert!(!text.contains(&home), "{text}");
        let project = Project::open(&path).unwrap();
        for run in project.summary().unwrap().runs {
            assert!(!text.contains(&run.name), "{} in {text}", run.name);
            assert!(
                !text.contains(&run.biosample),
                "{} in {text}",
                run.biosample
            );
        }
    }

    #[test]
    fn a_broken_project_is_reported_not_fatal() {
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("broken.tuff");
        fs::write(&path, "not a project").unwrap();
        let text = report(None, Some(&path), None, "installed");
        assert!(text.contains("project: can't be read"), "{text}");
        assert!(text.contains("panic logs: 0"), "{text}");
        assert!(!text.contains("broken.tuff"), "{text}");
    }

    #[test]
    fn panics_in_the_same_second_keep_their_own_logs() {
        let dir = tempfile::tempdir().unwrap();
        let first = write_crash(dir.path(), "first").unwrap();
        let second = write_crash(dir.path(), "second").unwrap();
        assert_ne!(first, second);
        assert_eq!(fs::read_to_string(first).unwrap(), "first");
    }

    #[test]
    fn the_home_folder_is_hidden_only_as_a_whole_path_component() {
        let dir = tempfile::tempdir().unwrap();
        write_crash(dir.path(), "at /srv/al/x and /srv/alice/y").unwrap();
        let text = report(Some(dir.path()), None, Some("/srv/al/"), "installed");
        assert!(text.contains("at ~/x and /srv/alice/y"), "{text}");
        assert!(report(None, None, Some("/"), "installed").contains("Read this"));
    }

    #[test]
    fn a_report_leaves_an_older_project_as_it_found_it() {
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("demo.tuff");
        crate::demo::create_demo_project(&path).unwrap();
        let before = fs::read(&path).unwrap();
        report(None, Some(&path), None, "installed");
        assert_eq!(fs::read(&path).unwrap(), before);
    }
}
