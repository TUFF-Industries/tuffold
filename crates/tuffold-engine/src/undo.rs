//! One level of undo, and redo, for the changes that can't be recomputed:
//! imports, searches, and sample-group edits. (Validation and thresholds are
//! derived from the data and the saved settings, so they are redone, not
//! undone.)
//!
//! Before such a change the whole project is copied to `<project>.undo`, next
//! to it and visible, so the disk it takes is plain to see. The copy is only
//! kept when the change actually changed something, so a failed import can't
//! overwrite the real undo point. Undo swaps the two files, so undoing again
//! redoes.

use std::path::{Path, PathBuf};
use std::time::{SystemTime, UNIX_EPOCH};

use tuffold_store::{Project, StoreError, UndoNote};

#[derive(Debug, thiserror::Error)]
pub enum UndoError {
    #[error("there is no change to undo in {0}")]
    NothingToUndo(String),
    #[error(
        "{0} was interrupted in the middle of a change (its -journal file is still there); \
         open the project once to let it recover, then undo"
    )]
    Interrupted(String),
    #[error(
        "{undo} belongs to a different project than {project} (the project was replaced or \
         recreated); delete {undo} if you don't need it"
    )]
    Unrelated { project: String, undo: String },
    #[error(
        "an earlier undo of {project} was interrupted: {swap} holds the project as it was \
         before that undo. Keep whichever of the two you want as the project, and delete \
         the other"
    )]
    Unfinished { project: String, swap: String },
    #[error(transparent)]
    Store(#[from] StoreError),
    #[error("{path}: {source}")]
    Io {
        path: String,
        source: std::io::Error,
    },
}

/// `study.tuff` -> `study.tuff.undo`.
pub fn undo_path(project: &Path) -> PathBuf {
    with_suffix(project, ".undo")
}

fn with_suffix(path: &Path, suffix: &str) -> PathBuf {
    let mut name = path.as_os_str().to_owned();
    name.push(suffix);
    PathBuf::from(name)
}

fn io(path: &Path) -> impl FnOnce(std::io::Error) -> UndoError + '_ {
    move |source| UndoError::Io {
        path: path.display().to_string(),
        source,
    }
}

fn now() -> i64 {
    SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .map_or(0, |elapsed| {
            i64::try_from(elapsed.as_secs()).unwrap_or(i64::MAX)
        })
}

/// A copy of the project taken before a change. [`Checkpoint::keep`] makes
/// it the undo point; dropped without keeping, it is deleted.
pub struct Checkpoint {
    copy: PathBuf,
    target: PathBuf,
    label: String,
}

impl Checkpoint {
    /// Copy `project` aside. `None` for a project with no file (in memory).
    pub fn take(project: &mut Project, label: &str) -> Result<Option<Self>, UndoError> {
        // Before the copy, so the copy carries it too.
        project.identity()?;
        let Some(path) = project.path().map(Path::to_path_buf) else {
            return Ok(None);
        };
        let copy = with_suffix(&path, ".undo.new");
        // A leftover from a crash is never anything to keep.
        if copy.exists() {
            std::fs::remove_file(&copy).map_err(io(&copy))?;
        }
        if let Err(error) = project.save_copy(&copy) {
            // A partial copy (a full disk, say) must not linger.
            let _ = std::fs::remove_file(&copy);
            return Err(error.into());
        }
        Ok(Some(Self {
            copy,
            target: undo_path(&path),
            label: label.to_string(),
        }))
    }

    /// The change happened: this copy is now the undo point.
    pub fn keep(self) -> Result<(), UndoError> {
        let note = UndoNote {
            label: self.label.clone(),
            redo: false,
            saved_at: now(),
        };
        Project::open(&self.copy)?.set_undo_note(Some(&note))?;
        replace(&self.copy, &self.target)
    }
}

impl Drop for Checkpoint {
    fn drop(&mut self) {
        // Kept copies were renamed away; this only removes unkept ones.
        let _ = std::fs::remove_file(&self.copy);
    }
}

/// Move `from` over `to`, replacing it atomically.
fn replace(from: &Path, to: &Path) -> Result<(), UndoError> {
    std::fs::rename(from, to).map_err(io(from))
}

/// What undo would do now: the change it undoes (or redoes, once undone).
pub fn last_change(project: &Path) -> Result<Option<UndoNote>, UndoError> {
    let copy = undo_path(project);
    if !copy.is_file() {
        return Ok(None);
    }
    Ok(Project::open(&copy)?.undo_note()?)
}

/// Undo the last change (or redo it, if it was just undone). Returns what
/// was undone or redone.
///
/// The swap never leaves the project missing: the current file is first
/// linked (or copied) aside, the copy then atomically replaces it, and only
/// then does the old state become the new copy. If that is interrupted, the
/// aside file is still there and the next undo refuses, saying so.
pub fn undo(project: &Path) -> Result<UndoNote, UndoError> {
    let copy = undo_path(project);
    let aside = with_suffix(project, ".undo.swap");
    if aside.exists() {
        return Err(UndoError::Unfinished {
            project: project.display().to_string(),
            swap: aside.display().to_string(),
        });
    }
    let Some(note) = last_change(project)? else {
        return Err(UndoError::NothingToUndo(project.display().to_string()));
    };
    for file in [project, copy.as_path()] {
        if with_suffix(file, "-journal").exists() {
            return Err(UndoError::Interrupted(file.display().to_string()));
        }
    }
    let mut current = Project::open(project)?;
    if current.identity()? != Project::open(&copy)?.identity()? {
        return Err(UndoError::Unrelated {
            project: project.display().to_string(),
            undo: copy.display().to_string(),
        });
    }
    // No other connection may write while the files move.
    current.exclusive(|| swap(project, &copy, &aside))??;
    drop(current);
    // The project carries no note; the old state becomes the copy that
    // reverses this undo.
    Project::open(project)?.set_undo_note(None)?;
    let reverse = UndoNote {
        label: note.label.clone(),
        redo: !note.redo,
        saved_at: now(),
    };
    Project::open(&copy)?.set_undo_note(Some(&reverse))?;
    Ok(note)
}

fn swap(project: &Path, copy: &Path, aside: &Path) -> Result<(), UndoError> {
    // A hard link is instant; a filesystem without them gets a copy.
    if std::fs::hard_link(project, aside).is_err() {
        std::fs::copy(project, aside).map_err(io(project))?;
    }
    if let Err(error) = replace(copy, project) {
        // Nothing moved: drop the extra link so the next undo isn't refused.
        let _ = std::fs::remove_file(aside);
        return Err(error);
    }
    replace(aside, copy)
}

/// How many imports a project has: imports and searches that changed
/// nothing leave it alone, so it tells whether a checkpoint is worth keeping.
pub(crate) fn import_count(project: &Project) -> usize {
    project.imports().map_or(0, |imports| imports.len())
}
