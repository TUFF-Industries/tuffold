#![allow(dead_code)] // each integration-test binary uses a subset

use std::fs;
use std::path::{Path, PathBuf};

pub fn fixture_bundle() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("../../fixtures/fragpipe_synthetic")
}

/// Copy the synthetic bundle into `dir/name` and return its path.
pub fn copy_bundle(dir: &Path, name: &str) -> PathBuf {
    let target = dir.join(name);
    copy_dir(&fixture_bundle(), &target);
    target
}

fn copy_dir(from: &Path, to: &Path) {
    fs::create_dir_all(to).unwrap();
    for entry in fs::read_dir(from).unwrap() {
        let path = entry.unwrap().path();
        let target = to.join(path.file_name().unwrap());
        if path.is_dir() {
            copy_dir(&path, &target);
        } else {
            fs::copy(&path, &target).unwrap();
        }
    }
}

pub fn write_tsv(path: &Path, header: &[&str], rows: &[&[&str]]) {
    let mut text = header.join("\t");
    text.push('\n');
    for row in rows {
        text.push_str(&row.join("\t"));
        text.push('\n');
    }
    fs::write(path, text).unwrap();
}

/// Remove one named column from every line of a TSV.
pub fn drop_column(path: &Path, column: &str) {
    let text = fs::read_to_string(path).unwrap();
    let mut lines = text.lines();
    let header: Vec<&str> = lines.next().unwrap().split('\t').collect();
    let dropped = header.iter().position(|name| *name == column).unwrap();
    let rebuilt: Vec<String> = text
        .lines()
        .map(|line| {
            line.split('\t')
                .enumerate()
                .filter(|(index, _)| *index != dropped)
                .map(|(_, field)| field)
                .collect::<Vec<_>>()
                .join("\t")
        })
        .collect();
    fs::write(path, rebuilt.join("\n") + "\n").unwrap();
}

/// Append a row to a TSV, filling named columns and leaving the rest empty.
pub fn append_row(path: &Path, values: &[(&str, &str)]) {
    let text = fs::read_to_string(path).unwrap();
    let header: Vec<&str> = text.lines().next().unwrap().split('\t').collect();
    let mut row = vec![""; header.len()];
    for (column, value) in values {
        row[header.iter().position(|name| name == column).unwrap()] = value;
    }
    fs::write(path, format!("{}{}\n", text, row.join("\t"))).unwrap();
}

/// A minimal bundle: one `control-1` folder, stub optional reports.
pub fn minimal_bundle(dir: &Path) -> PathBuf {
    let bundle = dir.join("bundle");
    fs::create_dir_all(bundle.join("control-1")).unwrap();
    fs::write(
        bundle.join("modified_peptide.tsv"),
        "Peptide Sequence\tModified Sequence\n",
    )
    .unwrap();
    fs::write(bundle.join("ion.tsv"), "Peptide Sequence\tCharge\n").unwrap();
    bundle
}
