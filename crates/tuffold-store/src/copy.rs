//! Whole-project copies, for undo: a consistent snapshot written with
//! `VACUUM INTO`, and a note in `meta` saying which change the copy undoes.

use std::path::Path;

use rusqlite::{Connection, OptionalExtension};

use crate::StoreError;

/// What a saved copy is for: restoring it undoes `label` (or, once undone,
/// redoes it).
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct UndoNote {
    pub label: String,
    /// Restoring this copy redoes the change instead of undoing it.
    pub redo: bool,
    /// Unix seconds when the copy was saved.
    pub saved_at: i64,
}

pub(crate) fn save_copy(conn: &Connection, to: &Path) -> Result<(), StoreError> {
    if to.exists() {
        return Err(StoreError::AlreadyExists(to.display().to_string()));
    }
    let target = to
        .to_str()
        .ok_or_else(|| StoreError::Damaged(format!("{} is not a UTF-8 path", to.display())))?;
    conn.execute("VACUUM INTO ?1", [target])?;
    Ok(())
}

/// This project's identity: a random-enough id stored in `meta` the first
/// time it is asked for. Copies share it, so a copy can be checked against
/// the project it was made from.
pub(crate) fn identity(conn: &Connection) -> Result<String, StoreError> {
    let existing: Option<String> = conn
        .query_row(
            "SELECT value FROM meta WHERE key = 'project_id'",
            [],
            |row| row.get(0),
        )
        .optional()?;
    if let Some(id) = existing {
        return Ok(id);
    }
    let nanos = std::time::SystemTime::now()
        .duration_since(std::time::UNIX_EPOCH)
        .map_or(0, |elapsed| elapsed.as_nanos());
    let id = format!("{nanos:x}-{:x}", std::process::id());
    conn.execute(
        "INSERT INTO meta (key, value) VALUES ('project_id', ?1)",
        [&id],
    )?;
    Ok(id)
}

/// Run `work` while holding SQLite's exclusive lock on the project, so no
/// other connection can be writing to it meanwhile.
pub(crate) fn exclusive<T>(conn: &Connection, work: impl FnOnce() -> T) -> Result<T, StoreError> {
    conn.execute_batch("BEGIN EXCLUSIVE")?;
    let result = work();
    // Nothing was written through this connection; just release the lock.
    conn.execute_batch("ROLLBACK")?;
    Ok(result)
}

pub(crate) fn undo_note(conn: &Connection) -> Result<Option<UndoNote>, StoreError> {
    let value = |key: &str| -> Result<Option<String>, StoreError> {
        Ok(conn
            .query_row("SELECT value FROM meta WHERE key = ?1", [key], |row| {
                row.get(0)
            })
            .optional()?)
    };
    let Some(label) = value("undo_label")? else {
        return Ok(None);
    };
    Ok(Some(UndoNote {
        label,
        redo: value("undo_redo")?.as_deref() == Some("1"),
        saved_at: value("undo_saved_at")?
            .and_then(|at| at.parse().ok())
            .unwrap_or_default(),
    }))
}

pub(crate) fn set_undo_note(
    conn: &mut Connection,
    note: Option<&UndoNote>,
) -> Result<(), StoreError> {
    let tx = conn.transaction()?;
    tx.execute(
        "DELETE FROM meta WHERE key IN ('undo_label', 'undo_redo', 'undo_saved_at')",
        [],
    )?;
    if let Some(note) = note {
        let redo = if note.redo { "1" } else { "0" };
        for (key, value) in [
            ("undo_label", note.label.clone()),
            ("undo_redo", redo.to_string()),
            ("undo_saved_at", note.saved_at.to_string()),
        ] {
            tx.execute(
                "INSERT INTO meta (key, value) VALUES (?1, ?2)",
                [key, &value],
            )?;
        }
    }
    tx.commit()?;
    Ok(())
}
