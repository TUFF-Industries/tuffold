//! The server binary over stdio: initialize, list tools, call one.

use std::io::{BufRead, BufReader, Write};
use std::process::{Command, Stdio};

use serde_json::{Value, json};

#[test]
fn the_server_speaks_mcp_over_stdio() {
    let dir = tempfile::tempdir().unwrap();
    let project = dir.path().join("demo.tuff");
    let mut child = Command::new(env!("CARGO_BIN_EXE_tuffold-mcp"))
        .stdin(Stdio::piped())
        .stdout(Stdio::piped())
        .stderr(Stdio::null())
        .spawn()
        .unwrap();
    let mut input = child.stdin.take().unwrap();
    let mut output = BufReader::new(child.stdout.take().unwrap());
    let mut send = |message: Value| {
        writeln!(input, "{message}").unwrap();
        input.flush().unwrap();
    };
    let mut receive = || {
        let mut line = String::new();
        output.read_line(&mut line).unwrap();
        serde_json::from_str::<Value>(&line).unwrap()
    };

    send(
        json!({"jsonrpc": "2.0", "id": 1, "method": "initialize", "params": {
            "protocolVersion": "2025-06-18", "capabilities": {},
            "clientInfo": {"name": "test", "version": "0"},
        }}),
    );
    let initialized = receive();
    assert!(
        initialized["result"]["capabilities"]["tools"].is_object(),
        "{initialized}"
    );
    let instructions = initialized["result"]["instructions"].as_str().unwrap();
    assert!(instructions.contains(".tuff"));
    send(json!({"jsonrpc": "2.0", "method": "notifications/initialized"}));

    send(json!({"jsonrpc": "2.0", "id": 2, "method": "tools/list"}));
    let listed = receive();
    let mut names: Vec<&str> = listed["result"]["tools"]
        .as_array()
        .unwrap()
        .iter()
        .map(|tool| tool["name"].as_str().unwrap())
        .collect();
    names.sort_unstable();
    assert_eq!(
        names,
        [
            "compare_all_categories",
            "compare_categories",
            "create_demo",
            "export_library",
            "export_mzid",
            "export_workbook",
            "import_data",
            "list_protein_groups",
            "match_between_runs",
            "methods_text",
            "peptide_spectra",
            "project_summary",
            "protein_group",
            "quantify_reporters",
            "remove_plex",
            "sample_qc",
            "search_spectra",
            "set_categories",
            "set_plex",
            "set_thresholds",
            "undo_last_change",
            "validate",
        ]
    );

    send(
        json!({"jsonrpc": "2.0", "id": 3, "method": "tools/call", "params": {
            "name": "create_demo", "arguments": {"project": project},
        }}),
    );
    let created = receive();
    assert_eq!(created["result"]["isError"], false, "{created}");
    // A failing call is a readable tool error, not a protocol error.
    send(
        json!({"jsonrpc": "2.0", "id": 4, "method": "tools/call", "params": {
            "name": "protein_group", "arguments": {"project": project, "group_id": -1},
        }}),
    );
    let failed = receive();
    assert_eq!(failed["result"]["isError"], true, "{failed}");
    let message = failed["result"]["content"][0]["text"].as_str().unwrap();
    assert!(message.contains("-1"), "{message}");

    drop(input);
    child.wait().unwrap();
}
