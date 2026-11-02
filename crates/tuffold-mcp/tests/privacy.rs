//! Privacy mode end to end: the server started with `--private=proteins`
//! never puts a sample, category, run, file, or protein name in a reply, and
//! aliases sent back reach the real names.

use std::io::{BufRead, BufReader, Write};
use std::path::Path;
use std::process::{ChildStdin, ChildStdout, Command, Stdio};

use serde_json::{Value, json};
use tuffold_store::Project;

struct Server {
    input: ChildStdin,
    output: BufReader<ChildStdout>,
    next_id: u64,
}

impl Server {
    fn start(args: &[&str]) -> (std::process::Child, Self) {
        let mut child = Command::new(env!("CARGO_BIN_EXE_tuffold-mcp"))
            .args(args)
            .env_remove("TUFFOLD_MCP_PRIVATE")
            .stdin(Stdio::piped())
            .stdout(Stdio::piped())
            .stderr(Stdio::null())
            .spawn()
            .unwrap();
        let mut server = Self {
            input: child.stdin.take().unwrap(),
            output: BufReader::new(child.stdout.take().unwrap()),
            next_id: 1,
        };
        let initialized = server.request(
            "initialize",
            json!({
                "protocolVersion": "2025-06-18", "capabilities": {},
                "clientInfo": {"name": "test", "version": "0"},
            }),
        );
        server.send(json!({"jsonrpc": "2.0", "method": "notifications/initialized"}));
        let instructions = initialized["result"]["instructions"].as_str().unwrap();
        assert!(instructions.contains("Privacy mode"), "{instructions}");
        (child, server)
    }

    fn send(&mut self, message: Value) {
        writeln!(self.input, "{message}").unwrap();
        self.input.flush().unwrap();
    }

    fn request(&mut self, method: &str, params: Value) -> Value {
        let id = self.next_id;
        self.next_id += 1;
        self.send(json!({"jsonrpc": "2.0", "id": id, "method": method, "params": params}));
        let mut line = String::new();
        self.output.read_line(&mut line).unwrap();
        serde_json::from_str(&line).unwrap()
    }

    /// The tool's reply text (JSON or an error message).
    fn call(&mut self, tool: &str, arguments: Value) -> String {
        let reply = self.request("tools/call", json!({"name": tool, "arguments": arguments}));
        reply["result"]["content"][0]["text"]
            .as_str()
            .unwrap_or_else(|| panic!("{reply}"))
            .to_string()
    }
}

/// Every name in the project a reply could leak.
fn real_names(project: &Path) -> Vec<String> {
    let names = Project::open(project).unwrap().names().unwrap();
    [
        names.biosamples,
        names.categories,
        names.plexes,
        names.runs,
        names.spectra_files,
        names.sources,
        names.accessions,
    ]
    .concat()
}

#[test]
fn replies_carry_aliases_and_aliases_reach_the_real_names() {
    let dir = tempfile::tempdir().unwrap();
    let project = dir.path().join("demo.tuff");
    let (mut child, mut server) = Server::start(&["--private=proteins"]);
    let p = json!(project);

    server.call("create_demo", json!({"project": p}));
    let names = real_names(&project);
    // The demo: four samples in two categories, 150 proteins and decoys.
    assert!(names.len() > 150, "{}", names.len());

    let summary: Value =
        serde_json::from_str(&server.call("project_summary", json!({"project": p}))).unwrap();
    let samples: Vec<(&str, &str)> = summary["biosamples"]
        .as_array()
        .unwrap()
        .iter()
        .map(|b| (b["name"].as_str().unwrap(), b["category"].as_str().unwrap()))
        .collect();
    assert_eq!(
        samples,
        [
            ("sample-1", "group-1"),
            ("sample-2", "group-1"),
            ("sample-3", "group-2"),
            ("sample-4", "group-2")
        ]
    );

    let mut replies = vec![summary.to_string()];
    let groups: Value = serde_json::from_str(
        &server.call("list_protein_groups", json!({"project": p, "limit": 500})),
    )
    .unwrap();
    let first = &groups["groups"][0];
    assert!(first["accession"].as_str().unwrap().starts_with("protein-"));
    assert!(first["description"].is_null());
    // An alias filters exactly.
    let alias = first["accession"].as_str().unwrap();
    let filtered: Value = serde_json::from_str(&server.call(
        "list_protein_groups",
        json!({"project": p, "filter": alias}),
    ))
    .unwrap();
    assert_eq!(filtered["total"], 1, "{filtered}");
    replies.push(groups.to_string());
    replies.push(server.call(
        "protein_group",
        json!({"project": p, "group_id": first["group_id"]}),
    ));
    replies.push(server.call(
        "compare_categories",
        json!({"project": p, "baseline": "group-1", "versus": "group-2"}),
    ));
    replies.push(server.call("compare_all_categories", json!({"project": p})));
    let qc = server.call("sample_qc", json!({"project": p}));
    assert!(qc.contains("\"name\":\"sample-1\""), "{qc}");
    replies.push(server.call("methods_text", json!({"project": p})));
    // An error quoting a name is rewritten too.
    replies.push(server.call(
        "compare_categories",
        json!({"project": p, "baseline": "group-1", "versus": "nobody"}),
    ));
    replies.push(server.call(
        "set_plex",
        json!({
            "project": p,
            "name": "secret_plex",
            "label": "tmt6",
            "runs": ["sample-1"],
            "channels": {"126": "sample-1", "127": "patient_channel"},
        }),
    ));
    replies.push(server.call(
        "quantify_reporters",
        json!({"project": p, "plex": "plex-1"}),
    ));
    replies.push(server.call("remove_plex", json!({"project": p, "name": "plex-1"})));

    let names = real_names(&project);
    for reply in replies.iter().chain([&qc]) {
        for name in names.iter().filter(|name| name.len() >= 3) {
            assert!(!reply.contains(name.as_str()), "{name} leaked in {reply}");
        }
    }
    assert!(
        replies[3].contains("\"baseline\":\"group-1\""),
        "{}",
        replies[3]
    );
    assert!(
        replies[7].contains("\"name\":\"plex-1\"")
            && replies[7].contains("\"biosample\":\"sample-5\""),
        "{}",
        replies[7]
    );
    assert!(
        replies[8].contains("\"plexes\":[{") && replies[8].contains("\"name\":\"plex-1\""),
        "{}",
        replies[8]
    );
    assert!(
        replies[9].contains("\"plexes\":[]") && replies[9].contains("\"name\":\"sample-5\""),
        "{}",
        replies[9]
    );
    assert!(
        Project::open(&project)
            .unwrap()
            .plexes()
            .unwrap()
            .is_empty()
    );

    // An alias sent back groups the real sample under the real category.
    let before = Project::open(&project).unwrap().names().unwrap();
    server.call(
        "set_categories",
        json!({"project": p, "assignments": {"sample-1": "group-2"}}),
    );
    let list = Project::open(&project).unwrap().protein_list().unwrap();
    let moved = list
        .biosamples
        .iter()
        .find(|b| b.name == before.biosamples[0])
        .unwrap();
    assert_eq!(
        moved.category.as_deref(),
        Some(before.categories[1].as_str())
    );

    drop(server);
    child.wait().unwrap();
}

#[test]
fn an_unknown_privacy_level_stops_the_server() {
    let status = Command::new(env!("CARGO_BIN_EXE_tuffold-mcp"))
        .arg("--private=some")
        .stdin(Stdio::null())
        .stdout(Stdio::null())
        .stderr(Stdio::null())
        .status()
        .unwrap();
    assert_eq!(status.code(), Some(2));
}
