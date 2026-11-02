"""Chat with a local Ollama model that can use tuffold-mcp's tools, so no
project data leaves the computer. Python standard library only.

    ollama pull qwen2.5          # any model with tool calling
    cargo build --release -p tuffold-mcp
    python3 tools/local-model/ollama_mcp.py --model qwen2.5 \\
        --server target/release/tuffold-mcp --private names \\
        "Create the demo project at /tmp/demo.tuff and summarize it."

With no question it asks for one line at a time (empty line or Ctrl-D
quits). --private names (or proteins) is passed to tuffold-mcp, which
then shows the model aliases instead of sample and file names; with a
local model that is belt and braces. Every tool call is printed, so you see
what the model asked for and what it got.
"""

import argparse
import json
import subprocess
import sys
import urllib.request

MAX_ROUNDS = 12  # tool rounds per question, so a confused model stops


class McpServer:
    """tuffold-mcp over stdio: one JSON-RPC message per line."""

    def __init__(self, command):
        self.process = subprocess.Popen(
            command, stdin=subprocess.PIPE, stdout=subprocess.PIPE, text=True
        )
        self.next_id = 1
        initialized = self.request("initialize", {
            "protocolVersion": "2025-06-18", "capabilities": {},
            "clientInfo": {"name": "ollama_mcp", "version": "1"},
        })
        self.instructions = initialized.get("instructions", "")
        self.send({"jsonrpc": "2.0", "method": "notifications/initialized"})
        self.tools = self.request("tools/list", {})["tools"]

    def send(self, message):
        self.process.stdin.write(json.dumps(message) + "\n")
        self.process.stdin.flush()

    def request(self, method, params):
        request_id = self.next_id
        self.next_id += 1
        self.send({"jsonrpc": "2.0", "id": request_id, "method": method, "params": params})
        while True:
            line = self.process.stdout.readline()
            if not line:
                raise RuntimeError("tuffold-mcp stopped")
            reply = json.loads(line)
            if reply.get("id") == request_id:
                if "error" in reply:
                    raise RuntimeError(reply["error"].get("message", reply["error"]))
                return reply["result"]

    def call(self, name, arguments):
        result = self.request("tools/call", {"name": name, "arguments": arguments})
        text = "".join(block.get("text", "") for block in result.get("content", []))
        return ("error: " if result.get("isError") else "") + text

    def close(self):
        self.process.stdin.close()
        self.process.wait()


def chat(host, model, messages, tools):
    body = json.dumps({"model": model, "messages": messages, "tools": tools, "stream": False})
    request = urllib.request.Request(
        f"{host}/api/chat", data=body.encode(), headers={"Content-Type": "application/json"}
    )
    with urllib.request.urlopen(request) as response:
        return json.load(response)["message"]


def answer(server, host, model, messages, tools):
    for _ in range(MAX_ROUNDS):
        message = chat(host, model, messages, tools)
        messages.append(message)
        calls = message.get("tool_calls") or []
        if not calls:
            return message.get("content", "")
        for call in calls:
            name = call["function"]["name"]
            arguments = call["function"].get("arguments") or {}
            if isinstance(arguments, str):
                arguments = json.loads(arguments)
            print(f"  [tool] {name} {json.dumps(arguments)}", file=sys.stderr)
            result = server.call(name, arguments)
            print(f"  [tool] -> {result[:300]}{'…' if len(result) > 300 else ''}", file=sys.stderr)
            messages.append({"role": "tool", "content": result, "tool_name": name})
    return f"(stopped after {MAX_ROUNDS} tool rounds)"


def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--model", default="qwen2.5")
    parser.add_argument("--server", default="target/release/tuffold-mcp")
    parser.add_argument("--host", default="http://127.0.0.1:11434")
    parser.add_argument("--private", choices=["names", "proteins"],
                        help="pass --private=names or --private=proteins to tuffold-mcp")
    parser.add_argument("question", nargs="*")
    args = parser.parse_args()
    if not args.host.startswith(("http://127.0.0.1", "http://localhost", "http://[::1]")):
        sys.exit("--host must be this computer: the point is that nothing leaves it")

    command = [args.server] + ([f"--private={args.private}"] if args.private else [])
    server = McpServer(command)
    tools = [{"type": "function", "function": {
        "name": tool["name"], "description": tool.get("description", ""),
        "parameters": tool["inputSchema"],
    }} for tool in server.tools]
    messages = [{"role": "system", "content": server.instructions}]
    try:
        questions = [" ".join(args.question)] if args.question else None
        while True:
            if questions is not None:
                if not questions:
                    break
                question = questions.pop()
            else:
                try:
                    question = input("> ").strip()
                except EOFError:
                    break
                if not question:
                    break
            messages.append({"role": "user", "content": question})
            print(answer(server, args.host, args.model, messages, tools))
    finally:
        server.close()


if __name__ == "__main__":
    main()
