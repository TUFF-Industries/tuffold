# Installing Tuffold

Tuffold is one desktop app plus two command-line programs: `tuffold` (every
workflow, scriptable) and `tuffold-mcp` (the server AI assistants connect
to; see [mcp.md](mcp.md)). Installers are built from this repository by
`tools/package/` and by the *Package* workflow in GitHub Actions; there is
no signed release yet.

## macOS (13 Ventura or later)

1. Open `Tuffold-<version>-macos-<arch>.dmg` and drag **Tuffold** to
   **Applications**. The DMG runs on the kind of Mac that built it (`arm64`,
   Apple silicon, from GitHub Actions; build on an Intel Mac for `x86_64`).
2. The app is signed ad hoc, not with an Apple Developer ID, and not
   notarized, so the first time macOS says it can't check it for malicious
   software and won't open it. Allow it once: on macOS 15 and later, try to
   open it, then go to System Settings › Privacy & Security and click
   *Open Anyway*; on 13 and 14, **right-click Tuffold › Open** and confirm.
   (Or, in a terminal, `xattr -dr com.apple.quarantine
   /Applications/Tuffold.app`.) It opens normally after that.
3. The command-line programs are inside the app. To use them from a
   terminal:

   ```sh
   sudo ln -s /Applications/Tuffold.app/Contents/MacOS/tuffold /usr/local/bin/tuffold
   sudo ln -s /Applications/Tuffold.app/Contents/MacOS/tuffold-mcp /usr/local/bin/tuffold-mcp
   ```

Build the DMG yourself on a Mac with Rust and Xcode's command-line tools:
`tools/package/macos.sh` (writes `target/package/`).

## Linux (x86_64 or arm64)

```sh
tar -xzf tuffold-<version>-linux-<arch>.tar.gz
tuffold-<version>-linux-<arch>/install.sh          # into ~/.local; or give a prefix
```

The installer copies `tuffold-app`, `tuffold`, and `tuffold-mcp` to
`~/.local/bin`, and a desktop entry and icon under `~/.local/share`, so
*Tuffold* appears in the application menu. Add `~/.local/bin` to `PATH`
for the commands. Remove it by deleting those three programs,
`~/.local/share/applications/tuffold.desktop`, and
`~/.local/share/icons/hicolor/scalable/apps/tuffold.svg`.

The desktop app draws with Vulkan and needs, on Debian or Ubuntu:

```sh
sudo apt install libvulkan1 mesa-vulkan-drivers libxkbcommon-x11-0 libwayland-client0 libfontconfig1
```

(`vulkan-loader`, `mesa-vulkan-drivers`, `libxkbcommon-x11`, `wayland`, and
`fontconfig` on Fedora.) The command-line programs need none of these and
run on servers without a display.

Build the tarball yourself with the GPUI libraries listed in
`.github/workflows/ci.yml`: `tools/package/linux.sh`.

## Windows

Not packaged, and not planned. The command-line programs may build with
`cargo build --release -p tuffold-cli -p tuffold-mcp`, but nothing tests
that on Windows.

## Signing

Signed, notarized macOS builds need a paid Apple Developer account
(`codesign --sign "Developer ID Application: …"`, then `xcrun notarytool`).
That is parked, so the packaging script signs ad hoc.
