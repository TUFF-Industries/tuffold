#!/bin/sh
# Build a Linux tarball of Tuffold from this checkout.
#
#     tools/package/linux.sh            # -> target/package/tuffold-<version>-linux-<arch>.tar.gz
#
# The tarball holds the desktop app, the `tuffold` command, the
# `tuffold-mcp` server, a desktop entry and icon, and install.sh, which
# copies them under ~/.local (or a prefix given as its argument). The
# desktop app needs Vulkan, xkbcommon, and Wayland or X11 at run time;
# see docs/install.md.
set -eu

root=$(cd "$(dirname "$0")/../.." && pwd)
cd "$root"

version=$(sed -n 's/^version = "\(.*\)"/\1/p' Cargo.toml | head -n 1)
arch=$(uname -m)
name="tuffold-$version-linux-$arch"
out="$root/target/package"
stage="$out/$name"

echo "building Tuffold $version for $arch"
cargo build --release --locked -p tuffold-app -p tuffold-cli -p tuffold-mcp

rm -rf "$stage" "$out/$name.tar.gz"
mkdir -p "$stage/bin" "$stage/share/applications" \
    "$stage/share/icons/hicolor/scalable/apps" "$stage/share/doc/tuffold"
cp target/release/tuffold-app target/release/tuffold target/release/tuffold-mcp "$stage/bin/"
cp app/assets/app-icon.svg "$stage/share/icons/hicolor/scalable/apps/tuffold.svg"
cp LICENSE "$stage/share/doc/tuffold/LICENSE"
cp NOTICE "$stage/share/doc/tuffold/NOTICE"
python3 tools/package/third_party.py --check "$stage/share/doc/tuffold/THIRD-PARTY-LICENSES"
cp app/assets/fonts/OFL.txt "$stage/share/doc/tuffold/IBM-Plex-OFL.txt"

cat > "$stage/share/applications/tuffold.desktop" <<'DESKTOP'
[Desktop Entry]
Type=Application
Name=Tuffold
GenericName=Proteomics validation
Comment=Validate and view MS/MS proteomics results
Exec=tuffold-app
Icon=tuffold
Terminal=false
Categories=Science;Biology;Chemistry;
DESKTOP

cat > "$stage/install.sh" <<'INSTALL'
#!/bin/sh
# Install Tuffold under PREFIX (default ~/.local): programs in PREFIX/bin,
# the desktop entry and icon under PREFIX/share. No root needed.
set -eu
here=$(cd "$(dirname "$0")" && pwd)
prefix=${1:-"$HOME/.local"}
mkdir -p "$prefix/bin" "$prefix/share/applications" \
    "$prefix/share/icons/hicolor/scalable/apps" "$prefix/share/doc/tuffold"
cp "$here/bin/tuffold-app" "$here/bin/tuffold" "$here/bin/tuffold-mcp" "$prefix/bin/"
sed "s|^Exec=tuffold-app|Exec=$prefix/bin/tuffold-app|" \
    "$here/share/applications/tuffold.desktop" > "$prefix/share/applications/tuffold.desktop"
cp "$here/share/icons/hicolor/scalable/apps/tuffold.svg" \
    "$prefix/share/icons/hicolor/scalable/apps/tuffold.svg"
cp "$here/share/doc/tuffold/"* "$prefix/share/doc/tuffold/"
command -v update-desktop-database >/dev/null 2>&1 &&
    update-desktop-database "$prefix/share/applications" 2>/dev/null || true
echo "installed into $prefix (add $prefix/bin to PATH for the tuffold command)"
INSTALL
chmod +x "$stage/install.sh"

tar -C "$out" -czf "$out/$name.tar.gz" "$name"
rm -rf "$stage"
echo "wrote $out/$name.tar.gz ($(du -h "$out/$name.tar.gz" | cut -f1))"
