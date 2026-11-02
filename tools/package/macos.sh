#!/bin/sh
# Build Tuffold.app and a DMG holding it, from this checkout, on a Mac.
#
#     tools/package/macos.sh            # -> target/package/Tuffold-<version>-macos-<arch>.dmg
#
# The app carries the desktop app and, beside it in Contents/MacOS, the
# `tuffold` command and the `tuffold-mcp` server. Only Apple's own tools
# are used (swift, sips, iconutil, codesign, hdiutil).
#
# The bundle is signed ad hoc, not with a Developer ID, and not notarized:
# macOS will say it can't check the app for malicious software. Open it the
# first time with right-click > Open (or System Settings > Privacy &
# Security > Open Anyway). Signing for distribution needs an Apple
# Developer account; see docs/install.md.
set -eu

root=$(cd "$(dirname "$0")/../.." && pwd)
cd "$root"

version=$(sed -n 's/^version = "\(.*\)"/\1/p' Cargo.toml | head -n 1)
arch=$(uname -m)
out="$root/target/package"
app="$out/Tuffold.app"
dmg="$out/Tuffold-$version-macos-$arch.dmg"

echo "building Tuffold $version for $arch"
cargo build --release --locked -p tuffold-app -p tuffold-cli -p tuffold-mcp

rm -rf "$app" "$dmg" "$out/dmg" "$out/icon.iconset"
mkdir -p "$app/Contents/MacOS" "$app/Contents/Resources" "$out/icon.iconset"
cp target/release/tuffold-app target/release/tuffold target/release/tuffold-mcp \
    "$app/Contents/MacOS/"
cp LICENSE "$app/Contents/Resources/LICENSE.txt"
cp NOTICE "$app/Contents/Resources/NOTICE.txt"
python3 tools/package/third_party.py --check "$app/Contents/Resources/THIRD-PARTY-LICENSES.txt"
cp app/assets/fonts/OFL.txt "$app/Contents/Resources/IBM-Plex-OFL.txt"

# The icon: the SVG at every size an .icns holds.
swift tools/package/render_icon.swift app/assets/app-icon.svg "$out/icon-1024.png" 1024
for size in 16 32 128 256 512; do
    sips -z "$size" "$size" "$out/icon-1024.png" \
        --out "$out/icon.iconset/icon_${size}x${size}.png" >/dev/null
    double=$((size * 2))
    sips -z "$double" "$double" "$out/icon-1024.png" \
        --out "$out/icon.iconset/icon_${size}x${size}@2x.png" >/dev/null
done
iconutil -c icns "$out/icon.iconset" -o "$app/Contents/Resources/Tuffold.icns"
rm -rf "$out/icon.iconset" "$out/icon-1024.png"

cat > "$app/Contents/Info.plist" <<PLIST
<?xml version="1.0" encoding="UTF-8"?>
<!DOCTYPE plist PUBLIC "-//Apple//DTD PLIST 1.0//EN" "http://www.apple.com/DTDs/PropertyList-1.0.dtd">
<plist version="1.0">
<dict>
    <key>CFBundleDevelopmentRegion</key><string>en</string>
    <key>CFBundleDisplayName</key><string>Tuffold</string>
    <key>CFBundleExecutable</key><string>tuffold-app</string>
    <key>CFBundleIconFile</key><string>Tuffold</string>
    <key>CFBundleIdentifier</key><string>io.github.tuff-industries.tuffold</string>
    <key>CFBundleInfoDictionaryVersion</key><string>6.0</string>
    <key>CFBundleName</key><string>Tuffold</string>
    <key>CFBundlePackageType</key><string>APPL</string>
    <key>CFBundleShortVersionString</key><string>$version</string>
    <key>CFBundleVersion</key><string>$version</string>
    <key>LSApplicationCategoryType</key><string>public.app-category.education</string>
    <key>LSMinimumSystemVersion</key><string>13.0</string>
    <key>NSHighResolutionCapable</key><true/>
    <key>NSHumanReadableCopyright</key><string>Copyright 2026 Ansar Salken (TUFF Industries). Apache License 2.0</string>
</dict>
</plist>
PLIST
plutil -lint "$app/Contents/Info.plist" >/dev/null

codesign --force --deep --sign - "$app"
codesign --verify --strict "$app"

mkdir -p "$out/dmg"
cp -R "$app" "$out/dmg/"
ln -s /Applications "$out/dmg/Applications"
hdiutil create -quiet -volname "Tuffold $version" -srcfolder "$out/dmg" -fs HFS+ \
    -format UDZO -ov "$dmg"
rm -rf "$out/dmg"

echo "wrote $dmg ($(du -h "$dmg" | cut -f1))"
