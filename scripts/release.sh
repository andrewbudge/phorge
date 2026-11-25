#!/usr/bin/env bash
# Build a versioned, installable phorge tarball.
#
# Native build:
#   bash scripts/release.sh
#
# Cross build (requires a matching toolchain and zlib/minizip for the target):
#   TARGET_OS=windows TARGET_ARCH=x86_64 \
#   TARGET_CC=x86_64-w64-mingw32-gcc \
#   TARGET_CFLAGS="-std=c11 -O2 -D_POSIX_C_SOURCE=200809L" \
#   TARGET_LDLIBS="-lz -lminizip" bash scripts/release.sh
set -euo pipefail

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
DIST="$ROOT/dist"
VERSION="$(sed -n 's/^#define PHORGE_VERSION "\([^"]*\)"$/\1/p' "$ROOT/src/version.h" | head -1)"
VERSION="${VERSION:-dev}"

OS="${TARGET_OS:-$(uname -s | tr '[:upper:]' '[:lower:]')}"
ARCH="${TARGET_ARCH:-$(uname -m)}"
NAME="phorge-${VERSION}-${OS}-${ARCH}"
BUILD="${TARGET_BUILD_DIR:-$ROOT/build/release-${OS}-${ARCH}}"
CC_BIN="${TARGET_CC:-${CC:-cc}}"
CFLAGS="${TARGET_CFLAGS:--std=c11 -O2 -Wall -Wextra -D_POSIX_C_SOURCE=200809L}"
LDLIBS="${TARGET_LDLIBS:--lm -lminizip -lz}"

# macOS ships libz but not libminizip headers; Homebrew provides them.
if [ "$(uname -s)" = "Darwin" ] && [ -z "${TARGET_CFLAGS:-}" ] && [ -z "${TARGET_LDLIBS:-}" ]; then
    for prefix in /opt/homebrew /usr/local; do
        if [ -f "$prefix/include/minizip/unzip.h" ]; then
            CFLAGS="$CFLAGS -I$prefix/include"
            LDLIBS="$LDLIBS -L$prefix/lib -Wl,-rpath,$prefix/lib"
            break
        fi
    done
fi

DEPS="Runtime dependencies: zlib (libz) and minizip (libminizip)."
if [ "$OS" = "darwin" ]; then
    DEPS="$DEPS On macOS: brew install minizip zlib."
fi

mkdir -p "$BUILD" "$DIST"
echo "Building $NAME with $CC_BIN"
# shellcheck disable=SC2086
"$CC_BIN" $CFLAGS -I"$ROOT/src" "$ROOT"/src/*.c $LDLIBS -o "$BUILD/phorge"

STAGE="$DIST/$NAME"
rm -rf "$STAGE"
mkdir -p "$STAGE/bin" "$STAGE/share/doc/phorge"
install -m 755 "$BUILD/phorge" "$STAGE/bin/phorge"
cp "$ROOT/README.md" "$STAGE/share/doc/phorge/README.md"

LICENSE_SRC=""
for candidate in "$ROOT/LICENSE" "$ROOT/rust-legacy/LICENSE" "$ROOT/../phorge/LICENSE"; do
    if [ -f "$candidate" ]; then LICENSE_SRC="$candidate"; break; fi
done
if [ -n "$LICENSE_SRC" ]; then
    cp "$LICENSE_SRC" "$STAGE/share/doc/phorge/LICENSE"
else
    echo "warning: no LICENSE found; tarball will not contain one" >&2
fi

cat > "$STAGE/share/doc/phorge/INSTALL.md" <<EOF
# Install phorge $VERSION

Extract the archive and copy the binary onto your PATH:

    tar -xzf $NAME.tar.gz
    cp $NAME/bin/phorge /usr/local/bin/

Or install to a user prefix:

    cp $NAME/bin/phorge ~/.local/bin/

${DEPS} The align and extract commands additionally need MAFFT/MUSCLE and
MMseqs2 on PATH.
EOF

TARBALL="$DIST/$NAME.tar.gz"
tar -czf "$TARBALL" -C "$DIST" "$NAME"
(cd "$DIST" && sha256sum "$(basename "$TARBALL")" > "$NAME.tar.gz.sha256")

echo "Wrote $TARBALL"
echo "Wrote $TARBALL.sha256"
