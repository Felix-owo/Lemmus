#!/usr/bin/env bash
set -euo pipefail
[[ $# == 1 && "$1" == /* && ! -e "$1" ]] || { echo 'Usage: build.sh /absolute/nonexistent-build-directory' >&2; exit 2; }
bundle="$(cd "$(dirname "$0")" && pwd -P)"
work="$1"
mkdir "$work"
python="${PYTHON:-python3}"
curl --fail --location --proto '=https' --tlsv1.2 \
  'https://codeload.github.com/COMBINE-lab/alevin-fry/tar.gz/refs/tags/v0.18.3' \
  --output "$work/upstream.tar.gz"
"$python" - "$bundle" "$work" <<'PY'
import hashlib,json,sys,tarfile
from pathlib import Path
bundle,work=map(Path,sys.argv[1:]);data=json.loads((bundle/'distribution.json').read_text())
for p,expected in [(work/'upstream.tar.gz',data['source_archive_sha256']),(bundle/'canonical-gene-ec.patch',data['patch_sha256'])]:
    if hashlib.sha256(p.read_bytes()).hexdigest()!=expected:raise ValueError(f'SHA256 mismatch: {p}')
with tarfile.open(work/'upstream.tar.gz') as archive:archive.extractall(work,filter='data')
lock=work/'alevin-fry-0.18.3/Cargo.lock'
if hashlib.sha256(lock.read_bytes()).hexdigest()!=data['cargo_lock_sha256']:raise ValueError('Cargo.lock SHA256 mismatch')
PY
patch --directory "$work/alevin-fry-0.18.3" -p1 < "$bundle/canonical-gene-ec.patch"
[[ "$(rustc --version)" == 'rustc 1.97.1 '* ]] || { echo 'Rust 1.97.1 required' >&2; exit 2; }
export CARGO_HOME="$work/cargo-home" CARGO_TARGET_DIR="$work/target"
export CC=/usr/bin/cc
export RUSTFLAGS="--remap-path-prefix=$work=/build/alevin-fry"
export CFLAGS="-ffile-prefix-map=$work=/build/alevin-fry"
if [[ "$(uname -s)" == Linux && "$(uname -m)" == x86_64 ]]; then
  export CARGO_TARGET_X86_64_UNKNOWN_LINUX_GNU_LINKER=/usr/bin/cc
  RUSTFLAGS="-C target-cpu=x86-64 $RUSTFLAGS"
  CFLAGS="-DLIBDEFLATE_ASSEMBLER_DOES_NOT_SUPPORT_AVX_VNNI $CFLAGS"
fi
unset LD_LIBRARY_PATH LD_RUN_PATH
args=(--release --locked)
if [[ -n "${CARGO_VENDOR_DIR:-}" ]]; then
  "$python" - "$work/vendor.toml" "$CARGO_VENDOR_DIR" <<'PY'
import json,sys
from pathlib import Path
Path(sys.argv[1]).write_text('[source.crates-io]\nreplace-with="vendored-sources"\n[source.vendored-sources]\ndirectory='+json.dumps(str(Path(sys.argv[2]).resolve(strict=True)))+'\n')
PY
  args+=(--offline --config "$work/vendor.toml")
fi
cd "$work/alevin-fry-0.18.3"
cargo build "${args[@]}"
cargo test "${args[@]}"
"$work/target/release/alevin-fry" lemmus-contract
"$python" - "$work/target/release/alevin-fry" <<'PY'
from pathlib import Path
import hashlib,sys
p=Path(sys.argv[1]);print(hashlib.sha256(p.read_bytes()).hexdigest(),p)
PY
