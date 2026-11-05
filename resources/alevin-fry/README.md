# Lemmus fixed alevin-fry distribution

This is **Lemmus-patched alevin-fry 0.18.3**, based on upstream commit
`ad05742d274230f9141b2eabeebd1c1b31692199`, under the included BSD-3-Clause
[license](LICENSE). It is not an unmodified upstream binary.

[distribution.json](distribution.json) fixes the upstream archive, Cargo.lock,
patch and platform binary SHA256 values. [canonical-gene-ec.patch](canonical-gene-ec.patch)
contains the complete source changes. macOS arm64 and Linux x86_64 binaries are
provided; Linux was built on glibc 2.28 with a generic x86-64 target and system C
linker. Other platforms are rejected rather than silently using unpatched code.

The patch numbers projected gene equivalence classes lexicographically by gene
labels. This removes read-order dependence from PUG vertex identities while
retaining UMI frequencies, gene sets, barcode-seeded hashing, the PG algorithm and
residual ambiguity rejection. `alevin-fry lemmus-contract` prints
`lexicographic_gene_labels_v1`; quant.json carries the same contract. The native
version remains 0.18.3 because upstream collate metadata requires that version
format. Deterministic resolution does not establish the biological source of
intrinsically ambiguous molecules.

Doctor verifies the selected binary SHA256, performs a real five-read
forward/reverse PG smoke, and atomically installs into a candidate rules prefix.
An existing rules environment is cloned before publication; active jobs retain
their old physical prefix. Repeated checks execute only the lightweight contract
query. Runtime preparation never compiles Rust or downloads native source.

## Source build

Use Python 3.12+, curl, patch, a C compiler and Rust/Cargo 1.97.1:

```sh
bash resources/alevin-fry/build.sh /absolute/path/to/new-build-directory
```

The destination must not exist. The script verifies the upstream archive,
Cargo.lock and patch, applies the patch, builds with locked dependencies and runs
the native tests. For an already verified vendor cache, set `CARGO_VENDOR_DIR` to
its absolute directory; the script then builds offline. All mutable Cargo state
and build outputs stay under the destination. Review the printed artifact hash
before replacing any distributed artifact or updating the manifest; the build
script never installs into a Pipeline environment. Compiler/platform differences
can produce different executable hashes; scientific and portability validation
are required for a replacement artifact.
