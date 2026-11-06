# DARLIN reference contract

DARLIN does **not** use an external genome FASTA or STAR index. The lineage-template reference sequences and primer boundaries are wet-lab validated scientific constants embedded in `core/lemmus/darlin.py`.

Supported templates and loci:

| Template | Locus | Internal unedited reference |
| --- | --- | --- |
| `cCARLIN` | `CA` | `CA_CARLIN` |
| `Rosa`, `Rosa_v2` | `RA` | `RA_CARLIN` |
| `Tigre`, `Tigre_2022`, `Tigre_2022_v2` | `TA` | `TA_CARLIN` |

The project `00_config/Barcode_Map.csv` is the experimental barcode definition. The Pipeline generates the RNA/DARLIN cell-barcode whitelist from this file; there is no separate DARLIN reference file to download.

Do not replace the embedded CARLIN sequences with genome FASTA. Changing these constants is a scientific-method change and requires regression validation.

## Related documents

- [Embedded implementation](../core/lemmus/darlin.py)
- [Extraction and optional RNA rescue](https://github.com/Felix-owo/Lemmus_dev/blob/main/docs/contracts/darlin-design.md)
- [Scientific and validation contract](https://github.com/Felix-owo/Lemmus_dev/blob/main/docs/contracts/darlin.md#5-darlin)
