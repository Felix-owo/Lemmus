# 固定45S竞争序列

| RNA reference | RefSeq来源 | 长度 | 固定FASTA SHA256 |
| --- | --- | ---: | --- |
| hg38 | [人 RNA45SN2 NR_146144.1](https://www.ncbi.nlm.nih.gov/nuccore/NR_146144.1) | 13,315 nt | `d5a3d20f8f84fe5770096d55644d6b7fc8d5116e64343b1259805df90c17ebcd` |
| mm10 | [小鼠 Rn45s NR_046233.2](https://www.ncbi.nlm.nih.gov/nuccore/NR_046233.2) | 13,400 nt | `1bea30ddff317d9a76f5e3827082dc4446d8f365115d9dedd1d6fdbb274875ee` |

来源为NCBI RefSeq，2026-11-02固定；序列保持原方向、转为大写DNA字母，FASTA行宽80，唯一header为 `RRNA45S_COMPETITOR`。Doctor校验固定哈希，并以 `simpleaf index --unspliced` 加入expanded reference。这里的U是参考建模标签，不能据此证明所有匹配reads都来自生物学未剪接pre-rRNA。

该feature代表45S家族竞争序列，不是单拷贝基因。生产PG拒绝仍存在的多gene歧义，不强行EM分配。技术计数保存在每个barcode的 `rna45s_competitor_counts` 和 `uns/rna_quantification`，从生物表达矩阵、emptyDrops、HVG及归一化分母排除。0技术计数不证明没有rRNA污染：无法唯一归属的分子可能已被PG拒绝。
