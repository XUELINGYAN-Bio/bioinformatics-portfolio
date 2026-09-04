# RNA-seq Demo Analysis Summary

## Scope

This summary reports the reproducible outputs from a small simulated RNA-seq count matrix. The dataset is designed for workflow practice, not for making biological claims.

## Input design

- Genes tested: 36
- Samples: 3 control and 3 stress replicates
- Model: DESeq2 design `~ condition`, contrast `stress / control`
- DEG threshold: adjusted p-value < 0.05 and absolute log2 fold change >= 1

## Key results

- Total teaching DEGs: 15
- Up-regulated under stress: 8
- Down-regulated under stress: 7
- Significant simulated gene sets: 2
- Mean within-condition sample correlation: 0.997
- Mean between-condition sample correlation: 0.00807

## Top up-regulated genes

| Gene | log2FC | adjusted p-value | Status |
|---|---:|---:|---|
| Gene003 | 3.52 | 2.28e-295 | Up |
| Gene005 | 2.95 | 7.98e-255 | Up |
| Gene002 | 3.85 | 9.91e-210 | Up |
| Gene007 | 3.02 | 1.11e-157 | Up |
| Gene006 | 3.57 | 7.88e-136 | Up |
| Gene001 | 4.26 | 3.23e-113 | Up |
| Gene008 | 3.69 | 2.28e-99 | Up |
| Gene004 | 3.97 | 1.85e-77 | Up |

## Top down-regulated genes

| Gene | log2FC | adjusted p-value | Status |
|---|---:|---:|---|
| Gene009 | -4.35 | 2.59e-236 | Down |
| Gene012 | -3.13 | 2.11e-209 | Down |
| Gene010 | -3.92 | 2.10e-160 | Down |
| Gene014 | -3.02 | 2.89e-140 | Down |
| Gene011 | -3.52 | 9.92e-99 | Down |
| Gene015 | -3.15 | 3.42e-86 | Down |
| Gene013 | -3.19 | 5.94e-61 | Down |

## Gene-set enrichment

| Gene set | DEG overlap | Set size | adjusted p-value |
|---|---:|---:|---:|
| Stress_response | 8 | 8 | 1.06e-03 |
| Photosynthesis | 7 | 8 | 1.17e-02 |
| Hormone_signaling | 3 | 6 | 8.21e-01 |
| Cell_cycle | 0 | 7 | 1.00e+00 |
| Primary_metabolism | 0 | 8 | 1.00e+00 |

## Interpretation

- PCA separates control and stress samples mainly along PC1, showing that the simulated treatment effect dominates sample variation.
- The volcano plot shows a balanced set of strongly up- and down-regulated teaching genes.
- The top-gene heatmap separates the two conditions and demonstrates how expression patterns can be reviewed across biological replicates.
- The enrichment demo recovers the expected simulated stress-response and photosynthesis gene sets, illustrating over-representation analysis logic.

## Quality-control checks added

- `library_size_barplot.png` checks whether samples have comparable total counts before normalization.
- `sample_correlation_heatmap.png` checks whether biological replicates are more similar to each other than to the other condition.

## Limitations and next steps

- Replace simulated genes with a public plant RNA-seq dataset from GEO or SRA.
- Add raw-read QC with FastQC/MultiQC and transcript quantification with Salmon or STAR + featureCounts.
- Use real gene identifiers and organism-specific GO/KEGG annotation.
- Connect candidate genes back to literature evidence before making biological conclusions.
