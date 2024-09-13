# stage_2-task-
**Gene Expression Analysis of a Glioblastoma dataset**

Here, we are visualizing and interpreting a gene expression dataset of glioblastoma samples generating heatmaps and performing downstream functional enrichment analysis of a set of deregulated genes.

**Software used for analysis**

R studio 4.4.0

**Workflow**

* download the following dataset of top 500+ differentially expressed genes under different conditions of glioblastoma samples: [https://raw.githubusercontent.com/HackBio-Internship/public\_datasets/main/Cancer2024/glioblastoma.csv](https://raw.githubusercontent.com/HackBio-Internship/public_datasets/main/Cancer2024/glioblastoma.csv)

* generation of heatmaps of the entire dataset using a sequential and diverging color palette

* generation of heatmaps clustering by (1) genes alone; (2) samples alone; (3) both genes and samples together

* subsetting of genes that are up- and down-regulated

* functional enrichment analysis using ShinyGO: [http://bioinformatics.sdstate.edu/go/](http://bioinformatics.sdstate.edu/go/)

