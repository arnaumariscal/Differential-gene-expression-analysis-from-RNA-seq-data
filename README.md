# **Differential-gene-expression-analysis-from-RNA-seq-data**

## Overview:

This project involves a **differential gene expression analysis** of blood samples from COVID-19 patients, bacterial pneumonia patients, and healthy individuals.
It was developed as part of the *Omics Data Analysis* module in my Master’s in Bioinformatics and Biostatistics.

The dataset was downloaded from the Gene Expression Omnibus database (GSE161731) and originates from the study: McClain et al., Nature Communications, 2021 (DOI: DOI: 10.1038/s41467-021-21289-y)
The analysis was conducted entirely in **R** using **Bioconductor** packages such as *SummarizedExperiment, edgeR* and *limma*.

Results show that bacterial infections triggered a much stronger and more diverse transcriptomic response than COVID-19 infections. In COVID-19 patients, overexpressed genes were primarily linked to cell cycle regulation and oligodendrocyte differentiation.

## Methodology:

 - RNA-seq data import and **SummarizedExperiment** creation
 - Metadata cleaning and cohort selection
 - Data preprocessing and transformation: Filtering low-expression genes and applying TMM normalization
 - Exploratory analysis (heatmaps, PCA, MDS)
 - Identification of **confounding variables** (age, group)
 - Design and contrast matrices
 - Differential expression analysis with **voom + limma**
 - Over-representation analysis using Gene Ontology

## Key results:

Differential expression profiles
 - COVID vs Healthy. Overexpression of genes in healthy individuals, reduced activation in COVID-19 patients
 - Bacterial vs Healthy. Two clusters: one overexpressed and one underexpressed in bacterial samples

   ***
 - 7 times higher number of differentially expressed genes in the bacterial comparison (more intense and diverse transcriptomic response??)

Biological significance analysis (over-representation analysis)
 - overexpressed genes in COVID-19 patients are involved in cell proliferation and oligodendrocyte differentiation
 
