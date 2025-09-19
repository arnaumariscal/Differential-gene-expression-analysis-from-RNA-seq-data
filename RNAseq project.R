## ----setup, include=FALSE----------------------------------------------------------------
# Document configuration
knitr::opts_chunk$set(
  echo = TRUE, 
  message = FALSE, 
  warning = FALSE)


## ----paquets, include=FALSE--------------------------------------------------------------
# Installation of BiocManager if not already installed
if(!require(BiocManager)){
      install.packages("BiocManager", dep=TRUE)
}

installifnot <- function (paquet, BioC=TRUE){
  if(BioC){
    if(!require(paquet, character.only=TRUE)){
      BiocManager::install(paquet)
    }
  }else{
    if(!require(paquet, character.only=TRUE)){
      install.packages(paquet, dep=TRUE)
    }
  }
}

# If needed, install the following Bioconductor packages
installifnot("SummarizedExperiment")
installifnot("EnsDb.Hsapiens.v86")
installifnot("org.Hs.eg.db")
installifnot("clusterProfiler")

# If needed, install the following CRAN packages
installifnot("dplyr", BioC = FALSE)
installifnot("tidyr", BioC = FALSE)
installifnot("stringr", BioC = FALSE)
installifnot("ggplot2", BioC = FALSE)  

# Other useful packages
installifnot("limma")
installifnot("edgeR")
installifnot("ggrepel", BioC = FALSE)
installifnot("pheatmap", BioC = FALSE)


## ----llibreries, include=FALSE-----------------------------------------------------------
# Library loading
library(SummarizedExperiment)
library(EnsDb.Hsapiens.v86)
library(org.Hs.eg.db)
library(clusterProfiler)
library(dplyr)
library(tidyr)
library(readr)
library(stringr)
library(ggplot2)
library(limma)


## ----include=FALSE-----------------------------------------------------------------------
#Directory creation
dir.create("dades", showWarnings = FALSE)

#Import the count matrix (called "comptatges") and the metadata (named "metadades") from the "dades" folder
matriu_comptatges <- "dades/GSE161731_counts.csv.gz"
metadades <- "dades/GSE161731_counts_key.csv.gz"

comptatges <- read.csv(matriu_comptatges, row.names = 1)
metadades <- read.csv(metadades, row.names = 1)


## ----echo=FALSE, results='hold'----------------------------------------------------------
cat("Output 1. Count matrix and metadadata: number of genes and samples\n")
cat("\n• Count matrix (genes x samples): \n")
dim(comptatges)
cat("\n• Metadata (genes x samples): \n")
dim(metadades)



## ----include=FALSE-----------------------------------------------------------------------
# For the count matrix: remove the X from numeric sample names and replace dots with dashes in letter-based sample names

colnames(comptatges) <- sub("^X", "", colnames(comptatges))
colnames(comptatges) <- sub("\\.", "-", colnames(comptatges))

# Find common samples in both the count matrix and the metadata
mostres_comunes <- intersect(colnames(comptatges), rownames(metadades))

# Keep only the common samples in both matrices
comptatges_filtrat <- comptatges[, mostres_comunes]
metadades_filtrat <- metadades[mostres_comunes, ]

dim(comptatges_filtrat)
dim(metadades_filtrat)

# Check that the number of samples is the same in the count matrix and the metadata
all.equal(colnames(comptatges_filtrat), rownames(metadades_filtrat))


## ----include=FALSE-----------------------------------------------------------------------
#Shortname is created: cohort name + sequential number
metadades_filtrat$shortname <- paste(metadades_filtrat$cohort, seq_along(metadades_filtrat$cohort), sep = "_")

#A  color is assigned  to each of the cohorts
metadades_filtrat$color <- NA

for (i in 1:nrow(metadades_filtrat)) {
  if (metadades_filtrat$cohort[i] == "COVID-19") {metadades_filtrat$color[i] <- "pink"}
  else if (metadades_filtrat$cohort[i] == "Bacterial") {metadades_filtrat$color[i] <- "blue"}
  else if (metadades_filtrat$cohort[i] == "healthy") {metadades_filtrat$color[i] <- "green"}
  else {metadades_filtrat$color[i] <- "gray"}}

# Name of the samples = Shortname 
metadades_filtrat$nom_mostra <- rownames(metadades_filtrat)
rownames(metadades_filtrat) <- metadades_filtrat$shortname
metadades_filtrat$shortname <- NULL
colnames(comptatges_filtrat) <- rownames(metadades_filtrat)


## ----include=FALSE-----------------------------------------------------------------------
library(EnsDb.Hsapiens.v86) # (4)
library(GenomicRanges)

# Get the genomic coordinates of all genes
anotacio <- genes(EnsDb.Hsapiens.v86)

# Keep only genes from the count matrix that have genomic coordinates
gens_comuns <- intersect(rownames(comptatges_filtrat), names(anotacio))

comptatges_filtrat <- comptatges_filtrat[gens_comuns, ]
rangs_gen <- anotacio[gens_comuns]

rangs_gen <- anotacio[match(gens_comuns, names(anotacio))]
# Ensure the order and names of the samples of the annotation match the ones in 'comptatges_filtrat'


## ----echo=FALSE, results='hold'----------------------------------------------------------
my_se <- SummarizedExperiment(
  assays = list(counts = as.matrix(comptatges_filtrat)),
  colData = metadades_filtrat,
  rowRanges = rangs_gen
)

cat("Output 2. SummarizedExperiment object: counts + metadata + genomic coordinates \n\n")
show(my_se)


## ----include=FALSE-----------------------------------------------------------------------
metadades_finals <- metadades_filtrat %>%
  # only COVID-19, Bacterial and healthy samples are kept
  dplyr::filter(cohort %in% c("COVID-19", "Bacterial", "healthy")) %>%
  
  # Duplicate samples are eliminated (taking into account the original name, no the given shortname)
  distinct(nom_mostra, .keep_all = TRUE) %>%
  
  # Standarization
  mutate(age = str_replace_all(as.character(age), "[^0-9.]", "")) %>%
  mutate(age = as.numeric(age)) %>%
  mutate(gender = as.factor(gender)) %>%
  mutate(race = as.factor(race)) %>%
  mutate(batch = as.factor(batch)) %>%
  mutate(across(where(is.character), ~ str_replace_all(., c(" " = "_", "-" = "_", "/" = "_")))) %>%
  mutate(across(where(is.factor), ~ factor(str_replace_all(as.character(.), c(" " = "_", "-" = "_", "/" = "_")))))

# Random seed
myseed <- sum(utf8ToInt("arnaumariscalpuig"))
set.seed(myseed)

# Only 75 samples are kept
mostres_aleatories <- sample(rownames(metadades_finals), size = 75)
metadades_finals <- metadades_finals[mostres_aleatories, ]

# update SummarizedExperiment
my_se <- my_se[, colnames(my_se) %in% rownames(metadades_finals)]
metadades_ordenades <- metadades_finals[match(colnames(my_se), rownames(metadades_finals)), ]
colData(my_se) <- DataFrame(metadades_ordenades)
dim(my_se)


## ----echo=FALSE, results='hold'----------------------------------------------------------

cat("Output 3. Number of counts per sample (first 5 samples)\n\n")

colSums(assay(my_se)[, 1:5])


## ----include=FALSE-----------------------------------------------------------------------
library(edgeR)

# Calculate CPMs
comptatges_cpm <- cpm(assay(my_se))

# Gene filtering with low expression
llindar<- comptatges_cpm >0.5
mantenim<- rowSums(llindar) >= 2

# Update SE
my_se <- my_se[mantenim, ]
dim(my_se)


## ----include=FALSE-----------------------------------------------------------------------
dge1 <- DGEList(
  counts = assay(my_se),
  samples = as.data.frame(colData(my_se)),
  group = colData(my_se)$cohort,
  genes = data.frame(nom_gens = rownames(assay(my_se))))


## ----echo=FALSE--------------------------------------------------------------------------
# Log transformation of CPMs
log_comptatges <- cpm(dge1,log=TRUE)

# Sample distribution using boxplot
boxplot(log_comptatges, 
        ylab="Log2-CPM",
        las=2, 
        xlab="", 
        cex.axis=0.8, 
        main="Figure 1. Boxplots of log2-CPM (*unnormalized)")
abline(h=median(log_comptatges), col="green") 


## ----echo=FALSE--------------------------------------------------------------------------
dge_normal1<- calcNormFactors(dge1)

# Boxplot of normalized data
log_comptatges_norm <- cpm(dge_normal1,log=TRUE)

boxplot(log_comptatges_norm, 
        ylab="Log2-CPM",
        las=2, 
        xlab="", 
        cex.axis=0.7, 
        main="Figure 2. Boxplots of log2-CPM (*normalized)")
abline(h=median(log_comptatges_norm), col="green")


## ----include=FALSE-----------------------------------------------------------------------
matriu_dist<- dist(t(log_comptatges_norm))

# Heatmap
library(factoextra)
plot_heatmap<- fviz_dist(matriu_dist) + ggtitle("Figure A1. Heatmap") + theme(plot.title = element_text(size = 16, face = "bold"))

#Dendogram
dendograma <- hclust(matriu_dist)
plot(dendograma, labels = colnames(log_comptatges_norm), main = "Figure A2. Hierarchical clustering of the samples", cex = 0.8)

#MDS plot
color <- dge_normal1$samples$color
plotMDS(log_comptatges_norm, col = color, main = "Figure A3. MDS plot", cex = 0.7)

#PCA
pca <- prcomp(t(log_comptatges_norm), scale. = TRUE)  # log_comptatges_norm needs to be transposed, since columns must represent samples

loads <- round((pca$sdev^2 / sum(pca$sdev^2)) * 100, 1) # percentage of variance explained

# Create a dataframe that combines sample coordinates with metadata
pc_df <- as.data.frame(pca$x)
pc_df$mostra_id <- colnames(log_comptatges_norm)
info_metadades <- as.data.frame(colData(my_se)) 
pc_df <- cbind(pc_df, info_metadades)

# PCA x cohort
pca_cohort<- ggplot(pc_df, aes(x = PC1, y = PC2, color = cohort)) +
  geom_point(size = 2) +
  labs(title = "Figure A4 (A). PCA x cohort",
       x = paste0("PC1 (", loads[1], "%)"),
       y = paste0("PC2 (", loads[2], "%)")) +
  theme_minimal() +
  theme(plot.title = element_text(size = 16, face = "bold"))

# PCA x batch
pca_batch<- ggplot(pc_df, aes(x = PC1, y = PC2, color = batch)) +
  geom_point(size = 2) +
  labs(title = "Figure A4 (B). PCA x batch") +
  theme_minimal() +
  theme(plot.title = element_text(size = 16, face = "bold"))

# PCA x sexe
pca_sex<- ggplot(pc_df, aes(x = PC1, y = PC2, color = gender)) +
  geom_point(size = 2) +
  labs(title = "Figure A4 (C). PCA x sex")  +
  theme_minimal() +
  theme(plot.title = element_text(size = 16, face = "bold"))

# PCA x race
pca_race<- ggplot(pc_df, aes(x = PC1, y = PC2, color = race)) +
  geom_point(size = 2) +
  labs(title = "Figure A4 (D). PCA x race")  +
  theme_minimal() +
  theme(plot.title = element_text(size = 16, face = "bold"))

# PCA x edat
pca_age<- ggplot(pc_df, aes(x = PC1, y = PC2, color = age)) +
  geom_point(size = 2) +
  scale_color_gradient(low = "blue", high = "red") +
  labs(title = "Figure A4 (E). PCA x age")  +
  theme_minimal() +
  theme(plot.title = element_text(size = 16, face = "bold"))



## ----include=FALSE-----------------------------------------------------------------------
outlier <- c("COVID-19_50")

# Outliers are removed
my_se <- my_se[, !colnames(my_se) %in% outlier]

# Update DGElist
dge2 <- DGEList(
  counts = assay(my_se),
  samples = as.data.frame(colData(my_se)),
  group = colData(my_se)$cohort,
  genes = data.frame(nom_gens = rownames(assay(my_se)))
)


## ----include=FALSE-----------------------------------------------------------------------
table1<- table(colData(my_se)$race, colData(my_se)$cohort)
chi1<- chisq.test(table(colData(my_se)$race, colData(my_se)$cohort))


## ----include=FALSE-----------------------------------------------------------------------
table2<- table(colData(my_se)$gender, colData(my_se)$cohort)
chi2<- chisq.test(table(colData(my_se)$gender, colData(my_se)$cohort))


## ----include=FALSE-----------------------------------------------------------------------
table3<- table(colData(my_se)$batch, colData(my_se)$cohort)
chi3<- chisq.test(table(colData(my_se)$batch, colData(my_se)$cohort))


## ----include=FALSE-----------------------------------------------------------------------
boxplot(age ~ cohort, data = as.data.frame(colData(my_se)), main="Figure A5. Age distribution x cohort")
table4<- anova(lm(age ~ cohort, data = as.data.frame(colData(my_se))))


## ----include=FALSE-----------------------------------------------------------------------
# Update DGE list again
dge_normal2<- calcNormFactors(dge2)
log_comptatges_norm2 <- cpm(dge_normal2,log=TRUE)

#Make sure right format of variables
dge_normal2$samples$group <- factor(dge_normal2$samples$group, levels = c("COVID_19", "Bacterial", "healthy"))
dge_normal2$samples$age <- as.numeric(dge_normal2$samples$age)

# Building of the Design matrix taking into account the chosen confounding variables.
matriu_disseny <- model.matrix(~0 + age + group, data = dge_normal2$samples)

# Rewriting of column names (otherwise problems come up later in the code)
colnames(matriu_disseny) <- gsub("group", "", colnames(matriu_disseny))
colnames(matriu_disseny) <- gsub("-", "_", colnames(matriu_disseny))


## ----echo=FALSE, results='hold'----------------------------------------------------------
# Contrast matrix
matriu_contrast <- makeContrasts(
  COVIDvsHealthy = COVID_19 - healthy,
  BacterialvsHealthy = Bacterial - healthy,
  levels = colnames(matriu_disseny)
)

cat("Output 4. Contrast matrix\n\n")
print(matriu_contrast)


## ----include=FALSE-----------------------------------------------------------------------
voom_arnau <- voom(dge_normal2, matriu_disseny)

ml <- lmFit(voom_arnau, matriu_disseny)
ml2 <- contrasts.fit(ml, matriu_contrast)
ml2 <- eBayes(ml2)


## ----include=FALSE-----------------------------------------------------------------------
# Bacterial vs Healthy
bact_vs_health <- topTable(ml2, coef = "BacterialvsHealthy", sort.by = "p", number = nrow(ml2))

# Covid vs Healthy
covid_vs_health <- topTable(ml2, coef = "COVIDvsHealthy", sort.by = "p", number = nrow(ml2))


## ----include=FALSE-----------------------------------------------------------------------
## COVID vs Healthy analysis

# Volcanoplot
volcanoplot(ml2, 
            coef = "COVIDvsHealthy", 
            highlight = 10,
            names = rownames(ml2$coefficients),
            main = "Figure A6. Volcanoplot: COVID-19 vs Healthy")
abline(v = c(-1.5, 1.5), col = "red", lty = 2)  # log2FC thresholds

# Heatmap
gens_top1 <- rownames(subset(covid_vs_health, (abs(logFC) > 1.5) & (adj.P.Val < 0.05))) 

library(pheatmap)
map <- log_comptatges_norm2[gens_top1, ]
map <- map - rowMeans(map)
pheatmap(map, fontsize_col = 14, fontsize_row = 1, main = "Figure A7. Heatmap: COVID-19 vs Healthy")



## ----include=FALSE-----------------------------------------------------------------------
## Bacterial vs Healthy analysis

# Volcanoplot
volcanoplot(ml2, 
            coef = "BacterialvsHealthy", 
            highlight = 10,
            names = rownames(ml2$coefficients),
            main = "Figure A8. Volcanoplot: Bacterial vs Healthy")
abline(v = c(-1.5, 1.5), col = "red", lty = 2)  # log2FC thresholds

# Heatmap
gens_top2 <- rownames(subset(bact_vs_health, (abs(logFC) > 1.5) & (adj.P.Val < 0.05)))

library(pheatmap)
map2 <- log_comptatges_norm2[gens_top2, ]
map2 <- map2 - rowMeans(map2)
pheatmap(map2, fontsize_col = 14, fontsize_row = 1, main =  "Figure A9. Heatmap: Bacterial vs Healthy")


## ----echo = FALSE------------------------------------------------------------------------
summa.fit <- decideTests(ml2, p.value = 0.05, lfc = 1.5)

vc<- vennCounts(summa.fit)

vennDiagram(vc,
            include=c("up", "down"),
            counts.col=c("red", "blue"),
            circle.col = c("red", "blue", "green3"), 
            cex=c(1,1,1),
            main = "Figure 3. Venn Diagram of differentially expressed genes")


## ----include=FALSE-----------------------------------------------------------------------
#Obtain selected gens and genes to be analyzed 
covid_vs_health <- topTable(ml2, coef = "COVIDvsHealthy", sort.by = "p", number = nrow(ml2))

gens_top1 <- subset(covid_vs_health, (logFC > 1.5) & (adj.P.Val < 0.05))

gens_seleccionats<- rownames(gens_top1)
gens_analitzats<- rownames(covid_vs_health)

#Selected genes (upregulated)
gens_seleccionats_entrez <- bitr(gens_seleccionats, 
                            fromType = "ENSEMBL", 
                            toType = "ENTREZID", 
                            OrgDb = org.Hs.eg.db)

#Genes to be analyzed
gens_analitzats_entrez <- bitr(gens_analitzats, 
                            fromType = "ENSEMBL", 
                            toType = "ENTREZID", 
                            OrgDb = org.Hs.eg.db)

#enrichGO analysis
egoda <- enrichGO(gene = gens_seleccionats_entrez$ENTREZID,
                universe = gens_analitzats_entrez$ENTREZID,
                OrgDb = org.Hs.eg.db,
                keyType = "ENTREZID",
                ont = "BP",
                pAdjustMethod = "BH",
                qvalueCutoff = 0.05,
                readable = TRUE)

#Biologic significance analysis plots
dotplot(egoda, showCategory=9, font.size=6)
cnetplot(egoda, font.size=5)
