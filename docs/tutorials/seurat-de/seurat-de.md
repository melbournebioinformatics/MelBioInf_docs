# scRNA-seq Integration and Differential Expression 

## Overview

**Topic**

* [ ] Genomics
* [x] Transcriptomics
* [ ] Proteomics
* [ ] Metabolomics
* [ ] Statistics and visualisation
* [ ] Structural Modelling
* [ ] Basic skills


**Skill level**

* [ ] Beginner  
* [ ] Intermediate  
* [ ] Advanced  

**Data:** IFNB-Stimulated and Control PBMCs.

**Tools:** R >=4.4.0 and associated packages:

* Seurat  
* SeuratData  
* tidyverse  
* DESeq2  
* patchwork  
* pheatmap  
* grid  
* metap  

**Pipeline:**  
*Section 1:* Setup and Intro.  
*Section 2:* Placeholder.  
*Section 3:* Differential Gene Expression when dealing with two treatment conditions.  
*Section 4:* Differential Expression using a pseudobulk approach and DESeq2.  

**Learning objectives:** Placeholder.

<!-- Remove this if not applicable.>
!!! warning "Disclaimer"
    This tutorial is partially based on some existing material.
-->

<!-- Use this divider between sections -->
---
<!-- Use this divider between sections -->

## Setup
**Install something**  

TO-DO
Tell users what to install. Subsection headers using bold, not the hash character.

!!! success "Well done!"
    Confirm that users are ready to start.


<!-- Use this divider between sections -->
---
<!-- Use this divider between sections -->

## Introduction

```r
knitr::opts_chunk$set(echo = TRUE)
```

### Step 1. Load the packages and data

Today we'll be working with Seurat (a popular scRNA-seq analysis package). SeuratData will be used to load in the experimental data we're analysing. Tidyverse is a fundamental and very popularly used set of tools to wrangle and visualise data.

We'll need to load the DESeq2 R package for when we explore pseudobulk DE approaches

pheatmap and grid are two really useful packages for creating custom heatmaps with our scRNA-seq data and exporting figures, respectively.

```r
library(Seurat)
library(SeuratData)
library(tidyverse)
library(DESeq2)
library(patchwork)
library(pheatmap)
library(grid)
library(metap)

set.seed(4242) #Set Seed for Reproducibility
setwd("/Users/vsalazar/Bio/MB/Seurat_DE_Workshop") ## Change this
```

We're using the ifnb public dataset provided by Seurat. This dataset contains PBMC data from 8 lupus patients before and after interferon beta treatment.

I strongly encourage you to explore the other datasets offered by the SeuratData package, it can be really good practice in your spare time.

The ifnb Seurat object we're loading in here was originally made in Seurat v4, there have since been a lot of changes from Seurat v4 to v5 so we'll use the UpdateSeuratObject() function to update the Seurat object so that it is compatible for today.

QUESTION: Looking at the output from the str() function on the ifnb seurat object, can you tell whether this seurat object is processed or unprocessed?

ANSWER: When loading in seurat objects, we can have a look at what processing steps have been performed on it by using the str() function. In the output we can tell that the ifnb Seurat object is unprocessed because:

1)  the scale.data slot is empty

2)  no variable features have been identified

3)  no dimensionality reduction functions have been performed.

```r
AvailableData() # if you want to see the available datasets use this function
InstallData("ifnb") # install our treatment vs control dataset for today

data("ifnb") # Load the dataset into our current R script
ifnb <- UpdateSeuratObject(ifnb) # Make sure the seurat object is in the format of Seurat v5


str(ifnb) # we can use this to take a look at the information in our Seurat Object, we know that this is an unprocessed Seurat object because the scale.data slot is empty, no variable features identified, no reductions performed..etc
```

Note that the `echo = FALSE` parameter was added to the code chunk to prevent printing of the R code that generated the plot.

### Step 2: Run QC, filter out low quality cells

Lets start by processing our data (run the standard seurat workflow steps including preprocessing and filtering).

First we need to take a look at QC metrics, then decide on the thresholds for filtering.

QUESTION: Looking

```r
# Step 2a: QC and filtering
ifnb$percent.mt <- PercentageFeatureSet(object = ifnb, pattern = "^MT-") # First let's annotate the mitochondrial percentage for each cell

View(ifnb@meta.data) # we can take a look mitochondrial percentages for the seurat object by viewing the seurat objects metadata

# Step 2b: Visualise QC metrics and identify filtering thresholds
qc.metric.plts <- VlnPlot(ifnb, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3) + 
  ggtitle("Before Filtering")

association.plt.raw <- FeatureScatter(ifnb, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") + geom_smooth(method = "lm") +
  ggtitle("Before Filtering")

qc.metric.plts
association.plt.raw
```

After visualising QC metrics, we'll move on to the actual filtering

```r
# Step 2c: filter out low-quality cells + visualise the metrics for our filtered seurat object
ifnb.filtered <- subset(ifnb, subset = nCount_RNA > 800 & 
                          nCount_RNA < 5000 &
                          nFeature_RNA > 200 &
                          nFeature_RNA < 1200 &
                          percent.mt < 5)

qc.metric.plts.filtered <- VlnPlot(ifnb.filtered, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3) +
  ggtitle("After Filtering")
association.plt.filtered <- FeatureScatter(ifnb.filtered, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") + geom_smooth(method = "lm") +
  ggtitle("After Filtering")

qc.metric.plts.filtered
association.plt.filtered
```

Let's check how many cells we've filtered out (looks like \~400 cells were removed):

```r
## Defining a couple helper functions to standardise x and y axis for two plots
get_plot_range <- function(plot) {
  data <- layer_data(plot)
  list(
    x = range(data$x, na.rm = TRUE),
    y = range(data$y, na.rm = TRUE)
  )
}
standardise_plt_scale <- function(plt1, plt2){

  # Get ranges for both plots
  range_raw <- get_plot_range(plt1)
  range_filtered <- get_plot_range(plt2)
  
  # Calculate overall range
  x_range <- range(c(range_raw$x, range_filtered$x))
  y_range <- range(c(range_raw$y, range_filtered$y))
  
  suppressMessages({
  # Update both plots with the same x and y scales
  association.plt.raw <- association.plt.raw +
    scale_x_continuous(limits = x_range) +
    scale_y_continuous(limits = y_range)
  
  association.plt.filtered <- association.plt.filtered +
    scale_x_continuous(limits = x_range) +
    scale_y_continuous(limits = y_range)
  })
  
  # Wrap the plots
  wrapped_plots <- wrap_plots(list(association.plt.raw, association.plt.filtered), 
                              ncol = 2)

  return(wrapped_plots)
}

wrap_plots(list(qc.metric.plts, qc.metric.plts.filtered), 
           ncol = 1)
association.plts <- standardise_plt_scale(association.plt.raw,
                                          association.plt.filtered)
association.plts


ifnb
ifnb.filtered

ifnb.filtered[["RNA"]] <- split(ifnb.filtered[["RNA"]], f = ifnb.filtered$stim) # Lets split our count matrices based on conditions (stored within different layers) -> needed for integration steps in Seurat v5
```

### Step 3: Before performing differential expression between the two conditions, let's assess whether we need to integrate our data

After filtering out low quality cells, we want to visualise our data to see how cells group by condition and if we need to perform batch-effect correction (integration)

```r
ifnb.filtered <- NormalizeData(ifnb.filtered)
ifnb.filtered <- FindVariableFeatures(ifnb.filtered)
ifnb.filtered <- ScaleData(ifnb.filtered)
ifnb.filtered <- RunPCA(ifnb.filtered)

ElbowPlot(ifnb.filtered) # Visualise the dimensionality of the data, looks like 15 PCs is adequate to capture the majority of the variation in the data, but we'll air on the higher side and consider all 20 dimensions.
ifnb.filtered <- RunUMAP(ifnb.filtered, dims = 1:20, reduction = 'pca')

DimPlot(ifnb.filtered, reduction = 'umap', group.by = 'stim') # lets see how our cells separate by condition and whether integration is necessary

DimPlot(ifnb.filtered, reduction = 'pca', group.by = 'stim') # lets see how our cells separate by condition and whether integration is necessary

# These are PBMCs before and after treatment, there should be cells that are similar between both conditions, it looks like we'll have to run some batch effect correction to overlay similar cell-types from both conditions to perform downstream analysis
```

QUESTION: Do you think we need to integrate our data? Hint: Look at the UMAP and PC1/PC2 plots we made above

QUESTION: What do you think would happen if we were to perform unsupervised clustering right now, without integrating our data (or overlaying similar cells ontop of each other from both conditions)?

### Step 4: Integrating our data using the harmony method

Seurat v5 has made it really easy to test different integration methods quickly, let's use a really popular approach (harmony) first.

```r
# code adapted from: https://satijalab.org/seurat/articles/seurat5_integration
ifnb.filtered <- IntegrateLayers(object = ifnb.filtered,
                                 method = HarmonyIntegration,
                                 orig.reduction = "pca", 
                                 new.reduction = "harmony")

ifnb.filtered <- RunUMAP(ifnb.filtered, reduction = "harmony", dims = 1:20, reduction.name = "umap.harmony")

after.harmony <- DimPlot(ifnb.filtered, reduction = "umap.harmony", group.by = "stim") + 
  ggtitle("After Harmony Integration")
before.integration <- DimPlot(ifnb.filtered, reduction = "umap", group.by = "stim") +
  ggtitle("Before Integration")

before.integration | after.harmony
```

Question: Looking at the UMAPs above, do you think integration was successful? Have a slide on what if its just different cell types

#### Trying a different integration method (Seurat CCA), lets see if our integration improves

Question: Try looking at the PC1 and PC2 plots for harmony and seurat as well

```r
ifnb.filtered <- IntegrateLayers(object = ifnb.filtered,
                                 method = CCAIntegration,
                                 orig.reduction = "pca", 
                                 new.reduction = "integrated.cca")

ifnb.filtered <- RunUMAP(ifnb.filtered, reduction = "integrated.cca", dims = 1:20, reduction.name = "umap.cca")

after.seuratCCA <- DimPlot(ifnb.filtered, reduction = "umap.cca", group.by = "stim") +
  ggtitle("After Seurat CCA Integration")

before.integration | after.seuratCCA
after.harmony | after.seuratCCA

## Show example slide of integration 'failing' but due to different cell types in each sample ***
```

Question: What do you think of the integration results now?

Challenge: look at the pc1/2 plots for each integration method -\> explain different pca slots here

\*\*Slide on the caveats of integration, you can lose biological signal when you regress out technical batch effects\*\*

### Step 5: Once we're happy with integration we need to rejoin layers before moving on to DE

This step collapses individual control and treatment datasets together and needs to be done before differential expression analysis

```r
ifnb.filtered <- FindNeighbors(ifnb.filtered, reduction = "integrated.cca", dims = 1:20)
ifnb.filtered <- FindClusters(ifnb.filtered, resolution = 0.5)

ifnb.filtered <- JoinLayers(ifnb.filtered)
```

## Part 2: Differential Gene Expression when dealing with two treatment conditions

```r
DimPlot(ifnb.filtered, reduction = "umap.cca", label = T)
DimPlot(ifnb.filtered, reduction = "umap.cca", group.by = "stim")
```

### Step 1. Find Conserved Markers to label our celltypes

```r
## Let's look at conserved markers in cluster 4 across our two conditions (compared to all other clusters)
markers.cluster.4 <- FindConservedMarkers(ifnb.filtered, ident.1 = 4,
                     grouping.var = 'stim')

head(markers.cluster.4)

## **Explain p_val adjusted and p_val and avg_log2FC** in slide
```

Let's visualise the top features identified in this cluster

(Will need to explain the min.cutoff value here)

Question: try running this function without defining a min.cutoff, or changing the value

```r
# Try looking up some of these markers here: https://www.proteinatlas.org/
FeaturePlot(ifnb.filtered, reduction = "umap.cca", 
            features = c('VMO1', 'FCGR3A', 'MS4A7', 'CXCL16'), min.cutoff = 'q10')

# I happen to know that the cells in cluster 3 are CD16 monocytes - lets rename this cluster
Idents(ifnb.filtered) # Let's look at the identities of our cells at the moment
ifnb.filtered <- RenameIdents(ifnb.filtered, '4' = 'CD16 Mono') # Let's rename cells in cluster 3 with a new cell type label
Idents(ifnb.filtered) # we can take a look at the cell identities again

DimPlot(ifnb.filtered, reduction = "umap.cca", label = T) +
  ggtitle("After changing the identity of cluster 4")
```

### Step 2: Set the identity of our clusters to the annotations provided

```r
Idents(ifnb.filtered) <- ifnb.filtered@meta.data$seurat_annotations
Idents(ifnb.filtered)

DimPlot(ifnb.filtered, reduction = "umap.cca", label = T)

## **refer to automatic cell type annotation in slide**
```

### Step 3: Find DEGs between our two conditions (using CD16 Mono cells as an example)

```r
# Make another column in metadata showing what cells belong to each treatment group (This will make more sense in a bit)
ifnb.filtered$celltype.and.stim <- paste0(ifnb.filtered$seurat_annotations, '_', ifnb.filtered$stim)
View(ifnb.filtered@meta.data)

Idents(ifnb.filtered) <- ifnb.filtered$celltype.and.stim

DimPlot(ifnb.filtered, reduction = "umap.cca", label = T) # each cluster is now made up of two labels (control or stimulated)
DimPlot(ifnb.filtered, reduction = "umap.cca", 
        label = T, split.by = "stim") # Lets separate by condition to see what we've done a bit more clearly

```

We'll now leverage these new identities to compare DEGs between our treatment groups

```r
treatment.response.CD16 <- FindMarkers(ifnb.filtered, ident.1 = 'CD16 Mono_STIM', 
                                       ident.2 = 'CD16 Mono_CTRL')
head(treatment.response.CD16) # These are the genes that are upregulated in the stimulated versus control group
```

### Step 5: Lets plot conserved features vs DEGs between conditions

```r
FeaturePlot(ifnb.filtered, reduction = 'umap.cca', 
            features = c('VMO1', 'FCGR3A', 'IFIT1', 'ISG15'),
            split.by = 'stim', min.cutoff = 'q10')
```

### Step 6: Create a Heatmap to visualise DEGs between our two conditions + cell types

```r
# # Find upregulated genes in each group (cell type and condition)
# ifnb.treatVsCtrl.markers <- FindAllMarkers(ifnb.filtered,
#                                           only.pos = TRUE)
# saveRDS(ifnb.treatVsCtrl.markers, "ifnb_stimVsCtrl_markers.rds")
```
Uncomment the top code block to run from scratch -\> We'll load in the saved output since this can take a while

Seurat's in built heatmap function can be quite messy and hard to interpret sometimes (we'll learn how to make better and clearer custom heatmaps from our Seurat scRNA-seq data in a bit).

```r
ifnb.treatVsCtrl.markers <- readRDS("ifnb_stimVsCtrl_markers.rds")

top5 <- ifnb.treatVsCtrl.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 5) %>%
  ungroup()

DEG.heatmap <- DoHeatmap(ifnb.filtered, features = top5$gene,
          label = FALSE)

DEG.heatmap
```

## Part 3 : Differential Expression using a pseudobulk approach and DESeq2

### Step 1: We need to import sample information for each cell from the original paper

Question: have a look at the ifnb.filtered seurat metadata, can you spot what we've done here?

```r
# defining a function here to retrieve that information (code from https://satijalab.org/seurat/articles/de_vignette)
loadDonorMetadata <- function(seu.obj){
  # load the inferred sample IDs of each cell
  ctrl <- read.table(url("https://raw.githubusercontent.com/yelabucsf/demuxlet_paper_code/master/fig3/ye1.ctrl.8.10.sm.best"), head = T, stringsAsFactors = F)
  stim <- read.table(url("https://raw.githubusercontent.com/yelabucsf/demuxlet_paper_code/master/fig3/ye2.stim.8.10.sm.best"), head = T, stringsAsFactors = F)
  info <- rbind(ctrl, stim)
  
  # rename the cell IDs by substituting the '-' into '.'
  info$BARCODE <- gsub(pattern = "\\-", replacement = "\\.", info$BARCODE)
  
  # only keep the cells with high-confidence sample ID
  info <- info[grep(pattern = "SNG", x = info$BEST), ]
  
  # remove cells with duplicated IDs in both ctrl and stim groups
  info <- info[!duplicated(info$BARCODE) & !duplicated(info$BARCODE, fromLast = T), ]
  
  # now add the sample IDs to ifnb 
  rownames(info) <- info$BARCODE
  info <- info[, c("BEST"), drop = F]
  names(info) <- c("donor_id")
  seu.obj <- AddMetaData(seu.obj, metadata = info)
  
  # remove cells without donor IDs
  seu.obj$donor_id[is.na(seu.obj$donor_id)] <- "unknown"
  seu.obj <- subset(seu.obj, subset = donor_id != "unknown")
  
  return(seu.obj)
}

ifnb.filtered <- loadDonorMetadata(ifnb.filtered)
```

### Step 2: Aggregate our counts based on treatment group, cell-type, and donor id

Condensing our single-cell matrix in this manner is referred to 'pseudobulking'. We're making a condensed count matrix that looks more like a bulk matrix so that we can use bulk differential expression algorithms like DESeq2. We can see this clearly when we have a look at the ifnb.pseudbulk.df we make in the following code block.

```r
ifnb.pseudobulk <- AggregateExpression(ifnb.filtered, assays = "RNA",
                                   group.by = c("stim", "donor_id", "seurat_annotations"),
                                   return.seurat = TRUE)
View(ifnb.pseudobulk@meta.data)

# If you want the pseudobulk matrix as a dataframe you can do this:
ifnb.pseudobulk.df <- AggregateExpression(ifnb.filtered, assays = "RNA",
                                          group.by = c("stim", "donor_id", "seurat_annotations")) %>% 
  as.data.frame()

head(ifnb.pseudobulk.df)
```

### Step 3: Perform Differential Expression using DESeq2

Just like before, lets make a new column containing the cell type and treatment group before DE with DESeq2. We do this so that we can leverage the Idents() function combined with the FindMarkers function.

To use DESeq2 in older versions of Seurat (prior to v5) we would have to perform a myriad of convoluted data wrangling steps to convert our Seurat pseudobulk matrix into a format compatible with DESeq2 (or whatever bulk DE package being used). Fortunately in Seurat v5, these different bulk DE algorithms/packages have been wrapped into the FindMarker() functions - simplifying pseudobulk DE approaches immensely.

It is worth taking a look at the Seurat DE vignette on your own to see the other pseudobulk DE methods you can use. Do keep in mind that to use some of the DE tests, you need to install that package and load it as a library separately to Seurat, as we've done here with DESeq2 at the start.

```r
ifnb.pseudobulk$celltype.and.stim <- paste(ifnb.pseudobulk$seurat_annotations, ifnb.pseudobulk$stim, sep = "_")
Idents(ifnb.pseudobulk) <- "celltype.and.stim"

View(ifnb.pseudobulk@meta.data)

# Lets run a DEG test between treated and control CD 16 monocytes using the same FindMarkers function but with DESeq2
treatment.response.CD16.pseudo <- FindMarkers(object = ifnb.pseudobulk, 
                                      ident.1 = 'CD16 Mono_STIM', 
                                      ident.2 = 'CD16 Mono_CTRL',
                                      test.use = "DESeq2")

head(treatment.response.CD16.pseudo)

# How are we able to use the same findmarkers function here?
head(Cells(ifnb.pseudobulk)) # our 'cells' are no longer barcodes, but have been renamed according to stim-donor-annotation when we aggregated our data earlier
```

### Step 4: Assessing differences between our pseudbulk DEGs and single-cell DEGs

First lets take a look at the DEG dataframes we made for both.

QUESTION: Having a look at them, can you identify any differences? Can you think of any underlying reasons behind those differences?

Hint: Look at the p_val and p_val_adj columns

```r
head(treatment.response.CD16)
head(treatment.response.CD16.pseudo)
```

Next let's take a look at the degree of overlap between the actual DEGs in both approaches.

To help with this, I'm just going to define a couple helper functions below.

```r
Merge_DEG_dataframes <- function(pseudobulk.de,
                                 singlecell.de){
  names(pseudobulk.de) <- paste0(names(pseudobulk.de), ".bulk")
  pseudobulk.de$gene <- rownames(pseudobulk.de)
  
  names(singlecell.de) <- paste0(names(singlecell.de), ".sc")
  singlecell.de$gene <- rownames(singlecell.de)
  
  merge_dat <- merge(singlecell.de, pseudobulk.de, by = "gene")
  merge_dat <- merge_dat[order(merge_dat$p_val.bulk), ]
  
  return(merge_dat)
}
Visualise_Overlapping_DEGs <- function(pseudobulk.de,
                                       singlecell.de){
  
  merge_dat <- Merge_DEG_dataframes(pseudobulk.de,
                                    singlecell.de)
  
  # Number of genes that are marginally significant in both; marginally significant only in bulk; and marginally     significant only in single-cell
  common <- merge_dat$gene[which(merge_dat$p_val_adj.bulk < 0.05 & 
                                   merge_dat$p_val_adj.sc < 0.05)]
  only_sc <- merge_dat$gene[which(merge_dat$p_val_adj.bulk > 0.05 & 
                                    merge_dat$p_val_adj.sc < 0.05)]
  only_pseudobulk <- merge_dat$gene[which(merge_dat$p_val_adj.bulk < 0.05 & 
                                      merge_dat$p_val_adj.sc > 0.05)]
  
  # Create a dataframe to plot overlaps
  overlap.info <- data.frame(
    category = c("Common", 
                 "Only in single-cell", 
                 "Only in pseudobulk"),
    count = c(length(common), length(only_sc), length(only_pseudobulk))
  )
  
  overlap.info$category <- factor(overlap.info$category, 
                                  levels = c("Common", 
                                             "Only in single-cell", 
                                             "Only in pseudobulk"))
  # Create the bar plot
  overlap.bar.plt <- ggplot(overlap.info, aes(x = category, y = count, fill = category)) +
    geom_bar(stat = "identity") +
    geom_text(aes(label = count), vjust = -0.5, size = 4) +
    theme_minimal() +
    labs(title = "Number of Overlapping and Unique DEGs from Single-Cell and Pseudobulk Tests",
         x = "",
         y = "Number of Genes") +
    theme(legend.position = "none",
          axis.text.x = element_text(angle = 45, hjust = 1))
  
  return(overlap.bar.plt)
}
```

Let's use the helper functions above to plot and view the overlap/agreement between our pseudobulk versus single-cell approaches

QUESTION: Can you explain why we're seeing the discrepencies we see between the two methods here?

```r
merged_deg_data <- Merge_DEG_dataframes(pseudobulk.de = treatment.response.CD16.pseudo,
                                        singlecell.de = treatment.response.CD16)
merged_deg_data %>% 
  dplyr::select(gene, 
                p_val.sc, p_val.bulk,
                p_val_adj.sc, p_val_adj.bulk) %>%
  head(10)

# How many DEGs overlap between our two methods? Is there anything in the merged_deg_data frame that stands out to you?
overlap.bar.plt <- Visualise_Overlapping_DEGs(pseudobulk.de = treatment.response.CD16.pseudo,
                                              singlecell.de = treatment.response.CD16)

overlap.bar.plt
```

### Step 5: Investigate the differences between pseudobulk DE and single-cell DE closer

Let's create lists of genes in each of our categories first

```r
common <- merged_deg_data$gene[which(merged_deg_data$p_val.bulk < 0.05 & 
                                       merged_deg_data$p_val.sc < 0.05)]
only_sc <- merged_deg_data$gene[which(merged_deg_data$p_val.bulk > 0.05 & 
                                        merged_deg_data$p_val.sc < 0.05)]
only_pseudobulk <- merged_deg_data$gene[which(merged_deg_data$p_val.bulk < 0.05 & 
                                          merged_deg_data$p_val.sc > 0.05)]
```

Now I want to look at the expression of genes that only appear in our sc deg test and not in our pseudobulk deg test I've picked two genes from the 'only_sc' variable we just defined

```r
# create a new column to annotate sample-condition-celltype in the single-cell dataset
ifnb.filtered$donor_id.and.stim <- paste0(ifnb.filtered$stim, "-", ifnb.filtered$donor_id)
Idents(ifnb.filtered) <- "celltype.and.stim"

# Explore some genes that only appear in the sc deg test---
VlnPlot(ifnb.filtered, features = c("PABPC1", "SRGN"), 
        idents = c("CD16 Mono_CTRL", "CD16 Mono_STIM"), 
        group.by = "stim") 

VlnPlot(ifnb.filtered, features = c("PABPC1", "SRGN"), 
        idents = c("CD16 Mono_CTRL", "CD16 Mono_STIM"), 
        group.by = "donor_id.and.stim", ncol = 1)
```

QUESTION: How would you interpret the plots above? What does this tell you about some of the pitfalls of single-cell DE approaches?

Now lets take a look at some genes that show agreement across both sc and pseudobulk deg tests

```r
VlnPlot(ifnb.filtered, features = c("IFIT2", "PSMA4"), 
        idents = c("CD16 Mono_CTRL", "CD16 Mono_STIM"), 
        group.by = "stim") 

VlnPlot(ifnb.filtered, features = c("IFIT2", "PSMA4"), 
        idents = c("CD16 Mono_CTRL", "CD16 Mono_STIM"), 
        group.by = "donor_id.and.stim", ncol = 1) 
```

### Step 6: Creating our own custom visualisations for DEG analysis between cell-types in two different experimental groups

So far, we've been using functions wrapped within Seurat to plot and visualise our data. But we're not limited to those functions. Let's use the pheatmap ('pretty heatmap') package to visualise our DEGs with adjusted p values \< 0.05.

First lets look at our significant DEGs as defined by our pseudobulk approach

```r
CD16.sig.markers <- treatment.response.CD16.pseudo %>% 
  dplyr::filter(p_val_adj < 0.05) %>%
  dplyr::mutate(gene = rownames(.))
```

This is how we can pull our average (scaled) pseudobulk expression values from our seurat obj:

```r
ifnb.filtered$celltype.stim.donor_id <- paste0(ifnb.filtered$seurat_annotations, "-",
                                               ifnb.filtered$stim, "-", ifnb.filtered$donor_id)
Idents(ifnb.filtered) <- "celltype.stim.donor_id"

all.sig.avg.Expression.mat <- AverageExpression(ifnb.filtered, 
                         features = CD16.sig.markers$gene, 
                         layer = 'scale.data')

View(all.sig.avg.Expression.mat %>%
  as.data.frame())
```

Now lets make sure we're only using data from the CD16 cell type

```r
CD16.sig.avg.Expression.mat <- all.sig.avg.Expression.mat$RNA %>%
  as.data.frame() %>%
  dplyr::select(starts_with("CD16 Mono"))

View(CD16.sig.avg.Expression.mat)

```

Let's finally view our heatmap with averaged pseudobulk scaled expression values for our signficant DEGs

```r
pheatmap::pheatmap(CD16.sig.avg.Expression.mat,
         cluster_rows = TRUE,
         show_rownames = FALSE, 
         border_color = NA, 
         fontsize = 10, 
         scale = "row", 
         fontsize_row = 10, 
         height = 20)
```

The cool thing about pheatmap is that we can customise our heatmap with added metadata

```r
cluster_metadata <- data.frame(
  row.names = colnames(CD16.sig.avg.Expression.mat)
) %>% 
  dplyr::mutate(
    Cell_Type = "CD16 Mono",
    Treatment_Group = ifelse(str_detect(row.names(.), "STIM|CTRL"), 
                      str_extract(row.names(.), "STIM|CTRL")))

sig.DEG.heatmap <- pheatmap::pheatmap(CD16.sig.avg.Expression.mat,
         cluster_rows = TRUE,
         show_rownames = FALSE,
         annotation = cluster_metadata[, c("Treatment_Group", "Cell_Type")], 
         border_color = NA, 
         fontsize = 10, 
         scale = "row", 
         fontsize_row = 10, 
         height = 20,
         annotation_names_col = FALSE)

sig.DEG.heatmap
```

Lets save our heatmap using the grid library (very useful package!) Code from: <https://stackoverflow.com/questions/43051525/how-to-draw-pheatmap-plot-to-screen-and-also-save-to-file>

```r
save_pheatmap_pdf <- function(x, filename, width=7, height=7) {
   stopifnot(!missing(x))
   stopifnot(!missing(filename))
   pdf(filename, width=width, height=height)
   grid::grid.newpage()
   grid::grid.draw(x$gtable)
   dev.off()
}

save_pheatmap_pdf(sig.DEG.heatmap, "sig_DEG_pseudo.pdf")
```

CHALLENGE QUESTION: do the same thing with significant DEGs from the sc approach. Do you see any differences?

Solution: run the same code as above but with the 'treatment.response.CD16' instead of 'treatment.response.CD16.pseudo' variable

CONCLUSION

QUESTION: What kind of analyses can we do next after having a list of DEGs for each cell type and treatment groups?


---

Add keypoints at the end of each section or subsection. Keypoints use the `note` admonition.

Read more about admonitions on the [Material for MkDocs documentation](https://squidfunk.github.io/mkdocs-material/reference/admonitions/).

!!! note "Keypoints"

    - Version control is like an unlimited ‘undo’.
    - Version control also allows many people to work in parallel.

!!! warning "TO-DO"
    - Add more admonitions to the template.

!!! question "Question"

    Questions have a nested, collapsed solution that uses the "example" admonition.

    This is because it's not uncommon to have "Keypoints" right after a question, and keypoints use the "note" admonition. So we prevent repetition.

    ??? example "Solution"
        Answer to the question above.


## Glossary

A glossary at the end of the page is always helpful. This way you can refer to [terms](./template.md#glossary) throughout the lesson.

**Terms:** are written in bold in the glossary.


Calculating cluster CD14 Mono_CTRL
  |++++++++++++++++++++++++++++++++++++++++++++++++++| 100% elapsed=22s  
Calculating cluster pDC_CTRL
  |++++++++++++++++++++++++++++++++++++++++++++++++++| 100% elapsed=01m 02s
Calculating cluster CD4 Memory T_CTRL
  |++++++++++++++++++++++++++++++++++++++++++++++++++| 100% elapsed=47s  
Calculating cluster T activated_CTRL
  |++++++++++++++++++++++++++++++++++++++++++++++++++| 100% elapsed=60s  
Calculating cluster CD4 Naive T_CTRL
  |++++++++++++++++++++++++++++++++++++++++++++++++++| 100% elapsed=49s  
Calculating cluster CD8 T_CTRL
  |++++++++++++++++++++++++++++++++++++++++++++++++++| 100% elapsed=58s  
Calculating cluster Mk_CTRL
  |++++++++++++++++++++++++++++++++++++++++++++++++++| 100% elapsed=38s  
Calculating cluster B Activated_CTRL
  |++++++++++++++++++++++++++++++++++++++++++++++++++| 100% elapsed=52s  
Calculating cluster B_CTRL
  |++++++++++++++++++++++++++++++++++++++++++++++++++| 100% elapsed=52s  
Calculating cluster CD16 Mono_CTRL
  |++++++++++++++++++++++++++++++++++++++++++++++++++| 100% elapsed=24s  
Calculating cluster NK_CTRL
  |++++++++++++++++++++++++++++++++++++++++++++++++++| 100% elapsed=01m 02s
Calculating cluster DC_CTRL
  |++++++++++++++++++++++++++++++++++++++++++++++++++| 100% elapsed=33s  
Calculating cluster Eryth_CTRL
  |++++++++++++++++++++++++++++++++++++++++++++++++++| 100% elapsed=13s  
Calculating cluster CD8 T_STIM
  |++++++++++++++++++++++++++++++++++++++++++++++++++| 100% elapsed=48s  
Calculating cluster pDC_STIM
  |++++++++++++++++++++++++++++++++++++++++++++++++++| 100% elapsed=56s  
Calculating cluster CD4 Naive T_STIM
  |++++++++++++++++++++++++++++++++++++++++++++++++++| 100% elapsed=41s  
Calculating cluster B_STIM
  |++++++++++++++++++++++++++++++++++++++++++++++++++| 100% elapsed=45s  
Calculating cluster CD14 Mono_STIM
  |++++++++++++++++++++++++++++++++++++++++++++++++++| 100% elapsed=17s  
Calculating cluster T activated_STIM
  |++++++++++++++++++++++++++++++++++++++++++++++++++| 100% elapsed=53s  
Calculating cluster CD4 Memory T_STIM
  |++++++++++++++++++++++++++++++++++++++++++++++++++| 100% elapsed=35s  
Calculating cluster B Activated_STIM
  |++++++++++++++++++++++++++++++++++++++++++++++++++| 100% elapsed=51s  
Calculating cluster NK_STIM
  |++++++++++++++++++++++++++++++++++++++++++++++++++| 100% elapsed=52s  
Calculating cluster CD16 Mono_STIM
  |++++++++++++++++++++++++++++++++++++++++++++++++++| 100% elapsed=21s  
Calculating cluster DC_STIM
  |++++++++++++++++++++++++++++++++++++++++++++++++++| 100% elapsed=23s  
Calculating cluster Mk_STIM
  |++++++++++++++++++++++++++++++++++++++++++++++++++| 100% elapsed=38s  
Calculating cluster Eryth_STIM
  |++++++++++++++++++++++++++++++++++++++++++++++++++| 100% elapsed=13s
