Aim 2 - Finding differentially expressed features (cluster biomarkers)
================
Yilin Qiu
07/04/2021



``` r
#load library
library(dplyr)
```

    ## 
    ## Attaching package: 'dplyr'

    ## The following objects are masked from 'package:stats':
    ## 
    ##     filter, lag

    ## The following objects are masked from 'package:base':
    ## 
    ##     intersect, setdiff, setequal, union

``` r
library(Seurat)
```

    ## Attaching SeuratObject

``` r
library(patchwork)
library(data.table)
```

    ## 
    ## Attaching package: 'data.table'

    ## The following objects are masked from 'package:dplyr':
    ## 
    ##     between, first, last

``` r
dt <- readRDS("Aim_3/All_Samples_Merged_Filtered.txt.gz")


# find all markers of cluster 1
cluster1.markers <- FindMarkers(dt, ident.1 = 2, min.pct = 0.25)
head(cluster1.markers, n = 5)
```

    ##       p_val avg_log2FC pct.1 pct.2 p_val_adj
    ## KIF1B     0  0.8289893 0.385 0.149         0
    ## PGD       0  1.2912575 0.596 0.210         0
    ## PADI2     0  0.5232170 0.267 0.121         0
    ## PADI4     0  1.2878607 0.321 0.063         0
    ## CDA       0  1.4002589 0.637 0.222         0

``` r
# find all markers distinguishing cluster 5 from clusters 0 and 3
cluster5.markers <- FindMarkers(dt, ident.1 = 5, ident.2 = c(0, 3), min.pct = 0.25)
head(cluster5.markers, n = 5)
```

    ##       p_val avg_log2FC pct.1 pct.2 p_val_adj
    ## ALPL      0 -1.2051144 0.172 0.461         0
    ## CSF3R     0 -0.5045555 0.857 0.961         0
    ## SMAP2     0 -1.2165508 0.269 0.570         0
    ## JUN       0  1.7010061 0.338 0.114         0
    ## CTSS      0 -0.9231722 0.447 0.691         0

``` r
# find markers for every cluster compared to all remaining cells, report only the positive ones
dt.markers <- FindAllMarkers(dt, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
```

    ## Calculating cluster 0

    ## Calculating cluster 1

    ## Calculating cluster 2

    ## Calculating cluster 3

    ## Calculating cluster 4

    ## Calculating cluster 5

    ## Calculating cluster 6

    ## Calculating cluster 7

    ## Calculating cluster 8

    ## Calculating cluster 9

    ## Calculating cluster 10

    ## Calculating cluster 11

    ## Calculating cluster 12

    ## Calculating cluster 13

``` r
dt.markers_2 <- dt.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC)

# VlnPlot() (shows expression probability distributions across clusters)
VlnPlot(dt, features = c("HBB", "SMAP2"))
```

![](Aim_2_files/figure-gfm/unnamed-chunk-1-1.png)<!-- -->

``` r
f1 <- FeaturePlot(dt, features = c("HBB", "SMAP2", "IFIT3", "IFIT2"))
f1
```

![](Aim_2_files/figure-gfm/unnamed-chunk-1-2.png)<!-- -->

``` r
f2 <- FeaturePlot(dt, features = c("S100A12", "MMP9", "IL1B", "IL1RN"))
f2
```

![](Aim_2_files/figure-gfm/unnamed-chunk-2-1.png)<!-- -->

``` r
#DoHeatmap() generates an expression heatmap for given cells and features. In this case, we are plotting the top 20 markers (or all markers if less than 20) for each cluster.
top10 <- dt.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
DoHeatmap(dt, features = top10$gene) + NoLegend()
```

    ## Warning in DoHeatmap(dt, features = top10$gene): The following features were
    ## omitted as they were not found in the scale.data slot for the RNA assay:
    ## ATP6V1B2, PHACTR1, PYGL, S100A9, TSPO, C10orf54, PTPRE, CXCR2, TXNIP, CELF2,
    ## MIDN, CPPED1, SORL1, SMAP2

![](Aim_2_files/figure-gfm/unnamed-chunk-3-1.png)<!-- -->

``` r
new.cluster.ids <- c("Neutrophils (1)", "Neutrophils (2)", "Monocytes (1)", "Macrophages (1)", "Neutrophils (3)", "Neutrophils (4)", "T cells, NK-cells", "Macrophages (2)", "Monocytes (2)", "Neutrophils (5)", "Neutrophils (6)", "Neutrophils (7)","Epithelial_cell", "B cells")

names(new.cluster.ids) <- levels(dt)
dt_2 <- RenameIdents(dt, new.cluster.ids)
DimPlot(dt_2, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()
```

![](Aim_2_files/figure-gfm/unnamed-chunk-4-1.png)<!-- -->

\`\`\`
