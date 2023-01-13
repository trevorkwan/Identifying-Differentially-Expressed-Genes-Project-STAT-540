GO Analysis
================
Madeline Iseminger (38390134)
09/04/2021

This script finds the genes that are differentially expressed between
male and female individuals and runs GO analysis on those genes.
Neutrophils across disease severities and sample types are included
here, but any combination could be used with get\_DE\_genes() from
Trevor’s script (copied here).

## Setup and load expression data

Run using Ubuntu 20.04, R 4.0.5, RStudio Version 1.4.1103.

Load the binary datafile of filtered and integrated data.

``` r
dt <- readRDS("~/Desktop/AllSamplesMergedFiltered/AllSamplesMergedFiltered")
```

## Get DE genes

This section is from Trevor’s script.

# list cluster ID to cell type conversions

``` r
cluster_ID_list <- list(Neutrophils = c("0", "1", "4", "5", "9", "10", "11"),
                        Monocytes = c("2", "8"), 
                        Macrophages = c("3", "7"), 
                        T_cells_NK_cells = "6", 
                        Epithelial_cells = "12", 
                        B_cells = "13")
```

# create function to identify top DE genes

``` r
get_DE_genes <- function(COVID_severity, tissue_type, cell_type){
  
  if (is.na(COVID_severity) == FALSE) {
  # subset by clinic status
  dt.sub <- subset(dt, clinic_status == COVID_severity)
  }
  
  if (is.na(tissue_type) == FALSE) {
  # subset by tissue type
  dt.sub <- subset(dt, tissue == tissue_type)
  }
  # subset by cell type
  if (cell_type == "Neutrophils") {
    cell_ID <- 1
  }
  if (cell_type == "Monocytes") {
    cell_ID <- 2
  }
  if (cell_type == "Macrophages") {
    cell_ID <- 3
  }
  if (cell_type == "T_cells_NK_cells") {
    cell_ID <- 4
  }
  if (cell_type == "Epithelial_cells") {
    cell_ID <- 5
  }
  if (cell_type == "B_cells") {
    cell_ID <- 6
  }
  dt.sub <- subset(dt.sub, seurat_clusters == cluster_ID_list[[cell_ID]])
  # create df
  Idents(dt.sub) <- c("male", "female")
  avg.dt.sub <- as.data.frame(log1p(AverageExpression(dt.sub, verbose = FALSE)$RNA))
  avg.dt.sub$gene <- rownames(avg.dt.sub)
  # find DE genes
  dt.response <- FindMarkers(dt.sub, ident.1 = "male", ident.2 = "female", verbose = FALSE)
  # dt.response <- head(dt.response, n = 15)
  return(dt.response)
}
```

## Identify DE genes for GO enrichment: Neutrophils

# get all significant DE genes for mild COVID and neutrophils

``` r
neut_mild_df <- get_DE_genes(COVID_severity = "Mild COVID", tissue_type = NA, cell_type = "Neutrophils")
```

    ## Warning in `==.default`(seurat_clusters, cluster_ID_list[[cell_ID]]): longer
    ## object length is not a multiple of shorter object length

    ## Warning in is.na(e1) | is.na(e2): longer object length is not a multiple of
    ## shorter object length

``` r
neut_mild_genes <- neut_mild_df %>% 
  rownames() %>% 
  as.vector()
```

# get all significant DE genes for healthy controls and neutrophils

``` r
neut_healthy_df <- get_DE_genes(COVID_severity = "Healthy control", tissue_type = NA, cell_type = "Neutrophils")
```

    ## Warning in `==.default`(seurat_clusters, cluster_ID_list[[cell_ID]]): longer
    ## object length is not a multiple of shorter object length

    ## Warning in is.na(e1) | is.na(e2): longer object length is not a multiple of
    ## shorter object length

``` r
neut_healthy_genes <- neut_healthy_df %>% 
  rownames() %>% 
  as.vector()
```

## GO enrichment analysis from a list of genes

Fetch GO terms for selected genes and plot them with their functions.

``` r
# select database
setEnrichrSite("Enrichr") # human genes
```

    ## Connection changed to https://maayanlab.cloud/Enrichr/

    ## Connection is Live!

``` r
websiteLive <- TRUE
# list available databases
dbs <- listEnrichrDbs()
```

Get the list of genes to query:

``` r
# list of genes to query
#geneList <- neut_healthy_genes # healthy, neutrophils
#geneList <- neut_mild_genes # mild, neutrophils
#geneList <- c("MT2A", "SERPINB9", "HBB", "HBA2") # severe, neutrophils
geneList <- c("IFNGR2", "INSIG1", "SERPINB9", "HBB") # BAL, neutrophils
```

Query GO and make a plot:

``` r
# query database
dbs <- c("GO_Molecular_Function_2015", "GO_Cellular_Component_2015", "GO_Biological_Process_2015")
if (websiteLive) {
    enriched <- enrichr(geneList, dbs)
}
```

    ## Uploading data to Enrichr... Done.
    ##   Querying GO_Molecular_Function_2015... Done.
    ##   Querying GO_Cellular_Component_2015... Done.
    ##   Querying GO_Biological_Process_2015... Done.
    ## Parsing results... Done.

``` r
# view results table
if (websiteLive) enriched[["GO_Biological_Process_2015"]]
```

    ##                                                                                                       Term
    ## 1                                                                cranial suture morphogenesis (GO:0060363)
    ## 2                                        regulation of ER to Golgi vesicle-mediated transport (GO:0060628)
    ## 3                                      negative regulation of fatty acid biosynthetic process (GO:0045717)
    ## 4                                                                            renal absorption (GO:0070293)
    ## 5                                                           craniofacial suture morphogenesis (GO:0097094)
    ## 6                                                         myeloid leukocyte mediated immunity (GO:0002444)
    ## 7                                                                            oxygen transport (GO:0015671)
    ## 8                                         negative regulation of fatty acid metabolic process (GO:0045922)
    ## 9                                                                               gas transport (GO:0015669)
    ## 10                                                        hydrogen peroxide catabolic process (GO:0042744)
    ## 11                                        negative regulation of steroid biosynthetic process (GO:0010894)
    ## 12                                                                   middle ear morphogenesis (GO:0042474)
    ## 13                                           negative regulation of steroid metabolic process (GO:0045939)
    ## 14                                  regulation of interferon-gamma-mediated signaling pathway (GO:0060334)
    ## 15                                                 regulation of response to interferon-gamma (GO:0060330)
    ## 16                                              regulation of fatty acid biosynthetic process (GO:0042304)
    ## 17                                                     cellular response to estrogen stimulus (GO:0071391)
    ## 18                                                                      bicarbonate transport (GO:0015701)
    ## 19                                   positive regulation of nitric oxide biosynthetic process (GO:0045429)
    ## 20                                                        hydrogen peroxide metabolic process (GO:0042743)
    ## 21                                                           cholesterol biosynthetic process (GO:0006695)
    ## 22                                          negative regulation of lipid biosynthetic process (GO:0051055)
    ## 23                                                                sterol biosynthetic process (GO:0016126)
    ## 24                                            negative regulation of fat cell differentiation (GO:0045599)
    ## 25                                                                       platelet aggregation (GO:0070527)
    ## 26                                            regulation of nitric oxide biosynthetic process (GO:0045428)
    ## 27                                                 regulation of steroid biosynthetic process (GO:0050810)
    ## 28                                                               homotypic cell-cell adhesion (GO:0034109)
    ## 29                                                                    inner ear morphogenesis (GO:0042472)
    ## 30                                                     cellular response to hydrogen peroxide (GO:0070301)
    ## 31                                                            regulation of blood vessel size (GO:0050880)
    ## 32                                                                    regulation of tube size (GO:0035150)
    ## 33                                             negative regulation of lipid metabolic process (GO:0045833)
    ## 34                                                interferon-gamma-mediated signaling pathway (GO:0060333)
    ## 35  negative regulation of cysteine-type endopeptidase activity involved in apoptotic process (GO:0043154)
    ## 36                                                                leukocyte mediated immunity (GO:0002443)
    ## 37                                negative regulation of cysteine-type endopeptidase activity (GO:2000117)
    ## 38                                                                       renal system process (GO:0003014)
    ## 39                                                 regulation of fatty acid metabolic process (GO:0019217)
    ## 40                                                    regulation of steroid metabolic process (GO:0019218)
    ## 41                                                                         palate development (GO:0060021)
    ## 42                                                              protein heterooligomerization (GO:0051291)
    ## 43                                                  reactive oxygen species metabolic process (GO:0072593)
    ## 44                                                     vascular process in circulatory system (GO:0003018)
    ## 45                                                             triglyceride metabolic process (GO:0006641)
    ## 46                                                     regulation of fat cell differentiation (GO:0045598)
    ## 47                                                      cellular response to interferon-gamma (GO:0071346)
    ## 48                                               cellular response to reactive oxygen species (GO:0034614)
    ## 49                                            negative regulation of protein complex assembly (GO:0031333)
    ## 50                                                             acylglycerol metabolic process (GO:0006639)
    ## 51                                                            neutral lipid metabolic process (GO:0006638)
    ## 52                                                               ER-nucleus signaling pathway (GO:0006984)
    ## 53                                          regulation of cytokine-mediated signaling pathway (GO:0001959)
    ## 54                                                              response to hydrogen peroxide (GO:0042542)
    ## 55                                                regulation of response to cytokine stimulus (GO:0060759)
    ## 56                                             negative regulation of intracellular transport (GO:0032387)
    ## 57                                                               alcohol biosynthetic process (GO:0046165)
    ## 58                                                              cholesterol metabolic process (GO:0008203)
    ## 59                                                               steroid biosynthetic process (GO:0006694)
    ## 60                                                               response to interferon-gamma (GO:0034341)
    ## 61                                                   regulation of lipid biosynthetic process (GO:0046890)
    ## 62                                                                   sterol metabolic process (GO:0016125)
    ## 63                                              cellular response to steroid hormone stimulus (GO:0071383)
    ## 64                                                                 circulatory system process (GO:0003013)
    ## 65                                                               regulation of blood pressure (GO:0008217)
    ## 66                                                        response to reactive oxygen species (GO:0000302)
    ## 67                                                      cellular response to oxidative stress (GO:0034599)
    ## 68                                                                       response to estrogen (GO:0043627)
    ## 69                                              organic hydroxy compound biosynthetic process (GO:1901617)
    ## 70                                            regulation of cellular ketone metabolic process (GO:0010565)
    ## 71           regulation of cysteine-type endopeptidase activity involved in apoptotic process (GO:0043281)
    ## 72                                         regulation of cysteine-type endopeptidase activity (GO:2000116)
    ## 73                                              negative regulation of endopeptidase activity (GO:0010951)
    ## 74                                                  negative regulation of peptidase activity (GO:0010466)
    ## 75                                                    regulation of anatomical structure size (GO:0090066)
    ## 76                                                       single organismal cell-cell adhesion (GO:0016337)
    ## 77                                              negative regulation of organelle organization (GO:0010639)
    ## 78                                                                  steroid metabolic process (GO:0008202)
    ## 79                                                     regulation of protein complex assembly (GO:0043254)
    ## 80                                                      regulation of lipid metabolic process (GO:0019216)
    ## 81                                                                          response to virus (GO:0009615)
    ## 82                                                       regulation of innate immune response (GO:0045088)
    ## 83                                                              single organism cell adhesion (GO:0098602)
    ## 84                                                               response to oxidative stress (GO:0006979)
    ## 85                                               cellular response to organic cyclic compound (GO:0071407)
    ## 86                                                             glycerolipid metabolic process (GO:0046486)
    ## 87                                                         negative regulation of proteolysis (GO:0045861)
    ## 88                                                   regulation of vesicle-mediated transport (GO:0060627)
    ## 89                                                                 cellular response to lipid (GO:0071396)
    ## 90                                                  negative regulation of protein processing (GO:0010955)
    ## 91                                                  negative regulation of protein maturation (GO:1903318)
    ## 92                                                                    organic anion transport (GO:0015711)
    ## 93                                                                  alcohol metabolic process (GO:0006066)
    ## 94                                                        cytokine-mediated signaling pathway (GO:0019221)
    ## 95                                                  negative regulation of hydrolase activity (GO:0051346)
    ## 96                                                       regulation of endopeptidase activity (GO:0052548)
    ## 97                                                           negative regulation of transport (GO:0051051)
    ## 98                                                                    protein oligomerization (GO:0051259)
    ## 99                                                                response to steroid hormone (GO:0048545)
    ## 100                                                           response to inorganic substance (GO:0010035)
    ## 101                                                          regulation of peptidase activity (GO:0052547)
    ## 102                                                       small molecule biosynthetic process (GO:0044283)
    ## 103                                                     regulation of intracellular transport (GO:0032386)
    ## 104                                                                   embryonic morphogenesis (GO:0048598)
    ## 105                                                                           anion transport (GO:0006820)
    ## 106                                                     cellular response to hormone stimulus (GO:0032870)
    ## 107                                                                response to other organism (GO:0051707)
    ## 108                                                               nitrogen compound transport (GO:0071705)
    ## 109                                                    cellular response to cytokine stimulus (GO:0071345)
    ## 110                                                                               coagulation (GO:0050817)
    ## 111                                                                         blood coagulation (GO:0007596)
    ## 112                                    negative regulation of cellular component organization (GO:0051129)
    ## 113                                                organic hydroxy compound metabolic process (GO:1901615)
    ## 114                                                                                hemostasis (GO:0007599)
    ## 115                                                                lipid biosynthetic process (GO:0008610)
    ##     Overlap     P.value Adjusted.P.value Old.P.value Old.Adjusted.P.value
    ## 1      1/10 0.001998613       0.03669208           0                    0
    ## 2      1/10 0.001998613       0.03669208           0                    0
    ## 3      1/10 0.001998613       0.03669208           0                    0
    ## 4      1/13 0.002597617       0.03669208           0                    0
    ## 5      1/15 0.002996803       0.03669208           0                    0
    ## 6      1/16 0.003196351       0.03669208           0                    0
    ## 7      1/16 0.003196351       0.03669208           0                    0
    ## 8      1/18 0.003595358       0.03669208           0                    0
    ## 9      1/20 0.003994245       0.03669208           0                    0
    ## 10     1/20 0.003994245       0.03669208           0                    0
    ## 11     1/21 0.004193643       0.03669208           0                    0
    ## 12     1/23 0.004592351       0.03669208           0                    0
    ## 13     1/23 0.004592351       0.03669208           0                    0
    ## 14     1/24 0.004791660       0.03669208           0                    0
    ## 15     1/25 0.004990939       0.03669208           0                    0
    ## 16     1/29 0.005787756       0.03669208           0                    0
    ## 17     1/30 0.005986886       0.03669208           0                    0
    ## 18     1/31 0.006185986       0.03669208           0                    0
    ## 19     1/32 0.006385056       0.03669208           0                    0
    ## 20     1/33 0.006584096       0.03669208           0                    0
    ## 21     1/34 0.006783106       0.03669208           0                    0
    ## 22     1/38 0.007578847       0.03669208           0                    0
    ## 23     1/39 0.007777708       0.03669208           0                    0
    ## 24     1/39 0.007777708       0.03669208           0                    0
    ## 25     1/40 0.007976539       0.03669208           0                    0
    ## 26     1/44 0.008771563       0.03879730           0                    0
    ## 27     1/53 0.010558622       0.04017826           0                    0
    ## 28     1/57 0.011352095       0.04017826           0                    0
    ## 29     1/60 0.011946886       0.04017826           0                    0
    ## 30     1/62 0.012343265       0.04017826           0                    0
    ## 31     1/63 0.012541409       0.04017826           0                    0
    ## 32     1/64 0.012739524       0.04017826           0                    0
    ## 33     1/65 0.012937609       0.04017826           0                    0
    ## 34     1/65 0.012937609       0.04017826           0                    0
    ## 35     1/69 0.013729650       0.04017826           0                    0
    ## 36     1/70 0.013927586       0.04017826           0                    0
    ## 37     1/71 0.014125493       0.04017826           0                    0
    ## 38     1/71 0.014125493       0.04017826           0                    0
    ## 39     1/72 0.014323369       0.04017826           0                    0
    ## 40     1/77 0.015312304       0.04017826           0                    0
    ## 41     1/78 0.015510002       0.04017826           0                    0
    ## 42     1/81 0.016102916       0.04017826           0                    0
    ## 43     1/84 0.016695563       0.04017826           0                    0
    ## 44     1/86 0.017090512       0.04017826           0                    0
    ## 45     1/87 0.017287942       0.04017826           0                    0
    ## 46     1/87 0.017287942       0.04017826           0                    0
    ## 47     1/88 0.017485343       0.04017826           0                    0
    ## 48     1/89 0.017682713       0.04017826           0                    0
    ## 49     1/91 0.018077365       0.04017826           0                    0
    ## 50     1/93 0.018471898       0.04017826           0                    0
    ## 51     1/94 0.018669120       0.04017826           0                    0
    ## 52     1/94 0.018669120       0.04017826           0                    0
    ## 53     1/95 0.018866312       0.04017826           0                    0
    ## 54     1/95 0.018866312       0.04017826           0                    0
    ## 55    1/101 0.020048841       0.04106862           0                    0
    ## 56    1/104 0.020639704       0.04106862           0                    0
    ## 57    1/104 0.020639704       0.04106862           0                    0
    ## 58    1/107 0.021230300       0.04106862           0                    0
    ## 59    1/108 0.021427107       0.04106862           0                    0
    ## 60    1/108 0.021427107       0.04106862           0                    0
    ## 61    1/113 0.022410692       0.04224966           0                    0
    ## 62    1/119 0.023590015       0.04375567           0                    0
    ## 63    1/123 0.024375638       0.04449521           0                    0
    ## 64    1/130 0.025749337       0.04626834           0                    0
    ## 65    1/133 0.026337621       0.04659733           0                    0
    ## 66    1/141 0.027905076       0.04862248           0                    0
    ## 67    1/148 0.029275046       0.05024821           0                    0
    ## 68    1/166 0.032791178       0.05465196           0                    0
    ## 69    1/166 0.032791178       0.05465196           0                    0
    ## 70    1/178 0.035129953       0.05771349           0                    0
    ## 71    1/186 0.036686779       0.05942225           0                    0
    ## 72    1/197 0.038824337       0.06201109           0                    0
    ## 73    1/222 0.043669195       0.06879394           0                    0
    ## 74    1/230 0.045215674       0.06915851           0                    0
    ## 75    1/232 0.045602001       0.06915851           0                    0
    ## 76    1/235 0.046181271       0.06915851           0                    0
    ## 77    1/237 0.046567304       0.06915851           0                    0
    ## 78    1/240 0.047146135       0.06915851           0                    0
    ## 79    1/244 0.047917499       0.06915851           0                    0
    ## 80    1/245 0.048110266       0.06915851           0                    0
    ## 81    1/250 0.049073667       0.06967249           0                    0
    ## 82    1/254 0.049843860       0.06990297           0                    0
    ## 83    1/259 0.050805944       0.07039378           0                    0
    ## 84    1/290 0.056754566       0.07704458           0                    0
    ## 85    1/291 0.056945990       0.07704458           0                    0
    ## 86    1/296 0.057902675       0.07742800           0                    0
    ## 87    1/302 0.059049736       0.07800050           0                    0
    ## 88    1/306 0.059813861       0.07800050           0                    0
    ## 89    1/315 0.061531440       0.07800050           0                    0
    ## 90    1/316 0.061722137       0.07800050           0                    0
    ## 91    1/316 0.061722137       0.07800050           0                    0
    ## 92    1/328 0.064008233       0.08001029           0                    0
    ## 93    1/340 0.066290148       0.08156444           0                    0
    ## 94    1/342 0.066670062       0.08156444           0                    0
    ## 95    1/354 0.068947107       0.08238425           0                    0
    ## 96    1/354 0.068947107       0.08238425           0                    0
    ## 97    1/361 0.070273458       0.08238425           0                    0
    ## 98    1/366 0.071219984       0.08238425           0                    0
    ## 99    1/369 0.071787552       0.08238425           0                    0
    ## 100   1/370 0.071976684       0.08238425           0                    0
    ## 101   1/372 0.072354860       0.08238425           0                    0
    ## 102   1/394 0.076507178       0.08542064           0                    0
    ## 103   1/394 0.076507178       0.08542064           0                    0
    ## 104   1/403 0.078201830       0.08647318           0                    0
    ## 105   1/443 0.085705414       0.09304245           0                    0
    ## 106   1/462 0.089253524       0.09304245           0                    0
    ## 107   1/462 0.089253524       0.09304245           0                    0
    ## 108   1/464 0.089626408       0.09304245           0                    0
    ## 109   1/471 0.090930599       0.09304245           0                    0
    ## 110   1/472 0.091116798       0.09304245           0                    0
    ## 111   1/472 0.091116798       0.09304245           0                    0
    ## 112   1/474 0.091489109       0.09304245           0                    0
    ## 113   1/476 0.091861306       0.09304245           0                    0
    ## 114   1/478 0.092233388       0.09304245           0                    0
    ## 115   1/491 0.094649139       0.09464914           0                    0
    ##     Odds.Ratio Combined.Score    Genes
    ## 1    740.25926     4600.93460   INSIG1
    ## 2    740.25926     4600.93460   INSIG1
    ## 3    740.25926     4600.93460   INSIG1
    ## 4    555.11111     3304.66577      HBB
    ## 5    475.76190     2764.27625   INSIG1
    ## 6    444.02222     2551.23866 SERPINB9
    ## 7    444.02222     2551.23866      HBB
    ## 8    391.74510     2204.78521   INSIG1
    ## 9    350.47368     1935.63138      HBB
    ## 10   350.47368     1935.63138      HBB
    ## 11   332.93333     1822.53879   INSIG1
    ## 12   302.63636     1629.20146   INSIG1
    ## 13   302.63636     1629.20146   INSIG1
    ## 14   289.46377     1545.99078   IFNGR2
    ## 15   277.38889     1470.19750   IFNGR2
    ## 16   237.71429     1224.70651   INSIG1
    ## 17   229.50575     1174.65261 SERPINB9
    ## 18   221.84444     1128.18302      HBB
    ## 19   214.67742     1084.93569      HBB
    ## 20   207.95833     1044.59515      HBB
    ## 21   201.64646     1006.88537   INSIG1
    ## 22   179.81081      877.90726   INSIG1
    ## 23   175.07018      850.22719   INSIG1
    ## 24   175.07018      850.22719   INSIG1
    ## 25   170.57265      824.07924      HBB
    ## 26   154.67442      732.57521      HBB
    ## 27   127.84615      581.80387   INSIG1
    ## 28   118.69048      531.53785      HBB
    ## 29   112.63842      498.68234   INSIG1
    ## 30   108.93443      478.72810      HBB
    ## 31   107.17204      469.27630      HBB
    ## 32   105.46561      460.15130      HBB
    ## 33   103.81250      451.33697   INSIG1
    ## 34   103.81250      451.33697   IFNGR2
    ## 35    97.68627      418.89804 SERPINB9
    ## 36    96.26570      411.42842 SERPINB9
    ## 37    94.88571      404.19171 SERPINB9
    ## 38    94.88571      404.19171      HBB
    ## 39    93.54460      397.17755   INSIG1
    ## 40    87.36842      365.12124   INSIG1
    ## 41    86.22944      359.25513   INSIG1
    ## 42    82.98333      342.61784      HBB
    ## 43    79.97189      327.29393      HBB
    ## 44    78.08235      317.73519      HBB
    ## 45    77.17054      313.13846   INSIG1
    ## 46    77.17054      313.13846   INSIG1
    ## 47    76.27969      308.65757   IFNGR2
    ## 48    75.40909      304.28833      HBB
    ## 49    73.72593      295.86912   INSIG1
    ## 50    72.11594      287.85112   INSIG1
    ## 51    71.33692      283.98403   INSIG1
    ## 52    71.33692      283.98403   INSIG1
    ## 53    70.57447      280.20727   IFNGR2
    ## 54    70.57447      280.20727      HBB
    ## 55    66.32000      259.28361   IFNGR2
    ## 56    64.37864      249.82381   INSIG1
    ## 57    64.37864      249.82381   INSIG1
    ## 58    62.54717      240.95208   INSIG1
    ## 59    61.95950      238.11647   INSIG1
    ## 60    61.95950      238.11647   IFNGR2
    ## 61    59.17857      224.77306   INSIG1
    ## 62    56.15254      210.39974   INSIG1
    ## 63    54.30055      201.68152 SERPINB9
    ## 64    51.33592      187.85590      HBB
    ## 65    50.16162      182.42560      HBB
    ## 66    47.27619      169.19896      HBB
    ## 67    45.00907      158.92792      HBB
    ## 68    40.06263      136.91786 SERPINB9
    ## 69    40.06263      136.91786   INSIG1
    ## 70    37.32392      124.98664   INSIG1
    ## 71    35.69550      117.98571 SERPINB9
    ## 72    33.67347      109.39527 SERPINB9
    ## 73    29.82655       93.39027 SERPINB9
    ## 74    28.77293       89.08994 SERPINB9
    ## 75    28.52092       88.06701      HBB
    ## 76    28.15100       86.56941      HBB
    ## 77    27.90960       85.59476   INSIG1
    ## 78    27.55509       84.16711   INSIG1
    ## 79    27.09602       82.32515   INSIG1
    ## 80    26.98361       81.87527   INSIG1
    ## 81    26.43507       79.68675   IFNGR2
    ## 82    26.01186       78.00592   IFNGR2
    ## 83    25.50129       75.98727      HBB
    ## 84    22.73010       65.21310      HBB
    ## 85    22.65057       64.90866 SERPINB9
    ## 86    22.26102       63.42145   INSIG1
    ## 87    21.81063       61.71046 SERPINB9
    ## 88    21.52022       60.61208   INSIG1
    ## 89    20.89384       58.25636 SERPINB9
    ## 90    20.82646       58.00402 SERPINB9
    ## 91    20.82646       58.00402 SERPINB9
    ## 92    20.04995       55.11217      HBB
    ## 93    19.32842       52.45180   INSIG1
    ## 94    19.21310       52.02906   IFNGR2
    ## 95    18.54863       49.60675 SERPINB9
    ## 96    18.54863       49.60675 SERPINB9
    ## 97    18.18148       48.27840   INSIG1
    ## 98    17.92785       47.36506      HBB
    ## 99    17.77899       46.83063 SERPINB9
    ## 100   17.72990       46.65469      HBB
    ## 101   17.63252       46.30605 SERPINB9
    ## 102   16.62680       42.73705   INSIG1
    ## 103   16.62680       42.73705   INSIG1
    ## 104   16.24710       41.40512   INSIG1
    ## 105   14.74661       36.23004      HBB
    ## 106   14.12509       34.13009 SERPINB9
    ## 107   14.12509       34.13009   IFNGR2
    ## 108   14.06263       33.92056      HBB
    ## 109   13.84823       33.20332   IFNGR2
    ## 110   13.81812       33.10286      HBB
    ## 111   13.81812       33.10286      HBB
    ## 112   13.75828       32.90341   INSIG1
    ## 113   13.69895       32.70590   INSIG1
    ## 114   13.64011       32.51029      HBB
    ## 115   13.26939       31.28362   INSIG1

``` r
# make plot
if (websiteLive) plotEnrich(enriched[[3]], showTerms = 20, numChar = 40, y = "Count", orderBy = "P.value")
```

![](GO_enrichment_files/figure-gfm/unnamed-chunk-8-1.png)<!-- -->
