This is the team repository of Team Quaranteam.

Our team proposal can be found [here.](https://github.com/trevorkwan/Identifying-Differentially-Expressed-Genes-Project-STAT-540/blob/main/reports/Proposal/project_proposal.md)

Our team progress report can be found [here.](https://github.com/trevorkwan/Identifying-Differentially-Expressed-Genes-Project-STAT-540/blob/main/reports/Progress%20Report/progressreport.md)

Our team presentation slides can be found [here.](https://github.com/trevorkwan/Identifying-Differentially-Expressed-Genes-Project-STAT-540/blob/main/final_presentation/Quaranteam.pdf)

Project Summary:

In this project, we analyzed a recent scRNA-seq dataset producted by [Bost et al. 2021. Nature Communications, "Deciphering the State of Immune Silence in Fatal COVID-19 Patients"](https://doi.org/10.1038/s41467-021-21702-6)

In short, we sought to investigate potential differential immunological and transcriptional responses between male and female patients infected with SARS-CoV-2. With the global spread of this virus, it has become apparent that a number of factors are associated with increased infection and mortality rates of COVID-19. One of these factors is the sex of the patient, with many studies reporting that males have higher morbidity and mortality rates of the disease. We hypothesized that a reason for this lies within the immunological response of male individuals. There are a number of sex-linked immunological traits related to the X and Y chromosome, and specific immune related genes have been implicated as well (including Interleukin genes, ACE2 receptor - SARS-CoV-2 entry receptor, and endrogen receptor activity). 

To investigate our hypothesis, we analyzed 54 samples from patients infected with SARS-CoV-2 (minus control patients) that have had scRNA-seq done on blood and lung tissues. The patients were characterized by their degree of severity of disease (healthy control, mild, and severe). To analyze this dataset, we first merged all of the scRNA data together and filtered to remove low quality cells. Count data was also normalized, scaled, and cells were clustered using UMAP/tSNE. See the script [here](https://github.com/trevorkwan/Identifying-Differentially-Expressed-Genes-Project-STAT-540/blob/main/results/Preprocessing-and-Filtering-scRNA-Data.md).

After filtering, we identified differential cell markers in each of the cell clusters to identify the specific cell type. See the pdf [here](https://github.com/trevorkwan/Identifying-Differentially-Expressed-Genes-Project-STAT-540/blob/main/results/Aim%202%20-%20Finding%20differentially%20expressed%20features%20(cluster%20biomarkers).pdf),and the Rmd [here](https://github.com/trevorkwan/Identifying-Differentially-Expressed-Genes-Project-STAT-540/blob/main/src/Aim%202%20-%20Finding%20differentially%20expressed%20features%20(cluster%20biomarkers).Rmd)

Once we identified our specific cell clusters, we performed differential gene expression analysis. We subset by various things, such as disease severity, tissue type, and cell cluster. See the script [here](https://github.com/trevorkwan/Identifying-Differentially-Expressed-Genes-Project-STAT-540/blob/main/results/Differential_Expression_Analysis.md), and pdf [here](https://github.com/trevorkwan/Identifying-Differentially-Expressed-Genes-Project-STAT-540/blob/main/results/Differential_Expression_Analysis.pdf)

Lastly, we performed GO analysis on some of the DE genes to identify specific pathways that may be biologically relevant as to why males have increased morbidity/mortality. See Rmd [here](https://github.com/trevorkwan/Identifying-Differentially-Expressed-Genes-Project-STAT-540/blob/main/src/DE_plots_GO_enrichment.Rmd) and md [here](https://github.com/trevorkwan/Identifying-Differentially-Expressed-Genes-Project-STAT-540/blob/main/results/GO_enrichment.md)

We were able to identify differentially expressed genes between male and female patients across different cell types. We as of yet need to undergo more research as to the biological relevance of these differentially expressed genes, and whether they offer an interesting pathway for why male patients are more prone to infection and mortality.
