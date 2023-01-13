## **Team 8: Quaranteam**


### **Motivation and Background**

SARS-CoV-2, the virus that causes COVID-19, has infected 110 million people worldwide to date and led to 2.4 million deaths1. Therefore, it is important to keep close watch on the strategies that can effectively and quickly predict the disease progression. 


SARS-CoV-2, as a single stranded RNA virus, has approximate 30 kb genomic RNA (gRNA) which is under frequent mutations2,3. According to previous research, both epithelial cells and immune cells can be infected with SARS-CoV-24 [4]. These immune cells include macrophages, natural killers (NK) plasma cells, T cells and many more. In response to SARS-CoV-2 infection, infected cells promote the secretion of large amounts of chemokines and cytokines4. Clinally, approximately 20% of laboratory-confirmed cases appear to be severe symptoms with acute respiratory distress syndrome, organ failure, while the remaining 80% cases have moderate symptoms5,6. The immune patterns can be closely related to disease progression of SARS-CoV-2 positive patients5. Therefore, the investigation on the changes in immune-related gene expression levels during SARS-CoV-2 infection can be of help to find out markers for predicting disease progression and death outcome. 




### **Division of Labour**

| Name | Background | Affiliation | Responsibilities |
| --- | --- | --- | --- |
| Madeline Iseminger | Biophysics | Bioinformatics | Data collection, Aim #4:DE genes of SARS-CoV-2 post-infected and negative patients |
| Jackson Moore | Microbiology | Genome Science and Technology | Data validation, Aim #5: Interaction effects between metadata |
| Trevor Kwan | Data Science | Biostatistics | Data visualization, Aim #3: DE genes of SARS-CoV-2 positive and negative patients |
| Yilin Qiu | Metagenomics | Genome Science and Technology | Data cleaning, Aim #6: GO analysis & Immunological pathways |
	

### **Dataset**


We will be using a single cell gene expression dataset called “Large-scale single-cell analysis reveals critical immune characteristics of COVID-19 patients” (GSE158055)7.  This dataset was accepted to GEO on February 4, 2021, and is associated with a preprint8.  It includes single cell RNAseq data from 171 individuals with mild/moderate, severe or recovered COVID-19 disease (22, 54, and 95 individuals respectively) and 25 controls.  There are three tissue types represented in the dataset: pleural effusion including sputum, human peripheral blood mononuclear cells, and bronchoalveolar lavage fluid, making a total of 284 samples.  The scRNA-seq data was generated using a 5’ single cell sequencing platform from 10X Genomics, with a protocol that included generating T cell receptor (TCR) or B cell receptor (BCR) data8.


The counts are presented as a sparse matrix; raw data is not publicly available because of privacy policies.  The metadata contains 74 variables: patient ID, number of days after symptom onset the sample was taken, sample type (eg. "fresh Sputum"), sampling time (eg. progression), patient outcome (eg. discharge from hospital), sex and age of individual, geographic location, cell type, SARS-CoV-2 status (positive or negative) and severity, single cell sequencing platform ID, immune cell concentration, known comorbidities, and medication taken by each individual.


The data is provided in the form of several files: a csv of barcodes, another csv of features, a count matrix, and a metadata file.


### **Aims and Methodology**

* **Processing of data for differential gene analysis.**

Packages such as dplyr and tidyr will be used to process, filter, and transform our data to a format that incorporates patient metadata and gene expression data ready for analysis.


* **Plotting gene expression for all the genes as a sanity check.**

The package ggplot2 will be used to compare and plot gene expression for all the genes between SARS-CoV-2 positive and negative (control) patients, amongst different immune cell types to see whether the gene expression distributions line up. 


* **Identify differentially expressed genes (DE) between SARS-CoV-2 positive and negative (control) patients, within different immune cell types.**


Bioconductor packages such as limma will be used to identify differentially expressed genes in the SARS-CoV-2 positive patients relative to control patients (that test negative for the virus). The Welch 2 Sample t-test will be used to test if SARS-CoV-2 positive patients are significantly different in their expression values compared to SARS-CoV-2 negative patients. We will perform differential gene expression subsetting on each of the following cell-types: CD4+ T cells, CD8+ T cells, B cells, natural killer cells, epithelial cells, and myeloid cells.


* **Compare the transcriptional profiles of patients sampled post-infection to transcriptional profiles of non-infected patients in their expression levels.**


The limma package will be used to identify DE genes between patients recovered from SARS-CoV-2 infection (convalescence) and control patients that are negative for the virus. The Welch 2 Sample t-test will be used to test if the expression values of recovered patients are significantly different from the expression values of infected patients. We will examine whether Covid-19 infection results in a long-term transcriptional response that continues post-infection.


* **Examine interaction effects between metadata factors: (patient age, sex, symptom severity, comorbidities, medication taken, and immune cell concentration), on expression values of DE genes in SARS-CoV-2 positive and negative patients.**

The patient metadata includes a number of factors that may interact with each other on gene expression during SARS-CoV-2 infection. Using the limma package, we will examine whether these factors (listed above) interact with each other on DE genes between infected and non-infected control patients. For example, we will see if some DE genes are up-regulated over factor A for one level of factor B, but down-regulated over factor A for another level of factor B.


* **Determine the immunological pathways that are implicated in the above analysis, and how they contribute to Covid-19 disease response and progression.**


We will perform GO analysis on the genes identified in the above aims to investigate what immunological pathways are involved during SARS-CoV-2 infection. We will also examine how these immunological pathways may differ depending on interaction effects, and whether these pathways remain altered post-infection.

### **References**

1. [Dong, E., Du, H. & Gardner, L. An interactive web-based dashboard to track COVID-19 in real time. Lancet Infect. Dis. 20, 533–534 (2020).](https://pubmed.ncbi.nlm.nih.gov/32087114/)
2. [Zhao, Z. et al. Moderate mutation rate in the SARS coronavirus genome and its implications. BMC Evol. Biol. 4, 21 (2004).](https://bmcecolevol.biomedcentral.com/articles/10.1186/1471-2148-4-21)
3. [Pal, M., Berhanu, G., Desalegn, C. & Kandi, V. Severe acute respiratory syndrome Coronavirus-2 (SARS-CoV-2): An update. Cureus 12, e7423 (2020).](https://www.cureus.com/articles/29589-severe-acute-respiratory-syndrome-coronavirus-2-sars-cov-2-an-update)
4. [Huang, C. et al. Clinical features of patients infected with 2019 novel coronavirus in Wuhan, China. Lancet 395, 497–506 (2020).](https://pubmed.ncbi.nlm.nih.gov/31986264/)
5. [Liu, T., Jia, P., Fang, B. & Zhao, Z. Differential Expression of Viral Transcripts From Single-Cell RNA Sequencing of Moderate and Severe COVID-19 Patients and Its Implications for Case Severity. Front. Microbiol. 11, 603509 (2020).](https://pubmed.ncbi.nlm.nih.gov/33178176/)
6. [Chen, N. et al. Epidemiological and clinical characteristics of 99 cases of 2019 novel coronavirus pneumonia in Wuhan, China: a descriptive study. Lancet 395, 507–513 (2020).](https://pubmed.ncbi.nlm.nih.gov/32007143/)
7. [GEO Accession viewer.](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE158055)
8. [Ren, X. et al. Large-scale single-cell analysis reveals critical immune characteristics of COVID-19 patients. Cold Spring Harbor Laboratory 2020.10.29.360479 (2020) doi:10.1101/2020.10.29.360479.](https://www.biorxiv.org/content/10.1101/2020.10.29.360479v1)
