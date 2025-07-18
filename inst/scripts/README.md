**Data Generation for MultiModalGraphics DatasetsData Generation for MultiModalGraphics Datasets**

This document describes how each dataset in inst/extdata/ was generated, including data sources, processing, and relevant context. The MultiModalGraphics package leverages publicly available biomedical data from major consortia and databases. 

***Informative Heatmap dataset***

- Pan_cancer_CESC_Mutated_CNV_STR_VART.csv
- Pan_Cancer_DNA_methylation.csv
- Pan_cancer_miRNA-mRNA_interaction_data.csv
- Pan_cancer_protein-mRNA_combined_data.csv

Datasets were obtained from [The Cancer Genome Atlas (TCGA)](https://www.google.com/url?sa=t&rct=j&q=&esrc=s&source=web&cd=&cad=rja&uact=8&ved=2ahUKEwj6j7SEidKDAxUJ_bsIHS0LAJcQFnoECAYQAQ&url=https%3A%2F%2Fwww.cancer.gov%2Fccg%2Fresearch%2Fgenome-sequencing%2Ftcga&usg=AOvVaw2DAuLgtHHkBJrTRkFONrpi&opi=89978449https://www.google.com/url?sa=t&rct=j&q=&esrc=s&source=web&cd=&cad=rja&uact=8&ved=2ahUKEwj6j7SEidKDAxUJ_bsIHS0LAJcQFnoECAYQAQ&url=https%3A%2F%2Fwww.cancer.gov%2Fccg%2Fresearch%2Fgenome-sequencing%2Ftcga&usg=AOvVaw2DAuLgtHHkBJrTRkFONrpi&opi=89978449), [International Cancer Genome Consortium](https://dcc.icgc.org/) and [cBioPortal for cancer genomics](https://www.google.com/url?sa=t&rct=j&q=&esrc=s&source=web&cd=&cad=rja&uact=8&ved=2ahUKEwiX4dvoidKDAxXegv0HHcDRBmsQFnoECAgQAQ&url=https%3A%2F%2Fwww.cbioportal.org%2F&usg=AOvVaw099aWQoFJloiAjhgTScY7c&opi=89978449). Highly expressed messenger RNAs, proteins, microRNAs and probes of methylation sites across the six-cancer types were significantly associated with cellular proliferations and antiapoptotic signaling pathways. The genetic datasets (copy number and structural variations) showed mutation of these genes though to a different extent of mutational frequencies and significances. ‘The six cancer types include CESC: Cervical squamous cell carcinoma and endocervical adenocarcinoma; OV: Ovarian serous cystadenocarcinoma; PRAD: Prostate adenocarcinoma; TGCT: Testicular Germ Cell Tumors; UCEC: Uterine Corpus Endometrial Carcinoma; UCS: Uterine Carcinosarcoma.’

***CompositeFeatureHeatmap dataset***

- scatterexample.csv

These datasets were downloaded from [Gene Expression Omnibus GSE45035](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE45035). Brain regions include AY: amygdala, HC: hippocampus, MPFC: medial prefrontal cortex, SE: septal region, ST: corpus striatum, and VS: ventral striatum.’


***ThresholdedScatterplot dataset***

- twoD_graphics.csv

Differentially expressed genes (DEGs) across spleen, heart, blood and brain regions at the five timepoints (T5R1, T10R1, T5R10, T10R28 and T10R42). The plots show both negative log10 p-values and log2 -fold changes at five timepoints (across) and from spleen, heart, blood and seven brain regions (down); the numbers with red and blue fonts show the numbers of transcripts that were up- and down-regulated in each group. ‘Keys – AY: amygdala, HC: hippocampus, MPFC: medial prefrontal cortex, SE: septal region, ST: corpus striatum, and VS: ventral striatum, T: number of trauma exposure days, R: post-trauma sample collection days.’ The method used to generate this scatter plot can be used to generate any scatter plot with any number of multi-datasets (which can only be limited by graphical canvas size and/or legibility). Showing the number of features that passed the given cutoff makes this tool unique, since this is an important piece of information that is not available with other scatter plot making tools (also it comes with a lot of flexibility for the user for adjusting the cutoff, the colors, and whole lot of arguments can be supplied by the user to customize the plot according to their need or taste). It is easier to see which have the largest changes across tissues or time points (both visually and the actual values). These datasets were downloaded from [Gene Expression Omnibus GSE68077](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE68077).
