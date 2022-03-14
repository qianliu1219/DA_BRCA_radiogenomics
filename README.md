# MRI-genomic mapping of invasive breast carcinoma using stacked deep convolutional denoise autoencoder

This repository is hosting the scource code of the Keras-based denoisie autoencoder and the following radiogenomic assciation analyses.

The raw data was dowloaded from [The Cancer Image Archive (TCIA)](https://wiki.cancerimagingarchive.net/display/Public/TCGA-BRCA) and [The Cancer Genome Atlas (TCGA)](https://portal.gdc.cancer.gov/) BRCA project. Please note that we do not have the authority to share the raw data here. The sample size and data exclusion criteria are shown in below figure.


![image1](https://github.com/qianliu1219/DA_BRCA_radiogenomics/blob/main/Figures/Figure1.png)

The raw images have different resolutions. The resolution distribution of the raw images are shown in below bar plot. To keep the size comparable and to save the computational costs, we compressed the images to 64X64 (the lowest resolution among the raw images).

![image2](https://github.com/qianliu1219/DA_BRCA_radiogenomics/blob/main/Figures/imageSize.png)

A Karas-based stacked convolutional denoise autoencoder was build to extract 4,096 deep radiomic features (DRFs) from each image. The code of the model is in `DA.py`. An illustration is shown below.

![image3](https://github.com/qianliu1219/DA_BRCA_radiogenomics/blob/main/Figures/Figure2.jpg)

Then we did sample-wise quantile normalization and visulization for the DRFs. 

![image4](https://github.com/qianliu1219/DA_BRCA_radiogenomics/blob/main/Figures/density.png)

![image5](https://github.com/qianliu1219/DA_BRCA_radiogenomics/blob/main/Figures/Figure3.png)

We also did unsupervised custering (hierarchical clustering and t-SNE clustering) on the DRFs. 

![image6](https://github.com/qianliu1219/DA_BRCA_radiogenomics/blob/main/Figures/Figure4.png)


We built LASSO model to classify 5 clinical features using the normalized DRFs and compared the classifcation performance with the baseline semi-auto radiomic features(RFs) which were provided by the TCGA Breast Phenotype Research Group. The R code could be found in `biglasso.R`. The performance could be found in below figure.

![image7](https://github.com/qianliu1219/DA_BRCA_radiogenomics/blob/main/Figures/Figure5.png)


We then did the radiogenomic analysis using linear mixed effect (LME) model to evaluate the association between genomic profiles and DRFs. We focused on three-level of gene expression features: the expression of  288 breast cancer risk genes, the value of 6 commonly used breast cancer gene signatures, and 182 KEGG pathway activity scores. The R code could be found in `lme.R`. The results were visulized as below. a is the association results of 288 risk genes and 4,096 DRFs. b is the association results of 288 risk genes and 36 traditional semi-auto RFs. c is the results of 6 gene signatures with DRFs. d is the results of 182 pathway activity scores with DRFs.

![image8](https://github.com/qianliu1219/DA_BRCA_radiogenomics/blob/main/Figures/Figure6.png)
