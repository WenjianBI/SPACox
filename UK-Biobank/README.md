
## This folder is to store results of UK-Biobank analysis for 21 cancers

* NoSmoking: pheno~PC1+PC2+PC3+PC4+strat(Gender). Note: for some gender-specific cancers, only the top 4 PCs were adjusted for.
* Smoking:   add smoking as a time-varying covariate

* Either directory includes 
1. A CSV file to summarize all SNPs with p values less than 5E-8
2. QQ plots of p values from both Score and ESPA methods
3. Manhattan plots of p values from both Score and ESPA methods
