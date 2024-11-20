# Burn pit smoke condensate-mediated toxicity in human airway epithelial cells 

> Exploratory analyses using biomarker data (mRNAs and cytokines) to explore differences in modulations after exposure to various biomasss burn conditions

> All analyses in this respository are designated by their figure number or table number in the manuscript in parantheses.

In the instance that the files are unable to rendered the NBViewer link can be clicked below.

[Link to NBViewer](https://nbviewer.org/github/UNC-CEMALB/Burn-pit-smoke-condensate-mediated-toxicity-in-human-airway-epithelial-cells-/tree/main/)

This script was generated to support the manuscript titled 'Burn pit smoke condensate-mediated toxicity in human airway epithelial cells', published in 2024 by Ghosh et. al (PMID:38652897). doi: 10.1021/acs.chemrestox.4c00064

![image](https://github.com/UNC-CEMALB/Burn-pit-smoke-condensate-mediated-toxicity-in-human-airway-epithelial-cells-/assets/69641855/fd33457d-0822-4788-91fc-ffbe044d2a6b)

# 1. Data Processing
- Cytokine Data Processing: Used variance stabalizing normalization (VSN) to normalize the data
- mRNA Data Processing: Averaged any replicates, applied VSN, and imputed missing data using the Quantile Regression Imputation of Left-Censored (QRILC) technique

# 2. Biomarker Friedman Analysis (Table S2-S5)
- Testing for statistical differences across burn conditions (control, smoldering, and flaming) for each individual biomarker (mRNAs or cytokines) within low (1ug/cm2) and high (25ug/cm2) doses after various exposure durations (4, 24, or 72 hours)
- Used a Friedman's test (non-parametric test for 3+ group comparisons) followed by Nemenyi's post hoc tests

# 3. Biomarker Cluster Analysis (Figures 4,5,7)
- Using Principal Component Analysis (PCA) to determine if and how exposure groups cluster based upon gene expression
- Using hierarchical clustering in heatmap visualizations to determine if and how genes cluster based upon abundance across exposure groups

# 4. Significant Biomarkers Visualization (Figure 6)
- Venn diagrams to visualize significant (p value < 0.1) gene expression modulations across exposure groups at the higher dose (25 ug/cm^2)

# 5. Correlation Analysis (Figure 8)
- Determining correlations between cytokine and PAHs using Spearman's rank correlation tests
