# Perturb-uCoTECH
A single-cell multiomics method for high-throughput joint profiling of sgRNAs, transcriptome and histone modifications/TFs.

# Single-cell read deconvolution
Read deconvolution is processed with UMI-tools.

# Pseudo-cell method to calculate perturbation effect
For accurate quantification of perturbations, our pipeline generated pseudo cells by clustering 50 neighboring single cells with the highest similarity in the transcriptome and classifying them into CRE j+ (with sgRNAs targeting to CRE j) and CRE j- cells, respectively. Perturbation effects of CRE j were expressed as fold change (FC) of gene expression and H3K27ac levels between CRE j+ and CRE j- pseudo cells.

