# Perturb-uCoTECH
New technology for Jointly profiling transcriptome, epigenome and perturbations in one single cell.

# Single-cell read decovolution
Read deconvolution is processed with UMI-tools.

# Pseudo-cell method to calculate perturbation effect
Given a heterogeneity in hESCs infected at high MOI, for accurate quantification of perturbations per cell, we developed a computational pipeline to generate pseudo cells consisting of a certain number of neighboring single cells with a highest similarity in the transcriptome. Three sgRNAs targeting one given element were summed in sgRNA(j). 50 neighbors on WNN-UMAP were selected and classified into sgRNA(j)+ and sgRNA(j)- cells, respectively. Perturbation effects of sgRNA(j) were expressed as fold change (FC) in RNA expression or H3K27ac of sgRNA(j)+/sgRNA(j)- pseudo cells.

