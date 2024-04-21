# Phages-in-the-Gut-Bugs-Trial

In this repository, we provide the custom codes used to characterize the gut phageome of the participants in the Gut Bugs Trial, as outlined in the manuscript Fecal Microbiota Transplantation alters gut phage communities in a clinical trial for obesity. 

Clustering of UViGs into vOTUs

The UViGs were clustered into vOTUs according to MIUViG standards, _i.e._ 95% ANI over 85% AF of the shortest contig ([Roux et al, 2019](https://www.nature.com/articles/nbt.4306)). Clustering of the UViGs was performed based on length using a [custom R script](https://github.com/miczuppi/Phages-in-the-Gut-Bugs-Trial/blob/main/Clustering%20based%20on%20length.R), and the longest UViG of the vOTU was selected as the representative. 

For the identification of engrafted vOTUs, the same [custom R script](https://github.com/miczuppi/Phages-in-the-Gut-Bugs-Trial/blob/main/Clustering%20based%20on%20length.R) was used, with the addition of the "engraftment clustering" snippet, which divided the contigs by sex and treatment group. 

