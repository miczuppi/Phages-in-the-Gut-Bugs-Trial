
library(tidyverse)
library(stringr)

output_files <- list.files("path_to_folder_containing_abundance.txt_files") 
 gene_abundances <- tibble()
 for (file in output_files) {
   sampleID = str_extract(file, "[:alnum:]+_[:alnum:]+")
   gene_abundances <- 
     read_tsv(str_c("path_to_folder_containing_abundance.txt_files", file),
              col_names = c("geneID", "RPK"),
              skip = 1) %>% 
     mutate(sampleID = sampleID) %>% 
     bind_rows(gene_abundances)
 }

#Remove the genes with no mapping reads
  gene_abundances <-
      gene_abundances %>%
      filter(RPK != 0)


#Import the file containing the gene_ID clustering with the gene centroids




file_path <- "path_to_gene_cluster.clstr" # Replace with the path to your .clstr file

# Read the file into a character vector
clstr_output <- readLines(file_path, warn = FALSE)

if (length(clstr_output) == 0) {
  stop("The file is empty or not read correctly.")
}

# Initialize variables
cluster_number <- NA
rep_gene <- ""
output <- data.frame(ClusterNumber=character(), GeneName=character(), RepresentativeGeneName=character(), stringsAsFactors=FALSE)

# Process each line
for (line in clstr_output) {
  if (grepl("^>Cluster", line)) {
    # Extract cluster number
    cluster_number <- as.integer(sub("^>Cluster ", "", line))
  } else if (!is.na(cluster_number) && nchar(line) > 0) {
    # Extract gene name and check if it's the representative gene
    gene_name <- sub(".*>,", "", line)
    gene_name <- sub("\\....$", "", gene_name)
    if (grepl("\\*$", line)) {
      rep_gene <- gene_name
    }
    # Add to output
    output <- rbind(output, data.frame(ClusterNumber=cluster_number, GeneName=gene_name, RepresentativeGeneName=rep_gene))
  }
}

# Assuming 'output' is your data frame and 'GeneName' is the column
# Modify 'GeneName' to match the actual column name in your data frame
output$GeneName <- sub(".*>", "", output$GeneName)
output$RepresentativeGeneName <- sub(".*>", "", output$RepresentativeGeneName)
output$GeneName <- sub("\\..*", "", output$GeneName)
output$RepresentativeGeneName <- sub("\\..*", "", output$RepresentativeGeneName)


gene_clusters <- 
  output %>%
  rename(Seq.Name = GeneName) %>% 
  rename(Representative = RepresentativeGeneName) %>% 
  select(-ClusterNumber)


gene_abundances <- 
   gene_clusters %>% 
   rename(geneID = "Representative") %>% 
   left_join(gene_abundances, ., by = "geneID") %>% 
   select(-geneID) %>% 
   rename(geneID = "Seq.Name")
 
total_rpk_per_sample <-
   gene_abundances %>%
   group_by(sampleID) %>%
   summarise(total_rpk = sum(RPK))

gene_abundances_cpm <-
   left_join(gene_abundances, total_rpk_per_sample) %>%
   mutate(CPM = RPK / total_rpk * 1000000) %>%
   select(geneID, sampleID, CPM)

gene_abundances_cpm$contig_id <- str_extract(gene_abundances_cpm$geneID, "^.*(?=(\\_))") 

contig_median_cpm<- 
  gene_abundances_cpm %>% 
  group_by(contig_id) %>% 
  summarise(median= median(CPM)) %>% 
  distinct(contig_id, .keep_all = T) %>% 
  ungroup()

  #This (contig_median_cpm) is the file used for the analysis


