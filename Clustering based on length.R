

#Cluster UViGs into vOTUs based on lenght.
#We cluster the contigs based on length. For each group of clustering contig, we select the longest one as representative of the vOTU.

#Import FastANI output####

fastani <- read_tsv("fastani_output", header = F)

#Rename the columns
fastani <- 
  fastani %>%
  rename(query = "V1") %>% 
  rename(reference = "V2") %>% 
  rename(ANI = "V3") %>% 
  mutate(AF = V4/V5) %>% 
  select(c(query, reference, ANI, AF))

#Polish the output
fastani$query <- str_sub(fastani$query, start = 3, end = -7)
fastani$reference <- str_sub(fastani$reference, start = 3, end = -7)

#The contigs were already clustered on 85% AF, select the ones with ANI >= 95%
fastani <- 
  fastani %>% 
  filter(ANI >= 95)

#Add the length
#This is necessary to cluster based on length

#Add the length to the query
length_query <- 
  UViG_length %>% 
  select(c(contig_id, length)) %>% 
  rename(query = "contig_id") %>% 
  rename(length_query = "length")

#Add the length to the reference
length_reference <- 
  UViG_length %>% 
  select(c(contig_id, length,)) %>% 
  rename(reference = "contig_id") %>% 
  rename(length_reference = "length")

#Add the length information to the FastANI output
Length_vOTU <- 
  fastani %>% 
  left_join(., length_query, by = "query") %>% 
  left_join(., length_reference, by = "reference") 

#Sort the clusters by length#####
clusters_sorted_by_length <- 
  Length_vOTU %>% 
  #Remove contigs aligning with themselves
  filter(query != reference) %>% 
  #Select the query contigs longer than or equal to the reference ones
  filter(length_query >= length_reference)

#Select the longer contigs in each vOTU

#Select the contigs that were never the longest in the cluster, either shorter or equal to in size

shorter_contigs <- 
  clusters_sorted_by_length %>% 
  select(reference) %>% 
  rename(query = "reference") %>% #rename them for anti_join() the next line
  #Just one observation per contig
  filter(!duplicated(query))  

#Contigs that are only present in the reference column. This are the shortest contigs that are never longer than others

only_short <- anti_join(shorter_contigs, clusters_sorted_by_length, by = "query") 

#These are the contigs that are shorter than some but also longer than others. 
ref_query <- anti_join(shorter_contigs, only_short, by = "query")

## Identify the representative contigs based on length!
#(These are the contigs that are only longer than other)

representative <- anti_join(clusters_sorted_by_length, ref_query, by ="query") 

#Now, we were working only with the longest ones. Now some contigs might have the same length as others

same_length <- 
  clusters_sorted_by_length %>% 
  filter(length_query == length_reference) 

#Among the contigs with the same length, select just one, based on alphabetical order
representative_same_length <- 
  same_length %>%
  select(c(query, reference)) %>% 
  #Create a combined column
  mutate(combo = pmap_chr(across(), ~paste(sort(c(...)), collapse = "|"))) %>% 
  select(combo) %>%
  distinct() %>% 
  #Separate combined columns into individual columns
  separate(combo, into = c("contig_id"), sep = "\\|") %>% 
  #Distinct names
  distinct() %>% 
  #Rename for downstream analysis
  rename(query = "contig_id")

#Find the contigs with the same length as the ones just selected

same_length_contigs <- 
  representative_same_length %>% 
  rename(reference = "query") %>% 
  anti_join(fastani, ., by = "reference")

#Create the table reporting the representative contigs and all the other contigs clustering with the with the same length 
representative_same_length <-
  representative_same_length 
  left_join(., same_length_contigs, by = "query") %>%
 #Remove the contigs aligning with themselves
  filter(query != reference)

#Remove, among the representative with same length, those ones that were clustering with longer representative        
representative_same_length <-  
  representative_same_length %>%
  rename(temp_ref = "reference") %>% 
  rename(reference = "query") %>% 
  anti_join(., representative, by = "reference")

#Rename for downstream analysis
representative_same_length <- 
  representative_same_length %>% 
  rename(query = "reference") %>%
  rename(reference = "temp_ref") 

#Create the final representative table

representative <- bind_rows(representative, representative_same_length)

#modify the representative table by adding the query matching against themselves

representative <- 
  representative %>% 
  filter(query != reference)


same_cluster_matching <- 
  representative %>% 
  select(query) %>% 
  distinct() %>% 
  mutate(reference = query)

representative <- 
  bind_rows(representative, same_cluster_matching)

vOTUs <- 
  representative %>% 
  select(-c(length_query, length_reference, AF)) %>% 
  mutate(ANI = replace_na(ANI, 100))

#create the table with all the contigs that did not cluster in any votu (were not included in the FastANI output)

singleton_vOTUs <- 
   UViG_length %>% 
  select(contig_id) %>% 
  mutate(reference = contig_id) %>% 
  mutate(ANI = 100) %>% 
  anti_join(., representative, by = c("contig_id" = "reference"))

singleton_vOTUs <- 
  singleton_vOTUs %>% 
  rename(query = contig_id)

vOTUs_GBT <- bind_rows(vOTUs, singleton_vOTUs)

#This (vOTUs_GBT) is the final table reporting the representative contigs of each vOTU (query) and the UViGs (reference) clustering with them


