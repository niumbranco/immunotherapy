# create phyloseq objects and rarefy counts for OTU and ITS analysis
library(tidyverse)
library(phyloseq)
library(readxl)
#library(DESeq2) # TODO to incorporate later

root_dir <- "~/Repositories/nium/metasoil"
out_dir <-  "~/Repositories/nium/metasoil/results"

# load data table
sampleinfo = read_tsv(file.path(root_dir, "processed_tables","sampleIDs_info_conversion.tsv"),
                       show_col_types = FALSE)
## read metadata
metadat = read_excel(file.path(root_dir, "data", "MASTER_meta_20220413.xlsx"))

samplemetadata = metadat %>% 
  select(-subsampleID, -Malte_order, -Kate_OTU_order, -clean_name, -RNA_DNA, -`16SNum_reads`) %>% 
  # redo column Lumgmt combination of landuse and management
  mutate(Lumgmt = paste0(gsub("^([A-Z]).+","\\1", LandUse), Mgmt)) %>% 
  distinct() %>% 
  group_by(sampleID) %>% 
  slice_max(pH_CaCl2) # some columns contain pH_CaCl2 and some dont

### read other tables
taxinfo_otus = read_tsv(file.path(root_dir, "processed_tables","taxonomy_table_OTUs.tsv"),
                        show_col_types = FALSE)
otutab_wide_D_f = read_tsv(file.path(root_dir, "processed_tables", "16S_OTUs_rawcounts_DNA.tsv"),
                           show_col_types = FALSE)
otutab_wide_R_f = read_tsv(file.path(root_dir, "processed_tables", "16S_OTUs_rawcounts_RNA.tsv"),
                           show_col_types = FALSE)

taxinfo_its = read_tsv(file.path(root_dir, "processed_tables", "ITS_taxonomy_table.tsv"),
                       show_col_types = FALSE)
itstab = read_tsv(file.path(root_dir, "processed_tables", "ITS_zOTUs_rawcounts.tsv"),
                   show_col_types = FALSE)

###################
#function to "normalize" with phyloseq function, input are tsv, first column is the ID column, output is a list of objects
# samplinfo used as a proxy

rarefied_physeq = function(otutab, taxtab, sampleinfo, cutoff) {
  removed = colnames(otutab)[-1][colSums(otutab[-1], na.rm=T)< cutoff]
  # create phyloseq obj
  TAX = taxtab %>% 
    column_to_rownames( colnames(taxtab)[1]) %>% 
    as.matrix() %>% 
    tax_table()
  OTU =  otutab %>% 
    select(!ends_with(removed )) %>%  # filter out samples
    column_to_rownames(colnames(otutab)[1]) %>% 
    as.matrix() %>% otu_table( taxa_are_rows = T)
  OTU_all =  otutab %>% 
    #select(!ends_with(removed )) %>%  # filter out samples
    column_to_rownames(colnames(otutab)[1]) %>% 
    as.matrix() %>% otu_table( taxa_are_rows = T)
  SAMPLE = sampleinfo %>% 
    column_to_rownames(colnames(sampleinfo)[1]) %>% 
    sample_data()
  phyobj = phyloseq(OTU, TAX,SAMPLE)
  
  # rarefaction
  phyobj_rarefied = rarefy_even_depth(phyobj, rngseed = 1, sample.size = cutoff, replace=T) 
  phyobj_all = phyloseq(OTU_all, TAX, SAMPLE)
  #deseq norm
  #phyobj_deseq = phyloseq_to_deseq2(phyobj, ~val)
  
  # return(list(removed, phyobj, phyobj_rarefied, phyobj_deseq))
  return(list(removed, phyobj, phyobj_rarefied, phyobj_all))
}

######################## ITS
its_outlist = rarefied_physeq(itstab, taxinfo_its, samplemetadata, 10000)
removed10k =  its_outlist[[1]]
its_phyobj = its_outlist[[2]]
its_phyobj_10k = its_outlist[[3]]
its_phyobj_all = its_outlist[[4]]

####################### OTUs DNA
otus_D_outlist = rarefied_physeq(otutab= otutab_wide_D_f, taxtab=  taxinfo_otus, sampleinfo= samplemetadata, cutoff= 30000)
removed30k_D =  otus_D_outlist[[1]]
D_phyobj = otus_D_outlist[[2]]
D_phyobj_30k = otus_D_outlist[[3]]
D_phyobj_30k_all = otus_D_outlist[[4]]

###################### OTUs RNA
otus_R_outlist = rarefied_physeq(otutab= otutab_wide_R_f, taxtab=  taxinfo_otus, sampleinfo= samplemetadata, cutoff= 30000)
removed30k_R =  otus_R_outlist[[1]]
R_phyobj = otus_R_outlist[[2]]
R_phyobj_30k = otus_R_outlist[[3]]
R_phyobj_30k_all = otus_R_outlist[[4]]

###################### make tsvs tables and write out
its_counts_out = otu_table(its_phyobj_10k) %>% 
  as.data.frame() %>% 
  rownames_to_column("zOTU")
otuD_counts_out = otu_table(D_phyobj_30k) %>% 
  as.data.frame() %>% 
  rownames_to_column("OTU")
otuR_counts_out = otu_table(R_phyobj_30k) %>% 
  as.data.frame() %>% 
  rownames_to_column("OTU")

write_tsv(its_counts_out, file.path(root_dir, "results","ITS_counts_rarefied10k.tsv"))
write_tsv(otuD_counts_out, file.path(root_dir, "results","OTU_DNA_counts_rarefied30k.tsv"))
write_tsv(otuR_counts_out, file.path(root_dir, "results","OTU_RNA_counts_rarefied30k.tsv"))      

# add removed samples to samplesinfo
sampleinfo_f = sampleinfo %>% 
  mutate(ITS_removed_rarefy10k = ifelse(sampleID %in% removed10k, "remove", NA)) %>% 
  mutate(OTU_DNA_rarefy30k = ifelse(grepl("_D($|_)",subsampleID) &  sampleID %in% removed30k_D, "remove", NA)) %>% 
  mutate(OTU_RNA_rarefy30k = ifelse(grepl("_R($|_)",subsampleID) &  sampleID %in% removed30k_R, "remove", NA))  

write_tsv(sampleinfo_f, file.path(root_dir, "processed_tables","sampleIDs_info_after_rarefaction.tsv"))

write_tsv(samplemetadata, file.path(root_dir, "processed_tables", "samplemetadata.tsv"))

## write out phyloseq objects with low count samples retained
dir.create(file.path(root_dir, "results/RDS_files"))

saveRDS(its_phyobj_all, file.path(root_dir,"results/RDS_files", "ITS_all_phyobj.rds"))
saveRDS(D_phyobj_30k_all, file.path(root_dir,"results/RDS_files","OTU_DNA_all_phyobj.rds"))
saveRDS(R_phyobj_30k_all, file.path(root_dir,"results/RDS_files", "OTU_RNA_all_phyobj.rds"))

saveRDS(its_phyobj_10k, file.path(root_dir,"results/RDS_files", "ITS_rarefied10k_phyobj.rds"))
saveRDS(D_phyobj_30k, file.path(root_dir,"results/RDS_files","OTU_DNA_rarefied30k_phyobj.rds"))
saveRDS(R_phyobj_30k, file.path(root_dir,"results/RDS_files", "OTU_RNA_rarefied30k_phyobj.rds"))
