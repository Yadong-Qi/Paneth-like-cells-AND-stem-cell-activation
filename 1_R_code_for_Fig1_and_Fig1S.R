# curatedMetagenomicData Database

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

#1-----install packages
BiocManager::install("ExperimentHubData")
BiocManager::install("curatedMetagenomicData")

library(curatedMetagenomicData)


#2------view datasets
View(combined_metadata)



table(combined_metadata$antibiotics_current_use)

table(combined_metadata$disease)
table(combined_metadata$body_site)
table(combined_metadata$antibiotics_current_use)
table(combined_metadata$age)



#3-filter microbiome data from database
curatedMetagenomicData("*metaphlan_bugs_list.stool*", dryrun = TRUE)

library(dplyr)
Health_project = filter(combined_metadata,
                        body_site %in% c("stool")& 
                        age != "NA" & 
                        antibiotics_current_use %in% c("no")& 
                        disease %in% c("healthy") 
                               )

Health_project_name = table(Health_project$dataset_name)

allprojectname = as.data.frame(curatedMetagenomicData("*metaphlan_bugs_list.stool*", dryrun = TRUE),row.names = )


Health_project_downloadname =paste0( names(Health_project_name),".metaphlan_bugs_list.stool")


loman <- curatedMetagenomicData(Health_project_downloadname, dryrun = FALSE)




eset <- mergeData(loman)



#4--cconvert data subjects to plot 
library(phyloseq)
loman.pseq = ExpressionSet2phyloseq( eset, phylogenetictree = TRUE)
loman.pseq


plot_bar(loman.pseq, "age_category",fill = "Species")
dev.off()


#5----define the age group
Young = subset_samples(loman.pseq, age %in% c(18:40))
Adult = subset_samples(loman.pseq, age %in% c(41:60))
Old = subset_samples(loman.pseq, age %in% c(61:120))

par(mar = c(15, 2, 2, 2) + 5) # make more room on bottom margin

sum(Young@otu_table@.Data["s__Bifidobacterium_adolescentis",])
sum(Young@otu_table)

barplot(sort(taxa_sums(Young), TRUE)[1:50]/nsamples(Young), las=2)
barplot(sort(taxa_sums(Adult), TRUE)[1:50]/nsamples(Adult), las=2)
barplot(sort(taxa_sums(Old), TRUE)[1:50]/nsamples(Old), las=2)

dev.off()


write.csv(sort(taxa_sums(Young), TRUE)[1:500]/nsamples(Young),file = "Young18-40_normoliz_TOP500.csv",quote = F,col.names = T,row.names = T)
write.csv(sort(taxa_sums(Adult), TRUE)[1:500]/nsamples(Adult),file = "Adult41-60_normoliz_TOP500.csv",quote = F,col.names = T,row.names = T)
write.csv(sort(taxa_sums(Old), TRUE)[1:500]/nsamples(Old),file = "Old61-120_normoliz_TOP500.csv",quote = F,col.names = T,row.names = T)





#6-get BA relative abundence
Speices = subset_taxa(loman.pseq,  Species =="Bifidobacterium_adolescentis")
pdf("demo23.pdf",width = 10,height = 10)
plot_bar(Speices, "age_category", fill="Species")
dev.off()

S = as.data.frame((Speices@otu_table))
S1 = as.data.frame(t(Speices@sam_data))
BA = as.data.frame(t(rbind(S,S1)))
write.csv(BA,file = "B.a abundenceOUT.csv",quote = F,col.names = T,row.names = T)


SpeicesJohn = subset_taxa(loman.pseq,  Species =="Lactobacillus_johnsonii")
pdf("demoJN.pdf",width = 10,height = 10)
plot_bar(SpeicesJohn, "age_category", fill="Species")
dev.off()

S = as.data.frame((SpeicesJohn@otu_table))
S1 = as.data.frame(t(SpeicesJohn@sam_data))
JN = as.data.frame(t(rbind(S,S1)))
write.csv(JN,file = "JN abundenceOUT.csv",quote = F,col.names = T,row.names = T)

#6---save and quit
save.image("JN_Aging.RData")


# sessionInfo()
# R version 3.6.3 (2020-02-29)
# Platform: x86_64-apple-darwin15.6.0 (64-bit)
# Running under: macOS Catalina 10.15.6
# 
# Matrix products: default
# BLAS:   /System/Library/Frameworks/Accelerate.framework/Versions/A/Frameworks/vecLib.framework/Versions/A/libBLAS.dylib
# LAPACK: /Library/Frameworks/R.framework/Versions/3.6/Resources/lib/libRlapack.dylib
# 
# locale:
#   [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
# 
# attached base packages:
#   [1] stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
#   [1] phyloseq_1.30.0
# 
# loaded via a namespace (and not attached):
#   [1] progress_1.2.2      tinytex_0.40        tidyselect_1.1.0   
# [4] xfun_0.29           purrr_0.3.4         reshape2_1.4.4     
# [7] splines_3.6.3       lattice_0.20-38     rhdf5_2.30.1       
# [10] colorspace_2.0-3    vctrs_0.4.1         generics_0.1.2     
# [13] stats4_3.6.3        mgcv_1.8-31         survival_3.1-8     
# [16] utf8_1.2.2          rlang_1.0.2         pillar_1.7.0       
# [19] glue_1.6.2          BiocGenerics_0.32.0 foreach_1.5.1      
# [22] lifecycle_1.0.1     plyr_1.8.6          stringr_1.4.0      
# [25] zlibbioc_1.32.0     Biostrings_2.54.0   munsell_0.5.0      
# [28] gtable_0.3.0        codetools_0.2-16    Biobase_2.46.0     
# [31] permute_0.9-5       IRanges_2.20.2      biomformat_1.14.0  
# [34] parallel_3.6.3      fansi_1.0.2         Rcpp_1.0.8.3       
# [37] scales_1.1.1        BiocManager_1.30.10 vegan_2.5-6        
# [40] S4Vectors_0.24.4    jsonlite_1.8.0      XVector_0.26.0     
# [43] ggplot2_3.3.5       hms_0.5.3           stringi_1.5.3      
# [46] dplyr_1.0.2         grid_3.6.3          ade4_1.7-16        
# [49] cli_3.2.0           tools_3.6.3         magrittr_2.0.3     
# [52] tibble_3.1.7        cluster_2.1.0       crayon_1.5.1       
# [55] ape_5.4-1           pkgconfig_2.0.3     MASS_7.3-51.5      
# [58] ellipsis_0.3.2      Matrix_1.2-18       data.table_1.14.2  
# [61] prettyunits_1.1.1   iterators_1.0.13    Rhdf5lib_1.8.0     
# [64] R6_2.5.1            multtest_2.42.0     igraph_1.2.5       
# [67] nlme_3.1-144        compiler_3.6.3     

