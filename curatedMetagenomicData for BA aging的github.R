curatedMetagenomicData 

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("ExperimentHubData")
BiocManager::install("curatedMetagenomicData")

library(curatedMetagenomicData)

View(combined_metadata)



table(combined_metadata$antibiotics_current_use)

table(combined_metadata$disease)
table(combined_metadata$body_site)
table(combined_metadata$antibiotics_current_use)
table(combined_metadata$age)


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


library(phyloseq)
loman.pseq = ExpressionSet2phyloseq( eset, phylogenetictree = TRUE)
loman.pseq


plot_bar(loman.pseq, "age_category",fill = "Species")
dev.off()



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


save.image("JN_Aging.RData")
