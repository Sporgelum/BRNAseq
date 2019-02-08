#Control analysis for the samples analyzing 
#BETA prior T for the CONTROL analysis

#Control approach where we are considering the following comparisson.


#RFP vs BFP


#betaprior is left to T.

#Clear and clean all references and saved data.
rm(list=ls(all=T))

#Establish where to work.
getwd()
setwd("/home/marius/Desktop/Ordered/")

#Set the pathway to the files and read recursively
all_files<-list.files(path = "/haus/marius/marcos_data/all_lanes/results", recursive = T, pattern = "*.h5", full.names = T)
head(all_files)



#Get the data and add names to localize it.
#Order the values numerically
all_files<-all_files[order(nchar(all_files),all_files)]
head(all_files)
name_files<-list.files("/haus/marius/marcos_data/all_lanes/all_lanes/", pattern = "*.fastq")
name_files
#Subsittute all the string for a shorter name.
name_files<-gsub(pattern = "*HVJ22BGX5_sox10_18s002081-1-1_",replacement = "", name_files)
head(name_files)
name_files<-gsub(pattern = ".fastq",replacement = "_seq",name_files)
head(name_files)
#Change the names of the files, select each lane with a name representative.
names(all_files)<-name_files
head(all_files)


#BIOMART CREATE the TX2GENE filename for TXIMPORT, importing kallisto dato to DESeq2.

#Fisrt create a table of genes/transcript or what is in the column interested to switch to.
filterValues=as.data.frame(read.table("/haus/marius/marcos_data/all_lanes/results/output1/abundance.tsv",header = T, sep = "\t"))
head(filterValues)

#List all the available datasets
biomaRt::listMarts()

#Biomart oneliner.
#Conect to ensembl and use the specific gene ensembl wanted.
ensembl=useEnsembl(biomart = "ensembl", dataset = "drerio_gene_ensembl")
#List all the attributes to select for the dat aif neecssay, right now just with necessary for tx2gnee will be enough.
listAttributes(ensembl)


#get the combination of attributes names for ensembl filtered out by your filtertype.
ensembl_annotation_transcripts=getBM(attributes = c('ensembl_gene_id',"ensembl_transcript_id_version"),filters = "ensembl_transcript_id_version",values =filterValues$target_id , mart = ensembl)#,'external_gene_name','hgnc_symbol','description'), 
head(ensembl_annotation_transcripts);dim(ensembl_annotation_transcripts)

#Create the geneIDs and TXNAME lists.
tx2gene<-merge(filterValues,ensembl_annotation_transcripts, by.x="target_id", by.y="ensembl_transcript_id_version", all.x=T,all.y=F)

#Take only the necessary data for tx2gene file.
head(tx2gene)
tx2gene<-data.frame(tx2gene$target_id, tx2gene$ensembl_gene_id)
colnames(tx2gene)<-c("TXNAME","GENEID")
head(tx2gene)
dim(tx2gene)

#Select the data interested in from all files in interested in separatring comparrisons.
name_files
#Create the data groups separately.
control_names<-c("Sande_lane1105_seq","Sande_lane1113_seq","Sande_lane1121_seq","Sande_lane1209_seq","Sande_lane1217_seq","Sande_lane1225_seq",
                 "Sande_lane1233_seq","Sande_lane1241_seq","Sande_lane1249_seq","Sande_lane1257_seq","Sande_lane1297_seq","Sande_lane1313_seq",
                 "Sande_lane1317_seq","Sande_lane1325_seq","Sande_lane1333_seq","Sande_lane1341_seq","Sande_lane1349_seq","Sande_lane1357_seq",
                 "Sande_lane1365_seq","Sande_lane181_seq","Sande_lane189_seq","Sande_lane197_seq" )

control_names<-sort(control_names,decreasing = F)
control_names
control_names_all_files<-all_files[control_names]

txi_control<-tximport(files=control_names_all_files,type="kallisto",tx2gene=tx2gene)

names(txi_control)
head(txi_control$counts)
head(txi_control$infReps$Sande_lane1105_seq)



#Create a folder for all the data related with this analysis.
getwd()
dir.create("./control")
#Set the analysis results into the control folder.
setwd("/home/marius/Desktop/Ordered/control/")


#Go for DESeq2 analysis.
#As DESeq2 uses raw counts the analysis can go straight froward since here!

head(txi_control$counts)



#Create the sampletable for later comparissons. Creating the cases of the condition.
sampleTable_control<-data.frame(CELLTYPE=factor(c(rep("RFP",1),rep("BFP",1),rep("RFP",1),rep("BFP",1),rep("RFP",2),
                                                  rep("BFP",1),rep("RFP",2),rep("BFP",2),rep("RFP",3),rep("BFP",1),rep("RFP",2),
                                                  rep("BFP",1),rep("RFP",2),rep("BFP",1),rep("RFP",1))),row.names = colnames(txi_control$counts))


#See and save the data
sampleTable_control



write.table(sampleTable_control,file = "./sampleTable_control_Cryoinjured_vs_Control.txt",sep = "\t",row.names = T, col.names = T)
#Create the DESeq2 DATA SET
dds_control<-DESeqDataSetFromTximport(txi_control,sampleTable_control, ~CELLTYPE)


dim(dds_control)
head(dds_control)
dds_control@colData
#Raw counts and average transcript length.
head(dds_control@assays$data$counts)


#Remove all the samples under 1 count of the ROWSUM.
dim(dds_control)
dds_control <- dds_control[ rowSums(counts(dds_control)) > 3, ]
dim(dds_control)
#Check the differences.



#DIfferential EXpression ANalysis
dds_control=DESeq(dds_control,betaPrior = T,parallel = TRUE)
dds_control@colData
head(dds_control$CELLTYPE)
#using pre-existing size factors
#estimating dispersions
#gene-wise dispersion estimates
#mean-dispersion relationship
#final dispersion estimates
#fitting model and testing
#See the dispersion of the data.
pdf("plot_dispersion.pdf", width = 6, height = 10)
plotDispEsts(dds_control)
dev.off()

#For columns colors.
#help("plotDispEsts")

#············································································
head(dds_control@assays$data$counts)
head(dds_control@assays$data$normalizationFactors)

#············································································
head(sampleTable_control)
#USING BH also known as FDR
#tidy=T the first column will be the rownames
#Results of the RNA-seq
#Specify which contrast to do, comparisson against.
dds_res_control=results(dds_control,contrast=c("CELLTYPE","RFP","BFP"), parallel=TRUE,pAdjustMethod = "BH", tidy=F)
head(dds_res_control)
dim(dds_res_control)

#Check what has been done.
#Only works if tidy is established to F.
mcols(dds_res_control)$description
#Annotate results.
#After Selecting the STATISTICAL THRESHOLDS we can start annotating.


#Sort the data frame annotated by the padj value and the -log2foldchange
dds_res_control<-dds_res_control[order(dds_res_control[,6], -dds_res_control[,2]),]
head(dds_res_control)


#Start annotating the data with the Human and Mouse homologous for later dat obtention.

#get the comboination of attributnames for ensembl filtered out by your filtertype.
listAttributes(ensembl)
ensembl_dds_res_control=getBM(attributes = c("ensembl_gene_id","description","hsapiens_homolog_associated_gene_name","hsapiens_homolog_ensembl_gene","mmusculus_homolog_associated_gene_name","mmusculus_homolog_ensembl_gene"),filters = "ensembl_gene_id",values =rownames(dds_res_control) , mart = ensembl)#,'external_gene_name','hgnc_symbol','description'), 
head(ensembl_dds_res_control);dim(ensembl_dds_res_control)


dim(subset(dds_res_control,padj < .05))

#Merge the annotated data with the different species ENTREZIDs.
#Select the data from the following packages --> NCBI  after ensembl.



#####################################################
#                                                   #
#           ZEBRAFISH                               #
#                                                   #
#####################################################

#ZEBRAFISH.
#Zebrafish annotating with ZEBRAFISH SPECIFIC ENTREZID
columns(org.Dr.eg.db)
anno_zebrafish_control<-AnnotationDbi::select(org.Dr.eg.db,
                                                 keys = ensembl_dds_res_control$ensembl_gene_id,
                                                 columns = c("SYMBOL","GENENAME","ENTREZID"),
                                                 keytype = "ENSEMBL") 

dim(anno_zebrafish_control)
head(anno_zebrafish_control)

#Merge the statistical information with the annotations realized
anno_zebrafish_control_stats<-merge(as.data.frame(dds_res_control),anno_zebrafish_control,by.x="row.names",by.y="ENSEMBL")
#Rename the ensemble gene id and wit the species.
colnames(anno_zebrafish_control_stats)[1]<-"drerio_ensembl_gene_id"
dim(anno_zebrafish_control_stats)
head(anno_zebrafish_control_stats)




#####################################################
#                                                   #
#           HUMAN                                   #
#                                                   #
#####################################################

#HUMAN
#Homo_sapiens
#get the Human orthologues ENTREZIDs
#USING Human specific ENTREZIDs
#NCBI package.
columns(org.Hs.eg.db)
anno_human_control<-AnnotationDbi::select(org.Hs.eg.db,
                                             keys=ensembl_dds_res_control$hsapiens_homolog_ensembl_gene,
                                             columns = c("SYMBOL","GENENAME","ENTREZID"),
                                             keytype = "ENSEMBL")


head(anno_human_control)
head(ensembl_dds_res_control)

#Merge the human data with the homologous data and select later the columns interesting.
anno_human_pre_stats_control<-merge(ensembl_dds_res_control[,c("ensembl_gene_id","hsapiens_homolog_ensembl_gene")],anno_human_control,by.x="hsapiens_homolog_ensembl_gene",by.y="ENSEMBL")
head(anno_human_pre_stats_control)
dim(anno_human_pre_stats_control)


#Merge the statistical information with the annotations realized
anno_human_control_stats<-merge(as.data.frame(dds_res_control),anno_human_pre_stats_control,by.x="row.names",by.y="ensembl_gene_id")
dim(anno_human_control_stats)
head(anno_human_control_stats)

#Create the same data frame as in zebrafish for easier working.
anno_human_control_stats<-anno_human_control_stats[,c("hsapiens_homolog_ensembl_gene","baseMean","log2FoldChange","lfcSE","stat","pvalue","padj","SYMBOL","GENENAME","ENTREZID")]
head(anno_human_control_stats)




#####################################################
#                                                   #
#           MOUSE                                   #
#                                                   #
#####################################################
#MOUSE
#Mus_musculus
#MOUSE entrez data 
#MOUSE specific ENTREZID data...
columns(org.Mm.eg.db)
anno_mouse_control<-AnnotationDbi::select(org.Mm.eg.db,
                                             keys=ensembl_dds_res_control$mmusculus_homolog_ensembl_gene,
                                             columns = c("SYMBOL","GENENAME","ENTREZID"),
                                             keytype = "ENSEMBL")

head(anno_mouse_control)
dim(anno_mouse_control)
head(ensembl_dds_res_control)
#Merge the human data with the homologous data and select later the columns interesting.
anno_mouse_pre_stats_chery<-merge(ensembl_dds_res_control[,c("ensembl_gene_id","mmusculus_homolog_ensembl_gene")],anno_mouse_control,by.x="mmusculus_homolog_ensembl_gene",by.y="ENSEMBL")
head(anno_mouse_pre_stats_chery)
dim(anno_mouse_pre_stats_chery)



#Merge the statistical information with the annotations realized
anno_mouse_control_stats<-merge(as.data.frame(dds_res_control),anno_mouse_pre_stats_chery,by.x="row.names",by.y="ensembl_gene_id")
dim(anno_mouse_control_stats)
head(anno_mouse_control_stats)

#Create the same data frame as in zebrafish for easier working.
anno_mouse_control_stats<-anno_mouse_control_stats[,c("mmusculus_homolog_ensembl_gene","baseMean","log2FoldChange","lfcSE","stat","pvalue","padj","SYMBOL","GENENAME","ENTREZID")]
head(anno_mouse_control_stats)





#####################################################
#             START WITH FDR ENRICHMENT             #
#####################################################
#         Using for this approach:                  #
#                                                   #
#           -1 → List DEgenes padj<0.05            #
#                                                   #
#           -2 → ENTREZID                          # 
#                                                   #
#       PERFORM ENRICHMET IN DIFFERENT DB.          #
#                                                   #
#####################################################








#Startat the enrichment.
getwd()


#Create the directories for the files using significant filtered genes and enrichment of that list, using HUMAN converted ID, and also mouse,zebrafish
#Start the analysis of the data with the significative genes.
dir.create(recursive=T,"./Enrichment/FDR/")

#Load the NCBI data packages for the annotations.
#ENRICHMENT GO with the adjusted values.

#Using enrichGO from CLUSTERPROFILER.
#Create a folder for the Go annotations 
dir.create("./Enrichment/FDR/GO")




#####################################################
#                                                   #
#           ZEBRAFISH                               #
#                                                   #
#####################################################
#                                                   #
#           GO ENRICHMENT with  FDR                 #
#                                                   #
#   Using GO→ BP                                    #  
#   Using GO→ MF                                    #
#   Using GO→ CC                                    #
######################################################

#Start with Zebrafish annotations with zebrafish
dir.create("./Enrichment/FDR/GO/DanioRerio/")
DrOrgDb<-org.Dr.eg.db


#Get the differential expressed genes, and use only the unique seuqences no repeated..
dr_control_diff_expressed_genes<-unique(anno_zebrafish_control_stats$ENTREZID[which(anno_zebrafish_control_stats$padj<.05)])


length(dr_control_diff_expressed_genes)
head(dr_control_diff_expressed_genes)


#Biologicall process enriched for zebrafish 
dr_control_BP_genelist_enrichGO <- enrichGO(gene = dr_control_diff_expressed_genes ,OrgDb  =DrOrgDb, ont = "BP", pAdjustMethod = "BH",pvalueCutoff  = 0.05,readable = T)


#Molecular Function enriched for zebrafish
dr_control_MF_genelist_enrichGO <- enrichGO(gene = dr_control_diff_expressed_genes ,OrgDb  = DrOrgDb, ont = "MF", pAdjustMethod = "BH",pvalueCutoff  = 0.05,readable = T)

#Cellular component enriched zebrafish
dr_control_CC_genelist_enrichGO <- enrichGO(gene = dr_control_diff_expressed_genes ,OrgDb  = DrOrgDb, ont = "CC", pAdjustMethod = "BH",pvalueCutoff  = 0.05,readable = T)


#Save the lists for Zebrafish enrichment.
#Save better all the results,head(dr_control_BP_genelist_enrichGO@result)
#Save BP list as dataframe.
write.table(x = data.frame(dr_control_BP_genelist_enrichGO@result),"./Enrichment/FDR/GO/DanioRerio/dr_control_BP_enrichGO.txt",sep="\t", quote = F, row.names = F, col.names = T)


#Save the MF list as dataframe
write.table(x = data.frame(dr_control_MF_genelist_enrichGO@result),"./Enrichment/FDR/GO/DanioRerio/dr_control_MF_enrichGO.txt",sep="\t", quote = F, row.names = F, col.names = T)

#Save the CC list as a dataframe
write.table(x = data.frame(dr_control_CC_genelist_enrichGO@result),"./Enrichment/FDR/GO/DanioRerio/dr_control_CC_enrichGO.txt",sep="\t", quote = F, row.names = F, col.names = T)







#Visualize the data obtained for Danio rerio Enrichment.

#Plot for the BProcess data obtained from the enrichment analysis in GO fpr Danio rerio control Cryoinjury vs Control.

#Biological processes plots.
#BARPLOT
pdf("./Enrichment/FDR/GO/DanioRerio/dr_control_BP_barplot.pdf", height = 20, width = 16)
barplot(dr_control_BP_genelist_enrichGO,drop=T,showCategory = 25,title = "dr BP enrichment Cherry \n genelist Selected by padj<0.05 ")
dev.off()

#DOTPLOT of BP in Dr
pdf("./Enrichment/FDR/GO/DanioRerio/dr_control_BP_dotplot.pdf",height = 20,width = 16)
dotplot(dr_control_BP_genelist_enrichGO, showCategory=25)
dev.off()


#Molecular Functions plots.
#BARPLOT
pdf("./Enrichment/FDR/GO/DanioRerio/dr_control_MF_barplot.pdf", height = 20, width = 16)
barplot(dr_control_MF_genelist_enrichGO,drop=T,showCategory = 25,title = "dr MF enrichment Cherry \n genelist Selected by padj<0.05 ")
dev.off()

#DOTPLOT of MF in Danio rerio
pdf("./Enrichment/FDR/GO/DanioRerio/dr_control_MF_dotplot.pdf",height = 20,width = 16)
dotplot(dr_control_MF_genelist_enrichGO, showCategory=25)
dev.off()


#Celullar Components plots.
#BARPLOT
pdf("./Enrichment/FDR/GO/DanioRerio/dr_control_CC_barplot.pdf", height = 20, width = 16)
barplot(dr_control_CC_genelist_enrichGO,drop=T,showCategory = 25,title = "dr CC enrichment Cherry \n genelist Selected by padj<0.05 ")
dev.off()

#DOTPLOT of CC in Danio rerio
pdf("./Enrichment/FDR/GO/DanioRerio/dr_control_CC_dotplot.pdf",height = 20,width = 16)
dotplot(dr_control_CC_genelist_enrichGO, showCategory=25)
dev.off()








#####################################################
#                                                   #
#           HUMAN                                   #
#                                                   #
#####################################################
#                                                   #
#           GO ENRICHMENT with  FDR                 #
#                                                   #
#   Using GO→ BP                                    #  
#   Using GO→ MF                                    #
#   Using GO→ CC                                    #
######################################################



#Second continue with human
dir.create("./Enrichment/FDR/GO/Homosapiens/")
HsOrgDb<-org.Hs.eg.db


#Get the differential expressed genes
hs_control_diff_expressed_genes<-unique(anno_human_control_stats$ENTREZID[which(anno_human_control_stats$padj<.05)])
length(hs_control_diff_expressed_genes)


#Biologicall process enriched for human 
hs_control_BP_genelist_enrichGO <- enrichGO(gene = hs_control_diff_expressed_genes ,OrgDb  =HsOrgDb, ont = "BP", pAdjustMethod = "BH",pvalueCutoff  = 0.05,readable = T)


#Molecular Function enriched for human
hs_control_MF_genelist_enrichGO <- enrichGO(gene = hs_control_diff_expressed_genes ,OrgDb  = HsOrgDb, ont = "MF", pAdjustMethod = "BH",pvalueCutoff  = 0.05,readable = T)

#Cellular component enriched human
hs_control_CC_genelist_enrichGO <- enrichGO(gene = hs_control_diff_expressed_genes ,OrgDb  = HsOrgDb, ont = "CC", pAdjustMethod = "BH",pvalueCutoff  = 0.05,readable = T)



#Save the lists for Human enrichment.
#Save better all the results,head(hs_control_BP_genelist_enrichGO@result)
#Save BP list as dataframe.
write.table(x = data.frame(hs_control_BP_genelist_enrichGO@result),"./Enrichment/FDR/GO/Homosapiens/hs_control_BP_enrichGO.txt",sep="\t", quote = F, row.names = F, col.names = T)


#Save the MF list as dataframe
write.table(x = data.frame(hs_control_MF_genelist_enrichGO@result),"./Enrichment/FDR/GO/Homosapiens/hs_control_MF_enrichGO.txt",sep="\t", quote = F, row.names = F, col.names = T)

#Save the CC list as a dataframe
write.table(x = data.frame(hs_control_CC_genelist_enrichGO@result),"./Enrichment/FDR/GO/Homosapiens/hs_control_CC_enrichGO.txt",sep="\t", quote = F, row.names = F, col.names = T)





#Visualize the data obtained for the human.
#Plot for the BProcess data obtained from the enrichment analysis in GO for Human in homologues for control Cryoinjury vs Control.
#Biological processes plots.
#BARPLOT
pdf("./Enrichment/FDR/GO/Homosapiens/hs_control_BP_barplot.pdf", height = 20, width = 16)
barplot(hs_control_BP_genelist_enrichGO,drop=T,showCategory = 25,title = "hs BP enrichment Cherry \n genelist Selected by padj<0.05 ")
dev.off()

#DOTPLOT of BP in Hs
pdf("./Enrichment/FDR/GO/Homosapiens/hs_control_BP_dotplot.pdf",height = 20,width = 16)
dotplot(hs_control_BP_genelist_enrichGO, showCategory=25)
dev.off()


#Molecular Functions plots.
#BARPLOT
pdf("./Enrichment/FDR/GO/Homosapiens/hs_control_MF_barplot.pdf", height = 20, width = 16)
barplot(hs_control_MF_genelist_enrichGO,drop=T,showCategory = 25,title = "hs MF enrichment Cherry \n genelist Selected by padj<0.05 ")
dev.off()

#DOTPLOT of MF in Human
pdf("./Enrichment/FDR/GO/Homosapiens/hs_control_MF_dotplot.pdf",height = 20,width = 16)
dotplot(hs_control_MF_genelist_enrichGO, showCategory=25)
dev.off()


#Celullar Components plots.
#BARPLOT
pdf("./Enrichment/FDR/GO/Homosapiens/hs_control_CC_barplot.pdf", height = 20, width = 16)
barplot(hs_control_CC_genelist_enrichGO,drop=T,showCategory = 25,title = "hs CC enrichment Cherry \n genelist Selected by padj<0.05 ")
dev.off()

#DOTPLOT of CC in Human
pdf("./Enrichment/FDR/GO/Homosapiens/hs_control_CC_dotplot.pdf",height = 20,width = 16)
dotplot(hs_control_CC_genelist_enrichGO, showCategory=25)
dev.off()








#####################################################
#                                                   #
#           MOUSE                                   #
#                                                   #
#####################################################
#                                                   #
#           GO ENRICHMENT with  FDR                 #
#                                                   #
#   Using GO→ BP                                    #  
#   Using GO→ MF                                    #
#   Using GO→ CC                                    #
######################################################


#Do the same analysis for Mouse

MmOrgDb<- org.Mm.eg.db
#Using Human Annotation Package
dir.create("./Enrichment/FDR/GO/Mmusculus/")



#Get the differential expressed genes
mm_control_diff_expressed_genes<-unique(anno_mouse_control_stats$ENTREZID[which(anno_mouse_control_stats$padj<.05)])
length(mm_control_diff_expressed_genes)

#Biologicall process enriched for Mouse 
mm_control_BP_genelist_enrichGO <- enrichGO(gene = mm_control_diff_expressed_genes ,OrgDb  =MmOrgDb, ont = "BP", pAdjustMethod = "BH",pvalueCutoff  = 0.05,readable = T)


#Molecular Function enriched for mouse
mm_control_MF_genelist_enrichGO <- enrichGO(gene = mm_control_diff_expressed_genes ,OrgDb  = MmOrgDb, ont = "MF", pAdjustMethod = "BH",pvalueCutoff  = 0.05,readable = T)

#Cellular component enriched mouse
mm_control_CC_genelist_enrichGO <- enrichGO(gene = mm_control_diff_expressed_genes ,OrgDb  = MmOrgDb, ont = "CC", pAdjustMethod = "BH",pvalueCutoff  = 0.05,readable = T)



#Save the lists for mouse homoogues enrichments.
#Save better all the results,head(mm_control_BP_genelist_enrichGO@result)
#Save BP list as dataframe.
write.table(x = data.frame(mm_control_BP_genelist_enrichGO@result),"./Enrichment/FDR/GO/Mmusculus/mm_control_BP_enrichGO.txt",sep="\t", quote = F, row.names = F, col.names = T)


#Save the MF list as dataframe
write.table(x = data.frame(mm_control_MF_genelist_enrichGO@result),"./Enrichment/FDR/GO/Mmusculus/mm_control_MF_enrichGO.txt",sep="\t", quote = F, row.names = F, col.names = T)

#Save the CC list as a dataframe
write.table(x = data.frame(mm_control_CC_genelist_enrichGO@result),"./Enrichment/FDR/GO/Mmusculus/mm_control_CC_enrichGO.txt",sep="\t", quote = F, row.names = F, col.names = T)


#Visualize the data obtained.

#Plot for the BProcess data obtained from the enrichment analysis in GO fpr Mouse homologues control Cryoinjury vs Control.

#Biological processes plots.
#BARPLOT
pdf("./Enrichment/FDR/GO/Mmusculus/mm_control_BP_barplot.pdf", height = 20, width = 16)
barplot(mm_control_BP_genelist_enrichGO,drop=T,showCategory = 12,title = "mm BP enrichment Cherry \n genelist Selected by padj<0.05 ")
dev.off()

#DOTPLOT of BP in Danio rerio
pdf("./Enrichment/FDR/GO/Mmusculus/mm_control_BP_dotplot.pdf",height = 20,width = 16)
dotplot(mm_control_BP_genelist_enrichGO, showCategory=12)
dev.off()


#Molecular Functions plots.
#BARPLOT
pdf("./Enrichment/FDR/GO/Mmusculus/mm_control_MF_barplot.pdf", height = 20, width = 16)
barplot(mm_control_MF_genelist_enrichGO,drop=T,showCategory = 12,title = "mm MF enrichment Cherry \n genelist Selected by padj<0.05 ")
dev.off()

#DOTPLOT of MF in mus musculus
pdf("./Enrichment/FDR/GO/Mmusculus/mm_control_MF_dotplot.pdf",height = 20,width = 16)
dotplot(mm_control_MF_genelist_enrichGO, showCategory=12)
dev.off()


#Celullar Components plots.
#BARPLOT
pdf("./Enrichment/FDR/GO/Mmusculus/mm_control_CC_barplot.pdf", height = 20, width = 16)
barplot(mm_control_CC_genelist_enrichGO,drop=T,showCategory = 12,title = "mm CC enrichment Cherry \n genelist Selected by padj<0.05 ")
dev.off()

#DOTPLOT of CC in Danio rerio
pdf("./Enrichment/FDR/GO/Mmusculus/mm_control_CC_dotplot.pdf",height = 20,width = 16)
dotplot(mm_control_CC_genelist_enrichGO, showCategory=12)
dev.off()







#####################################################
#                                                   #
#           ZEBRAFISH                               #
#                                                   #
#####################################################
#                                                   #
#           KEGG ENRICHMENT with  FDR               #
#                                                   #
#           PATHWAYS ANALYSIS                       #
#####################################################


#Significant data for pathways enrichment, padj<0.05


#Move to KEGG enrichment.
#Create directory for FDR --> KEGG
dir.create("./Enrichment/FDR/KEGG/")
#Obtain the same results pfor significant pathways.
#Following the same process
#Check for supported organisms.

search_kegg_organism(str = "Danio",ignore.case = T,by = "scientific_name")

#Create directory for Danio rerior enrichment
dir.create("./Enrichment/FDR/KEGG/Drerio/")

#Start with the zebrafish process again.
dr_control_genelist_enrichKEGG<-enrichKEGG(gene = dr_control_diff_expressed_genes,organism = "dre",keyType = "kegg",pvalueCutoff = 0.05,pAdjustMethod = "BH",minGSSize = 10)
dr_control_genelist_enrichKEGG@result[,8]
head(dr_control_genelist_enrichKEGG@result)

#Write the results of the table enriched with the KEGG significant data.
write.table(data.frame(dr_control_genelist_enrichKEGG@result),"./Enrichment/FDR/KEGG/Drerio/dr_control_enrichmentKEGG.txt",sep = "\t",quote = F, row.names = F,col.names = T)


head(dr_control_genelist_enrichKEGG@result)
#Plot the decided figures as required.
#Those figures got downloaded in the working directory
dre04371<-pathview(gene.data = dr_control_diff_expressed_genes,
                   pathway.id = "dre04371",
                   species="dre",
                   limit = list(gene=5,cpd=1))






#####################################################
#                                                   #
#           HUMAN                                   #
#                                                   #
#####################################################
#                                                   #
#           KEGG ENRICHMENT with  FDR               #
#                                                   #
#           PATHWAYS ANALYSIS                       #
#####################################################


###HUMAN KEGG

#Human Enrichment using KEGG
#Start by creating the folder into the FDR/KEGG directory
dir.create("./Enrichment/FDR/KEGG/Hsapiens/")

#check for the species name supported
search_kegg_organism("Sapiens", by="scientific_name",ignore.case = T)


#KEGG enrichment using significative info.
hs_control_genelist_enrichKEGG<-enrichKEGG(gene = hs_control_diff_expressed_genes,organism = "hsa",keyType = "kegg",pvalueCutoff = 0.05,pAdjustMethod = "BH",minGSSize = 10)
head(hs_control_genelist_enrichKEGG@result)


#Write the results of the table enriched with the KEGG significant data.
write.table(data.frame(hs_control_genelist_enrichKEGG@result),"./Enrichment/FDR/KEGG/Hsapiens/hs_control_enrichmentKEGG.txt",sep = "\t",quote = F, row.names = F,col.names = T)
head(hs_control_genelist_enrichKEGG@result)


#Plot the decided figures as required.
#Those figures got downloaded in the working directory
hsa00601<-pathview(gene.data = hs_control_diff_expressed_genes,
                   pathway.id = "hsa00601",
                   species="hsa",
                   limit = list(gene=5,cpd=1))








#####################################################
#                                                   #
#           MOUSE                                   #
#                                                   #
#####################################################
#                                                   #
#           KEGG ENRICHMENT with  FDR               #
#                                                   #
#           PATHWAYS ANALYSIS                       #
#####################################################


###MOUSE KEGG

#Start by creating the folder into the FDR/KEGG directory
dir.create("./Enrichment/FDR/KEGG/Mmusculus/")

#check for the species name supported
search_kegg_organism("musculus", by="scientific_name",ignore.case = T)

length(mm_control_diff_expressed_genes)
#KEGG enrichment using significative info.
mm_control_genelist_enrichKEGG<-enrichKEGG(gene = mm_control_diff_expressed_genes,organism = "mmu",keyType = "kegg",pvalueCutoff = 0.05,pAdjustMethod = "BH",minGSSize = 10)
head(mm_control_genelist_enrichKEGG@result)


#Write the results of the table enriched with the KEGG significant data.
write.table(data.frame(mm_control_genelist_enrichKEGG@result),"./Enrichment/FDR/KEGG/Mmusculus/mm_control_enrichmentKEGG.txt",sep = "\t",quote = F, row.names = F,col.names = T)


head(mm_control_genelist_enrichKEGG@result)
#Plot the decided figures as required.
#Those figures got downloaded in the working directory
mmu04020<-pathview(gene.data = mm_control_diff_expressed_genes,
                   pathway.id = "mmu04020",
                   species="mmu",
                   limit = list(gene=5,cpd=1))








#####################################################
#               START WITH GSE                      #
#####################################################
#         Using for this approach:                  #
#                                                   #
#           -1 → Log2FoldChange  (DE Genes)        #
#                                                   #
#           -2 → ENTREZID                          # 
#                                                   #
#       PERFORM ENRICHMET IN DIFFERENT DB.          #
#                                                   #
#####################################################


#Gene set enrichment METHOD.
#Now go to ENRICHMENT with GSE method, using ranked list of genes.
#Start with the enrichment creating the filedirectory

#In clusterProfiler some function need a genelist that has specific features as described below and this “genelist” is use as input to GeneSetenrichment analysis in GO 
#or in Pathways:

#geneList contains three features:
#numeric vector: fold change or other type of numerical variable
#named vector: every number has a name, the corresponding gene ID
#sorted vector: number should be sorted in decreasing order

#First start with the directory.
getwd()
dir.create("./Enrichment/GSE/")





#####################################################
#                                                   #
#           GO ENRICHMENT with  GSE                 #
#                                                   #
#   Using GO→ BP                                    #  
#   Using GO→ MF                                    #
#    Using GO→ CC                                   #
#####################################################



#Start doing the enrichment analysis with gseGO
#Set the respective folder for that.

dir.create("./Enrichment/GSE/GO/")

#####################################################
#                                                   #
#               ZEBRAFISH                           #
#                                                   #
#####################################################
#Start with zebrafish
dir.create("./Enrichment/GSE/GO/Drerio")
#Crete the ranked_genelist
#First clean the data a little bit.

#Remove all duplicated data from the stats table, leaving all the unique ENTREZID
anno_zebrafish_control_stats<-anno_zebrafish_control_stats[-which(duplicated(anno_zebrafish_control_stats$ENTREZID)),]
dim(anno_zebrafish_control_stats)
#first numeric vector using the log2FoldChange
dr_control_ranked_genes<-anno_zebrafish_control_stats$log2FoldChange

#second step name the vector
names(dr_control_ranked_genes)<-as.character(anno_zebrafish_control_stats$ENTREZID)

#Sort the list in decreasing order.
dr_control_ranked_genes<-sort(dr_control_ranked_genes,decreasing = T)



##############################
#        GSE                 #
#   Using GO→ BP            # 
#                            #
##############################
#Start with the GSE in GO databases with BP.
dr_control_GSE_GO_BP<-gseGO(geneList = dr_control_ranked_genes,ont="BP",OrgDb = DrOrgDb,
                               keyType = "ENTREZID",
                               nPerm=1000,
                               minGSSize = 10,
                               pvalueCutoff = 0.05,
                               pAdjustMethod = "BH")


head(dr_control_GSE_GO_BP@result)

##############################
#        GSE                 #
#   Using GO→ MF            # 
#                            #
##############################

#GSE in GO database with MF.
dr_control_GSE_GO_MF<-gseGO(geneList = dr_control_ranked_genes,ont="MF",OrgDb = DrOrgDb,
                               keyType = "ENTREZID",
                               nPerm=1000,
                               minGSSize = 10,
                               pvalueCutoff = 0.05,
                               pAdjustMethod = "BH")


head(dr_control_GSE_GO_MF@result)



##############################
#        GSE                 #
#   Using GO→ CC            # 
#                            #
##############################


#GSE in GO database with CC.
dr_control_GSE_GO_CC<-gseGO(geneList = dr_control_ranked_genes,ont="CC",OrgDb = DrOrgDb,
                               keyType = "ENTREZID",
                               nPerm=1000,
                               minGSSize = 10,
                               pvalueCutoff = 0.05,
                               pAdjustMethod = "BH")


head(dr_control_GSE_GO_CC@result)




#SAVE THE TABLES
#Save the data of Dr into the specific file each table to pvalueCutOff=0.05.
#GSE for GO:BP 
write.table(data.frame(dr_control_GSE_GO_BP@result),"./Enrichment/GSE/GO/Drerio/dr_control_GSE_GO_BP.txt",sep="\t", row.names = F,col.names = T,quote = F)

#GSE for GO:MF
write.table(data.frame(dr_control_GSE_GO_MF@result),"./Enrichment/GSE/GO/Drerio/dr_control_GSE_GO_MF.txt",sep="\t",row.names = F,col.names = T,quote = F)

#GSE for CC
write.table(data.frame(dr_control_GSE_GO_CC@result),"./Enrichment/GSE/GO/Drerio/dr_control_GSE_GO_BP.txt",sep="\t",row.names = F,col.names = T,quote = F)


#Represent the figure.
colors_val<-c("#377EB8","#E41A1C")


#Represent the figures, start with the BP,MF,CC.

#GSE GO BP
pdf("./Enrichment/GSE/GO/Drerio/dr_gseGO_BP_GO.pdf",width = 13,height = 12)

#Separate the lists in Up and Down regulated. 
dr_control_GSE_GO_BP_topUp_pathways<-dr_control_GSE_GO_BP@result%>%
  filter(NES>0)%>%
  top_n(20,wt=-pvalue)

dr_control_GSE_GO_BP_topDown_pathways<-dr_control_GSE_GO_BP@result%>%
  filter(NES<0)%>%
  top_n(20,wt=-pvalue)

dr_control_GSE_GO_BP_topPathWays<-bind_rows(dr_control_GSE_GO_BP_topUp_pathways,rev(dr_control_GSE_GO_BP_topDown_pathways))
head(dr_control_GSE_GO_BP_topPathWays)

#Plot the table converted
ggplot(dr_control_GSE_GO_BP_topPathWays,aes(x=NES,y=Description,
                                               size=pvalue,
                                               col=ifelse(NES > 0,"Upregulated","Downregulated")))+
  #shape=ifelse(p.adjust < 0.25,"False Discovery\nRate < 25%","False Discovery\nRate > 25%")))+
  geom_point()+
  expand_limits(x=c(-5,5))+
  theme_classic()+
  scale_color_manual(name="Normalized Enrichment Score", values =colors_val)+
  #scale_shape_manual(name="False Discoveries\nAllowance for the\nGene Set Enrichment",values = c(2,19))+
  labs(x="Normalized Enrichment Score",y="GO:BP Gene Sets",subtitle = "Cherry Cryoinjury vs Control",title = "\tgseGO GO Gene Sets Enrichment",size="P-Value\nSignificance")

dev.off()


#GSE GO MF
pdf("./Enrichment/GSE/GO/Drerio/dr_gseGO_MF_GO.pdf",width = 13,height = 12)
#Separate the lists in Up and Down regulated. 
dr_control_GSE_GO_MF_topUp_pathways<-dr_control_GSE_GO_MF@result%>%
  filter(NES>0)%>%
  top_n(20,wt=-pvalue)

dr_control_GSE_GO_MF_topDown_pathways<-dr_control_GSE_GO_MF@result%>%
  filter(NES<0)%>%
  top_n(20,wt=-pvalue)

dr_control_GSE_GO_MF_topPathWays<-bind_rows(dr_control_GSE_GO_MF_topUp_pathways,rev(dr_control_GSE_GO_MF_topDown_pathways))
head(dr_control_GSE_GO_MF_topPathWays)

#Plot the table converted
ggplot(dr_control_GSE_GO_MF_topPathWays,aes(x=NES,y=Description,
                                               size=pvalue,
                                               col=ifelse(NES > 0,"Upregulated","Downregulated")))+
  #shape=ifelse(p.adjust < 0.25,"False Discovery\nRate < 25%","False Discovery\nRate > 25%")))+
  geom_point()+
  expand_limits(x=c(-5,5))+
  theme_classic()+
  scale_color_manual(name="Normalized Enrichment Score", values =colors_val)+
  #scale_shape_manual(name="False Discoveries\nAllowance for the\nGene Set Enrichment",values = c(2,19))+
  labs(x="Normalized Enrichment Score",y="GO:MF Gene Sets",subtitle = "dr Cherry Cryoinjury vs Control",title = "\tgseGO GO Gene Sets Enrichment",size="P-Value\nSignificance")

dev.off()

#GSE GO CC
pdf("./Enrichment/GSE/GO/Drerio/dr_gseGO_CC_GO.pdf",width = 13,height = 12)
#Separate the lists in Up and Down regulated. 
dr_control_GSE_GO_CC_topUp_pathways<-dr_control_GSE_GO_CC@result%>%
  filter(NES>0)%>%
  top_n(20,wt=-pvalue)

dr_control_GSE_GO_CC_topDown_pathways<-dr_control_GSE_GO_CC@result%>%
  filter(NES<0)%>%
  top_n(20,wt=-pvalue)

dr_control_GSE_GO_CC_topPathWays<-bind_rows(dr_control_GSE_GO_CC_topUp_pathways,rev(dr_control_GSE_GO_CC_topDown_pathways))
head(dr_control_GSE_GO_CC_topPathWays)

#Plot the table converted
ggplot(dr_control_GSE_GO_CC_topPathWays,aes(x=NES,y=Description,
                                               size=pvalue,
                                               col=ifelse(NES > 0,"Upregulated","Downregulated")))+
  #shape=ifelse(p.adjust < 0.25,"False Discovery\nRate < 25%","False Discovery\nRate > 25%")))+
  geom_point()+
  expand_limits(x=c(-5,5))+
  theme_classic()+
  scale_color_manual(name="Normalized Enrichment Score", values =colors_val)+
  #scale_shape_manual(name="False Discoveries\nAllowance for the\nGene Set Enrichment",values = c(2,19))+
  labs(x="Normalized Enrichment Score",y="GO:CC Gene Sets",subtitle = "drCherry Cryoinjury vs Control",title = "\tgseGO GO Gene Sets Enrichment",size="P-Value\nSignificance")

dev.off()


#####################################################
#                                                   #
#               HUMAN                               #
#                                                   #
#####################################################
#Create the folder for the HUMAN directrory
dir.create("./Enrichment/GSE/GO/Hsapiens")
#Crete the ranked_genelist
#First clean the data a little bit.

#Remove all duplicated data from the stats table, leaving all the unique ENTREZID
dim(anno_human_control_stats)
anno_human_control_stats<-anno_human_control_stats[-which(duplicated(anno_human_control_stats$ENTREZID)),]
dim(anno_human_control_stats)

#first numeric vector using the log2FoldChange
hs_control_ranked_genes<-anno_human_control_stats$log2FoldChange

#second step name the vector
names(hs_control_ranked_genes)<-as.character(anno_human_control_stats$ENTREZID)

#Sort the list in decreasing order.
hs_control_ranked_genes<-sort(hs_control_ranked_genes,decreasing = T)
length(hs_control_ranked_genes)


##############################
#        GSE                 #
#   Using GO→ BP            # 
#                            #
##############################
#Start with the GSE in GO databases with BP.
hs_control_GSE_GO_BP<-gseGO(geneList = hs_control_ranked_genes,ont="BP",OrgDb = HsOrgDb,
                               keyType = "ENTREZID",
                               nPerm=1000,
                               minGSSize = 10,
                               pvalueCutoff = 0.05,
                               pAdjustMethod = "BH")


head(hs_control_GSE_GO_BP@result)

##############################
#        GSE                 #
#   Using GO→ MF            # 
#                            #
##############################

#GSE in GO database with MF.
hs_control_GSE_GO_MF<-gseGO(geneList = hs_control_ranked_genes,ont="MF",OrgDb = HsOrgDb,
                               keyType = "ENTREZID",
                               nPerm=1000,
                               minGSSize = 10,
                               pvalueCutoff = 0.05,
                               pAdjustMethod = "BH")


head(subset(hs_control_GSE_GO_MF@result,NES<0))


##############################
#        GSE                 #
#   Using GO→ CC            # 
#                            #
##############################


#GSE in GO database with CC.
hs_control_GSE_GO_CC<-gseGO(geneList = hs_control_ranked_genes,ont="CC",OrgDb = HsOrgDb,
                               keyType = "ENTREZID",
                               nPerm=1000,
                               minGSSize = 10,
                               pvalueCutoff = 0.05,
                               pAdjustMethod = "BH")


head(hs_control_GSE_GO_CC@result)




#SAVE THE TABLES
#Save the data of Dr into the specific file each table to pvalueCutOff=0.05.
#GSE for GO:BP 
write.table(data.frame(hs_control_GSE_GO_BP@result),"./Enrichment/GSE/GO/Hsapiens/hs_control_GSE_GO_BP.txt",sep="\t", row.names = F,col.names = T,quote = F)

#GSE for GO:MF
write.table(data.frame(hs_control_GSE_GO_MF@result),"./Enrichment/GSE/GO/Hsapiens/hs_control_GSE_GO_MF.txt",sep="\t",row.names = F,col.names = T,quote = F)

#GSE for CC
write.table(data.frame(hs_control_GSE_GO_CC@result),"./Enrichment/GSE/GO/Hsapiens/hs_control_GSE_GO_BP.txt",sep="\t",row.names = F,col.names = T,quote = F)



#Represent the figure.
#Represent the figures, start with the BP,MF,CC.
colors_val<-c("#377EB8","#E41A1C")


#GSE GO BP
#Plot the HUMAN GSE in GO using BP
pdf("./Enrichment/GSE/GO/Hsapiens/hs_gseGO_BP_GO.pdf",width = 13,height = 12)

#Separate the lists in Up and Down regulated. 
hs_control_GSE_GO_BP_topUp_pathways<-hs_control_GSE_GO_BP@result%>%
  filter(NES>0)%>%
  top_n(30,wt=-pvalue)
head(hs_control_GSE_GO_BP_topUp_pathways,30)

hs_control_GSE_GO_BP_topDown_pathways<-hs_control_GSE_GO_BP@result%>%
  filter(NES<0)%>%
  top_n(30,wt=-pvalue)
head(hs_control_GSE_GO_BP_topDown_pathways)

hs_control_GSE_GO_BP_topPathWays<-bind_rows(hs_control_GSE_GO_BP_topUp_pathways,rev(hs_control_GSE_GO_BP_topDown_pathways))
head(hs_control_GSE_GO_BP_topPathWays)

#Plot the table converted

ggplot(hs_control_GSE_GO_BP_topPathWays,aes(x=NES,y=Description,
                                               size=pvalue,
                                               col=ifelse(NES > 0,"Upregulated","Downregulated")))+
  #shape=ifelse(p.adjust < 0.25,"False Discovery\nRate < 25%","False Discovery\nRate > 25%")))+
  geom_point()+
  expand_limits(x=c(-5,5))+
  theme_classic()+
  scale_color_manual(name="Normalized Enrichment Score", values =colors_val)+
  #scale_shape_manual(name="False Discoveries\nAllowance for the\nGene Set Enrichment",values = c(2,19))+
  labs(x="Normalized Enrichment Score",y="GO:BP Gene Sets",subtitle = "Hs Cherry Cryoinjury vs Control",title = "\tgseGO GO Gene Sets Enrichment",size="P-Value\nSignificance")

dev.off()



#GSE GO MF
pdf("./Enrichment/GSE/GO/Hsapiens/hs_gseGO_MF_GO.pdf",width = 13,height = 12)

#Separate the lists in Up and Down regulated. 
hs_control_GSE_GO_MF_topUp_pathways<-hs_control_GSE_GO_MF@result%>%
  filter(NES>0)%>%
  top_n(30,wt=-pvalue)
head(hs_control_GSE_GO_MF_topUp_pathways,30)

hs_control_GSE_GO_MF_topDown_pathways<-hs_control_GSE_GO_MF@result%>%
  filter(NES<0)%>%
  top_n(30,wt=-pvalue)
head(hs_control_GSE_GO_MF_topDown_pathways)

hs_control_GSE_GO_MF_topPathWays<-bind_rows(hs_control_GSE_GO_MF_topUp_pathways,rev(hs_control_GSE_GO_MF_topDown_pathways))
head(hs_control_GSE_GO_MF_topPathWays)

#Plot the table converted

ggplot(hs_control_GSE_GO_MF_topPathWays,aes(x=NES,y=Description,
                                               size=pvalue,
                                               col=ifelse(NES > 0,"Upregulated","Downregulated")))+
  #shape=ifelse(p.adjust < 0.25,"False Discovery\nRate < 25%","False Discovery\nRate > 25%")))+
  geom_point()+
  expand_limits(x=c(-5,5))+
  theme_classic()+
  scale_color_manual(name="Normalized Enrichment Score", values =colors_val)+
  #scale_shape_manual(name="False Discoveries\nAllowance for the\nGene Set Enrichment",values = c(2,19))+
  labs(x="Normalized Enrichment Score",y="GO:MF Gene Sets",subtitle = "Hs Cherry Cryoinjury vs Control",title = "\tgseGO GO Gene Sets Enrichment",size="P-Value\nSignificance")

dev.off()


#GSE GO CC
pdf("./Enrichment/GSE/GO/Hsapiens/hs_gseGO_CC_GO.pdf",width = 13,height = 12)

#Separate the lists in Up and Down regulated. 
hs_control_GSE_GO_CC_topUp_pathways<-hs_control_GSE_GO_CC@result%>%
  filter(NES>0)%>%
  top_n(30,wt=-pvalue)
head(hs_control_GSE_GO_CC_topUp_pathways,30)

hs_control_GSE_GO_CC_topDown_pathways<-hs_control_GSE_GO_CC@result%>%
  filter(NES<0)%>%
  top_n(30,wt=-pvalue)
head(hs_control_GSE_GO_CC_topDown_pathways)

hs_control_GSE_GO_CC_topPathWays<-bind_rows(hs_control_GSE_GO_CC_topUp_pathways,rev(hs_control_GSE_GO_CC_topDown_pathways))
head(hs_control_GSE_GO_CC_topPathWays)

#Plot the table converted

ggplot(hs_control_GSE_GO_CC_topPathWays,aes(x=NES,y=Description,
                                               size=pvalue,
                                               col=ifelse(NES > 0,"Upregulated","Downregulated")))+
  #shape=ifelse(p.adjust < 0.25,"False Discovery\nRate < 25%","False Discovery\nRate > 25%")))+
  geom_point()+
  expand_limits(x=c(-5,5))+
  theme_classic()+
  scale_color_manual(name="Normalized Enrichment Score", values =colors_val)+
  #scale_shape_manual(name="False Discoveries\nAllowance for the\nGene Set Enrichment",values = c(2,19))+
  labs(x="Normalized Enrichment Score",y="GO:CC Gene Sets",subtitle = "Hs Cherry Cryoinjury vs Control",title = "\tgseGO GO Gene Sets Enrichment",size="P-Value\nSignificance")

dev.off()






#####################################################
#                                                   #
#               MOUSE                               #
#                                                   #
#####################################################
#Create the folder for the HUMAN directrory
dir.create("./Enrichment/GSE/GO/Mmusculus")
#Crete the ranked_genelist
#First clean the data a little bit.

#Remove all duplicated data from the stats table, leaving all the unique ENTREZID
dim(anno_mouse_control_stats)
anno_mouse_control_stats<-anno_mouse_control_stats[-which(duplicated(anno_mouse_control_stats$ENTREZID)),]
dim(anno_mouse_control_stats)

#first numeric vector using the log2FoldChange
mm_control_ranked_genes<-anno_mouse_control_stats$log2FoldChange

#second step name the vector
names(mm_control_ranked_genes)<-as.character(anno_mouse_control_stats$ENTREZID)

#Sort the list in decreasing order.
mm_control_ranked_genes<-sort(mm_control_ranked_genes,decreasing = T)
length(mm_control_ranked_genes)


##############################
#        GSE                 #
#   Using GO→ BP            # 
#                            #
##############################
#Start with the GSE in GO databases with BP.
mm_control_GSE_GO_BP<-gseGO(geneList = mm_control_ranked_genes,ont="BP",OrgDb = MmOrgDb,
                               keyType = "ENTREZID",
                               nPerm=1000,
                               minGSSize = 10,
                               pvalueCutoff = 0.05,
                               pAdjustMethod = "BH")


head(mm_control_GSE_GO_BP@result)

##############################
#        GSE                 #
#   Using GO→ MF            # 
#                            #
##############################

#GSE in GO database with MF.
mm_control_GSE_GO_MF<-gseGO(geneList =mm_control_ranked_genes,ont="MF",OrgDb = MmOrgDb,
                               keyType = "ENTREZID",
                               nPerm=1000,
                               minGSSize = 10,
                               pvalueCutoff = 0.05,
                               pAdjustMethod = "BH")


head(mm_control_GSE_GO_MF@result)



##############################
#        GSE                 #
#   Using GO→ CC            # 
#                            #
##############################

#GSE in GO database with CC.
mm_control_GSE_GO_CC<-gseGO(geneList = mm_control_ranked_genes,ont="CC",OrgDb = MmOrgDb,
                               keyType = "ENTREZID",
                               nPerm=1000,
                               minGSSize = 10,
                               pvalueCutoff = 0.05,
                               pAdjustMethod = "BH")


head(mm_control_GSE_GO_CC@result)




#SAVE THE TABLES
#Save the data of Dr into the specific file each table to pvalueCutOff=0.05.
#GSE for GO:BP 
write.table(data.frame(mm_control_GSE_GO_BP@result),"./Enrichment/GSE/GO/Mmusculus/mm_control_GSE_GO_BP.txt",sep="\t", row.names = F,col.names = T,quote = F)


#GSE for GO:MF
write.table(data.frame(mm_control_GSE_GO_MF@result),"./Enrichment/GSE/GO/Mmusculus/mm_control_GSE_GO_MF.txt",sep="\t",row.names = F,col.names = T,quote = F)

#GSE for CC
write.table(data.frame(mm_control_GSE_GO_CC@result),"./Enrichment/GSE/GO/Mmusculus/mm_control_GSE_GO_BP.txt",sep="\t",row.names = F,col.names = T,quote = F)


#Represent the figure.
#Represent the figures, start with the BP,MF,CC.
colors_val<-c("#377EB8","#E41A1C")


#GSE GO BP
#Plot the HUMAN GSE in GO using BP
pdf("./Enrichment/GSE/GO/Mmusculus/mm_gseGO_BP_GO.pdf",width = 13,height = 12)

#Separate the lists in Up and Down regulated. 
mm_control_GSE_GO_BP_topUp_pathways<-mm_control_GSE_GO_BP@result%>%
  filter(NES>0)%>%
  top_n(30,wt=-pvalue)
head(mm_control_GSE_GO_BP_topUp_pathways,30)

mm_control_GSE_GO_BP_topDown_pathways<-mm_control_GSE_GO_BP@result%>%
  filter(NES<0)%>%
  top_n(30,wt=-pvalue)
head(mm_control_GSE_GO_BP_topDown_pathways)

mm_control_GSE_GO_BP_topPathWays<-bind_rows(mm_control_GSE_GO_BP_topUp_pathways,rev(mm_control_GSE_GO_BP_topDown_pathways))
head(mm_control_GSE_GO_BP_topPathWays)

#Plot the table converted

ggplot(mm_control_GSE_GO_BP_topPathWays,aes(x=NES,y=Description,
                                               size=pvalue,
                                               col=ifelse(NES > 0,"Upregulated","Downregulated")))+
  #shape=ifelse(p.adjust < 0.25,"False Discovery\nRate < 25%","False Discovery\nRate > 25%")))+
  geom_point()+
  expand_limits(x=c(-5,5))+
  theme_classic()+
  scale_color_manual(name="Normalized Enrichment Score", values =colors_val)+
  #scale_shape_manual(name="False Discoveries\nAllowance for the\nGene Set Enrichment",values = c(2,19))+
  labs(x="Normalized Enrichment Score",y="GO:BP Gene Sets",subtitle = "Mm Cherry Cryoinjury vs Control",title = "\tgseGO GO Gene Sets Enrichment",size="P-Value\nSignificance")

dev.off()



#GSE GO MF
pdf("./Enrichment/GSE/GO/Mmusculus/mm_gseGO_MF_GO.pdf",width = 23,height = 16)

#Separate the lists in Up and Down regulated. 
mm_control_GSE_GO_MF_topUp_pathways<-mm_control_GSE_GO_MF@result%>%
  filter(NES>0)%>%
  top_n(30,wt=-pvalue)
head(mm_control_GSE_GO_MF_topUp_pathways,30)

mm_control_GSE_GO_MF_topDown_pathways<-mm_control_GSE_GO_MF@result%>%
  filter(NES<0)%>%
  top_n(30,wt=-pvalue)
head(mm_control_GSE_GO_MF_topDown_pathways)

mm_control_GSE_GO_MF_topPathWays<-bind_rows(mm_control_GSE_GO_MF_topUp_pathways,rev(mm_control_GSE_GO_MF_topDown_pathways))
head(mm_control_GSE_GO_MF_topPathWays)

#Plot the table converted

ggplot(mm_control_GSE_GO_MF_topPathWays,aes(x=NES,y=Description,
                                               size=pvalue,
                                               col=ifelse(NES > 0,"Upregulated","Downregulated")))+
  #shape=ifelse(p.adjust < 0.25,"False Discovery\nRate < 25%","False Discovery\nRate > 25%")))+
  geom_point()+
  expand_limits(x=c(-5,5))+
  theme_classic()+
  scale_color_manual(name="Normalized Enrichment Score", values =colors_val)+
  #scale_shape_manual(name="False Discoveries\nAllowance for the\nGene Set Enrichment",values = c(2,19))+
  labs(x="Normalized Enrichment Score",y="GO:MF Gene Sets",subtitle = "Mm Cherry Cryoinjury vs Control",title = "\tgseGO GO Gene Sets Enrichment",size="P-Value\nSignificance")

dev.off()


#GSE GO CC
pdf("./Enrichment/GSE/GO/Mmusculus/mm_gseGO_CC_GO.pdf",width = 13,height = 12)

#Separate the lists in Up and Down regulated. 
mm_control_GSE_GO_CC_topUp_pathways<-mm_control_GSE_GO_CC@result%>%
  filter(NES>0)%>%
  top_n(30,wt=-pvalue)
head(mm_control_GSE_GO_CC_topUp_pathways,30)

mm_control_GSE_GO_CC_topDown_pathways<-mm_control_GSE_GO_CC@result%>%
  filter(NES<0)%>%
  top_n(30,wt=-pvalue)
head(mm_control_GSE_GO_CC_topDown_pathways)

mm_control_GSE_GO_CC_topPathWays<-bind_rows(mm_control_GSE_GO_CC_topUp_pathways,rev(mm_control_GSE_GO_CC_topDown_pathways))
head(mm_control_GSE_GO_CC_topPathWays)

#Plot the table converted

ggplot(mm_control_GSE_GO_CC_topPathWays,aes(x=NES,y=Description,
                                               size=pvalue,
                                               col=ifelse(NES > 0,"Upregulated","Downregulated")))+
  #shape=ifelse(p.adjust < 0.25,"False Discovery\nRate < 25%","False Discovery\nRate > 25%")))+
  geom_point()+
  expand_limits(x=c(-5,5))+
  theme_classic()+
  scale_color_manual(name="Normalized Enrichment Score", values =colors_val)+
  #scale_shape_manual(name="False Discoveries\nAllowance for the\nGene Set Enrichment",values = c(2,19))+
  labs(x="Normalized Enrichment Score",y="GO:CC Gene Sets",subtitle = "Mm Cherry Cryoinjury vs Control",title = "\tgseGO GO Gene Sets Enrichment",size="P-Value\nSignificance")

dev.off()











#####################################################
#                                                   #
#           KEGGG ENRICHMENT with  GSE              #
#                                                   #
#            USINGGENE SET FOR KEGG                 #  
#                   PATHWAYS                        #
#                                                   #
#####################################################



#Start doing the enrichment analysis with gseGO
#Set the respective folder for that.

dir.create("./Enrichment/GSE/KEGG/")

#Use the ranked lists created above as they are both the same for the enrichment process with gseKEGG.

#####################################################
#                                                   #
#               ZEBRAFISH                           #
#                                                   #
#####################################################
#Start with zebrafish
dir.create("./Enrichment/GSE/KEGG/Drerio")

#Refresh for the organism used for GSE KEGG
search_kegg_organism("danio", by = "scientific_name", ignore.case = T)

#Start the enrichment calling KEGG
dr_control_GSE_KEGG<-gseKEGG(geneList = dr_control_ranked_genes, organism = "dre",
                                nPerm = 1000,
                                minGSSize = 10,
                                pvalueCutoff = 0.05,
                                pAdjustMethod = "BH")


#Save the pathway table for the GSE with KEGG.
write.table(data.frame(dr_control_GSE_KEGG@result),"./Enrichment/GSE/KEGG/Drerio/dr_control_GSE_KEGG.txt",sep = "\t", quote = F, row.names = F,col.names = T)
head(dr_control_GSE_KEGG@result)
#Save the figures desired.


#####################################################
#                                                   #
#               HUMAN                               #
#                                                   #
#####################################################
#Start with Human
dir.create("./Enrichment/GSE/KEGG/Hsapiens")

#Refresh for the organism used for GSE KEGG
search_kegg_organism("sapiens", by = "scientific_name", ignore.case = T)

#Start the enrichment calling KEGG
hs_control_GSE_KEGG<-gseKEGG(geneList = hs_control_ranked_genes, organism = "hsa",
                                nPerm = 1000,
                                minGSSize = 10,
                                pvalueCutoff = 0.05,
                                pAdjustMethod = "BH")

#Save the pathway table for the GSE with KEGG.
write.table(data.frame(hs_control_GSE_KEGG@result), "./Enrichment/GSE/KEGG/Hsapiens/hs_control_GSE_KEGG.txt",sep = "\t", quote = F, row.names = F,col.names = T)



#Save the figures desired.


#####################################################
#                                                   #
#               MOUSE                               #
#                                                   #
#####################################################
#Start with mouse
dir.create("./Enrichment/GSE/KEGG/Mmusculus")

#Refresh for the organism used for GSE KEGG
search_kegg_organism("musculus", by = "scientific_name", ignore.case = T)

#Start the enrichment calling KEGG
mm_control_GSE_KEGG<-gseKEGG(geneList = mm_control_ranked_genes, organism = "mmu",
                                nPerm = 1000,
                                minGSSize = 10,
                                pvalueCutoff = 0.05,
                                pAdjustMethod = "BH")

head(mm_control_GSE_KEGG@result)
#Save the pathway table for the GSE with KEGG.
write.table(data.frame(mm_control_GSE_KEGG@result), "./Enrichment/GSE/KEGG/Mmusculus/mm_control_GSE_KEGG.txt",sep = "\t", quote = F, row.names = F,col.names = T)




#Start with fGSEA

#RANKED log2FoldChange + ENTREZID list.


#####################################################
#                                                   #
#           fGSEA ENRICHMENT with  HALLMARKS        #
#                                                   #
#            USING L2FC ranked list                 #  
#                                                   #
#             HALLMARKS                             #
#                                                   #
#####################################################

#Download data in Rdat format from http://bioinf.wehi.edu.au/software/MSigDB/

#These process is only doable with HUMAN and MOUSE yet.
#It can be done using SYMBOLS or ENTREZIDs→ using ENTEREZIDS for this analysis.



#####################################################
#                                                   #
#               HUMAN                               #
#                                                   #
#               HALLMARKS                           #
#                                                   #
#####################################################
#Start with Human

#Create the pathway fo the Hallmarks data.
dir.create(recursive = T, "./Enrichment/GSEA/Hallmarks/Hsapiens")

#Load the Hallmark data
load("../gsea_bundle/human_H_v5p2.rdata")

#Save the 50 pathways composing the hallmark GSEAinto an object
hs_hallmark_fgsea<-Hs.H


#Use the ranked ENTREZIDs used before.
length(hs_control_ranked_genes)

#Conduct analysis the GSEA enrichment for the Hallmarks!
hs_control_fgsea_hallmarks<-fgsea(pathways = hs_hallmark_fgsea,stats = hs_control_ranked_genes, minSize = 15, maxSize =500, nperm = 1000 )

#Now lets see the 30 top adjusted.
head(hs_control_fgsea_hallmarks, n=30)


#Map the leadingEdges from ENTREZIDs to SYMBOLS for being more human readable.
hs_control_fgsea_hallmarks[,leadingEdge:=lapply(leadingEdge,mapIds, x=org.Hs.eg.db, keytype="ENTREZID",column="SYMBOL")]
head(hs_control_fgsea_hallmarks)


#Order the data according to padj value
hs_control_fgsea_hallmarks<-hs_control_fgsea_hallmarks[order(hs_control_fgsea_hallmarks[,3])]

#Remove the data which padj=1
hs_control_fgsea_hallmarks<-hs_control_fgsea_hallmarks[hs_control_fgsea_hallmarks$padj!=1]

#Save the data of the GSEA HALLMARK PATHWAYS, using a special function for writting the data table data frame with data.table::fwrite..
data.table::fwrite(hs_control_fgsea_hallmarks, "./Enrichment/GSEA/Hallmarks/Hsapiens/hs_control_fgsea_hallmark.txt",sep = "\t",col.names = T, row.names = F, quote = F)
head(hs_control_fgsea_hallmarks)


#Choose from the list any of the interested gsea plot enrichment interested.
pdf("./Enrichment/GSEA/Hallmarks/Hsapiens/hs_fgsea_hallmarks.pdf",width = 16,height = 10)
plotEnrichment(hs_hallmark_fgsea[["HALLMARK_HYPOXIA"]], hs_control_ranked_genes) + labs(title="Hypoxia Hallmark")
dev.off()



#Plot the fgsea table for significance based on Enrichment Scores.
#GSEA table plot.

hs_control_fgsea_topUp_hallmarks<-hs_control_fgsea_hallmarks%>%
  filter(ES > 0) %>%
  top_n(10,wt=-padj)
hs_control_fgsea_topDown_hallmarks<-hs_control_fgsea_hallmarks%>%
  filter(ES<0)%>%
  top_n(10,wt=-padj)
head(hs_control_fgsea_topUp_hallmarks)
head(hs_control_fgsea_topDown_hallmarks)
hs_control_fgsea_topPathways_hallmarks<-bind_rows(hs_control_fgsea_topUp_hallmarks,hs_control_fgsea_topDown_hallmarks)%>%
  arrange(-ES)
head(hs_control_fgsea_topPathways_hallmarks)

pdf("./Enrichment/GSEA/Hallmarks/Hsapiens/hs_control_fgsea_topPathWays_hallmarks_TablePlot.pdf",width = 9,height = 6)
#png("./plots/gseaPlottableHallmarks.png",1000,1000, pointsize = 20)
plotGseaTable(pathways = hs_hallmark_fgsea[hs_control_fgsea_topPathways_hallmarks$pathway],
              hs_control_ranked_genes,
              hs_control_fgsea_hallmarks,
              gseaParam = 0.5)
dev.off()



#Plot in a more visual way the pathways up and down regulated.
#Create the plot of the gsea file using ggplplot
colors_val<-c("#377EB8","#E41A1C")

pdf("./Enrichment/GSEA/Hallmarks/Hsapiens/hs_control_fgsea_hallmark_ggplot.pdf",width = 8,height = 12)
hs_control_fgsea_hallmarks%>%
  top_n(30,wt =- sort(padj))%>%
  ggplot(aes(x=NES,y=pathway,
             size=pval,
             col=ifelse(NES > 0,"Upregulated/Cryoinjury","Downregulated/Control")))+
  #shape=ifelse(padj < 0.25,"False Discovery\nRate < 25%","False Discovery\nRate > 25%")))+
  geom_point()+
  expand_limits(x=c(-5,5))+
  theme_classic()+
  scale_color_manual(name="Normalized Enrichment Score", values =colors_val)+
  #scale_shape_manual(name="False Discoveries\nAllowance for the\nGene Set Enrichment",values = c(2,19))+
  labs(x="Normalized Enrichment Score",y="Hallmark Gene Sets",subtitle = "Hs BFP Cryoinjury vs Control",title = "\tFGSEA HallMkark Gene Sets Enrichment Analysis",size="P-Value\nSignificance")

dev.off()






#####################################################
#                                                   #
#               MOUSE                               #
#                                                   #
#               HALLMARKS                           #
#                                                   #
#####################################################
#Start with Mouse Hallmarks

#Create the pathway fo the Hallmarks data.
dir.create(recursive = T, "./Enrichment/GSEA/Hallmarks/Mmusculus")

#Load the Hallmark data
load("../gsea_bundle/mouse_H_v5p2.rdata")

#Save the 50 pathways composing the hallmark GSEAinto an object
mm_hallmark_fgsea<-Mm.H


#Use the ranked ENTREZIDs used before.
length(mm_control_ranked_genes)

#Conduct analysis the GSEA enrichment for the Hallmarks!
mm_control_fgsea_hallmarks<-fgsea(pathways = mm_hallmark_fgsea,stats = mm_control_ranked_genes, minSize = 15, maxSize =500, nperm = 1000 )

#Now lets see the 30 top adjusted.
head(mm_control_fgsea_hallmarks, n=30)


#Map the leadingEdges from ENTREZIDs to SYMBOLS for being more human readable.
mm_control_fgsea_hallmarks[,leadingEdge:=lapply(leadingEdge,mapIds, x=org.Mm.eg.db, keytype="ENTREZID",column="SYMBOL")]
head(mm_control_fgsea_hallmarks)


#Order the data according to padj value
mm_control_fgsea_hallmarks<-mm_control_fgsea_hallmarks[order(mm_control_fgsea_hallmarks[,3])]

#Remove the data which padj=1
mm_control_fgsea_hallmarks<-mm_control_fgsea_hallmarks[mm_control_fgsea_hallmarks$padj!=1]

#Save the data of the GSEA HALLMARK PATHWAYS, using a special function for writting the data table data frame with data.table::fwrite..
data.table::fwrite(mm_control_fgsea_hallmarks, "./Enrichment/GSEA/Hallmarks/Mmusculus/mm_control_fgsea_hallmark.txt",sep = "\t",col.names = T, row.names = F, quote = F)

head(mm_control_fgsea_hallmarks)

#Choose from the list any of the interested gsea plot enrichment interested.
pdf("./Enrichment/GSEA/Hallmarks/Mmusculus/mm_fgsea_hallmarks.pdf",width = 16,height = 10)
plotEnrichment(mm_hallmark_fgsea[["HALLMARK_G2M_CHECKPOINT"]], mm_control_ranked_genes) + labs(title="G2M_CHECKPOINT")
dev.off()



#Plot the fgsea table for significance based on Enrichment Scores.
#GSEA table plot.

mm_control_fgsea_topUp_hallmarks<-mm_control_fgsea_hallmarks%>%
  filter(ES > 0) %>%
  top_n(10,wt=-padj)
mm_control_fgsea_topDown_hallmarks<-mm_control_fgsea_hallmarks%>%
  filter(ES<0)%>%
  top_n(10,wt=-padj)
head(mm_control_fgsea_topUp_hallmarks)
head(mm_control_fgsea_topDown_hallmarks)
mm_control_fgsea_topPathways_hallmarks<-bind_rows(mm_control_fgsea_topUp_hallmarks,mm_control_fgsea_topDown_hallmarks)%>%
  arrange(-ES)
head(mm_control_fgsea_topPathways_hallmarks)

pdf("./Enrichment/GSEA/Hallmarks/Mmusculus/mm_control_fgsea_topPathWays_hallmarks_TablePlot.pdf",width = 9,height = 6)
#png("./plots/gseaPlottableHallmarks.png",1000,1000, pointsize = 20)
plotGseaTable(pathways = mm_hallmark_fgsea[mm_control_fgsea_topPathways_hallmarks$pathway],
              mm_control_ranked_genes,
              mm_control_fgsea_hallmarks,
              gseaParam = 0.5)
dev.off()



#Plot in a more visual way the pathways up and down regulated.
#Create the plot of the gsea file using ggplplot
colors_val<-c("#377EB8","#E41A1C")

pdf("./Enrichment/GSEA/Hallmarks/Mmusculus/mm_control_fgsea_hallmark_ggplot.pdf",width = 8,height = 12)
mm_control_fgsea_hallmarks%>%
  top_n(30,wt =- sort(padj))%>%
  ggplot(aes(x=NES,y=pathway,
             size=pval,
             col=ifelse(NES > 0,"Upregulated/Cryoinjury","Downregulated/Control")))+
  #shape=ifelse(padj < 0.25,"False Discovery\nRate < 25%","False Discovery\nRate > 25%")))+
  geom_point()+
  expand_limits(x=c(-5,5))+
  theme_classic()+
  scale_color_manual(name="Normalized Enrichment Score", values =colors_val)+
  #scale_shape_manual(name="False Discoveries\nAllowance for the\nGene Set Enrichment",values = c(2,19))+
  labs(x="Normalized Enrichment Score",y="Hallmark Gene Sets",subtitle = "Mm BFP Cryoinjury vs Control",title = "\tFGSEA HallMkark Gene Sets Enrichment Analysis",size="P-Value\nSignificance")

dev.off()





#####################################################
#                                                   #
#           fGSEA ENRICHMENT with  GO               #
#                                                   #
#            USING L2FC ranked list                 #  
#                                                   #
#             GO DATABASES                          #
#                                                   #
#####################################################

#Download data in Rdat format from http://bioinf.wehi.edu.au/software/MSigDB/

#These process is only doable with HUMAN and MOUSE yet.
#It can be done using SYMBOLS or ENTREZIDs→ using ENTEREZIDS for this analysis.



#####################################################
#                                                   #
#               HUMAN                               #
#                                                   #
#               GO DATABASES                        #
#                                                   #
#####################################################
#Start with Human

#Create the pathway fo the Hallmarks data.
dir.create(recursive = T, "./Enrichment/GSEA/GO/Hsapiens")

#Load the Hallmark data
load("../gsea_bundle/human_c5_v5p2.rdata")

#Save the 50 pathways composing the hallmark GSEAinto an object
hs_GO_fgsea<-Hs.c5


#Use the ranked ENTREZIDs used before.
length(hs_control_ranked_genes)

#Conduct analysis the GSEA enrichment for the Hallmarks!
hs_control_fgsea_GO<-fgsea(pathways = hs_GO_fgsea,stats = hs_control_ranked_genes, minSize = 15, maxSize =500, nperm = 1000 )

#Now lets see the 30 top adjusted.
head(hs_control_fgsea_GO, n=30)


#FOR GO ENRICHMENT REDUCE THE LIST UP TO padj=0.5
hs_control_fgsea_GO<-subset(hs_control_fgsea_GO,padj<0.05)


#Map the leadingEdges from ENTREZIDs to SYMBOLS for being more human readable.
hs_control_fgsea_GO[,leadingEdge:=lapply(leadingEdge,mapIds, x=org.Hs.eg.db, keytype="ENTREZID",column="SYMBOL")]
head(hs_control_fgsea_GO)

#Check how much significative data would be taken with different padj values.
dim(hs_control_fgsea_GO[hs_control_fgsea_GO$padj<.05])

#Order the data according to padj value
hs_control_fgsea_GO<-hs_control_fgsea_GO[order(hs_control_fgsea_GO[,3])]

#Remove the data which padj=1
hs_control_fgsea_GO<-hs_control_fgsea_GO[hs_control_fgsea_GO$padj!=1]
head(hs_control_fgsea_GO)

#Save the data of the GSEA HALLMARK PATHWAYS, using a special function for writting the data table data frame with data.table::fwrite..
data.table::fwrite(hs_control_fgsea_GO, "./Enrichment/GSEA/GO/Hsapiens/hs_control_fgsea_hallmark.txt",sep = "\t",col.names = T, row.names = F, quote = F)



#Choose from the list any of the interested gsea plot enrichment interested.
pdf("./Enrichment/GSEA/GO/Hsapiens/hs_fgsea_GO.pdf",width = 16,height =16)
plotEnrichment(hs_GO_fgsea[["GO_CARDIAC_CHAMBER_DEVELOPMENT"]], hs_control_ranked_genes) + labs(title="GO CARDIAC CHAMBER DEVELOPMENT")
dev.off()


#Plot the fgsea table for significance based on Enrichment Scores.
#GSEA table plot.
library(dplyr)
hs_control_fgsea_topUp_GO<-hs_control_fgsea_GO%>%
  filter(ES > 0) %>%
  top_n(n=30,wt=-pval)

hs_control_fgsea_topDown_GO<-hs_control_fgsea_GO %>%
  filter(ES<0) %>%
  top_n(n = 30,wt=-pval)
#Keep i nmind if the padj value are the same, then the top_n function will no work properly.

hs_control_fgsea_topPathways_GO<-bind_rows(hs_control_fgsea_topUp_GO,hs_control_fgsea_topDown_GO)%>%
  arrange(-ES)
head(hs_control_fgsea_topPathways_GO)

pdf("./Enrichment/GSEA/GO/Hsapiens/hs_control_fgsea_topPathWays_GO_TablePlot.pdf",width =30,height = 50)
#png("./plots/gseaPlottableHallmarks.png",1000,1000, pointsize = 20)
plotGseaTable(pathways = hs_GO_fgsea[hs_control_fgsea_topPathways_GO$pathway],
              hs_control_ranked_genes,
              hs_control_fgsea_GO,
              gseaParam = 0.5)
dev.off()

#Plot in a more visual way the pathways up and down regulated.
#Create the plot of the gsea file using ggplplot
colors_val<-c("#377EB8","#E41A1C")

pdf("./Enrichment/GSEA/GO/Hsapiens/hs_control_fgsea_GO_ggplot.pdf",width = 20,height = 25)
hs_control_fgsea_topPathways_GO%>%
  top_n(60,wt =- sort(pval))%>%
  ggplot(aes(x=NES,y=pathway,
             size=pval,
             col=ifelse(NES > 0,"Upregulated/Cryoinjury","Downregulated/Control")))+
  #shape=ifelse(padj < 0.25,"False Discovery\nRate < 25%","False Discovery\nRate > 25%")))+
  geom_point()+
  expand_limits(x=c(-5,5))+
  theme_classic()+
  scale_color_manual(name="Normalized Enrichment Score", values =colors_val)+
  #scale_shape_manual(name="False Discoveries\nAllowance for the\nGene Set Enrichment",values = c(2,19))+
  labs(x="Normalized Enrichment Score",y="Hallmark Gene Sets",subtitle = "Hs BFP Cryoinjury vs Control",title = "\tFGSEA GO Gene Sets Enrichment Analysis",size="P-Value\nSignificance")

dev.off()









#PROCEED with the same process but for MOUSE DATABASE.
#####################################################
#                                                   #
#               MOUSE                               #
#                                                   #
#               GO DATABASES                        #
#                                                   #
#####################################################
#Start with MOUSE

#Create the pathway fo the Hallmarks data.
dir.create(recursive = T, "./Enrichment/GSEA/GO/Mmusculus")

#Load the Hallmark data
load("../gsea_bundle/mouse_c5_v5p2.rdata")

#Save the 50 pathways composing the hallmark GSEAinto an object
mm_GO_fgsea<-Mm.c5


#Use the ranked ENTREZIDs used before.
length(mm_control_ranked_genes)

#Conduct analysis the GSEA enrichment for the Hallmarks!
mm_control_fgsea_GO<-fgsea(pathways = mm_GO_fgsea,stats = mm_control_ranked_genes, minSize = 15, maxSize =500, nperm = 1000 )

#Now lets see the 30 top adjusted.
head(mm_control_fgsea_GO, n=30)


#FOR GO ENRICHMENT REDUCE THE LIST UP TO padj=0.05
mm_control_fgsea_GO<-subset(mm_control_fgsea_GO,padj<0.05)


#Map the leadingEdges from ENTREZIDs to SYMBOLS for being more human readable.
mm_control_fgsea_GO[,leadingEdge:=lapply(leadingEdge,mapIds, x=org.Mm.eg.db, keytype="ENTREZID",column="SYMBOL")]
head(mm_control_fgsea_GO)

#Check how much significative data would be taken with different padj values.
dim(mm_control_fgsea_GO[mm_control_fgsea_GO$padj<.05])

#Order the data according to padj value
mm_control_fgsea_GO<-mm_control_fgsea_GO[order(mm_control_fgsea_GO[,3])]

#Remove the data which padj=1
mm_control_fgsea_GO<-mm_control_fgsea_GO[mm_control_fgsea_GO$padj!=1]
head(mm_control_fgsea_GO)

#Save the data of the GSEA HALLMARK PATHWAYS, using a special function for writting the data table data frame with data.table::fwrite..
data.table::fwrite(mm_control_fgsea_GO, "./Enrichment/GSEA/GO/Mmusculus/mm_control_fgsea_GO.txt",sep = "\t",col.names = T, row.names = F, quote = F)



#Choose from the list any of the interested gsea plot enrichment interested.
pdf("./Enrichment/GSEA/GO/Mmusculus/mm_fgsea_GO.pdf",width = 16,height =16)
plotEnrichment(mm_GO_fgsea[["GO_DNA_DEPENDENT_DNA_REPLICATION_MAINTENANCE_OF_FIDELITY"]], mm_control_ranked_genes) + labs(title="GO GO_DNA_DEPENDENT_DNA_REPLICATION_MAINTENANCE_OF_FIDELITY MOUSE")
dev.off()



#Plot the fgsea table for significance based on Enrichment Scores.
#GSEA table plot.
library(dplyr)
mm_control_fgsea_topUp_GO<-mm_control_fgsea_GO%>%
  filter(ES > 0) %>%
  top_n(n=30,wt=-pval)

mm_control_fgsea_topDown_GO<-mm_control_fgsea_GO %>%
  filter(ES<0) %>%
  top_n(n = 30,wt=-pval)
#Keep i nmind if the padj value are the same, then the top_n function will no work properly.

mm_control_fgsea_topPathways_GO<-bind_rows(mm_control_fgsea_topUp_GO,mm_control_fgsea_topDown_GO)%>%
  arrange(-ES)
head(mm_control_fgsea_topPathways_GO)

pdf("./Enrichment/GSEA/GO/Mmusculus/mm_control_fgsea_topPathWays_GO_TablePlot.pdf",width =30,height = 50)
#png("./plots/gseaPlottableHallmarks.png",1000,1000, pointsize = 20)
plotGseaTable(pathways = mm_GO_fgsea[mm_control_fgsea_topPathways_GO$pathway],
              mm_control_ranked_genes,
              mm_control_fgsea_GO,
              gseaParam = 0.5)
dev.off()

#Plot in a more visual way the pathways up and down regulated.
#Create the plot of the gsea file using ggplplot
colors_val<-c("#377EB8","#E41A1C")

pdf("./Enrichment/GSEA/GO/Mmusculus/mm_control_fgsea_GO_ggplot.pdf",width = 20,height = 25)
mm_control_fgsea_topPathways_GO%>%
  top_n(60,wt =- sort(pval))%>%
  ggplot(aes(x=NES,y=pathway,
             size=pval,
             col=ifelse(NES > 0,"Upregulated/Cryoinjury","Downregulated/Control")))+
  #shape=ifelse(padj < 0.25,"False Discovery\nRate < 25%","False Discovery\nRate > 25%")))+
  geom_point()+
  expand_limits(x=c(-5,5))+
  theme_classic()+
  scale_color_manual(name="Normalized Enrichment Score", values =colors_val)+
  #scale_shape_manual(name="False Discoveries\nAllowance for the\nGene Set Enrichment",values = c(2,19))+
  labs(x="Normalized Enrichment Score",y="Hallmark Gene Sets",subtitle = "Mm BFP Cryoinjury vs Control",title = "\tFGSEA GO Gene Sets Enrichment Analysis",size="P-Value\nSignificance")

dev.off()






#####################################################
#                                                   #
#                   PLOTS                           #
#                                                   #
#                                                   #  
#                                                   #
#                                                   #
#                                                   #
#####################################################



#####################################################
#                                                   #
#                   VOLCANO                         #
#                                                   #
#                                                   #
#####################################################

#####################################################
#                                                   #
#                   HUMAN                           #
#                                                   #
#                                                   #
#####################################################

#Create the directory for plots.
dir.create(recursive = T,"./Plots/Hsapiens/")


#Use the anotated table for human.
head(anno_human_control_stats)

pdf("./Plots/Hsapiens/VolcanoHuman.pdf", width = 12, height =18)
#png("../plots/enhancedVolcano.png", 1000,1000,pointsize = 20)
EnhancedVolcano(anno_human_control_stats,lab=anno_human_control_stats$SYMBOL, x = "log2FoldChange", y="pvalue",legend = c("NS","Log (base 2) fold-change","Adjusted p-value", "Adjusted p-value & Log (base 2) fold-change"),
                legendPosition = "bottom", title = "BFP Cryoinjury vs Control Hs", cutoffLineCol = "deeppink",FCcutoff = 2)#, selectLab = keep)

dev.off()


#####################################################
#                                                   #
#                   ZEBRAFISH                       #
#                                                   #
#                                                   #
#####################################################


#Create the directory for plots.
dir.create(recursive = T,"./Plots/Zebrafish/")


#Use the anotated table for zebrafish.
head(anno_zebrafish_control_stats)

pdf("./Plots/Zebrafish/VolcanoZebra.pdf", width = 12, height =18)
#png("../plots/enhancedVolcano.png", 1000,1000,pointsize = 20)
EnhancedVolcano(anno_zebrafish_control_stats,lab=anno_zebrafish_control_stats$SYMBOL, x = "log2FoldChange", y="pvalue",legend = c("NS","Log (base 2) fold-change","Adjusted p-value", "Adjusted p-value & Log (base 2) fold-change"),
                legendPosition = "bottom", title = "BFP Cryoinjury vs Control Zebra", cutoffLineCol = "deeppink",FCcutoff = 2)#, selectLab = keep)

dev.off()



#####################################################
#                                                   #
#                   MOUSE                           #
#                                                   #
#                                                   #
#####################################################


#Create the directory for plots.
dir.create(recursive = T,"./Plots/Mouse/")


#Use the anotated table for zebrafish.
head(anno_mouse_control_stats)

pdf("./Plots/Mouse/VolcanoMouse.pdf", width = 12, height =18)
#png("../plots/enhancedVolcano.png", 1000,1000,pointsize = 20)
EnhancedVolcano(anno_mouse_control_stats,lab=anno_mouse_control_stats$SYMBOL, x = "log2FoldChange", y="pvalue",legend = c("NS","Log (base 2) fold-change","Adjusted p-value", "Adjusted p-value & Log (base 2) fold-change"),
                legendPosition = "bottom", title = "BFP Cryoinjury vs Control Mouse", cutoffLineCol = "deeppink",FCcutoff = 2)#, selectLab = keep)

dev.off()








#####################################################
#                                                   #
#                                                   #
#                     PCA                           #
#####################################################

#Create the data used for the PCA and  Heatmap by normalizitaion of the samples.
#Use rlog from DESeq2 to regularized log tranformation ignoring information about experimental groups, blind=T
dds_rlog_control=rlog(dds_control,blind = T)
dim(dds_rlog_control)
head(assay(dds_rlog_control))
head(sampleTable_control)

#Create the PCA using the DESeq2 funciton
dir.create(recursive = T, path = "./PCA/")
#Plot PCA
pdf(file = "./PCA/control_PCA.pdf",width = 12,height = 12)
z<-plotPCA(dds_rlog_control,intgroup=c("CELLTYPE"))
z+geom_label(aes(label=rownames(sampleTable_control)))
dev.off()
sampleTable_control
#####################################################
#                                                   #
#                   HUMAN                           #
#                                                   #
#                                                   #
#####################################################

#Directory
dir.create(recursive = T,path = "./Heatmap/Hsapiens/")

#Merge to obtain the human ENSEMBL homologue.
hs_control_heatmap<-merge(anno_human_pre_stats_control,data.frame(assay(dds_rlog_control)),by.x="ensembl_gene_id",by.y="row.names")
head(hs_control_heatmap)

#Merge for the statiscal data.
hs_control_heatmap<-merge(anno_human_control_stats,hs_control_heatmap,by.y="hsapiens_homolog_ensembl_gene",by.x="hsapiens_homolog_ensembl_gene")
head(hs_control_heatmap)

#Filter the duplicated ENSEMBL identifiers and keep only the significant data.
#Filter out duplicates
hs_control_heatmap<-hs_control_heatmap[-which(duplicated(hs_control_heatmap$ENTREZID.x)),]

dim(hs_control_heatmap2)
#Save the relevant data for the heatmap top 30 and down top 30.

hs_control_heatmap_topUp<-hs_control_heatmap%>%
  filter(log2FoldChange > 2)%>%
  arrange(desc(log2FoldChange,padj))%>%
  top_n(50,wt=-padj)


#Top down 30.
hs_control_heatmap_topDown<-hs_control_heatmap%>%
  filter(log2FoldChange < -2)%>%
  arrange(log2FoldChange,padj)%>%
  top_n(50,wt=-padj)


#Merge and save the data for morpheus.
hs_control_heatmap_tops<-bind_rows(hs_control_heatmap_topUp,hs_control_heatmap_topUp)
dim(hs_control_heatmap_tops)


#select and reorder the data 

hs_control_heatmap_tops<-hs_control_heatmap_tops%>%
  dplyr::select("SYMBOL.x","GENENAME.x",starts_with("Sande"))


#Save the table ready to heatmap.
write.table(x = hs_control_heatmap_tops, "./Heatmap/Hsapiens/hs_control_heatmap.txt",sep = "\t",col.names = T,row.names = F,quote = F)



#####################################################
#                                                   #
#                   ZEBRAFISH                       #
#                                                   #
#                                                   #
#####################################################

#Directory
dir.create(recursive = T,path = "./Heatmap/Drerio/")
head(assay(dds_rlog_control))
#Merge to obtain the human ENSEMBL homologue.
dr_control_heatmap<-merge(anno_zebrafish_control_stats,data.frame(assay(dds_rlog_control)),by.x="drerio_ensembl_gene_id",by.y="row.names")


#Filter the duplicated ENSEMBL identifiers and keep only the significant data.
#Filter out duplicates
dr_control_heatmap<-dr_control_heatmap[-which(duplicated(dr_control_heatmap$SYMBOL)),]
dim(dr_control_heatmap)
#Save the relevant data for the heatmap top 50 and down top 50.

dr_control_heatmap_topUp<-dr_control_heatmap%>%
  filter(log2FoldChange > 2)%>%
  arrange(desc(log2FoldChange),padj)%>%
  top_n(50,wt=-padj)
head(dr_control_heatmap_topUp)


#Top down 50.
#Use arrange normally as it does it by default ascending.
dr_control_heatmap_topDown<-dr_control_heatmap%>%
  filter(log2FoldChange < -2)%>%
  arrange(log2FoldChange,padj)%>%
  top_n(50,wt=-padj)

#Merge and save the data for morpheus.
dr_control_heatmap_tops<-bind_rows(dr_control_heatmap_topUp,dr_control_heatmap_topUp)
colnames(dr_control_heatmap_tops)


#Observe and select MANUALLY the data for the heatmap.
dr_control_heatmap_tops<-dr_control_heatmap_tops%>%
  dplyr::select("SYMBOL","GENENAME",starts_with("Sande_"))


#Save the table ready to heatmap.
write.table(x = dr_control_heatmap_tops, "./Heatmap/Drerio/dr_control_heatmap.txt",sep = "\t",col.names = T,row.names = F,quote = F)






#####################################################
#                                                   #
#           Differentially expressed                #
#                   Genes                           #
#                    List                           #
#####################################################

#Report a list with all the genes differentially expressed and annotated in Zebrafish.

#Create a folder for this task.
dir.create("./DEgenes")

#Save full annotated data for; Danio, Human and Mouse.
write.table(x = anno_zebrafish_control_stats,"./DEgenes/all_annotated_DExpressed_genes_control_zebrafish.txt",sep="\t",col.names=T,row.names =F,quote=F)
write.table(anno_human_control_stats,"./DEgenes/all_annotated_DExpressed_genes_control_human.txt",sep = "\t",col.names = T,row.names = F,quote = T)
write.table(anno_mouse_control_stats,"./DEgenes/all_annotated_DExpressed_genes_control_mouse.txt",sep = "\t",col.names = T,row.names = F,quote = T)



#Filtered data table for absolute Log2FoldChange > 2 and padj < 0.05
write.table(anno_zebrafish_control_stats[which(abs(anno_zebrafish_control_stats$log2FoldChange) > 2 & anno_zebrafish_control_stats$padj < 0.05),],"./DEgenes/filtered_DExpressed_genes_zebrafish.txt",sep="\t",row.names=F,col.names=T,quote=F)

#Check the length of the file if filtered.
dim(anno_zebrafish_control_stats[which(abs(anno_zebrafish_control_stats$log2FoldChange>2) & anno_zebrafish_control_stats$padj<.05),])

#Write a report of the analysis, contrast used and what have been done.
description_text_st_diet<-c("The following analysis was done for the comparison of """""";
                            
                            Contrast used: .........
                            
                            The results were annotated and enriched using different approaches, using Danio and Human and Mouse homologues.")


#Save the descript file.
write.table(description_text_st_diet,"./Description_text.txt",sep = "\t",col.names = F,row.names = F,quote = F)


#Save session with all the data stored in R current session.
dir.create("./R")
save.image(file = "./R/st_diet_session.RData")


#####################################################
#                                                   #
#                   DONE                            #
#   add extra info such as the graph in all_4       #
#   do enrich_david and plot significant genes      #
#    also do reactome database.db                   #
#   for pathways enrichment                         #
#####################################################
