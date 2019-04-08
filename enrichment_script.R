#Super extended Enrichment for Marcos comparing in Cherry negative injured vs control.
#Clear and clean all references and saved data.
rm(list=ls(all=T))
#Load packages
require(pacman)
p_load("matrixStats","MASS","tximport","sleuth","dplyr","rhdf5","Rhdf5lib","biomaRt","org.Hs.eg.db","org.Mm.eg.db","org.Dr.eg.db","AnnotationFilter","genefilter","GenomicFeatures","AnnotationDbi","clusterProfiler","DOSE",
       "DESeq2","rjson","RDAVIDWebService","RcppParallel","RColorBrewer","stats","stats4","parallel","RcppArmadillo","EnhancedVolcano","fgsea","ReactomePA","pathview","GSEABase","GOplot","GenomicAlignments","rJava","GenomicRanges","ggplot2","data.table","GenomeInfoDb")
#Set pathways
 
getwd()
setwd("/mnt/haus/marius/Carol_data/carol_ordered/")
#Reading transcript levels from Kallistoto interpretate the further data to a Differential Gene expressed DESeq2 or edgeR analysis.
#Transcript level information is summarized to gene level, for DEG analysis.
#install.packages("readr")
#Set the pathway to the files and read recursively
all_files<-list.files(path = "../carol_analysis/results", recursive = T, pattern = "*.h5", full.names = T)
head(all_files)

#Get the data and add names to localize it.
#Order the values numerically
all_files<-all_files[order(nchar(all_files),all_files)]
head(all_files)
name_files<-list.files("/mnt/haus/marius/Carol_data/fastq/", pattern = "*.fastq.gz")
name_files
#Subsittute all the string for a shorter name.
name_files<-gsub(pattern = "*..CarolinaPoyatos_RNA_Seq_Directional_S.._L005_R1_001.fastq.gz",replacement = '', name_files)
(name_files)

#Change the names of the files.
names(all_files)<-name_files
(all_files)


#BIOMART DATA ANNOTATED.

#Fisrt create a table of genes/transcript or what is in the column interested to switch to.
filterValues=as.data.frame(read.table("../carol_analysis/results/output_1/abundance.tsv",header = T, sep = "\t"))
head(filterValues)

#List all the available datasets
biomaRt::listMarts()

#Biomart oneliner.
#Conect to ensembl and use the specific gene ensembl wanted.
ensembl=useEnsembl(biomart = "ensembl", dataset = "drerio_gene_ensembl")
head(listAttributes(ensembl))


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
all_files

#Create the data groups separately.
wt_phenotype<-c("S0WT_DD_2","S0WT_DD_3","S0WT_DD_4","S0WT_S_6","S0WT_S_7","S0WT_S_B","S0WT_S_C")

#Orden de los datos no altera.KALLISTo en este caso.
#IMPORT COUNTS WITH GENE NAMES


wt_phenotype<-sort(wt_phenotype,decreasing = F)
wt_phenotype
wt_phenotype_all_files<-all_files[wt_phenotype]
wt_phenotype_all_files
txi_wt_phenotype<-tximport(files=wt_phenotype_all_files,type="kallisto",tx2gene=tx2gene)

names(txi_wt_phenotype)
head(txi_wt_phenotype$counts)




#Create a folder for all the data related with this analysis.
getwd()

#Set the analysis results into the wt_phenotype folder.
setwd("/home/marius/Desktop/carol_data/carol_kallisto_deseq2/")
dir.create(recursive = T, "/home/marius/Desktop/carol_data/carol_kallisto_deseq2/wt_super_extended")
#Set the analysis results into the wtso_ddko_phenotype folder.
setwd("/home/marius/Desktop/carol_data/carol_kallisto_deseq2/wt_super_extended")



#Go for DESeq2 analysis.
#As DESeq2 uses raw counts the analysis can go straight froward since here!
head(txi_wt_phenotype$counts)

#Create the sampletable for later comparissons. Creating the cases of the condition.
sampleTable_wt_phenotype<-data.frame(treatment_type=rep("WT DD_vs_SD",7),CELLTYPE=factor(c(rep("WT_DD",3),rep("WT_SO",4))),
                                     row.names = colnames(txi_wt_phenotype$counts))


#Check for the sampleTable is ready.

#Control the sample table and the data are logically related
if (all(colnames(txi_wt_phenotype$counts)!=(rownames(sampleTable_wt_phenotype)))){
  cat("STOP\nSTOP\nSTOP\nDATA MUST BE REARRANGED\nFOR A CORRECT READING OF THE SAMPLE TABLE!")
  print("STOP\nSTOP\nSTOP\nDATA MUST BE REARRANGED\nFOR A CORRECT READING OF THE SAMPLE TABLE!")
}else{
  cat("All good, sampletable is correlated to the samples!")
  print("All good, sampletable is correlated to the samples!")
}
#Table is READY if results in TRUE.

#See and save the data
sampleTable_wt_phenotype


#Save the table
write.table(sampleTable_wt_phenotype,file = "./sampleTable_wt_phenotype.txt",sep = "\t",row.names = T, col.names = T)
#Create the DESeq2 DATA SET
dds_wt_phenotype<-DESeqDataSetFromTximport(txi_wt_phenotype,sampleTable_wt_phenotype, ~CELLTYPE)

dim(dds_wt_phenotype)
head(dds_wt_phenotype)
dds_wt_phenotype@colData
#Raw counts and average transcript length.
head(dds_wt_phenotype@assays$data$counts)


#Remove all the samples under 3 count of the ROWSUM.
dim(dds_wt_phenotype)
dds_wt_phenotype <- dds_wt_phenotype[ rowSums(dds_wt_phenotype@assays$data$counts) > 3, ]

#Removed rows with less than 3 counts for all.
dim(dds_wt_phenotype)
#Check the differences.




#DIfferential EXpression ANalysis
dds_wt_phenotype=DESeq(dds_wt_phenotype,betaPrior = T,parallel = TRUE)
dds_wt_phenotype@colData
head(dds_wt_phenotype$CELLTYPE)
#using pre-existing size factors
#estimating dispersions
#gene-wise dispersion estimates
#mean-dispersion relationship
#final dispersion estimates
#fitting model and testing
#See the dispersion of the data.
dir.create("./plots")
pdf("./plots/hist_counts_pdf.pdf",width = 12,height = 12)
plotDispEsts(dds_wt_phenotype)
dev.off()

#For columns colors.
#help("plotDispEsts")

#············································································
head(dds_wt_phenotype@assays$data$counts)
head(dds_wt_phenotype@assays$data$normalizationFactors)

#············································································
head(sampleTable_wt_phenotype)
#USING BH also known as FDR
#tidy=T the first column will be the rownames
#Results of the RNA-seq
#Specify which contrast to do, comparisson against.
dds_res_wt_phenotype=results(dds_wt_phenotype,contrast=c("CELLTYPE","WT_DD","WT_SO"), parallel=TRUE,pAdjustMethod = "BH", tidy=F)
head(dds_res_wt_phenotype)
dim(dds_res_wt_phenotype)

#Check what has been done.
#Only works if tidy is established to F.
mcols(dds_res_wt_phenotype)$description


#Annotate results.
#After Selecting the STATISTICAL THRESHOLDS we can start annotating.
#Sort the data frame annotated by the padj value and the -log2foldchange
dds_res_wt_phenotype<-dds_res_wt_phenotype[order(dds_res_wt_phenotype[,6], -dds_res_wt_phenotype[,3]),]
head(dds_res_wt_phenotype)
dim(dds_res_wt_phenotype)


#Start annotating the data with the Human and Mouse homologous for later dat obtention.

#Conect to ensembl and use the specific gene ensembl wanted.
ensembl=useEnsembl(biomart = "ensembl", dataset = "drerio_gene_ensembl")
head(listAttributes(ensembl))

#get the comboination of attributnames for ensembl filtered out by your filtertype.
ensembl_dds_res_wt_phenotype=getBM(attributes = c("ensembl_gene_id","description","hsapiens_homolog_associated_gene_name","hsapiens_homolog_ensembl_gene","mmusculus_homolog_associated_gene_name","mmusculus_homolog_ensembl_gene"),filters = "ensembl_gene_id",values =dds_res_wt_phenotype[,1] , mart = ensembl)#,'external_gene_name','hgnc_symbol','description'), 
head(ensembl_dds_res_wt_phenotype);dim(ensembl_dds_res_wt_phenotype)


dim(subset(dds_res_wt_phenotype,padj < .05))

#Merge the annotated data with the different species ENTREZIDs.
#Select the data from the following packages --> NCBI  after ensembl.


#####################################################
#                                                   #
#           ZEBRAFISH                               #
#                                                   #
#####################################################
cat("Drerio database ensembl annotation")
print("Drerio database ensembl annotation")
#ZEBRAFISH.
#Zebrafish annotating with ZEBRAFISH SPECIFIC ENTREZID
columns(org.Dr.eg.db)
anno_zebrafish_wt_phenotype<-AnnotationDbi::select(org.Dr.eg.db,
                                                          keys = as.character(ensembl_dds_res_wt_phenotype$ensembl_gene_id),
                                                          columns = c("ENTREZID","SYMBOL","GENENAME"),
                                                          keytype = "ENSEMBL") 

dim(anno_zebrafish_wt_phenotype)
head(anno_zebrafish_wt_phenotype)

#Merge the statistical information with the annotations realized

anno_zebrafish_wt_phenotype_stats<-merge(as.data.frame(dds_res_wt_phenotype),anno_zebrafish_wt_phenotype,by.x="row.names",by.y="ENSEMBL")
#Rename the ensemble gene id and wit the species.
colnames(anno_zebrafish_wt_phenotype_stats)[1]<-"drerio_ensembl_gene_id"
anno_zebrafish_wt_phenotype_stats<-anno_zebrafish_wt_phenotype_stats[-which(duplicated(anno_zebrafish_wt_phenotype_stats$drerio_ensembl_gene_id)),]

anno_zebrafish_wt_phenotype_stats<-anno_zebrafish_wt_phenotype_stats[,c(1,3,6,7,8,9,10)]
colnames(anno_zebrafish_wt_phenotype_stats)<-c("drerio_ensembl_gene_id","log2FoldChange","pvalue","padj","ENTREZID","SYMBOL","GENENAME")
dim(anno_zebrafish_wt_phenotype_stats)
head(anno_zebrafish_wt_phenotype_stats)

#####################################################
#                                                   #
#           HUMAN                                   #
#                                                   #
#####################################################
cat("Hsapiens database ensembl annotation")
print("Hsapiens database ensembl annotation")
#HUMAN
#Homo_sapiens
#get the Human orthologues ENTREZIDs
#USING Human specific ENTREZIDs
#NCBI package.
columns(org.Hs.eg.db)

head(ensembl_dds_res_wt_phenotype)
length(ensembl_dds_res_wt_phenotype$ensembl_gene_id)
length(ensembl_dds_res_wt_phenotype$hsapiens_homolog_associated_gene_name)
anno_human_wt_phenotype<-AnnotationDbi::select(org.Hs.eg.db,
                                                      keys=ensembl_dds_res_wt_phenotype$hsapiens_homolog_ensembl_gene,
                                                      columns = c("ENTREZID","SYMBOL","GENENAME"),
                                                      keytype = "ENSEMBL")



head(anno_human_wt_phenotype)
dim(anno_human_wt_phenotype)

#Control the number of duplicates..
#Keep the unique ENSEMBL only as not all the 
length(unique(anno_human_wt_phenotype$ENSEMBL))
length(unique(anno_human_wt_phenotype$ENTREZID))


anno_human_wt_phenotype<-anno_human_wt_phenotype[-which(duplicated(anno_human_wt_phenotype$ENSEMBL)),]
dim(anno_human_wt_phenotype)


#Merge the human data with the homologous data and select later the columns interesting.
anno_human_pre_stats_wt_phenotype<-merge(ensembl_dds_res_wt_phenotype[,c("ensembl_gene_id","hsapiens_homolog_ensembl_gene")],anno_human_wt_phenotype,by.x="hsapiens_homolog_ensembl_gene",by.y="ENSEMBL")
head(anno_human_pre_stats_wt_phenotype)
dim(anno_human_pre_stats_wt_phenotype)




#Merge the statistical information with the annotations realized
anno_human_wt_phenotype_stats<-merge(anno_human_pre_stats_wt_phenotype,as.data.frame(dds_res_wt_phenotype),by.y="row.names",by.x="ensembl_gene_id")
dim(anno_human_wt_phenotype_stats)
#head(dds_res_wt_phenotype)
head(anno_human_wt_phenotype_stats)
colnames(anno_human_wt_phenotype_stats)[2]<-"hsapiens_homolog_ensembl_gene_id"

anno_human_wt_phenotype_stats<-anno_human_wt_phenotype_stats[-which(duplicated(anno_human_wt_phenotype_stats$hsapiens_homolog_ensembl_gene_id)),]

#Create the same data frame as in zebrafish for easier working.
#anno_human_wt_phenotype_stats<-anno_human_wt_phenotype_stats[,c("hsapiens_homolog_ensembl_gene_id","log2FoldChange","pvalue","padj","ENTREZID","SYMBOL","GENENAME")]
anno_human_wt_phenotype_stats<-anno_human_wt_phenotype_stats[,c("hsapiens_homolog_ensembl_gene_id","log2FoldChange","pvalue","padj","ENTREZID","SYMBOL","GENENAME")]
#colnames(anno_human_wt_phenotype_stats)<-c("hsapiens_homolog_ensembl_gene_id","log2FoldChange","pvalue","padj","ENTREZID","SYMBOL","GENENAME")
head(anno_human_wt_phenotype_stats)
dim(anno_human_wt_phenotype_stats)



#####################################################
#                                                   #
#           MOUSE                                   #
#                                                   #
#####################################################
cat("Mmusculus database ensembl annotation")
print("Mmusculus database ensembl annotation")
#MOUSE
#Mus_musculus
#MOUSE entrez data 
#MOUSE specific ENTREZID data...
columns(org.Mm.eg.db)
anno_mouse_wt_phenotype<-AnnotationDbi::select(org.Mm.eg.db,
                                                      keys=ensembl_dds_res_wt_phenotype$mmusculus_homolog_ensembl_gene,
                                                      columns = c("ENTREZID","SYMBOL","GENENAME"),
                                                      keytype = "ENSEMBL")

head(anno_mouse_wt_phenotype)
dim(anno_mouse_wt_phenotype)

anno_mouse_wt_phenotype<-anno_mouse_wt_phenotype[-which(duplicated(anno_mouse_wt_phenotype$ENSEMBL)),]
dim(anno_mouse_wt_phenotype)


#Merge the human data with the homologous data and select later the columns interesting.
anno_mouse_pre_stats_wt_phenotype<-merge(ensembl_dds_res_wt_phenotype[,c("ensembl_gene_id","mmusculus_homolog_ensembl_gene")],anno_mouse_wt_phenotype,by.x="mmusculus_homolog_ensembl_gene",by.y="ENSEMBL")
head(anno_mouse_pre_stats_wt_phenotype,10)
dim(anno_mouse_pre_stats_wt_phenotype)



#Merge the statistical information with the annotations realized
anno_mouse_wt_phenotype_stats<-merge(as.data.frame(dds_res_wt_phenotype),anno_mouse_pre_stats_wt_phenotype,by.x="row.names",by.y="ensembl_gene_id")
dim(anno_mouse_wt_phenotype_stats)
head(anno_mouse_wt_phenotype_stats)
colnames(anno_mouse_wt_phenotype_stats)[8]<-"mmusculus_homolog_ensembl_gene_id"
#Create the same data frame as in zebrafish for easier working.
anno_mouse_wt_phenotype_stats<-anno_mouse_wt_phenotype_stats[-which(duplicated(anno_mouse_wt_phenotype_stats$mmusculus_homolog_ensembl_gene_id)),]
anno_mouse_wt_phenotype_stats<-anno_mouse_wt_phenotype_stats[,c("mmusculus_homolog_ensembl_gene_id","log2FoldChange","pvalue","padj","ENTREZID","SYMBOL","GENENAME")]
#colnames(anno_mouse_wt_phenotype_stats)<-c("mmusculus_homolog_ensembl_gene","log2FoldChange","pvalue","padj","ENTREZID","SYMBOL","GENENAME")
head(anno_mouse_wt_phenotype_stats)

#print the samples to verify
head(anno_human_wt_phenotype_stats);dim(anno_human_wt_phenotype_stats)
head(anno_mouse_wt_phenotype_stats);dim(anno_mouse_wt_phenotype_stats)
head(anno_zebrafish_wt_phenotype_stats);dim(anno_zebrafish_wt_phenotype_stats)


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

#Add the zscore for the pathways and select the top expressed for up and down pathways, retruning a topUPandDOWN object, ready for plotting.
top_paths_FDR<-function(FDR_pathways,specie_anno_stats,n=15){
  print(ifelse(is.null(specie_anno_stats),break,"Specie selected!!"))
  print(ifelse(is.null(n),"Number of up and down rows not selected,using default!!",paste("n=",n)))
  if(length(FDR_pathways@result$ID)==0){
    warning("ERROR → The pathway enrichment was not returning informtion with these parameters!!!\n")
    cat("ERROR → The pathway enrichment was not returning informtion with these parameters!!!\n" )
    break
  }else{
    cat("\nAdding zscore values to the tables...\n")
    #Create the zscore new column.
    symbols<-list()
    symbols_extracted<-list()
    symbols_extracted_double<-list()
    for (i in 1:length(FDR_pathways@result$ID)){
      #Extract from the table all the SYMBOLS involved in the pathways from the gnees.
      target<-paste0("target",i)
      #print(target)
      
      #Clean the genelist of each row and obtain the genes.
      symbols[[target]]<-unique(unlist(strsplit(as.character(FDR_pathways@result$geneID[i]),split = "\\/")))
      symbols_extracted[[target]]<-specie_anno_stats[specie_anno_stats$SYMBOL%in%symbols[[i]],]
      
      #Select only the SYMBOLS and the logFC for calculating the zscore.
      symbols_extracted_double[[target]]<-symbols_extracted[[i]][,c("log2FoldChange","SYMBOL")]
      
      #Get the zscore added
      FDR_pathways@result$zscore[i]<-sum(as.numeric(symbols_extracted_double[[i]][,1]))/sqrt(length(symbols[[i]]))}
    cat("....\nSelecting top UP and DOWN values...\n......\n")
    topUP<-FDR_pathways@result%>%
      filter(zscore>0)%>%
      top_n(n,wt = -p.adjust)
    topDOWN<-FDR_pathways@result%>%
      filter(zscore<0)%>%
      top_n(n,wt = -p.adjust)
    topPATHS<-bind_rows(topUP,topDOWN)
    cat("\nCompleted!")
    return(topPATHS)}
}

#a<-as.data.frame(anno_zebrafish_wt_phenotype_stats[anno_zebrafish_wt_phenotype_stats$SYMBOL=="cldn7a",])
#a<-tidyr::drop_na(a)
#a

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
cat("FDR GO Drerio ENRICHMENTS")
print("FDR GO Drerio ENRICHMENTS")
#Start with Zebrafish annotations with zebrafish
dir.create("./Enrichment/FDR/GO/DanioRerio/")
DrOrgDb<-org.Dr.eg.db


#Get the differential expressed genes, and use only the unique seuqences no repeated..
dr_wt_phenotype_diff_expressed_genes<-unique(anno_zebrafish_wt_phenotype_stats$ENTREZID[which(anno_zebrafish_wt_phenotype_stats$padj<.05)])


length(dr_wt_phenotype_diff_expressed_genes)
head(dr_wt_phenotype_diff_expressed_genes)


#Biologicall process enriched for zebrafish 
dr_wt_phenotype_BP_genelist_enrichGO <- enrichGO(gene = dr_wt_phenotype_diff_expressed_genes ,OrgDb  =DrOrgDb, ont = "BP", pAdjustMethod = "BH",pvalueCutoff  = 0.05,readable = T)


#Molecular Function enriched for zebrafish
dr_wt_phenotype_MF_genelist_enrichGO <- enrichGO(gene = dr_wt_phenotype_diff_expressed_genes ,OrgDb  = DrOrgDb, ont = "MF", pAdjustMethod = "BH",pvalueCutoff  = 0.05,readable = T)

#Cellular component enriched zebrafish
dr_wt_phenotype_CC_genelist_enrichGO <- enrichGO(gene = dr_wt_phenotype_diff_expressed_genes ,OrgDb  = DrOrgDb, ont = "CC", pAdjustMethod = "BH",pvalueCutoff  = 0.05,readable = T)

head(dr_wt_phenotype_CC_genelist_enrichGO@result)
#Save the lists for Zebrafish enrichment.
#Save better all the results,head(dr_wt_phenotype_BP_genelist_enrichGO@result)
#Save BP list as dataframe.
write.table(x = data.frame(dr_wt_phenotype_BP_genelist_enrichGO@result),"./Enrichment/FDR/GO/DanioRerio/dr_wt_phenotype_BP_enrichGO.txt",sep="\t", quote = F, row.names = F, col.names = T)


#Save the MF list as dataframe
write.table(x = data.frame(dr_wt_phenotype_MF_genelist_enrichGO@result),"./Enrichment/FDR/GO/DanioRerio/dr_wt_phenotype_MF_enrichGO.txt",sep="\t", quote = F, row.names = F, col.names = T)

#Save the CC list as a dataframe
write.table(x = data.frame(dr_wt_phenotype_CC_genelist_enrichGO@result),"./Enrichment/FDR/GO/DanioRerio/dr_wt_phenotype_CC_enrichGO.txt",sep="\t", quote = F, row.names = F, col.names = T)



head(dr_wt_phenotype_BP_genelist_enrichGO@result,20)



#Visualize the data obtained for Danio rerio Enrichment.

#Plot for the BProcess data obtained from the enrichment analysis in GO fpr Danio rerio.

#Biological processes plots.
#BARPLOT
pdf("./Enrichment/FDR/GO/DanioRerio/dr_wt_phenotype_BP_barplot.pdf", height = 20, width = 16)
barplot(dr_wt_phenotype_BP_genelist_enrichGO,drop=T,showCategory = 25,title = "dr BP enrichment wt_phenotype \n genelist Selected by padj<0.05 ")
dev.off()

#DOTPLOT of BP in Dr
pdf("./Enrichment/FDR/GO/DanioRerio/dr_wt_phenotype_BP_dotplot.pdf",height = 20,width = 16)
dotplot(dr_wt_phenotype_BP_genelist_enrichGO, showCategory=25)
dev.off()


#Molecular Functions plots.
#BARPLOT
pdf("./Enrichment/FDR/GO/DanioRerio/dr_wt_phenotype_MF_barplot.pdf", height = 20, width = 16)
barplot(dr_wt_phenotype_MF_genelist_enrichGO,drop=T,showCategory = 25,title = "dr MF enrichment wt_phenotype \n genelist Selected by padj<0.05 ")
dev.off()

#DOTPLOT of MF in Danio rerio
pdf("./Enrichment/FDR/GO/DanioRerio/dr_wt_phenotype_MF_dotplot.pdf",height = 20,width = 16)
dotplot(dr_wt_phenotype_MF_genelist_enrichGO, showCategory=25)
dev.off()


#Celullar Components plots.
#BARPLOT
pdf("./Enrichment/FDR/GO/DanioRerio/dr_wt_phenotype_CC_barplot.pdf", height = 20, width = 16)
barplot(dr_wt_phenotype_CC_genelist_enrichGO,drop=T,showCategory = 25,title = "dr CC enrichment wt_phenotype \n genelist Selected by padj<0.05 ")
dev.off()

#DOTPLOT of CC in Danio rerio
pdf("./Enrichment/FDR/GO/DanioRerio/dr_wt_phenotype_CC_dotplot.pdf",height = 20,width = 16)
dotplot(dr_wt_phenotype_CC_genelist_enrichGO, showCategory=25)
dev.off()






#################################
#                               #
#           ZSCORE PLOT         #
#                               #
#################################
#Add the zscore and select the top 15 of each enriched functions in FDR for each GO enrichment.
#GO BP

dr_wt_phenotype_BP_genelist_enrichGO_topPathways<-top_paths_FDR(dr_wt_phenotype_BP_genelist_enrichGO,anno_zebrafish_wt_phenotype_stats,20)
head(dr_wt_phenotype_BP_genelist_enrichGO_topPathways[,c(2:7,10)]);tail(dr_wt_phenotype_BP_genelist_enrichGO_topPathways[,c(2:7,10)])


pdf("./Enrichment/FDR/GO/DanioRerio/dr_wt_phenotype_BP_genelist_enrichGO_topPathways.pdf",height = 12,width = 18)
ggplot(dr_wt_phenotype_BP_genelist_enrichGO_topPathways,aes(x=reorder(Description,zscore),y=zscore))+
  geom_bar(stat = "identity",fill="darkorchid")+
  #geom_text(aes(label= round(p.adjust,3)),hjust=1.2,vjust=-.00001,size=3.5,color="white")+
  #labs(text="padjust")+
  geom_point(aes(label=p.adjust,size=rev(p.adjust)))+
  coord_flip()+
  labs(title = "Dr FDR GO:BP",subtitle = sampleTable_wt_phenotype$treatment_type,size ="p.adjusted",x="Biological Processes")

dev.off()
#Save table 
write.table(dr_wt_phenotype_BP_genelist_enrichGO_topPathways,"./Enrichment/FDR/GO/DanioRerio/dr_wt_phenotype_BP_genelist_enrichGO_topPathways.tsv",row.names = F,col.names = T,sep = "\t")


#GO MF

dr_wt_phenotype_MF_genelist_enrichGO_topPathways<-top_paths_FDR(dr_wt_phenotype_MF_genelist_enrichGO,anno_zebrafish_wt_phenotype_stats,20)
head(dr_wt_phenotype_MF_genelist_enrichGO_topPathways[,c(2:7,10)]);tail(dr_wt_phenotype_MF_genelist_enrichGO_topPathways[,c(2:7,10)])


pdf("./Enrichment/FDR/GO/DanioRerio/dr_wt_phenotype_MF_genelist_enrichGO_topPathways.pdf",height = 12,width = 18)
ggplot(dr_wt_phenotype_MF_genelist_enrichGO_topPathways,aes(x=reorder(Description,zscore),y=zscore))+
  geom_bar(stat = "identity",fill="darkorchid")+
  #geom_text(aes(label= round(p.adjust,3)),hjust=1.2,vjust=-.00001,size=3.5,color="white")+
  #labs(text="padjust")+
  geom_point(aes(label=p.adjust,size=rev(p.adjust)))+
  coord_flip()+
  labs(title = "Dr FDR GO:MF",subtitle = sampleTable_wt_phenotype$treatment_type,size ="p.adjusted",x="Molecular Function")

dev.off()
#Save table 
write.table(dr_wt_phenotype_MF_genelist_enrichGO_topPathways,"./Enrichment/FDR/GO/DanioRerio/dr_wt_phenotype_MF_genelist_enrichGO_topPathways.tsv",row.names = F,col.names = T,sep = "\t")


#GO CC
dr_wt_phenotype_CC_genelist_enrichGO_topPathways<-top_paths_FDR(dr_wt_phenotype_CC_genelist_enrichGO,anno_zebrafish_wt_phenotype_stats,20)
head(dr_wt_phenotype_CC_genelist_enrichGO_topPathways[,c(2:7,10)]);tail(dr_wt_phenotype_CC_genelist_enrichGO_topPathways[,c(2:7,10)])


pdf("./Enrichment/FDR/GO/DanioRerio/dr_wt_phenotype_CC_genelist_enrichGO_topPathways.pdf",height = 12,width = 18)
ggplot(dr_wt_phenotype_CC_genelist_enrichGO_topPathways,aes(x=reorder(Description,zscore),y=zscore))+
  geom_bar(stat = "identity",fill="darkorchid")+
  #geom_text(aes(label= round(p.adjust,3)),hjust=1.2,vjust=-.00001,size=3.5,color="white")+
  #labs(text="padjust")+
  geom_point(aes(label=p.adjust,size=rev(p.adjust)))+
  coord_flip()+
  labs(title = "Dr FDR GO:CC",subtitle = sampleTable_wt_phenotype$treatment_type,size ="p.adjusted",x="Cellular Component")

dev.off()
#Save table 
write.table(dr_wt_phenotype_CC_genelist_enrichGO_topPathways,"./Enrichment/FDR/GO/DanioRerio/dr_wt_phenotype_CC_genelist_enrichGO_topPathways.tsv",row.names = F,col.names = T,sep = "\t")












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
cat("FDR GO Hsapiens ENRICHMENTS")
print("FDR GO Hsapiens ENRICHMENTS")


#Second continue with human
dir.create("./Enrichment/FDR/GO/Homosapiens/")
HsOrgDb<-org.Hs.eg.db


#Get the differential expressed genes
hs_wt_phenotype_diff_expressed_genes<-unique(anno_human_wt_phenotype_stats$ENTREZID[which(anno_human_wt_phenotype_stats$padj<.05)])
length(hs_wt_phenotype_diff_expressed_genes)


#Biologicall process enriched for human 
hs_wt_phenotype_BP_genelist_enrichGO <- enrichGO(gene = hs_wt_phenotype_diff_expressed_genes ,OrgDb  =HsOrgDb, ont = "BP", pAdjustMethod = "BH",pvalueCutoff  = 0.05,readable = T)


#Molecular Function enriched for human
hs_wt_phenotype_MF_genelist_enrichGO <- enrichGO(gene = hs_wt_phenotype_diff_expressed_genes ,OrgDb  = HsOrgDb, ont = "MF", pAdjustMethod = "BH",pvalueCutoff  = 0.05,readable = T)

#Cellular component enriched human
hs_wt_phenotype_CC_genelist_enrichGO <- enrichGO(gene = hs_wt_phenotype_diff_expressed_genes ,OrgDb  = HsOrgDb, ont = "CC", pAdjustMethod = "BH",pvalueCutoff  = 0.05,readable = T)



#Save the lists for Human enrichment.
#Save better all the results,head(hs_wt_phenotype_BP_genelist_enrichGO@result)
#Save BP list as dataframe.
write.table(x = data.frame(hs_wt_phenotype_BP_genelist_enrichGO@result),"./Enrichment/FDR/GO/Homosapiens/hs_wt_phenotype_BP_enrichGO.txt",sep="\t", quote = F, row.names = F, col.names = T)


#Save the MF list as dataframe
write.table(x = data.frame(hs_wt_phenotype_MF_genelist_enrichGO@result),"./Enrichment/FDR/GO/Homosapiens/hs_wt_phenotype_MF_enrichGO.txt",sep="\t", quote = F, row.names = F, col.names = T)

#Save the CC list as a dataframe
write.table(x = data.frame(hs_wt_phenotype_CC_genelist_enrichGO@result),"./Enrichment/FDR/GO/Homosapiens/hs_wt_phenotype_CC_enrichGO.txt",sep="\t", quote = F, row.names = F, col.names = T)

head(hs_wt_phenotype_BP_genelist_enrichGO@result)



#Visualize the data obtained for the human.
#Plot for the BProcess data obtained from the enrichment analysis in GO for Human in homologues for wt_phenotype.
#Biological processes plots.
#BARPLOT
pdf("./Enrichment/FDR/GO/Homosapiens/hs_wt_phenotype_BP_barplot.pdf", height = 20, width = 16)
barplot(hs_wt_phenotype_BP_genelist_enrichGO,drop=T,showCategory = 25,title = "hs BP enrichment wt_phenotype \n genelist Selected by padj<0.05 ")
dev.off()

#DOTPLOT of BP in Hs
pdf("./Enrichment/FDR/GO/Homosapiens/hs_wt_phenotype_BP_dotplot.pdf",height = 20,width = 16)
dotplot(hs_wt_phenotype_BP_genelist_enrichGO, showCategory=25)
dev.off()


#Molecular Functions plots.
#BARPLOT
pdf("./Enrichment/FDR/GO/Homosapiens/hs_wt_phenotype_MF_barplot.pdf", height = 20, width = 16)
barplot(hs_wt_phenotype_MF_genelist_enrichGO,drop=T,showCategory = 25,title = "hs MF enrichment wt_phenotype \n genelist Selected by padj<0.05 ")
dev.off()

#DOTPLOT of MF in Human
pdf("./Enrichment/FDR/GO/Homosapiens/hs_wt_phenotype_MF_dotplot.pdf",height = 20,width = 16)
dotplot(hs_wt_phenotype_MF_genelist_enrichGO, showCategory=25)
dev.off()


#Celullar Components plots.
#BARPLOT
pdf("./Enrichment/FDR/GO/Homosapiens/hs_wt_phenotype_CC_barplot.pdf", height = 20, width = 16)
barplot(hs_wt_phenotype_CC_genelist_enrichGO,drop=T,showCategory = 25,title = "hs CC enrichment wt_phenotype \n genelist Selected by padj<0.05 ")
dev.off()

#DOTPLOT of CC in Human
pdf("./Enrichment/FDR/GO/Homosapiens/hs_wt_phenotype_CC_dotplot.pdf",height = 20,width = 16)
dotplot(hs_wt_phenotype_CC_genelist_enrichGO, showCategory=25)
dev.off()





#################################
#                               #
#           ZSCORE PLOT         #
#                               #
#################################


#Add the zscore and select the top 15 of each enriched functions in FDR for each GO enrichment.
#GO BP

hs_wt_phenotype_BP_genelist_enrichGO_topPathways<-top_paths_FDR(hs_wt_phenotype_BP_genelist_enrichGO,anno_human_wt_phenotype_stats,20)
head(hs_wt_phenotype_BP_genelist_enrichGO_topPathways[,c(2:7,10)]);tail(hs_wt_phenotype_BP_genelist_enrichGO_topPathways[,c(2:7,10)])


pdf("./Enrichment/FDR/GO/Homosapiens/hs_wt_phenotype_BP_genelist_enrichGO_topPathways.pdf",height = 12,width = 18)
ggplot(hs_wt_phenotype_BP_genelist_enrichGO_topPathways,aes(x=reorder(Description,zscore),y=zscore))+
  geom_bar(stat = "identity",fill="darkorchid")+
  #geom_text(aes(label= round(p.adjust,3)),hjust=1.2,vjust=-.00001,size=3.5,color="white")+
  #labs(text="padjust")+
  geom_point(aes(label=p.adjust,size=rev(p.adjust)))+
  coord_flip()+
  labs(title = "Hs FDR GO:BP",subtitle = sampleTable_wt_phenotype$treatment_type,size ="p.adjusted",x="Biological Processes")

dev.off()

#Save table 
write.table(hs_wt_phenotype_BP_genelist_enrichGO_topPathways,"./Enrichment/FDR/GO/Homosapiens/hs_wt_phenotype_BP_genelist_enrichGO_topPathways.tsv",row.names = F,col.names = T,sep = "\t")


#GO MF

hs_wt_phenotype_MF_genelist_enrichGO_topPathways<-top_paths_FDR(hs_wt_phenotype_MF_genelist_enrichGO,anno_human_wt_phenotype_stats,20)
head(hs_wt_phenotype_MF_genelist_enrichGO_topPathways[,c(2:7,10)]);tail(hs_wt_phenotype_MF_genelist_enrichGO_topPathways[,c(2:7,10)])


pdf("./Enrichment/FDR/GO/Homosapiens/hs_wt_phenotype_MF_genelist_enrichGO_topPathways.pdf",height = 12,width = 18)
ggplot(hs_wt_phenotype_MF_genelist_enrichGO_topPathways,aes(x=reorder(Description,zscore),y=zscore))+
  geom_bar(stat = "identity",fill="darkorchid")+
  #geom_text(aes(label= round(p.adjust,3)),hjust=1.2,vjust=-.00001,size=3.5,color="white")+
  #labs(text="padjust")+
  geom_point(aes(label=p.adjust,size=rev(p.adjust)))+
  coord_flip()+
  labs(title = "Hs FDR GO:MF",subtitle = sampleTable_wt_phenotype$treatment_type,size ="p.adjusted",x="Molecular Function")

dev.off()

#Save table 
write.table(hs_wt_phenotype_MF_genelist_enrichGO_topPathways,"./Enrichment/FDR/GO/Homosapiens/hs_wt_phenotype_MF_genelist_enrichGO_topPathways.tsv",row.names = F,col.names = T,sep = "\t")



#GO CC
hs_wt_phenotype_CC_genelist_enrichGO_topPathways<-top_paths_FDR(hs_wt_phenotype_CC_genelist_enrichGO,anno_human_wt_phenotype_stats,20)
head(hs_wt_phenotype_CC_genelist_enrichGO_topPathways[,c(2:7,10)]);tail(hs_wt_phenotype_CC_genelist_enrichGO_topPathways[,c(2:7,10)])


pdf("./Enrichment/FDR/GO/Homosapiens/hs_wt_phenotype_CC_genelist_enrichGO_topPathways.pdf",height = 12,width = 18)
ggplot(hs_wt_phenotype_CC_genelist_enrichGO_topPathways,aes(x=reorder(Description,zscore),y=zscore))+
  geom_bar(stat = "identity",fill="darkorchid")+
  #geom_text(aes(label= round(p.adjust,3)),hjust=1.2,vjust=-.00001,size=3.5,color="white")+
  #labs(text="padjust")+
  geom_point(aes(label=p.adjust,size=rev(p.adjust)))+
  coord_flip()+
  labs(title = "Hs FDR GO:CC",subtitle = sampleTable_wt_phenotype$treatment_type,size ="p.adjusted",x="Cellular Component")

dev.off()
#Save table 
write.table(hs_wt_phenotype_CC_genelist_enrichGO_topPathways,"./Enrichment/FDR/GO/Homosapiens/hs_wt_phenotype_CC_genelist_enrichGO_topPathways.tsv",row.names = F,col.names = T,sep = "\t")




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
cat("FDR GO Mmusculus ENRICHMENTS")
print("FDR GO Mmusculus ENRICHMENTS")

#Do the same analysis for Mouse

MmOrgDb<- org.Mm.eg.db
#Using Human Annotation Package
dir.create("./Enrichment/FDR/GO/Mmusculus/")



#Get the differential expressed genes
mm_wt_phenotype_diff_expressed_genes<-unique(anno_mouse_wt_phenotype_stats$ENTREZID[which(anno_mouse_wt_phenotype_stats$padj<.05)])
length(mm_wt_phenotype_diff_expressed_genes)

#Biologicall process enriched for Mouse 
mm_wt_phenotype_BP_genelist_enrichGO <- enrichGO(gene = mm_wt_phenotype_diff_expressed_genes ,OrgDb  =MmOrgDb, ont = "BP", pAdjustMethod = "BH",pvalueCutoff  = 0.05,readable = T)


#Molecular Function enriched for mouse
mm_wt_phenotype_MF_genelist_enrichGO <- enrichGO(gene = mm_wt_phenotype_diff_expressed_genes ,OrgDb  = MmOrgDb, ont = "MF", pAdjustMethod = "BH",pvalueCutoff  = 0.05,readable = T)

#Cellular component enriched mouse
mm_wt_phenotype_CC_genelist_enrichGO <- enrichGO(gene = mm_wt_phenotype_diff_expressed_genes ,OrgDb  = MmOrgDb, ont = "CC", pAdjustMethod = "BH",pvalueCutoff  = 0.05,readable = T)



#Save the lists for mouse homoogues enrichments.
#Save better all the results,head(mm_wt_phenotype_BP_genelist_enrichGO@result)
#Save BP list as dataframe.
write.table(x = data.frame(mm_wt_phenotype_BP_genelist_enrichGO@result),"./Enrichment/FDR/GO/Mmusculus/mm_wt_phenotype_BP_enrichGO.txt",sep="\t", quote = F, row.names = F, col.names = T)


#Save the MF list as dataframe
write.table(x = data.frame(mm_wt_phenotype_MF_genelist_enrichGO@result),"./Enrichment/FDR/GO/Mmusculus/mm_wt_phenotype_MF_enrichGO.txt",sep="\t", quote = F, row.names = F, col.names = T)

#Save the CC list as a dataframe
write.table(x = data.frame(mm_wt_phenotype_CC_genelist_enrichGO@result),"./Enrichment/FDR/GO/Mmusculus/mm_wt_phenotype_CC_enrichGO.txt",sep="\t", quote = F, row.names = F, col.names = T)


head(mm_wt_phenotype_BP_genelist_enrichGO@result)

#Visualize the data obtained.

#Plot for the BProcess data obtained from the enrichment analysis in GO fpr Mouse homologues wt_phenotype.

#Biological processes plots.
#BARPLOT
pdf("./Enrichment/FDR/GO/Mmusculus/mm_wt_phenotype_BP_barplot.pdf", height = 20, width = 16)
barplot(mm_wt_phenotype_BP_genelist_enrichGO,drop=T,showCategory = 12,title = "mm BP enrichment wt_phenotype \n genelist Selected by padj<0.05 ")
dev.off()
#DOTPLOT of BP
pdf("./Enrichment/FDR/GO/Mmusculus/mm_wt_phenotype_BP_dotplot.pdf",height = 20,width = 16)
dotplot(mm_wt_phenotype_BP_genelist_enrichGO, showCategory=12)
dev.off()


#Molecular Functions plots.
#BARPLOT
pdf("./Enrichment/FDR/GO/Mmusculus/mm_wt_phenotype_MF_barplot.pdf", height = 20, width = 16)
barplot(mm_wt_phenotype_MF_genelist_enrichGO,drop=T,showCategory = 12,title = "mm MF enrichment wt_phenotype \n genelist Selected by padj<0.05 ")
dev.off()

#DOTPLOT of MF in mus musculus
pdf("./Enrichment/FDR/GO/Mmusculus/mm_wt_phenotype_MF_dotplot.pdf",height = 20,width = 16)
dotplot(mm_wt_phenotype_MF_genelist_enrichGO, showCategory=12)
dev.off()


#Celullar Components plots.
#BARPLOT
pdf("./Enrichment/FDR/GO/Mmusculus/mm_wt_phenotype_CC_barplot.pdf", height = 20, width = 16)
barplot(mm_wt_phenotype_CC_genelist_enrichGO,drop=T,showCategory = 12,title = "mm CC enrichment wt_phenotype \n genelist Selected by padj<0.05 ")
dev.off()

#DOTPLOT of CC in Danio rerio
pdf("./Enrichment/FDR/GO/Mmusculus/mm_wt_phenotype_CC_dotplot.pdf",height = 20,width = 16)
dotplot(mm_wt_phenotype_CC_genelist_enrichGO, showCategory=12)
dev.off()




#################################
#                               #
#           ZSCORE PLOT         #
#                               #
#################################
#Add the zscore and select the top 15 of each enriched functions in FDR for each GO enrichment.


#GO BP
head(mm_wt_phenotype_BP_genelist_enrichGO@result)
dim(mm_wt_phenotype_BP_genelist_enrichGO@result)
head(mm_wt_phenotype_BP_genelist_enrichGO_topPathways)
tail(mm_wt_phenotype_BP_genelist_enrichGO_topPathways)
dim(mm_wt_phenotype_BP_genelist_enrichGO_topPathways)

mm_wt_phenotype_BP_genelist_enrichGO_topPathways<-top_paths_FDR(mm_wt_phenotype_BP_genelist_enrichGO,anno_mouse_wt_phenotype_stats,20)
head(mm_wt_phenotype_BP_genelist_enrichGO_topPathways[,c(2:7,10)]);tail(mm_wt_phenotype_BP_genelist_enrichGO_topPathways[,c(2:7,10)])


pdf("./Enrichment/FDR/GO/Mmusculus/mm_wt_phenotype_BP_genelist_enrichGO_topPathways.pdf",height = 12,width = 18)
ggplot(mm_wt_phenotype_BP_genelist_enrichGO_topPathways,aes(x=reorder(Description,zscore),y=zscore))+
  geom_bar(stat = "identity",fill="darkorchid")+
  #geom_text(aes(label= round(p.adjust,3)),hjust=1.2,vjust=-.00001,size=3.5,color="white")+
  #labs(text="padjust")+
  geom_point(aes(label=p.adjust,size=rev(p.adjust)))+
  coord_flip()+
  labs(title = "Mm FDR GO:BP",subtitle = sampleTable_wt_phenotype$treatment_type,size ="p.adjusted",x="Biological Processes")

dev.off()

#Save table 
write.table(mm_wt_phenotype_BP_genelist_enrichGO_topPathways,"./Enrichment/FDR/GO/Mmusculus/mm_wt_phenotype_BP_genelist_enrichGO_topPathways.tsv",row.names = F,col.names = T,sep = "\t")


#GO MF
mm_wt_phenotype_MF_genelist_enrichGO_topPathways<-top_paths_FDR(mm_wt_phenotype_MF_genelist_enrichGO,anno_mouse_wt_phenotype_stats,20)
head(mm_wt_phenotype_MF_genelist_enrichGO_topPathways[,c(2:7,10)]);tail(mm_wt_phenotype_MF_genelist_enrichGO_topPathways[,c(2:7,10)])


pdf("./Enrichment/FDR/GO/Mmusculus/mm_wt_phenotype_MF_genelist_enrichGO_topPathways.pdf",height = 12,width = 18)
ggplot(mm_wt_phenotype_MF_genelist_enrichGO_topPathways,aes(x=reorder(Description,zscore),y=zscore))+
  geom_bar(stat = "identity",fill="darkorchid")+
  #geom_text(aes(label= round(p.adjust,3)),hjust=1.2,vjust=-.00001,size=3.5,color="white")+
  #labs(text="padjust")+
  geom_point(aes(label=p.adjust,size=rev(p.adjust)))+
  coord_flip()+
  labs(title = "Mm FDR GO:MF",subtitle = sampleTable_wt_phenotype$treatment_type,size ="p.adjusted",x="Molecular Function")

dev.off()

#Save table 
write.table(mm_wt_phenotype_MF_genelist_enrichGO_topPathways,"./Enrichment/FDR/GO/Mmusculus/mm_wt_phenotype_MF_genelist_enrichGO_topPathways.tsv",row.names = F,col.names = T,sep = "\t")



#GO CC
mm_wt_phenotype_CC_genelist_enrichGO_topPathways<-top_paths_FDR(mm_wt_phenotype_CC_genelist_enrichGO,anno_mouse_wt_phenotype_stats,20)
head(mm_wt_phenotype_CC_genelist_enrichGO_topPathways[,c(2:7,10)]);tail(mm_wt_phenotype_CC_genelist_enrichGO_topPathways[,c(2:7,10)])


pdf("./Enrichment/FDR/GO/Mmusculus/mm_wt_phenotype_CC_genelist_enrichGO_topPathways.pdf",height = 12,width = 18)
ggplot(mm_wt_phenotype_CC_genelist_enrichGO_topPathways,aes(x=reorder(Description,zscore),y=zscore))+
  geom_bar(stat = "identity",fill="darkorchid")+
  #geom_text(aes(label= round(p.adjust,3)),hjust=1.2,vjust=-.00001,size=3.5,color="white")+
  #labs(text="padjust")+
  geom_point(aes(label=p.adjust,size=rev(p.adjust)))+
  coord_flip()+
  labs(title = "Mm FDR GO:CC",subtitle = sampleTable_wt_phenotype$treatment_type,size ="p.adjusted",x="Cellular Componnent")

dev.off()


#Save table 
write.table(mm_wt_phenotype_CC_genelist_enrichGO_topPathways,"./Enrichment/FDR/GO/Mmusculus/mm_wt_phenotype_CC_genelist_enrichGO_topPathways.tsv",row.names = F,col.names = T,sep = "\t")


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
cat("FDR KEGG Drerio ENRICHMENTS")
print("FDR KEGG Drerio ENRICHMENTS")

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
dr_wt_phenotype_genelist_enrichKEGG<-enrichKEGG(gene = dr_wt_phenotype_diff_expressed_genes,organism = "dre",keyType = "kegg",pvalueCutoff = 0.05,pAdjustMethod = "BH",minGSSize = 10)
dr_wt_phenotype_genelist_enrichKEGG@result[,8]
head(dr_wt_phenotype_genelist_enrichKEGG@result)

#Translate the ENTREZID
dr_wt_phenotype_genelist_enrichKEGG<-setReadable(dr_wt_phenotype_genelist_enrichKEGG,OrgDb = org.Dr.eg.db,keyType = "ENTREZID")


#Write the results of the table enriched with the KEGG significant data.
write.table(data.frame(dr_wt_phenotype_genelist_enrichKEGG@result),"./Enrichment/FDR/KEGG/Drerio/dr_wt_phenotype_enrichmentKEGG.txt",sep = "\t",quote = F, row.names = F,col.names = T)


head(dr_wt_phenotype_genelist_enrichKEGG@result)
#Plot the decided figures as required.
#Those figures got downloaded in the working directory
dre04068<-pathview(gene.data = dr_wt_phenotype_diff_expressed_genes,
                   pathway.id = "dre04068",
                   species="dre",
                   limit = list(gene=5,cpd=1))


#################################
#                               #
#           ZSCORE PLOT         #
#                               #
#################################

#Plot the pathways.
#Call the functions and plot it.
dr_wt_phenotype_genelist_enrichKEGG_topPathWays<-top_paths_FDR(dr_wt_phenotype_genelist_enrichKEGG,anno_zebrafish_wt_phenotype_stats,20)
head(dr_wt_phenotype_genelist_enrichKEGG_topPathWays[,c(2,5,6,10)]);tail(dr_wt_phenotype_genelist_enrichKEGG_topPathWays[,c(2,5,6,10)])


pdf("./Enrichment/FDR/KEGG/Drerio/dr_wt_phenotype_genelist_enrichKEGG_topPathWays.pdf",height = 12,width = 18)
ggplot(dr_wt_phenotype_genelist_enrichKEGG_topPathWays,aes(x=reorder(Description,zscore),y=zscore))+
  geom_bar(stat = "identity",fill="orange")+
  #geom_text(aes(label= round(p.adjust,3)),hjust=1.2,vjust=-.00001,size=3.5,color="white")+
  #labs(text="padjust")+
  geom_point(aes(label=p.adjust,size=rev(p.adjust)))+
  coord_flip()+
  labs(title = "Dr FDR KEGG",subtitle = sampleTable_wt_phenotype$treatment_type,size ="p.adjusted",x="Pathways")

dev.off()

write.table(dr_wt_phenotype_genelist_enrichKEGG_topPathWays,"./Enrichment/FDR/KEGG/Drerio/dr_wt_phenotype_genelist_enrichKEGG_topPathWays.tsv",sep = "\t",row.names = F,col.names = T)


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
cat("FDR KEGG Hsapiens ENRICHMENTS")
print("FDR KEGG Hsapiens ENRICHMENTS")

###HUMAN KEGG

#Human Enrichment using KEGG
#Start by creating the folder into the FDR/KEGG directory
dir.create("./Enrichment/FDR/KEGG/Hsapiens/")

#check for the species name supported
search_kegg_organism("Sapiens", by="scientific_name",ignore.case = T)


#KEGG enrichment using significative info.
hs_wt_phenotype_genelist_enrichKEGG<-enrichKEGG(gene = hs_wt_phenotype_diff_expressed_genes,organism = "hsa",keyType = "kegg",pvalueCutoff = 0.05,pAdjustMethod = "BH",minGSSize = 10)
head(hs_wt_phenotype_genelist_enrichKEGG@result)

#Translate the ENTREZID
hs_wt_phenotype_genelist_enrichKEGG<-setReadable(hs_wt_phenotype_genelist_enrichKEGG,OrgDb = org.Hs.eg.db,keyType = "ENTREZID")

#Write the results of the table enriched with the KEGG significant data.
write.table(data.frame(hs_wt_phenotype_genelist_enrichKEGG@result),"./Enrichment/FDR/KEGG/Hsapiens/hs_wt_phenotype_enrichmentKEGG.txt",sep = "\t",quote = F, row.names = F,col.names = T)
head(hs_wt_phenotype_genelist_enrichKEGG@result)


#Plot the decided figures as required.
#Those figures got downloaded in the working directory
hsa04512<-pathview(gene.data = hs_wt_phenotype_diff_expressed_genes,
                   pathway.id = "hsa04512",
                   species="hsa",
                   limit = list(gene=5,cpd=1))

#Plot the pathways.
#Call the functions and plot it.
hs_wt_phenotype_genelist_enrichKEGG_topPathWays<-top_paths_FDR(hs_wt_phenotype_genelist_enrichKEGG,anno_human_wt_phenotype_stats,10)
head(hs_wt_phenotype_genelist_enrichKEGG_topPathWays[,c(2,5,6,10)]);tail(hs_wt_phenotype_genelist_enrichKEGG_topPathWays[,c(2,5,6,10)])


pdf("./Enrichment/FDR/KEGG/Hsapiens/hs_wt_phenotype_genelist_enrichKEGG_topPathWays.pdf",height = 12,width = 18)
ggplot(hs_wt_phenotype_genelist_enrichKEGG_topPathWays,aes(x=reorder(Description,zscore),y=zscore))+
  geom_bar(stat = "identity",fill="orange")+
  #geom_text(aes(label= round(p.adjust,3)),hjust=1.2,vjust=-.00001,size=3.5,color="white")+
  #labs(text="padjust")+
  geom_point(aes(label=p.adjust,size=rev(p.adjust)))+
  coord_flip()+
  labs(title = "Hs FDR KEGG",subtitle = sampleTable_wt_phenotype$treatment_type,size ="p.adjusted",x="Pathways")

dev.off()

#Save
write.table(hs_wt_phenotype_genelist_enrichKEGG_topPathWays,"./Enrichment/FDR/KEGG/Hsapiens/hs_wt_phenotype_genelist_enrichKEGG_topPathWays.tsv",sep = "\t",row.names = F,col.names = T)








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
cat("FDR KEGG Mmusculus ENRICHMENTS")
print("FDR KEGG Mmusculus ENRICHMENTS")

###MOUSE KEGG

#Start by creating the folder into the FDR/KEGG directory
dir.create("./Enrichment/FDR/KEGG/Mmusculus/")

#check for the species name supported
search_kegg_organism("musculus", by="scientific_name",ignore.case = T)

length(mm_wt_phenotype_diff_expressed_genes)
#KEGG enrichment using significative info.
mm_wt_phenotype_genelist_enrichKEGG<-enrichKEGG(gene = mm_wt_phenotype_diff_expressed_genes,organism = "mmu",keyType = "kegg",pvalueCutoff = 0.05,pAdjustMethod = "BH",minGSSize = 10)
head(mm_wt_phenotype_genelist_enrichKEGG@result)


#Translate the ENTREZID
mm_wt_phenotype_genelist_enrichKEGG<-setReadable(mm_wt_phenotype_genelist_enrichKEGG,OrgDb = org.Mm.eg.db,keyType = "ENTREZID")


#Write the results of the table enriched with the KEGG significant data.
write.table(data.frame(mm_wt_phenotype_genelist_enrichKEGG@result),"./Enrichment/FDR/KEGG/Mmusculus/mm_wt_phenotype_enrichmentKEGG.txt",sep = "\t",quote = F, row.names = F,col.names = T)


head(mm_wt_phenotype_genelist_enrichKEGG@result)
#Plot the decided figures as required.
#Those figures got downloaded in the working directory
mmu04062<-pathview(gene.data = mm_wt_phenotype_diff_expressed_genes,
                   pathway.id = "mmu04062",
                   species="mmu",
                   limit = list(gene=5,cpd=1))

#Plot the pathways.
#Call the functions and plot it.
mm_wt_phenotype_genelist_enrichKEGG_topPathWays<-top_paths_FDR(mm_wt_phenotype_genelist_enrichKEGG,anno_mouse_wt_phenotype_stats,15)
head(mm_wt_phenotype_genelist_enrichKEGG_topPathWays[,c(2,5,6,10)]);tail(mm_wt_phenotype_genelist_enrichKEGG_topPathWays[,c(2,5,6,10)])



pdf("./Enrichment/FDR/KEGG/Mmusculus/mm_wt_phenotype_genelist_enrichKEGG_topPathWays.pdf",height = 12,width = 18)
ggplot(mm_wt_phenotype_genelist_enrichKEGG_topPathWays,aes(x=reorder(Description,zscore),y=zscore))+
  geom_bar(stat = "identity",fill="orange")+
  #geom_text(aes(label= round(p.adjust,3)),hjust=1.2,vjust=-.00001,size=3.5,color="white")+
  #labs(text="padjust")+
  geom_point(aes(label=p.adjust,size=rev(p.adjust)))+
  coord_flip()+
  labs(title = "Mm FDR KEGG",subtitle = sampleTable_wt_phenotype$treatment_type,size ="p.adjusted",x="Pathways")

dev.off()

#Save
write.table(mm_wt_phenotype_genelist_enrichKEGG_topPathWays,"./Enrichment/FDR/KEGG/Mmusculus/mm_wt_phenotype_genelist_enrichKEGG_topPathWays.tsv",sep = "\t",row.names = F,col.names = T)







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


#Represent the figure.
colors_val<-c("#377EB8","#E41A1C")


#Plot selection function

top_paths_GSE<-function(GSE_pathways,n=15){
  if(length(GSE_pathways@result$ID)==0){
    warning("ERROR → The pathway enrichment was not returning informtion with these parameters!!!\n")
    cat("ERROR → The pathway enrichment was not returning informtion with these parameters!!!\n" )
    break
    
  }else{
    topUP<-GSE_pathways@result%>%
      filter(NES>0)%>%
      top_n(n,wt = -pvalue)
    topDOWN<-GSE_pathways@result%>%
      filter(NES<0)%>%
      top_n(n,wt = -pvalue)
    topPATHS<-bind_rows(topUP,topDOWN)
    return(topPATHS)}
}






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
cat("GSE GO Drerio ENRICHMENTS")
print("GSE GO Drerio ENRICHMENTS")
#Start with zebrafish
dir.create("./Enrichment/GSE/GO/Drerio")
#Crete the ranked_genelist
#First clean the data a little bit.

#Remove all duplicated data from the stats table, leaving all the unique ENTREZID
#dr_wt_phenotype_ranked_genes<-anno_zebrafish_wt_phenotype_stats[-which(duplicated(anno_zebrafish_wt_phenotype_stats$ENTREZID)),]
dim(anno_zebrafish_wt_phenotype_stats)
#first numeric vector using the log2FoldChange
dr_wt_phenotype_ranked_genes<-anno_zebrafish_wt_phenotype_stats$log2FoldChange

#second step name the vector
names(dr_wt_phenotype_ranked_genes)<-as.character(anno_zebrafish_wt_phenotype_stats$ENTREZID)
#length(anno_zebrafish_wt_phenotype_stats$ENTREZID)
#Sort the list in decreasing order.
dr_wt_phenotype_ranked_genes<-sort(dr_wt_phenotype_ranked_genes,decreasing = T)

#View the genes
head(dr_wt_phenotype_ranked_genes)

##############################
#        GSE                 #
#   Using GO→ BP            # 
#                            #
##############################
#Start with the GSE in GO databases with BP.
dr_wt_phenotype_GSE_GO_BP<-gseGO(geneList = dr_wt_phenotype_ranked_genes,ont="BP",OrgDb = DrOrgDb,
                                        keyType = "ENTREZID",
                                        nPerm=1000,
                                        minGSSize = 10,
                                        pvalueCutoff = 0.05,
                                        pAdjustMethod = "BH")

#Translate the ENTREZID
dr_wt_phenotype_GSE_GO_BP<-setReadable(dr_wt_phenotype_GSE_GO_BP,OrgDb = org.Dr.eg.db,keyType = "ENTREZID")

head(dr_wt_phenotype_GSE_GO_BP@result)

##############################
#        GSE                 #
#   Using GO→ MF            # 
#                            #
##############################

#GSE in GO database with MF.
dr_wt_phenotype_GSE_GO_MF<-gseGO(geneList = dr_wt_phenotype_ranked_genes,ont="MF",OrgDb = DrOrgDb,
                                        keyType = "ENTREZID",
                                        nPerm=1000,
                                        minGSSize = 10,
                                        pvalueCutoff = 0.05,
                                        pAdjustMethod = "BH")



#Translate the ENTREZID
dr_wt_phenotype_GSE_GO_MF<-setReadable(dr_wt_phenotype_GSE_GO_MF,OrgDb = org.Dr.eg.db,keyType = "ENTREZID")
head(dr_wt_phenotype_GSE_GO_MF@result)

##############################
#        GSE                 #
#   Using GO→ CC            # 
#                            #
##############################


#GSE in GO database with CC.
dr_wt_phenotype_GSE_GO_CC<-gseGO(geneList = dr_wt_phenotype_ranked_genes,ont="CC",OrgDb = DrOrgDb,
                                        keyType = "ENTREZID",
                                        nPerm=1000,
                                        minGSSize = 10,
                                        pvalueCutoff = 0.05,
                                        pAdjustMethod = "BH")



#Translate the ENTREZID
dr_wt_phenotype_GSE_GO_CC<-setReadable(dr_wt_phenotype_GSE_GO_CC,OrgDb = org.Dr.eg.db,keyType = "ENTREZID")
head(dr_wt_phenotype_GSE_GO_CC@result)



#SAVE THE TABLES
#Save the data of Dr into the specific file each table to pvalueCutOff=0.05.
#GSE for GO:BP 
write.table(data.frame(dr_wt_phenotype_GSE_GO_BP@result),"./Enrichment/GSE/GO/Drerio/dr_wt_phenotype_GSE_GO_BP.txt",sep="\t", row.names = F,col.names = T,quote = F)

#GSE for GO:MF
write.table(data.frame(dr_wt_phenotype_GSE_GO_MF@result),"./Enrichment/GSE/GO/Drerio/dr_wt_phenotype_GSE_GO_MF.txt",sep="\t",row.names = F,col.names = T,quote = F)

#GSE for CC
write.table(data.frame(dr_wt_phenotype_GSE_GO_CC@result),"./Enrichment/GSE/GO/Drerio/dr_wt_phenotype_GSE_GO_BP.txt",sep="\t",row.names = F,col.names = T,quote = F)





#Represent the figures, start with the BP,MF,CC.


#Plot the gseGO.



#GSE GO BP
#Call the functions and plot it.
dr_wt_phenotype_GSE_GO_BP_topPathWays<-top_paths_GSE(dr_wt_phenotype_GSE_GO_BP,20)
head(dr_wt_phenotype_GSE_GO_BP_topPathWays[,c(1,2,5,6,7)]);tail(dr_wt_phenotype_GSE_GO_BP_topPathWays[,c(1,2,5,6,7)])



pdf("./Enrichment/GSE/GO/Drerio/dr_wt_phenotype_GSE_GO_BP_topPathWays.pdf",height = 12,width = 18)
ggplot(dr_wt_phenotype_GSE_GO_BP_topPathWays,aes(x=reorder(Description,NES),y=NES))+
  geom_bar(stat = "identity",fill="darkorchid4")+
  #geom_text(aes(label= round(p.adjust,3)),hjust=1.2,vjust=-.00001,size=3.5,color="white")+
  #labs(text="padjust")+
  geom_point(aes(label=p.adjust,size=rev(p.adjust)))+
  coord_flip()+
  labs(title = "Dr GSE GO:BP",subtitle = sampleTable_wt_phenotype$treatment_type,size ="p-adjust",x="Biological Process",y="Normalized Enrichment Score")

dev.off()

#Save
write.table(dr_wt_phenotype_GSE_GO_BP_topPathWays,"./Enrichment/GSE/GO/Drerio/dr_wt_phenotype_GSE_GO_BP_topPathWays.tsv",sep = "\t",row.names = F,col.names = T)


#GSE GO MF
#Call the functions and plot it.
dr_wt_phenotype_GSE_GO_MF_topPathWays<-top_paths_GSE(dr_wt_phenotype_GSE_GO_MF,20)
head(dr_wt_phenotype_GSE_GO_MF_topPathWays[,c(1,2,5,6,7)]);tail(dr_wt_phenotype_GSE_GO_MF_topPathWays[,c(1,2,5,6,7)])



pdf("./Enrichment/GSE/GO/Drerio/dr_wt_phenotype_GSE_GO_MF_topPathWays.pdf",height = 12,width = 18)
ggplot(dr_wt_phenotype_GSE_GO_MF_topPathWays,aes(x=reorder(Description,NES),y=NES))+
  geom_bar(stat = "identity",fill="darkorchid4")+
  #geom_text(aes(label= round(p.adjust,3)),hjust=1.2,vjust=-.00001,size=3.5,color="white")+
  #labs(text="padjust")+
  geom_point(aes(label=p.adjust,size=rev(p.adjust)))+
  coord_flip()+
  labs(title = "Dr GSE GO:MF",subtitle = sampleTable_wt_phenotype$treatment_type,size ="p-adjust",x="Molecular Function",y="Normalized Enrichment Score")

dev.off()

#Save
write.table(dr_wt_phenotype_GSE_GO_MF_topPathWays,"./Enrichment/GSE/GO/Drerio/dr_wt_phenotype_GSE_GO_MF_topPathWays.tsv",sep = "\t",row.names = F,col.names = T)

#GSE GO CC
#Call the functions and plot it.
dr_wt_phenotype_GSE_GO_CC_topPathWays<-top_paths_GSE(dr_wt_phenotype_GSE_GO_CC,20)
head(dr_wt_phenotype_GSE_GO_CC_topPathWays[,c(1,2,5,6,7)]);tail(dr_wt_phenotype_GSE_GO_CC_topPathWays[,c(1,2,5,6,7)])



pdf("./Enrichment/GSE/GO/Drerio/dr_wt_phenotype_GSE_GO_CC_topPathWays.pdf",height = 12,width = 18)
ggplot(dr_wt_phenotype_GSE_GO_CC_topPathWays,aes(x=reorder(Description,NES),y=NES))+
  geom_bar(stat = "identity",fill="darkorchid4")+
  #geom_text(aes(label= round(p.adjust,3)),hjust=1.2,vjust=-.00001,size=3.5,color="white")+
  #labs(text="padjust")+
  geom_point(aes(label=p.adjust,size=rev(p.adjust)))+
  coord_flip()+
  labs(title = "Dr GSE GO:CC",subtitle = sampleTable_wt_phenotype$treatment_type,size ="p-adjust",x="Cellular Component",y="Normalized Enrichment Score")

dev.off()

#Save
write.table(dr_wt_phenotype_GSE_GO_CC_topPathWays,"./Enrichment/GSE/GO/Drerio/dr_wt_phenotype_GSE_GO_CC_topPathWays.tsv",sep = "\t",row.names = F,col.names = T)


#####################################################
#                                                   #
#               HUMAN                               #
#                                                   #
#####################################################
cat("GSE GO Hsapiens ENRICHMENTS")
print("GSE GO Hsapiens ENRICHMENTS")
#Create the folder for the HUMAN directrory
dir.create("./Enrichment/GSE/GO/Hsapiens")
#Crete the ranked_genelist
#First clean the data a little bit.

#Remove all duplicated data from the stats table, leaving all the unique ENTREZID
dim(anno_human_wt_phenotype_stats)
#anno_human_wt_phenotype_stats<-anno_human_wt_phenotype_stats[-which(duplicated(anno_human_wt_phenotype_stats$ENTREZID)),]
dim(anno_human_wt_phenotype_stats)

#first numeric vector using the log2FoldChange
hs_wt_phenotype_ranked_genes<-anno_human_wt_phenotype_stats$log2FoldChange

#second step name the vector
names(hs_wt_phenotype_ranked_genes)<-as.character(anno_human_wt_phenotype_stats$ENTREZID)

#Sort the list in decreasing order.
hs_wt_phenotype_ranked_genes<-sort(hs_wt_phenotype_ranked_genes,decreasing = T)
length(hs_wt_phenotype_ranked_genes)


##############################
#        GSE                 #
#   Using GO→ BP            # 
#                            #
##############################
#Start with the GSE in GO databases with BP.
hs_wt_phenotype_GSE_GO_BP<-gseGO(geneList = hs_wt_phenotype_ranked_genes,ont="BP",OrgDb = HsOrgDb,
                                        keyType = "ENTREZID",
                                        nPerm=1000,
                                        minGSSize = 10,
                                        pvalueCutoff = 0.05,
                                        pAdjustMethod = "BH")


#Translate the ENTREZID
hs_wt_phenotype_GSE_GO_BP<-setReadable(hs_wt_phenotype_GSE_GO_BP,OrgDb = org.Hs.eg.db,keyType = "ENTREZID")

head(hs_wt_phenotype_GSE_GO_BP@result)


##############################
#        GSE                 #
#   Using GO→ MF            # 
#                            #
##############################

#GSE in GO database with MF.
hs_wt_phenotype_GSE_GO_MF<-gseGO(geneList = hs_wt_phenotype_ranked_genes,ont="MF",OrgDb = HsOrgDb,
                                        keyType = "ENTREZID",
                                        nPerm=1000,
                                        minGSSize = 10,
                                        pvalueCutoff = 0.05,
                                        pAdjustMethod = "BH")


head(subset(hs_wt_phenotype_GSE_GO_MF@result,NES<0))
#Translate the ENTREZID
hs_wt_phenotype_GSE_GO_MF<-setReadable(hs_wt_phenotype_GSE_GO_MF,OrgDb = org.Hs.eg.db,keyType = "ENTREZID")
head(hs_wt_phenotype_GSE_GO_MF@result)

##############################
#        GSE                 #
#   Using GO→ CC            # 
#                            #
##############################


#GSE in GO database with CC.
hs_wt_phenotype_GSE_GO_CC<-gseGO(geneList = hs_wt_phenotype_ranked_genes,ont="CC",OrgDb = HsOrgDb,
                                        keyType = "ENTREZID",
                                        nPerm=1000,
                                        minGSSize = 10,
                                        pvalueCutoff = 0.05,
                                        pAdjustMethod = "BH")




#Translate the ENTREZID
hs_wt_phenotype_GSE_GO_CC<-setReadable(hs_wt_phenotype_GSE_GO_CC,OrgDb = org.Hs.eg.db,keyType = "ENTREZID")


head(hs_wt_phenotype_GSE_GO_CC@result)



#SAVE THE TABLES
#Save the data of Dr into the specific file each table to pvalueCutOff=0.05.
#GSE for GO:BP 
write.table(data.frame(hs_wt_phenotype_GSE_GO_BP@result),"./Enrichment/GSE/GO/Hsapiens/hs_wt_phenotype_GSE_GO_BP.txt",sep="\t", row.names = F,col.names = T,quote = F)

#GSE for GO:MF
write.table(data.frame(hs_wt_phenotype_GSE_GO_MF@result),"./Enrichment/GSE/GO/Hsapiens/hs_wt_phenotype_GSE_GO_MF.txt",sep="\t",row.names = F,col.names = T,quote = F)

#GSE for CC
write.table(data.frame(hs_wt_phenotype_GSE_GO_CC@result),"./Enrichment/GSE/GO/Hsapiens/hs_wt_phenotype_GSE_GO_BP.txt",sep="\t",row.names = F,col.names = T,quote = F)



#Represent the figure.
#Represent the figures, start with the BP,MF,CC.


#GSE GO BP

#Call the functions and plot it.
hs_wt_phenotype_GSE_GO_BP_topPathWays<-top_paths_GSE(hs_wt_phenotype_GSE_GO_BP,20)
head(hs_wt_phenotype_GSE_GO_BP_topPathWays[,c(1,2,5,6,7)]);tail(hs_wt_phenotype_GSE_GO_BP_topPathWays[,c(1,2,5,6,7)])



pdf("./Enrichment/GSE/GO/Hsapiens/hs_wt_phenotype_GSE_GO_BP_topPathWays.pdf",height = 12,width = 18)
ggplot(hs_wt_phenotype_GSE_GO_BP_topPathWays,aes(x=reorder(Description,NES),y=NES))+
  geom_bar(stat = "identity",fill="darkorchid4")+
  #geom_text(aes(label= round(p.adjust,3)),hjust=1.2,vjust=-.00001,size=3.5,color="white")+
  #labs(text="padjust")+
  geom_point(aes(label=p.adjust,size=rev(p.adjust)))+
  coord_flip()+
  labs(title = "Hs GSE GO:BP",subtitle = sampleTable_wt_phenotype$treatment_type,size ="p-adjust",x="Biological Process",y="Normalized Enrichment Score")

dev.off()

#Save
write.table(hs_wt_phenotype_GSE_GO_BP_topPathWays,"./Enrichment/GSE/GO/Hsapiens/hs_wt_phenotype_GSE_GO_BP_topPathWays.tsv",sep = "\t",row.names = F,col.names = T)


#GSE GO MF
#Call the functions and plot it.
hs_wt_phenotype_GSE_GO_MF_topPathWays<-top_paths_GSE(hs_wt_phenotype_GSE_GO_MF,20)
head(hs_wt_phenotype_GSE_GO_MF_topPathWays[,c(1,2,5,6,7)]);tail(hs_wt_phenotype_GSE_GO_MF_topPathWays[,c(1,2,5,6,7)])



pdf("./Enrichment/GSE/GO/Hsapiens/hs_wt_phenotype_GSE_GO_MF_topPathWays.pdf",height = 12,width = 18)
ggplot(hs_wt_phenotype_GSE_GO_MF_topPathWays,aes(x=reorder(Description,NES),y=NES))+
  geom_bar(stat = "identity",fill="darkorchid4")+
  #geom_text(aes(label= round(p.adjust,3)),hjust=1.2,vjust=-.00001,size=3.5,color="white")+
  #labs(text="padjust")+
  geom_point(aes(label=p.adjust,size=rev(p.adjust)))+
  coord_flip()+
  labs(title = "Hs GSE GO:MF",subtitle = sampleTable_wt_phenotype$treatment_type,size ="p-adjust",x="Molecular Function",y="Normalized Enrichment Score")

dev.off()

#Save
write.table(hs_wt_phenotype_GSE_GO_MF_topPathWays,"./Enrichment/GSE/GO/Hsapiens/hs_wt_phenotype_GSE_GO_MF_topPathWays.tsv",sep = "\t",row.names = F,col.names = T)


#GSE GO CC
#Call the functions and plot it.
hs_wt_phenotype_GSE_GO_CC_topPathWays<-top_paths_GSE(hs_wt_phenotype_GSE_GO_CC,20)
head(hs_wt_phenotype_GSE_GO_CC_topPathWays[,c(1,2,5,6,7)]);tail(hs_wt_phenotype_GSE_GO_CC_topPathWays[,c(1,2,5,6,7)])



pdf("./Enrichment/GSE/GO/Hsapiens/hs_wt_phenotype_GSE_GO_CC_topPathWays.pdf",height = 12,width = 18)
ggplot(hs_wt_phenotype_GSE_GO_CC_topPathWays,aes(x=reorder(Description,NES),y=NES))+
  geom_bar(stat = "identity",fill="darkorchid4")+
  #geom_text(aes(label= round(p.adjust,3)),hjust=1.2,vjust=-.00001,size=3.5,color="white")+
  #labs(text="padjust")+
  geom_point(aes(label=p.adjust,size=rev(p.adjust)))+
  coord_flip()+
  labs(title = "Hs GSE GO:CC",subtitle = sampleTable_wt_phenotype$treatment_type,size ="p-adjust",x="Cellular Component",y="Normalized Enrichment Score")

dev.off()

#Save
write.table(hs_wt_phenotype_GSE_GO_CC_topPathWays,"./Enrichment/GSE/GO/Hsapiens/hs_wt_phenotype_GSE_GO_CC_topPathWays.tsv",sep = "\t",row.names = F,col.names = T)




#####################################################
#                                                   #
#               MOUSE                               #
#                                                   #
#####################################################
cat("GSE GO Mmusculus ENRICHMENTS")
print("GSE GO Mmusculus ENRICHMENTS")
#Create the folder for the HUMAN directrory
dir.create("./Enrichment/GSE/GO/Mmusculus")
#Crete the ranked_genelist
#First clean the data a little bit.

#Remove all duplicated data from the stats table, leaving all the unique ENTREZID
dim(anno_mouse_wt_phenotype_stats)
#anno_mouse_wt_phenotype_stats<-anno_mouse_wt_phenotype_stats[-which(duplicated(anno_mouse_wt_phenotype_stats$ENTREZID)),]
dim(anno_mouse_wt_phenotype_stats)

#first numeric vector using the log2FoldChange
mm_wt_phenotype_ranked_genes<-anno_mouse_wt_phenotype_stats$log2FoldChange

#second step name the vector
names(mm_wt_phenotype_ranked_genes)<-as.character(anno_mouse_wt_phenotype_stats$ENTREZID)

#Sort the list in decreasing order.
mm_wt_phenotype_ranked_genes<-sort(mm_wt_phenotype_ranked_genes,decreasing = T)
length(mm_wt_phenotype_ranked_genes)


##############################
#        GSE                 #
#   Using GO→ BP            # 
#                            #
##############################
#Start with the GSE in GO databases with BP.
mm_wt_phenotype_GSE_GO_BP<-gseGO(geneList = mm_wt_phenotype_ranked_genes,ont="BP",OrgDb = MmOrgDb,
                                        keyType = "ENTREZID",
                                        nPerm=1000,
                                        minGSSize = 10,
                                        pvalueCutoff = 0.05,
                                        pAdjustMethod = "BH")


#Translate the ENTREZID
mm_wt_phenotype_GSE_GO_BP<-setReadable(mm_wt_phenotype_GSE_GO_BP,OrgDb = org.Mm.eg.db,keyType = "ENTREZID")

head(mm_wt_phenotype_GSE_GO_BP@result)


##############################
#        GSE                 #
#   Using GO→ MF            # 
#                            #
##############################

#GSE in GO database with MF.
mm_wt_phenotype_GSE_GO_MF<-gseGO(geneList =mm_wt_phenotype_ranked_genes,ont="MF",OrgDb = MmOrgDb,
                                        keyType = "ENTREZID",
                                        nPerm=1000,
                                        minGSSize = 10,
                                        pvalueCutoff = 0.05,
                                        pAdjustMethod = "BH")


#Translate the ENTREZID
mm_wt_phenotype_GSE_GO_MF<-setReadable(mm_wt_phenotype_GSE_GO_MF,OrgDb = org.Mm.eg.db,keyType = "ENTREZID")
head(mm_wt_phenotype_GSE_GO_MF@result)



##############################
#        GSE                 #
#   Using GO→ CC            # 
#                            #
##############################

#GSE in GO database with CC.
mm_wt_phenotype_GSE_GO_CC<-gseGO(geneList = mm_wt_phenotype_ranked_genes,ont="CC",OrgDb = MmOrgDb,
                                        keyType = "ENTREZID",
                                        nPerm=1000,
                                        minGSSize = 10,
                                        maxGSSize = 500,
                                        pvalueCutoff = 0.05,
                                        pAdjustMethod = "BH")



#Translate the ENTREZID
mm_wt_phenotype_GSE_GO_CC<-setReadable(mm_wt_phenotype_GSE_GO_CC,OrgDb = org.Mm.eg.db,keyType = "ENTREZID")

head(mm_wt_phenotype_GSE_GO_CC@result)


#SAVE THE TABLES
#Save the data of Dr into the specific file each table to pvalueCutOff=0.05.
#GSE for GO:BP 
write.table(data.frame(mm_wt_phenotype_GSE_GO_BP@result),"./Enrichment/GSE/GO/Mmusculus/mm_wt_phenotype_GSE_GO_BP.txt",sep="\t", row.names = F,col.names = T,quote = F)


#GSE for GO:MF
write.table(data.frame(mm_wt_phenotype_GSE_GO_MF@result),"./Enrichment/GSE/GO/Mmusculus/mm_wt_phenotype_GSE_GO_MF.txt",sep="\t",row.names = F,col.names = T,quote = F)

#GSE for CC
write.table(data.frame(mm_wt_phenotype_GSE_GO_CC@result),"./Enrichment/GSE/GO/Mmusculus/mm_wt_phenotype_GSE_GO_BP.txt",sep="\t",row.names = F,col.names = T,quote = F)


#Represent the figure.
#Represent the figures, start with the BP,MF,CC.
#GSE GO BP

#Call the functions and plot it.
mm_wt_phenotype_GSE_GO_BP_topPathWays<-top_paths_GSE(mm_wt_phenotype_GSE_GO_BP,20)
head(mm_wt_phenotype_GSE_GO_BP_topPathWays[,c(1,2,5,6,7)]);tail(mm_wt_phenotype_GSE_GO_BP_topPathWays[,c(1,2,5,6,7)])


pdf("./Enrichment/GSE/GO/Mmusculus/mm_wt_phenotype_GSE_GO_BP_topPathWays.pdf",height = 12,width = 18)
ggplot(mm_wt_phenotype_GSE_GO_BP_topPathWays,aes(x=reorder(Description,NES),y=NES))+
  geom_bar(stat = "identity",fill="darkorchid4")+
  #geom_text(aes(label= round(p.adjust,3)),hjust=1.2,vjust=-.00001,size=3.5,color="white")+
  #labs(text="padjust")+
  geom_point(aes(label=p.adjust,size=rev(p.adjust)))+
  coord_flip()+
  labs(title = "Mm GSE:GO BP",subtitle = sampleTable_wt_phenotype$treatment_type,size ="p-adjust",x="Biological Process",y="Normalized Enrichment Score")

dev.off()

#Save
write.table(mm_wt_phenotype_GSE_GO_BP_topPathWays,"./Enrichment/GSE/GO/Mmusculus/mm_wt_phenotype_GSE_GO_BP_topPathWays.tsv",sep = "\t",row.names = F,col.names = T)

#GSE GO MF

#Call the functions and plot it.
mm_wt_phenotype_GSE_GO_MF_topPathWays<-top_paths_GSE(mm_wt_phenotype_GSE_GO_MF,20)
head(mm_wt_phenotype_GSE_GO_MF_topPathWays[,c(1,2,5,6,7)]);tail(mm_wt_phenotype_GSE_GO_MF_topPathWays[,c(1,2,5,6,7)])


pdf("./Enrichment/GSE/GO/Mmusculus/mm_wt_phenotype_GSE_GO_MF_topPathWays.pdf",height = 12,width = 18)
ggplot(mm_wt_phenotype_GSE_GO_MF_topPathWays,aes(x=reorder(Description,NES),y=NES))+
  geom_bar(stat = "identity",fill="darkorchid4")+
  #geom_text(aes(label= round(p.adjust,3)),hjust=1.2,vjust=-.00001,size=3.5,color="white")+
  #labs(text="padjust")+
  geom_point(aes(label=p.adjust,size=rev(p.adjust)))+
  coord_flip()+
  labs(title = "Mm GSE GO:MF",subtitle = sampleTable_wt_phenotype$treatment_type,size ="p-adjust",x="Molecular Function",y="Normalized Enrichment Score")

dev.off()

#Save
write.table(mm_wt_phenotype_GSE_GO_MF_topPathWays,"./Enrichment/GSE/GO/Mmusculus/mm_wt_phenotype_GSE_GO_MF_topPathWays.tsv",sep = "\t",row.names = F,col.names = T)


#GSE GO CC

#Call the functions and plot it.
mm_wt_phenotype_GSE_GO_CC_topPathWays<-top_paths_GSE(mm_wt_phenotype_GSE_GO_CC,20)
head(mm_wt_phenotype_GSE_GO_CC_topPathWays[,c(1,2,5,6,7)]);tail(mm_wt_phenotype_GSE_GO_CC_topPathWays[,c(1,2,5,6,7)])


pdf("./Enrichment/GSE/GO/Mmusculus/mm_wt_phenotype_GSE_GO_CC_topPathWays.pdf",height = 12,width = 18)
ggplot(mm_wt_phenotype_GSE_GO_CC_topPathWays,aes(x=reorder(Description,NES),y=NES))+
  geom_bar(stat = "identity",fill="darkorchid4")+
  #geom_text(aes(label= round(p.adjust,3)),hjust=1.2,vjust=-.00001,size=3.5,color="white")+
  #labs(text="padjust")+
  geom_point(aes(label=p.adjust,size=rev(p.adjust)))+
  coord_flip()+
  labs(title = "Mm GSE GO:CC",subtitle = sampleTable_wt_phenotype$treatment_type,size ="p-adjust",x="Cellular Component",y="Normalized Enrichment Score")

dev.off()

#Save
write.table(mm_wt_phenotype_GSE_GO_CC_topPathWays,"./Enrichment/GSE/GO/Mmusculus/mm_wt_phenotype_GSE_GO_CC_topPathWays.tsv",sep = "\t",row.names = F,col.names = T)


#####################################################
#                                                   #
#           KEGGG ENRICHMENT with  GSE              #
#                                                   #
#            USINGGENE SET FOR KEGG                 #  
#                   PATHWAYS                        #
#                                                   #
#####################################################
cat("GSE KEGG ENRICHMENTS")
print("GSE KEGG ENRICHMENTS")

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
cat("GSE KEGG Drerio ENRICHMENTS")
print("GSE KEGG Drerio ENRICHMENTS")
dir.create("./Enrichment/GSE/KEGG/Drerio")

#Refresh for the organism used for GSE KEGG
search_kegg_organism("danio", by = "scientific_name", ignore.case = T)

#Start the enrichment calling KEGG
dr_wt_phenotype_GSE_KEGG<-gseKEGG(geneList = dr_wt_phenotype_ranked_genes, organism = "dre",
                                         nPerm = 1000,
                                         minGSSize = 10,
                                         maxGSSize = 500,
                                         pvalueCutoff = 0.05,
                                         pAdjustMethod = "BH")



#Translate the ENTREZID
dr_wt_phenotype_GSE_KEGG<-setReadable(dr_wt_phenotype_GSE_KEGG,OrgDb = org.Dr.eg.db,keyType = "ENTREZID")

dim(dr_wt_phenotype_GSE_KEGG@result)

#Paint a pathway using the list max and min, for the colors.
library(pathview)
#pdf("./plots/Oxidative_phosphorylation.pdf",width = 15,height = 15)
dre00190<-pathview(gene.data = hs_wt_phenotype_ranked_genes,
                   pathway.id = "dre00190",
                   species="dre",
                   limit = list(gene=max(abs(hs_wt_phenotype_ranked_genes)),cpd=1))

#Save the pathway table for the GSE with KEGG.
write.table(data.frame(dr_wt_phenotype_GSE_KEGG@result),"./Enrichment/GSE/KEGG/Drerio/dr_wt_phenotype_GSE_KEGG.txt",sep = "\t", quote = F, row.names = F,col.names = T)
#Save the figures desired.




#Plot the pathways.

#Call the functions and plot it.
dr_wt_phenotype_GSE_KEGG_topPathWays<-top_paths_GSE(dr_wt_phenotype_GSE_KEGG,20)
head(dr_wt_phenotype_GSE_KEGG_topPathWays[,c(1,2,5,6,7)],10);tail(dr_wt_phenotype_GSE_KEGG_topPathWays[,c(1,2,5,6,7)],10)
pdf("./Enrichment/GSE/KEGG/Drerio/dr_wt_phenotype_GSE_KEGG_topPathWays.pdf",height = 12,width = 18)
ggplot(dr_wt_phenotype_GSE_KEGG_topPathWays,aes(x=reorder(Description,NES),y=NES))+
  geom_bar(stat = "identity",fill="orange")+
  #geom_text(aes(label= round(p.adjust,3)),hjust=1.2,vjust=-.00001,size=3.5,color="white")+
  #labs(text="padjust")+
  geom_point(aes(label=p.adjust,size=rev(p.adjust)))+
  coord_flip()+
  labs(title = "Dr KEGG GSE",subtitle = sampleTable_wt_phenotype$treatment_type,size ="p-adjusted",x="Pathways",y="Normalized Enrichment Score")

dev.off()

#Save table top path ways
write.table(dr_wt_phenotype_GSE_KEGG_topPathWays,"./Enrichment/GSE/KEGG/Drerio/dr_wt_phenotype_GSE_KEGG_topPathWays.tsv",sep = "\t",col.names = T,row.names = F)

#####################################################
#                                                   #
#               HUMAN                               #
#                                                   #
#####################################################
#Start with Human
cat("GSE KEGG Hsapiens ENRICHMENTS")
print("GSE KEGG Hsapiens ENRICHMENTS")
dir.create("./Enrichment/GSE/KEGG/Hsapiens")

#Refresh for the organism used for GSE KEGG
search_kegg_organism("sapiens", by = "scientific_name", ignore.case = T)

#Start the enrichment calling KEGG
hs_wt_phenotype_GSE_KEGG<-gseKEGG(geneList = hs_wt_phenotype_ranked_genes, organism = "hsa",
                                         nPerm = 1000,
                                         minGSSize = 10,
                                         pvalueCutoff = 0.05,
                                         pAdjustMethod = "BH")



#Translate the EntrezIDs to symbols.
#columns(org.Hs.eg.db)
hs_wt_phenotype_GSE_KEGG<-setReadable(hs_wt_phenotype_GSE_KEGG,OrgDb = org.Hs.eg.db,keyType = "ENTREZID")

head(hs_wt_phenotype_GSE_KEGG@result)



#Paint a pathway using the list max and min, for the colors.
library(pathview)
hsa04621<-pathview(gene.data = hs_wt_phenotype_ranked_genes,
                   pathway.id = "hsa04621",
                   species="hsa",
                   limit = list(gene=max(abs(hs_wt_phenotype_ranked_genes)),cpd=1))

#Save the pathway table for the GSE with KEGG.
write.table(data.frame(hs_wt_phenotype_GSE_KEGG@result), "./Enrichment/GSE/KEGG/Hsapiens/hs_wt_phenotype_GSE_KEGG.txt",sep = "\t", quote = F, row.names = F,col.names = T)


#Plot the pathways.

hs_wt_phenotype_GSE_KEGG_topPathWays<-top_paths_GSE(hs_wt_phenotype_GSE_KEGG,20)
head(hs_wt_phenotype_GSE_KEGG_topPathWays[,c(1,2,5,6,7)],10);tail(hs_wt_phenotype_GSE_KEGG_topPathWays[,c(1,2,5,6,7)],10)




pdf("./Enrichment/GSE/KEGG/Hsapiens/hs_wt_phenotype_GSE_KEGG_topPathWays.pdf",height = 12,width = 18)
ggplot(hs_wt_phenotype_GSE_KEGG_topPathWays,aes(x=reorder(Description,NES),y=NES))+
  geom_bar(stat = "identity",fill="orange")+
  #geom_text(aes(label= round(p.adjust,3)),hjust=1.2,vjust=-.00001,size=3.5,color="white")+
  #labs(text="padjust")+
  geom_point(aes(label=p.adjust,size=rev(p.adjust)))+
  coord_flip()+
  labs(title = "Hs KEGG GSE",subtitle = sampleTable_wt_phenotype$treatment_type,size ="p-adjusted",x="Pathways",y="Normalized Enrichment Score")
dev.off()

#Save table top path ways
write.table(hs_wt_phenotype_GSE_KEGG_topPathWays,"./Enrichment/GSE/KEGG/Hsapiens/hs_wt_phenotype_GSE_KEGG_topPathWays.tsv",sep = "\t",col.names = T,row.names = F)


#####################################################
#                                                   #
#               MOUSE                               #
#                                                   #
#####################################################
cat("GSE KEGG Mmusculus ENRICHMENTS")
print("GSE KEGG Mmusculus ENRICHMENTS")
#Start with mouse
dir.create("./Enrichment/GSE/KEGG/Mmusculus")

#Refresh for the organism used for GSE KEGG
search_kegg_organism("musculus", by = "scientific_name", ignore.case = T)

#Start the enrichment calling KEGG
mm_wt_phenotype_GSE_KEGG<-gseKEGG(geneList = mm_wt_phenotype_ranked_genes, organism = "mmu",
                                         nPerm = 1000,
                                         minGSSize = 10,
                                         pvalueCutoff = 0.05,
                                         pAdjustMethod = "BH")

#Translate the ENTREZID
columns(org.Mm.eg.db)
mm_wt_phenotype_GSE_KEGG<-setReadable(mm_wt_phenotype_GSE_KEGG,OrgDb = org.Mm.eg.db,keyType = "ENTREZID")

head(mm_wt_phenotype_GSE_KEGG@result,20)

mmu04740<-pathview(gene.data = mm_wt_phenotype_ranked_genes,
                   pathway.id = "mmu04740",
                   species="mmu",
                   limit = list(gene=max(abs(mm_wt_phenotype_ranked_genes)),cpd=1))

#Save the pathway table for the GSE with KEGG.
write.table(data.frame(mm_wt_phenotype_GSE_KEGG@result), "./Enrichment/GSE/KEGG/Mmusculus/mm_wt_phenotype_GSE_KEGG.txt",sep = "\t", quote = F, row.names = F,col.names = T)


#Plot the pathways.
mm_wt_phenotype_GSE_KEGG_topPathWays<-top_paths_GSE(mm_wt_phenotype_GSE_KEGG,20)
head(mm_wt_phenotype_GSE_KEGG_topPathWays[,c(1,2,5,6,7)],10);tail(mm_wt_phenotype_GSE_KEGG_topPathWays[,c(1,2,5,6,7)],10)


pdf("./Enrichment/GSE/KEGG/Mmusculus/mm_wt_phenotype_GSE_KEGG_topPathWays.pdf",height = 12,width = 18)
ggplot(mm_wt_phenotype_GSE_KEGG_topPathWays,aes(x=reorder(Description,NES),y=NES))+
  geom_bar(stat = "identity",fill="orange")+
  #geom_text(aes(label= round(p.adjust,3)),hjust=1.2,vjust=-.00001,size=3.5,color="white")+
  #labs(text="padjust")+
  geom_point(aes(label=p.adjust,size=rev(p.adjust)))+
  coord_flip()+
  labs(title = "Mm KEGG GSE",subtitle = sampleTable_wt_phenotype$treatment_type,size ="p-adjusted",x="Pathways",y="Normalized Enrichment Score")
dev.off()


#Save table top path ways
write.table(mm_wt_phenotype_GSE_KEGG_topPathWays,"./Enrichment/GSE/KEGG/Mmusculus/mm_wt_phenotype_GSE_KEGG_topPathWays.tsv",sep = "\t",col.names = T,row.names = F)



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


#Plot selection function
top_paths_fGSEA<-function(fGSEA_pathways,n=15){
  
  if(length(fGSEA_pathways$pathway)==0){
    
    warning("ERROR → The pathway enrichment was not returning informtion with these parameters!!!\n")
    
    cat("ERROR → The pathway enrichment was not returning informtion with these parameters!!!\n" )
    
    break
    
  }else{
    topUP<-fGSEA_pathways %>%
      dplyr::filter(NES>0 & padj <= 0.05) %>%
      dplyr::top_n(n,wt = -pval)
    
    topDOWN<-fGSEA_pathways%>%
      dplyr::filter(NES<0 & padj <= 0.05) %>%
      dplyr::top_n(n,wt = -pval)
    
    topPATHS<-dplyr::bind_rows(topUP,topDOWN)
    
    return(topPATHS)}
}


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
print("HUMAN HALLMARKS fGSEA")
cat("HUMAN HALLMARKS fGSEA")
#Create the pathway fo the Hallmarks data.
dir.create(recursive = T, "./Enrichment/GSEA/Hallmarks/Hsapiens")

#Load the Hallmark data
load("/mnt/haus/marius/Carol_data/gsea_bundle/human_H_v5p2.rdata")

#Save the 50 pathways composing the hallmark GSEAinto an object
hs_hallmark_fgsea<-Hs.H


#Use the ranked ENTREZIDs used before.
length(hs_wt_phenotype_ranked_genes)

#Conduct analysis the GSEA enrichment for the Hallmarks!
hs_wt_phenotype_fgsea_hallmarks<-fgsea(pathways = hs_hallmark_fgsea,stats = hs_wt_phenotype_ranked_genes, 
                                              minSize = 10, 
                                              maxSize =500, 
                                              nperm = 1000,
                                              gseaParam = 1)

#Now lets see the 30 top adjusted.
head(hs_wt_phenotype_fgsea_hallmarks, n=30)


#Map the leadingEdges from ENTREZIDs to SYMBOLS for being more human readable.
hs_wt_phenotype_fgsea_hallmarks[,leadingEdge:=lapply(leadingEdge,mapIds, x=org.Hs.eg.db, keytype="ENTREZID",column="SYMBOL")]
head(hs_wt_phenotype_fgsea_hallmarks)

#Order the data according to padj value
hs_wt_phenotype_fgsea_hallmarks<-hs_wt_phenotype_fgsea_hallmarks[order(hs_wt_phenotype_fgsea_hallmarks[,3])]

#Remove the data which padj=1
hs_wt_phenotype_fgsea_hallmarks<-hs_wt_phenotype_fgsea_hallmarks[hs_wt_phenotype_fgsea_hallmarks$padj!=1]


#head(hs_wt_phenotype_fgsea_hallmarks$leadingEdge)
#hs_fgsearestidy<-select(hs_wt_phenotype_fgsea_hallmarks,-leadingEdge,-ES,-nMoreExtreme)%>%
#  arrange(padj)
#head(hs_fgsearestidy)
#pdf("./Enrichment/GSEA/Hallmarks/Hsapiens/hs_plot_gsea_bars_bars.pdf",width = 20,height = 30)
#plot_gsea_bars_coord<-ggplot(hs_fgsearestidy, aes(reorder(pathway, NES), NES)) +
#  geom_col(aes(fill=ifelse(pval < 0.05,"Significative","No-significative"))) +
#  coord_flip()+
#  labs(x="Pathway", y="Normalized Enrichment Score",
#       title="Hallmark pathways NES from GSEA") + 
#  scale_fill_manual(name="p-value", values =colors_val)+
#  theme_minimal()
#dev.off()
#plot_gsea_bars_coord



#Save the data of the GSEA HALLMARK PATHWAYS, using a special function for writting the data table data frame with data.table::fwrite..
data.table::fwrite(hs_wt_phenotype_fgsea_hallmarks, "./Enrichment/GSEA/Hallmarks/Hsapiens/hs_wt_phenotype_fgsea_hallmark.txt",sep = "\t",col.names = T, row.names = F, quote = F)

#Choose from the list any of the interested gsea plot enrichment interested.
pdf("./Enrichment/GSEA/Hallmarks/Hsapiens/hs_fgsea_hallmarks_function1.pdf",width = 16,height = 10)
plotEnrichment(hs_hallmark_fgsea[[hs_wt_phenotype_fgsea_hallmarks[[1]][1]]],
               hs_wt_phenotype_ranked_genes) + labs(title=hs_wt_phenotype_fgsea_hallmarks[[1]][1])
dev.off()





#Plot the fgsea table for significance based on Enrichment Scores.
#GSEA table plot.
hs_wt_phenotype_fgsea_hallmarks_topPathWays<-top_paths_fGSEA(hs_wt_phenotype_fgsea_hallmarks,10)
pdf("./Enrichment/GSEA/Hallmarks/Hsapiens/hs_wt_phenotype_fgsea_topPathWays_hallmarks_TablePlot.pdf",width = 9,height = 6)
#png("./plots/gseaPlottableHallmarks.png",1000,1000, pointsize = 20)
plotGseaTable(pathways = hs_hallmark_fgsea[hs_wt_phenotype_fgsea_hallmarks_topPathWays$pathway],
              hs_wt_phenotype_ranked_genes,
              hs_wt_phenotype_fgsea_hallmarks,
              gseaParam = 1)
dev.off()


#Plot in a more visual way the pathways up and down regulated.
#Create the plot of the gsea file using ggplplot



#Plot the pathways.
head(hs_wt_phenotype_fgsea_hallmarks_topPathWays[,c(1,2,3,5)],10);tail(hs_wt_phenotype_fgsea_hallmarks_topPathWays[,c(1,2,3,5)],10)


pdf("./Enrichment/GSEA/Hallmarks/Hsapiens/hs_wt_phenotype_fgsea_hallmarks_topPathWays.pdf",height = 12,width = 18)
ggplot(hs_wt_phenotype_fgsea_hallmarks_topPathWays,aes(x=reorder(pathway,NES),y=NES))+
  geom_bar(stat = "identity",fill="deepskyblue1")+
  #geom_text(aes(label= round(p.adjust,3)),hjust=1.2,vjust=-.00001,size=3.5,color="white")+
  #labs(text="padjust")+
  geom_point(aes(label=padj,size=rev(padj)))+
  coord_flip()+
  labs(title = "Hs fGSEA Hallmarks",subtitle = sampleTable_wt_phenotype$treatment_type,size ="p-adjusted",x="HUMAN HALLMARKS",y="Normalized Enrichment Score")
dev.off()




#Save the filtered data table.
fwrite(hs_wt_phenotype_fgsea_hallmarks_topPathWays,"./Enrichment/GSEA/Hallmarks/Hsapiens/hs_wt_phenotype_fgsea_hallmarks_topPathWays.tsv",sep = "\t",col.names = T,row.names = F)
#####################################################
#                                                   #
#               MOUSE                               #
#                                                   #
#               HALLMARKS                           #
#                                                   #
#####################################################
#Start with Mouse Hallmarks
print("MOUSE HALLMARKS fGSEA")
cat("MOUSE HALLMARKS fGSEA")
#Create the pathway fo the Hallmarks data.
dir.create(recursive = T, "./Enrichment/GSEA/Hallmarks/Mmusculus")

#Load the Hallmark data
load("/mnt/haus/marius/Carol_data/gsea_bundle/mouse_H_v5p2.rdata")

#Save the 50 pathways composing the hallmark GSEAinto an object
mm_hallmark_fgsea<-Mm.H


#Use the ranked ENTREZIDs used before.
length(mm_wt_phenotype_ranked_genes)

#Conduct analysis the GSEA enrichment for the Hallmarks!
mm_wt_phenotype_fgsea_hallmarks<-fgsea(pathways = mm_hallmark_fgsea,
                                              stats = mm_wt_phenotype_ranked_genes,
                                              minSize = 15,
                                              maxSize =500,
                                              nperm = 1000,
                                              gseaParam = 1)

#Now lets see the 30 top adjusted.
head(mm_wt_phenotype_fgsea_hallmarks[order(mm_wt_phenotype_fgsea_hallmarks$padj),], n=30)


#Map the leadingEdges from ENTREZIDs to SYMBOLS for being more human readable.
mm_wt_phenotype_fgsea_hallmarks[,leadingEdge:=lapply(leadingEdge,mapIds, x=org.Mm.eg.db, keytype="ENTREZID",column="SYMBOL")]
head(mm_wt_phenotype_fgsea_hallmarks)


#Order the data according to padj value
mm_wt_phenotype_fgsea_hallmarks<-mm_wt_phenotype_fgsea_hallmarks[order(mm_wt_phenotype_fgsea_hallmarks[,3])]  

#Remove the data which padj=1
mm_wt_phenotype_fgsea_hallmarks<-mm_wt_phenotype_fgsea_hallmarks[mm_wt_phenotype_fgsea_hallmarks$padj!=1]

#Save the data of the GSEA HALLMARK PATHWAYS, using a special function for writting the data table data frame with data.table::fwrite..
data.table::fwrite(mm_wt_phenotype_fgsea_hallmarks, "./Enrichment/GSEA/Hallmarks/Mmusculus/mm_wt_phenotype_fgsea_hallmark.txt",sep = "\t",col.names = T, row.names = F, quote = F)



#Check the name of the Hallmarks to plot, if they represent the same data in Juman and Mouse, use the second one for mouse..
if (mm_wt_phenotype_fgsea_hallmarks[[1]][1]==hs_wt_phenotype_fgsea_hallmarks[[1]][1]){
  function_to_plot<-mm_wt_phenotype_fgsea_hallmarks[[1]][2]
}else{
  function_to_plot<-mm_wt_phenotype_fgsea_hallmarks[[1]][1]
}

#Choose from the list any of the interested gsea plot enrichment interested.
pdf("./Enrichment/GSEA/Hallmarks/Mmusculus/mm_fgsea_hallmarks_function1.pdf",width = 16,height = 10)
plotEnrichment(mm_hallmark_fgsea[[function_to_plot]], mm_wt_phenotype_ranked_genes) + labs(title=function_to_plot)
dev.off()





#Get the filterd and organized hallmarks from the enrichment with gsea. (NES >0 and NES <0 & padj<=0.05)
mm_wt_phenotype_fgsea_hallmarks_topPathWays<-top_paths_fGSEA(mm_wt_phenotype_fgsea_hallmarks)
head(mm_wt_phenotype_fgsea_hallmarks_topPathWays[,c(1,2,3,5)]);tail(mm_wt_phenotype_fgsea_hallmarks_topPathWays[,c(1,2,3,5)])

#Plot the fgsea table for significance based on Enrichment Scores.
#GSEA table plot.
pdf("./Enrichment/GSEA/Hallmarks/Mmusculus/mm_wt_phenotype_fgsea_topPathWays_hallmarks_TablePlot.pdf",width = 9,height = 6)
#png("./plots/gseaPlottableHallmarks.png",1000,1000, pointsize = 20)
plotGseaTable(pathways = mm_hallmark_fgsea[mm_wt_phenotype_fgsea_hallmarks_topPathWays$pathway],
              mm_wt_phenotype_ranked_genes,
              mm_wt_phenotype_fgsea_hallmarks,
              gseaParam = 1)
dev.off()



#Plot in a more visual way the pathways up and down regulated.
#Create the plot of the gsea file using ggplplot


pdf("./Enrichment/GSEA/Hallmarks/Mmusculus/mm_wt_phenotype_fgsea_hallmarks_topPathWays.pdf",height = 12,width = 18)
ggplot(mm_wt_phenotype_fgsea_hallmarks_topPathWays,aes(x=reorder(pathway,NES),y=NES))+
  geom_bar(stat = "identity",fill="deepskyblue1")+
  #geom_text(aes(label= round(p.adjust,3)),hjust=1.2,vjust=-.00001,size=3.5,color="white")+
  #labs(text="padjust")+
  geom_point(aes(label=padj,size=rev(padj)))+
  coord_flip()+
  labs(title = "Mm fGSEA Hallmarks",subtitle = sampleTable_wt_phenotype$treatment_type,size ="p-adjusted",x="MOUSE HALLMARKS",y="Normalized Enrichment Score")
dev.off()




#Save the filtered data table.
fwrite(mm_wt_phenotype_fgsea_hallmarks_topPathWays,"./Enrichment/GSEA/Hallmarks/Mmusculus/mm_wt_phenotype_fgsea_hallmarks_topPathWays.tsv",sep = "\t",col.names = T,row.names = F)



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
print("HUMAN GO fGSEA")
cat("HUMAN GO fGSEA")
#Create the pathway fo the Hallmarks data.
dir.create(recursive = T, "./Enrichment/GSEA/GO/Hsapiens")

#Load the Hallmark data
load("/mnt/haus/marius/Carol_data/gsea_bundle/human_c5_v5p2.rdata")

#Save the 50 pathways composing the hallmark GSEAinto an object
hs_GO_fgsea<-Hs.c5


#Use the ranked ENTREZIDs used before.
length(hs_wt_phenotype_ranked_genes)

#Conduct analysis the GSEA enrichment for the Hallmarks!
hs_wt_phenotype_fgsea_GO<-fgsea(pathways = hs_GO_fgsea,
                                       stats = hs_wt_phenotype_ranked_genes,
                                       minSize = 15,
                                       maxSize =500,
                                       nperm = 1000,
                                       gseaParam = 1)

#Now lets see the 30 top adjusted.
head(hs_wt_phenotype_fgsea_GO[order(hs_wt_phenotype_fgsea_GO$padj,decreasing = F),])


#Check how much significative data would be taken with different padj values.
hs_wt_phenotype_fgsea_GO<-hs_wt_phenotype_fgsea_GO[hs_wt_phenotype_fgsea_GO$padj<.5]

#Map the leadingEdges from ENTREZIDs to SYMBOLS for being more human readable.
hs_wt_phenotype_fgsea_GO[,leadingEdge:=lapply(leadingEdge,mapIds, x=org.Hs.eg.db, keytype="ENTREZID",column="SYMBOL")]
head(hs_wt_phenotype_fgsea_GO)

#Order the data according to padj value
hs_wt_phenotype_fgsea_GO<-hs_wt_phenotype_fgsea_GO[order(hs_wt_phenotype_fgsea_GO[,3])]

#Remove the data which padj=1
hs_wt_phenotype_fgsea_GO<-hs_wt_phenotype_fgsea_GO[hs_wt_phenotype_fgsea_GO$padj!=1]
head(hs_wt_phenotype_fgsea_GO)



#Save the data of the GSEA HALLMARK PATHWAYS, using a special function for writting the data table data frame with data.table::fwrite..
data.table::fwrite(hs_wt_phenotype_fgsea_GO, "./Enrichment/GSEA/GO/Hsapiens/hs_wt_phenotype_fgsea_GO.txt",sep = "\t",col.names = T, row.names = F, quote = F)


#Choose from the list any of the interested gsea plot enrichment interested.
pdf("./Enrichment/GSEA/GO/Hsapiens/hs_fgsea_GO_function1.pdf",width = 16,height =16)
plotEnrichment(hs_GO_fgsea[[hs_wt_phenotype_fgsea_GO[[1]][1]]],
               hs_wt_phenotype_ranked_genes) + labs(title=hs_wt_phenotype_fgsea_GO[[1]][1])
dev.off()


#Plot the fgsea table for significance based on Enrichment Scores.
#GSEA table plot.
hs_wt_phenotype_fgsea_GO_topPathways<-top_paths_fGSEA(hs_wt_phenotype_fgsea_GO)
head(hs_wt_phenotype_fgsea_GO_topPathways[,c(1,2,3,5)]);tail(hs_wt_phenotype_fgsea_GO_topPathways[,c(1,2,3,5)])


pdf("./Enrichment/GSEA/GO/Hsapiens/hs_wt_phenotype_fgsea_GO_topPathways_TablePlot.pdf",width =30,height = 50)
#png("./plots/gseaPlottableHallmarks.png",1000,1000, pointsize = 20)
plotGseaTable(pathways = hs_GO_fgsea[hs_wt_phenotype_fgsea_GO_topPathways$pathway],
              hs_wt_phenotype_ranked_genes,
              hs_wt_phenotype_fgsea_GO,
              gseaParam = 1)
dev.off()



#Plot in a more visual way the pathways up and down regulated.
#Create the plot of the gsea file using ggplplot

pdf("./Enrichment/GSEA/GO/Hsapiens/hs_wt_phenotype_fgsea_GO_topPathways.pdf",width = 20,height = 25)
ggplot(hs_wt_phenotype_fgsea_GO_topPathways,aes(x=reorder(pathway,NES),y=NES))+
  geom_bar(stat = "identity",fill="deepskyblue1")+
  #geom_text(aes(label= round(p.adjust,3)),hjust=1.2,vjust=-.00001,size=3.5,color="white")+
  #labs(text="padjust")+
  geom_point(aes(label=padj,size=rev(padj)))+
  coord_flip()+
  labs(title = "Hs fGSEA GO",subtitle = sampleTable_wt_phenotype$treatment_type,size ="p-adjusted",x="HUMAN fGSEA:GO",y="Normalized Enrichment Score")
dev.off()


#Save the filtered data table.
fwrite(hs_wt_phenotype_fgsea_GO_topPathways,"./Enrichment/GSEA/GO/Hsapiens/hs_wt_phenotype_fgsea_GO_topPathways.tsv",sep = "\t",col.names = T,row.names = F)

#PROCEED with the same process but for MOUSE DATABASE.
#####################################################
#                                                   #
#               MOUSE                               #
#                                                   #
#               GO DATABASES                        #
#                                                   #
#####################################################
print("MOUSE GO fGSEA")
cat("MOUSE GO fGSEA")
#Start with MOUSE

#Create the pathway fo the Hallmarks data.
dir.create(recursive = T, "./Enrichment/GSEA/GO/Mmusculus")

#Load the Hallmark data
load("/mnt/haus/marius/Carol_data/gsea_bundle/mouse_c5_v5p2.rdata")

#Save the 50 pathways composing the hallmark GSEAinto an object
mm_GO_fgsea<-Mm.c5

#Conduct analysis the GSEA enrichment for the Hallmarks!
mm_wt_phenotype_fgsea_GO<-fgsea(pathways = mm_GO_fgsea,
                                       stats = mm_wt_phenotype_ranked_genes,
                                       minSize = 15,
                                       maxSize =500,
                                       nperm = 1000,
                                       gseaParam = 1)

#Now lets see the 30 top adjusted.
head(mm_wt_phenotype_fgsea_GO[order(mm_wt_phenotype_fgsea_GO$padj),], n=30)


#FOR GO ENRICHMENT REDUCE THE LIST UP TO padj=0.05
mm_wt_phenotype_fgsea_GO<-mm_wt_phenotype_fgsea_GO[mm_wt_phenotype_fgsea_GO$padj<0.5]

#Map the leadingEdges from ENTREZIDs to SYMBOLS for being more human readable.
mm_wt_phenotype_fgsea_GO[,leadingEdge:=lapply(leadingEdge,mapIds, x=org.Mm.eg.db, keytype="ENTREZID",column="SYMBOL")]
head(mm_wt_phenotype_fgsea_GO)

#Check how much significative data would be taken with different padj values.
dim(mm_wt_phenotype_fgsea_GO[mm_wt_phenotype_fgsea_GO$padj < 0.05])

#Order the data according to padj value
mm_wt_phenotype_fgsea_GO<-mm_wt_phenotype_fgsea_GO[order(mm_wt_phenotype_fgsea_GO[,3])]

#Remove the data which padj=1
mm_wt_phenotype_fgsea_GO<-mm_wt_phenotype_fgsea_GO[mm_wt_phenotype_fgsea_GO$padj!=1]
head(mm_wt_phenotype_fgsea_GO)

#Save the data of the GSEA HALLMARK PATHWAYS, using a special function for writting the data table data frame with data.table::fwrite..
data.table::fwrite(mm_wt_phenotype_fgsea_GO, "./Enrichment/GSEA/GO/Mmusculus/mm_wt_phenotype_fgsea_GO.txt",sep = "\t",col.names = T, row.names = F, quote = F)




#Select the function from the go_data to plot a different from human, that will say they got the same results.
if (mm_wt_phenotype_fgsea_GO[[1]][1]==hs_wt_phenotype_fgsea_GO[[1]][1]){
  function_to_plot2<-mm_wt_phenotype_fgsea_GO[[1]][2]
}else{
  function_to_plot2<-mm_wt_phenotype_fgsea_GO[[1]][1]
}


#Choose from the list any of the interested gsea plot enrichment interested.
pdf("./Enrichment/GSEA/GO/Mmusculus/mm_fgsea_GO_function1.pdf",width = 16,height =16)
plotEnrichment(mm_GO_fgsea[[function_to_plot2]], mm_wt_phenotype_ranked_genes) + labs(title=function_to_plot2)
dev.off()




mm_wt_phenotype_fgsea_GO_topPathways<-top_paths_fGSEA(mm_wt_phenotype_fgsea_GO,20)
head(mm_wt_phenotype_fgsea_GO_topPathways[,c(1,2,3,5)]);tail(mm_wt_phenotype_fgsea_GO_topPathways[,c(1,2,3,5)])
pdf("./Enrichment/GSEA/GO/Mmusculus/mm_wt_phenotype_fgsea_topPathWays_GO_TablePlot.pdf",width =30,height = 50)
#png("./plots/gseaPlottableHallmarks.png",1000,1000, pointsize = 20)
plotGseaTable(pathways = mm_GO_fgsea[mm_wt_phenotype_fgsea_GO_topPathways$pathway],
              mm_wt_phenotype_ranked_genes,
              mm_wt_phenotype_fgsea_GO,
              gseaParam = 1)





#Plot in a more visual way the pathways up and down regulated.
#Create the plot of the gsea file using ggplplot

pdf("./Enrichment/GSEA/GO/Mmusculus/mm_wt_phenotype_fgsea_GO_topPathways.pdf",width = 30,height = 25)
ggplot(mm_wt_phenotype_fgsea_GO_topPathways,aes(x=reorder(pathway,NES),y=NES))+
  geom_bar(stat = "identity",fill="deepskyblue1")+
  #geom_text(aes(label= round(p.adjust,3)),hjust=1.2,vjust=-.00001,size=3.5,color="white")+
  #labs(text="padjust")+
  geom_point(aes(label=padj,size=rev(padj)))+
  coord_flip()+
  labs(title = "Mm fGSEA GO",subtitle = sampleTable_wt_phenotype$treatment_type,size ="p-adjusted",x="MOUSE HALLMARKS",y="Normalized Enrichment Score")
dev.off()


#Save the filtered data table.
fwrite(mm_wt_phenotype_fgsea_GO_topPathways,"./Enrichment/GSEA/GO/Mmusculus/mm_wt_phenotype_fgsea_GO_topPathways.tsv",sep = "\t",col.names = T,row.names = F)


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
cat("Painting the Volcanos")
print("Painting the Volcanos")

#####################################################
#                                                   #
#                   HUMAN                           #
#                                                   #
#                                                   #
#####################################################

#Create the directory for plots.
dir.create(recursive = T,"./Plots/Hsapiens/")


#Use the anotated table for human.
head(anno_human_wt_phenotype_stats)

pdf("./Plots/Hsapiens/VolcanoHuman.pdf", width = 12, height =18)
#png("../plots/enhancedVolcano.png", 1000,1000,pointsize = 20)
EnhancedVolcano(anno_human_wt_phenotype_stats,lab=anno_human_wt_phenotype_stats$SYMBOL, x = "log2FoldChange", y="pvalue",legend = c("NS","Log (base 2) fold-change","Adjusted p-value", "Adjusted p-value & Log (base 2) fold-change"),
                legendPosition = "bottom", title = "wt_phenotype Hs", cutoffLineCol = "deeppink",FCcutoff = 1)#, selectLab = keep)

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
head(anno_zebrafish_wt_phenotype_stats)

pdf("./Plots/Zebrafish/VolcanoZebra.pdf", width = 12, height =18)
#png("../plots/enhancedVolcano.png", 1000,1000,pointsize = 20)
EnhancedVolcano(anno_zebrafish_wt_phenotype_stats,lab=anno_zebrafish_wt_phenotype_stats$SYMBOL, x = "log2FoldChange", y="pvalue",legend = c("NS","Log (base 2) fold-change","Adjusted p-value", "Adjusted p-value & Log (base 2) fold-change"),
                legendPosition = "bottom", title = "wt_phenotype Dre", cutoffLineCol = "deeppink",FCcutoff = 1)#, selectLab = keep)

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
head(anno_mouse_wt_phenotype_stats)

pdf("./Plots/Mouse/VolcanoMouse.pdf", width = 12, height =18)
#png("../plots/enhancedVolcano.png", 1000,1000,pointsize = 20)
EnhancedVolcano(anno_mouse_wt_phenotype_stats,lab=anno_mouse_wt_phenotype_stats$SYMBOL, x = "log2FoldChange", y="pvalue",legend = c("NS","Log (base 2) fold-change","Adjusted p-value", "Adjusted p-value & Log (base 2) fold-change"),
                legendPosition = "bottom", title = "wt_phenotype Mmus", cutoffLineCol = "deeppink",FCcutoff = 1)#, selectLab = keep)

dev.off()




#####################################################
#                                                   #
#                                                   #
#                     PCA                           #
#####################################################
cat("Painting the PCA's")
print("Painting the PCA's")
#Create the data used for the PCA and  Heatmap by normalizitaion of the samples.
#Use rlog from DESeq2 to regularized log tranformation ignoring information about experimental groups, blind=T
dds_rlog_wt_phenotype=DESeq2::rlog(dds_wt_phenotype,blind = T)

dim(dds_rlog_wt_phenotype)
head(assay(dds_rlog_wt_phenotype))#_wt_phenotype))
head(sampleTable_wt_phenotype)
head(assay(dds_rlog_wt_phenotype))

#Create the PCA using the DESeq2 funciton
dir.create(recursive = T, path = "./PCA/")
#Plot PCA

dds_rlog_wt_phenotype@assays$data
#DESeq2::plotPCA(dds_rlog_wt_phenotype,"CELLTYPE")

pdf("./PCA/wt_phenotype_PCA.pdf",width = 12,height = 12)
z<-DESeq2::plotPCA(dds_rlog_wt_phenotype,intgroup = "CELLTYPE")
z+geom_label(aes(label=rownames(sampleTable_wt_phenotype)))
dev.off()


#####################################################
#                                                   #
#                   HUMAN                           #
#                                                   #
#                                                   #
#####################################################

#Directory
dir.create(recursive = T,path = "./Heatmap/Hsapiens/")

#Merge to obtain the human ENSEMBL homologue.
hs_wt_phenotype_heatmap<-merge(anno_human_pre_stats_wt_phenotype,data.frame(assay(dds_rlog_wt_phenotype)),by.x="ensembl_gene_id",by.y="row.names")
head(hs_wt_phenotype_heatmap)
dim(hs_wt_phenotype_heatmap)
hs_wt_phenotype_heatmap<-hs_wt_phenotype_heatmap[-which(duplicated(hs_wt_phenotype_heatmap$hsapiens_homolog_ensembl_gene)),]


#Clean the list before it gets too longer.
#anno_human_wt_phenotype_stats<-anno_human_wt_phenotype_stats[-which(duplicated(anno_human_wt_phenotype_stats)),]

#Merge for the statiscal data.
hs_wt_phenotype_heatmap<-merge(anno_human_wt_phenotype_stats,hs_wt_phenotype_heatmap,by.x="hsapiens_homolog_ensembl_gene_id",by.y="hsapiens_homolog_ensembl_gene")


#Filter the duplicated ENSEMBL identifiers and keep only the significant data.
#Filter out duplicates
#hs_wt_phenotype_heatmap<-hs_wt_phenotype_heatmap[-which(duplicated(hs_wt_phenotype_heatmap$ENTREZID.x)),]

dim(hs_wt_phenotype_heatmap)
#Save the relevant data for the heatmap top 30 and down top 30.

hs_wt_phenotype_heatmap_topUp<-hs_wt_phenotype_heatmap%>%
  filter(log2FoldChange > 1)%>%
  arrange(desc(log2FoldChange,padj))%>%
  top_n(50,wt=-padj)


#Top down 30.
hs_wt_phenotype_heatmap_topDown<-hs_wt_phenotype_heatmap%>%
  filter(log2FoldChange < -1)%>%
  arrange(log2FoldChange,padj)%>%
  top_n(50,wt=-padj)


#Merge and save the data for morpheus.
hs_wt_phenotype_heatmap_tops<-bind_rows(hs_wt_phenotype_heatmap_topUp,hs_wt_phenotype_heatmap_topUp)
dim(hs_wt_phenotype_heatmap_tops)
head(hs_wt_phenotype_heatmap_tops)

#select and reorder the data 
head(sampleTable_wt_phenotype)
#colnames(hs_wt_phenotype_heatmap_tops[ncol(hs_wt_phenotype_heatmap_tops)])
hs_wt_phenotype_heatmap_tops<-hs_wt_phenotype_heatmap_tops%>%
  dplyr::select("SYMBOL.x","GENENAME.x",starts_with(substr(x=colnames(hs_wt_phenotype_heatmap_tops[ncol(hs_wt_phenotype_heatmap_tops)]),1,1)))



colnames(hs_wt_phenotype_heatmap_tops)
head(hs_wt_phenotype_heatmap_tops)
#Save the table ready to heatmap.
write.table(x = hs_wt_phenotype_heatmap_tops, "./Heatmap/Hsapiens/hs_wt_phenotype_heatmap.txt",sep = "\t",col.names = T,row.names = F,quote = F)



#####################################################
#                                                   #
#                   ZEBRAFISH                       #
#                                                   #
#                                                   #
#####################################################

#Directory
dir.create(recursive = T,path = "./Heatmap/Drerio/")
head(dds_rlog_wt_phenotype)
#Merge to obtain the human ENSEMBL homologue.
dr_wt_phenotype_heatmap<-merge(anno_zebrafish_wt_phenotype_stats,data.frame(assay(dds_rlog_wt_phenotype)),by.x="drerio_ensembl_gene_id",by.y="row.names")


#Filter the duplicated ENSEMBL identifiers and keep only the significant data.
#Filter out duplicates
#dr_wt_phenotype_heatmap<-dr_wt_phenotype_heatmap[-which(duplicated(dr_wt_phenotype_heatmap$drerio_ensembl_gene_id)),]
dim(dr_wt_phenotype_heatmap)
#Save the relevant data for the heatmap top 50 and down top 50.

dr_wt_phenotype_heatmap_topUp<-dr_wt_phenotype_heatmap%>%
  filter(log2FoldChange > 1)%>%
  arrange(desc(log2FoldChange),padj)%>%
  top_n(50,wt=-padj)
head(dr_wt_phenotype_heatmap_topUp)


#Top down 50.
#Use arrange normally as it does it by default ascending.
dr_wt_phenotype_heatmap_topDown<-dr_wt_phenotype_heatmap%>%
  filter(log2FoldChange < -1)%>%
  arrange(log2FoldChange,padj)%>%
  top_n(50,wt=-padj)

#Merge and save the data for morpheus.
dr_wt_phenotype_heatmap_tops<-bind_rows(dr_wt_phenotype_heatmap_topUp,dr_wt_phenotype_heatmap_topUp)
dim(dr_wt_phenotype_heatmap_tops)
colnames(dr_wt_phenotype_heatmap_tops)


#Observe and select MANUALLY the data for the heatmap.
dr_wt_phenotype_heatmap_tops<-dr_wt_phenotype_heatmap_tops%>%
  dplyr::select("SYMBOL","GENENAME",starts_with(substr(x=colnames(hs_wt_phenotype_heatmap_tops[ncol(hs_wt_phenotype_heatmap_tops)]),1,1)))
head(dr_wt_phenotype_heatmap_tops)
#Save the table ready to heatmap.
write.table(x = dr_wt_phenotype_heatmap_tops, "./Heatmap/Drerio/dr_wt_phenotype_heatmap.txt",sep = "\t",col.names = T,row.names = F,quote = F)






#####################################################
#                                                   #
#           Differentially expressed                #
#                   Genes                           #
#                    List                           #
#####################################################
print("Saving the DEG tables.")
cat("Saving the DEG tables.")
#Report a list with all the genes differentially expressed and annotated in Zebrafish.

#Create a folder for this task.
dir.create("./DEgenes")

write.table(x = anno_zebrafish_wt_phenotype_stats,"./DEgenes/all_annotated_DExpressed_genes_wt_phenotype_zebrafish.txt",sep="\t",col.names=T,row.names =F,quote=F)
write.table(anno_human_wt_phenotype_stats,"./DEgenes/all_annotated_DExpressed_genes_wt_phenotype_human.txt",sep = "\t",col.names = T,row.names = F,quote = T)
write.table(anno_mouse_wt_phenotype_stats,"./DEgenes/all_annotated_DExpressed_genes_wt_phenotype_mouse.txt",sep = "\t",col.names = T,row.names = F,quote = T)


#Filtered data table for absolute Log2FoldChange > 1 and padj < 0.05
write.table(anno_zebrafish_wt_phenotype_stats[which(abs(anno_zebrafish_wt_phenotype_stats$log2FoldChange) > 1 & anno_zebrafish_wt_phenotype_stats$padj < 0.05),],"./DEgenes/filtered_DExpressed_genes_zebrafish.txt",sep="\t",row.names=F,col.names=T,quote=F)
dim(anno_zebrafish_wt_phenotype_stats[which(abs(anno_zebrafish_wt_phenotype_stats$log2FoldChange>1) & anno_zebrafish_wt_phenotype_stats$padj<.05),])



#Write a report of the analysis, contrast used and what have been done.
description_text_wt_phenotype<-c("The following analysis was done for the comparison of WT cell in two conditios;
                                        
                                        Contrast used:  Double Diet vs Single Diet cells.
                                        
                                        The results were annotated and enriched using different approaches, using Danio and Human and Mouse homologues.
                                        
                                        For a more detailed understanding of the enrichments go to Analysis Guide.")


write.table(description_text_wt_phenotype,"./Description_text.txt",sep = "\t",col.names = F,row.names = F,quote = F)



#####################################################
#                                                   #
#                                                   #
#                                                   #
#               DISEASE ENRICHMENT                  #
#                DE    AND    ALL                   #
#             padj<0.05      LFC >1                 #
#####################################################

dir.create("./Disease")
cat("Disease enrichment")
print("Disease enrichment")
#List of significnat differentiually expressed genes.
length(hs_wt_phenotype_diff_expressed_genes)
length(hs_wt_phenotype_ranked_genes)

###########################
#                         #
#                         #
#     DE padj < 0.05      #
#                         #
###########################

#List of all genes entrezids with its respective log2FC
hs_wt_phenotype_ranked_genes
hs_wt_phenotype_DO_DE<-enrichDO(gene = hs_wt_phenotype_diff_expressed_genes,
                                       minGSSize = 5,
                                       ont = "DO",
                                       pAdjustMethod = "BH",
                                       readable = T,
                                       pvalueCutoff = 0.05)


###########################
#                         #
#                         #
#     ALL LFC > 1         #
#                         #
###########################
hs_wt_phenotype_DO_ALL<-enrichDO(gene = names(hs_wt_phenotype_ranked_genes[abs(hs_wt_phenotype_ranked_genes)>1]),
                                        minGSSize = 5,
                                        ont = "DO",
                                        pAdjustMethod = "BH",
                                        readable = T,
                                        pvalueCutoff = 0.05)

#see the only padj<0.05
#hs_wt_phenotype_disease_enrichment<-hs_wt_phenotype_disease_enrichment@result[hs_wt_phenotype_disease_enrichment@result$p.adjust<0.05,]
head(hs_wt_phenotype_DO_ALL@result,10)
head(hs_wt_phenotype_DO_DE@result,10)



#save the table of all the diseases and their significance.
write.table(hs_wt_phenotype_DO_DE@result,"./Disease/human_disease_enrichment_wt_phenotype_DE.tsv",sep = "\t",col.names = T,row.names = F)
write.table(hs_wt_phenotype_DO_ALL@result,"./Disease/human_disease_enrichment_wt_phenotype_ALL.tsv",sep = "\t",col.names = T,row.names = F)




#Save some plots
#Disease plot.
#BARPLOT
pdf("./Disease/hs_wt_phenotype_DO_DE_barplot.pdf", height = 20, width = 16)
barplot(hs_wt_phenotype_DO_DE,drop=T,showCategory = 20,title = "hs Disease enrichment wt_phenotype \n genelist Selected by padj<0.05 ")
dev.off()

#All
pdf("./Disease/hs_wt_phenotype_DO_ALL_barplot.pdf", height = 20, width = 16)
barplot(hs_wt_phenotype_DO_ALL,drop=T,showCategory = 20,title = "hs Disease enrichment wt_phenotype \n genelist Selected by LFC > 1 ")
dev.off()

#DOTPLOT of the Disease enrichment.
pdf("./Disease/hs_star_compare_DO_DE_dotplot.pdf", height = 20, width = 16)
dotplot(hs_wt_phenotype_DO_DE, showCategory=20, title="hs Disease enrichment wt_phenotype \n genelist Selected by padj<0.05 ")
dev.off()

#All
pdf("./Disease/hs_star_compare_DO_ALL_dotplot.pdf", height = 20, width = 16)
dotplot(hs_wt_phenotype_DO_ALL, showCategory=20, title="hs Disease enrichment wt_phenotype \n genelist Selected by LFC > 1 ")
dev.off()







##################################
#                                #
#       ENRICH DAVID             #
#                                #
#       USING DIFFERENT          #
#                                #
#          DATASETS              #
#                                #
##################################
#Add the zscore for the pathways and select the top expressed for up and down pathways, retruning a topUPandDOWN object, ready for plotting.
rm(symbols);rm(symbols_extracted);rm(symbols_extracted_double)
top_paths_DAVID<-function(DAVID_dataset,specie_anno_stats,n=15){
  print(ifelse(is.null(specie_anno_stats),break,"Specie selected!!"))
  print(ifelse(is.null(n),"Number of up and down rows not selected,using default!!",paste("n=",n)))
  if(length(DAVID_dataset@result$ID)==0){
    warning("ERROR → The pathway enrichment was not returning informtion with these parameters!!!\n")
    cat("ERROR → The pathway enrichment was not returning informtion with these parameters!!!\n" )
    break
  }else{
    cat("\nAdding zscore values to the tables...\n")
    #Create the zscore new column.
    symbols<-list()
    symbols_extracted<-list()
    symbols_extracted_double<-list()
    for (i in 1:length(DAVID_dataset@result$ID)){
      #Extract from the table all the SYMBOLS involved in the pathways from the gnees.
      target<-paste0("target",i)
      #print(target)
      
      #Clean the genelist of each row and obtain the genes.
      symbols[[target]]<-unique(unlist(strsplit(as.character(toupper(DAVID_dataset@result$geneID[i])),split = "\\/")))
      symbols_extracted[[target]]<-specie_anno_stats[specie_anno_stats$SYMBOL%in%symbols[[i]],]
      
      #Select only the SYMBOLS and the logFC for calculating the zscore.
      symbols_extracted_double[[target]]<-symbols_extracted[[i]][,c("log2FoldChange","SYMBOL")]
      
      #Get the zscore added
      DAVID_dataset@result$zscore[i]<-sum(as.numeric(symbols_extracted_double[[i]][,1]))/sqrt(length(symbols[[i]]))
    }
    cat("....\nSelecting top UP and DOWN values...\n......\n")
    topUP<-DAVID_dataset@result%>%
      filter(zscore>0)%>%
      top_n(n,wt = -p.adjust)
    topDOWN<-DAVID_dataset@result%>%
      filter(zscore<0)%>%
      top_n(n,wt = -p.adjust)
    topPATHS<-bind_rows(topUP,topDOWN)
    cat("\nCompleted!")
    return(topPATHS)}
}

#Using the genes for enrichment in DAVID data.
#?enrichDAVID


head(anno_human_wt_phenotype_stats)
#FDR significant genes.
length(hs_wt_phenotype_diff_expressed_genes)
head(hs_wt_phenotype_diff_expressed_genes)

#ALL genes.
names(hs_wt_phenotype_ranked_genes)
#All gnees with a logFC > 1. Otherwise the function does not return results..too much data.
names(hs_wt_phenotype_ranked_genes[abs(hs_wt_phenotype_ranked_genes)>1])

##################################
#                                #
#           KEGG                 #
#                                #
##################################
#                                #
#           HUMAN                #
#                                #
##################################
cat("DAVID KEGG Enrichments HUMAN")
print("DAVID KEGG Enrichments HUMAN")

#Plot the KEGG from DAVID using all and DE genes...
dir.create(recursive = T,"./DAVID/KEGG/Hsapiens/")

#use the entrez id siginficantly selected.


###########################
#                         #
#                         #
#     ALL LFC  1          #
#                         #
###########################
hs_wt_phenotype_david_KEGG_ALL <- enrichDAVID(gene = names(hs_wt_phenotype_ranked_genes[abs(hs_wt_phenotype_ranked_genes)>1]),
                                                     idType = "ENTREZ_GENE_ID",
                                                     minGSSize = 10,
                                                     maxGSSize = 500,
                                                     annotation = "KEGG_PATHWAY",
                                                     pvalueCutoff = 0.05,
                                                     pAdjustMethod = "BH",
                                                     david.user = "clusterProfiler@hku.hk")

###########################
#                         #
#                         #
#     DE padj < 0.05      #
#                         #
###########################

hs_wt_phenotype_david_KEGG_DE <- enrichDAVID(gene = hs_wt_phenotype_diff_expressed_genes,
                                                    idType = "ENTREZ_GENE_ID",
                                                    minGSSize = 10,
                                                    maxGSSize = 500,
                                                    annotation = "KEGG_PATHWAY",
                                                    pvalueCutoff = 0.05,
                                                    pAdjustMethod = "BH",
                                                    david.user = "clusterProfiler@hku.hk")

# annotation = "KEGG_PATHWAY",
#Translate the ENTREZID to SYMBOL
columns(org.Hs.eg.db)
hs_wt_phenotype_david_KEGG_ALL<-setReadable(hs_wt_phenotype_david_KEGG_ALL,OrgDb = org.Hs.eg.db,keyType = "ENTREZID")
hs_wt_phenotype_david_KEGG_DE<-setReadable(hs_wt_phenotype_david_KEGG_DE,OrgDb = org.Hs.eg.db,keyType = "ENTREZID")

#grep(pattern = "wt1",x = hs_wt_phenotype_david_KEGG@result$geneID,ignore.case = T)


pdf("./DAVID/KEGG/Hsapiens/hs_wt_phenotype_DAVID_KEGG_DE_barplot.pdf", height = 20, width = 16)
barplot(hs_wt_phenotype_david_KEGG_DE,drop=T,showCategory = 12,title = "hs DAVID KEGG enrichment wt_phenotype \n genelist Selected by padj<0.05 ")
dev.off()


#All
pdf("./DAVID/KEGG/Hsapiens/hs_wt_phenotype_DAVID_KEGG_ALL_barplot.pdf", height = 20, width = 16)
barplot(hs_wt_phenotype_david_KEGG_ALL,drop=T,showCategory = 12,title = "hs DAVID KEGG enrichment wt_phenotype \n genelist Selected by LFC > 1 ")
dev.off()

#Save the table
write.table(hs_wt_phenotype_david_KEGG_DE@result,"./DAVID/KEGG/Hsapiens/hs_wt_phenotype_DAVID_KEGG_DE.tsv",sep = "\t",row.names = F,col.names = T)
write.table(hs_wt_phenotype_david_KEGG_ALL@result,"./DAVID/KEGG/Hsapiens/hs_wt_phenotype_DAVID_KEGG_ALL.tsv",sep = "\t",row.names = F,col.names = T)







############################
#                          #
#       ZSCORED            #
#                          #
############################

#ALL
#Select the top pathways based on the ZSCORE added to the table.
hs_wt_phenotype_david_KEGG_ALL_topPathWays<-top_paths_DAVID(hs_wt_phenotype_david_KEGG_ALL,anno_human_wt_phenotype_stats)
head(hs_wt_phenotype_david_KEGG_ALL_topPathWays[order(hs_wt_phenotype_david_KEGG_ALL_topPathWays$zscore,decreasing = T),][,c(1,2,5,6,10)])
tail(hs_wt_phenotype_david_KEGG_ALL_topPathWays[order(hs_wt_phenotype_david_KEGG_ALL_topPathWays$zscore,decreasing = T),][,c(1,2,5,6,10)])

#Plot in a more visual way the pathways up and down regulated.
#Create the plot of the gsea file using ggplplot

pdf("./DAVID/KEGG/Hsapiens/hs_wt_phenotype_david_KEGG_ALL_topPathWays.pdf",width = 30,height = 25)
ggplot(hs_wt_phenotype_david_KEGG_ALL_topPathWays,aes(x=reorder(Description,zscore),y=zscore))+
  geom_bar(stat = "identity",fill="orange")+
  #geom_text(aes(label= round(p.adjust,3)),hjust=1.2,vjust=-.00001,size=3.5,color="white")+
  #labs(text="padjust")+
  geom_point(aes(label=p.adjust,size=rev(p.adjust)))+
  coord_flip()+
  labs(title = "Hs DAVID KEGG ALL genes LFC > 1",subtitle = sampleTable_wt_phenotype$treatment_type,size ="p-adjusted",x="Human DAVID:KEGG PATHWAYS",y="Normalized Enrichment Score")
dev.off()


#Save the filtered data table.
fwrite(hs_wt_phenotype_david_KEGG_ALL_topPathWays,"./DAVID/KEGG/Hsapiens/hs_wt_phenotype_david_KEGG_ALL_topPathWays.tsv",sep = "\t",col.names = T,row.names = F)


#DE
#Select the top pathways based on the ZSCORE added to the table.
hs_wt_phenotype_david_KEGG_DE_topPathWays<-top_paths_DAVID(hs_wt_phenotype_david_KEGG_DE,anno_human_wt_phenotype_stats)
head(hs_wt_phenotype_david_KEGG_DE_topPathWays[order(hs_wt_phenotype_david_KEGG_DE_topPathWays$zscore,decreasing = T),][,c(1,2,5,6,10)])
tail(hs_wt_phenotype_david_KEGG_DE_topPathWays[order(hs_wt_phenotype_david_KEGG_DE_topPathWays$zscore,decreasing = T),][,c(1,2,5,6,10)])

#Plot in a more visual way the pathways up and down regulated.
#Create the plot of the gsea file using ggplplot
pdf("./DAVID/KEGG/Hsapiens/hs_wt_phenotype_david_KEGG_DE_topPathWays.pdf",width = 30,height = 25)
ggplot(hs_wt_phenotype_david_KEGG_DE_topPathWays,aes(x=reorder(Description,zscore),y=zscore))+
  geom_bar(stat = "identity",fill="orange")+
  #geom_text(aes(label= round(p.adjust,3)),hjust=1.2,vjust=-.00001,size=3.5,color="white")+
  #labs(text="padjust")+
  geom_point(aes(label=p.adjust,size=rev(p.adjust)))+
  coord_flip()+
  labs(title = "Hs DAVID KEGG  DE genes selected padj < 0.05",subtitle = sampleTable_wt_phenotype$treatment_type,size ="p.adjusted",x="Human DAVID:KEGG PATHWAYS",y="Normalized Enrichment Score")
dev.off()


#Save the filtered data table.
fwrite(hs_wt_phenotype_david_KEGG_DE_topPathWays,"./DAVID/KEGG/Hsapiens/hs_wt_phenotype_david_KEGG_DE_topPathWays.tsv",sep = "\t",col.names = T,row.names = F)

##################################
#                                #
#           KEGG                 #
#                                #
##################################
#                                #
#           MOUSE                #
#                                #
##################################
dir.create(recursive = T,"./DAVID/KEGG/Mmusculus/")
cat("DAVID KEGG Enrichments MOUSE")
print("DAVID KEGG Enrichments MOUSE")


###########################
#                         #
#                         #
#     ALL LFC  1          #
#                         #
###########################
mm_wt_phenotype_david_KEGG_ALL <- enrichDAVID(gene = names(mm_wt_phenotype_ranked_genes[abs(mm_wt_phenotype_ranked_genes)>1]),
                                                     idType = "ENTREZ_GENE_ID",
                                                     minGSSize = 10,
                                                     maxGSSize = 500,
                                                     annotation = "KEGG_PATHWAY",
                                                     pvalueCutoff = 0.05,
                                                     pAdjustMethod = "BH",
                                                     david.user = "clusterProfiler@hku.hk")

###########################
#                         #
#                         #
#     DE padj < 0.05      #
#                         #
###########################

mm_wt_phenotype_david_KEGG_DE <- enrichDAVID(gene = mm_wt_phenotype_diff_expressed_genes,
                                                    idType = "ENTREZ_GENE_ID",
                                                    minGSSize = 10,
                                                    maxGSSize = 500,
                                                    annotation = "KEGG_PATHWAY",
                                                    pvalueCutoff = 0.05,
                                                    pAdjustMethod = "BH",
                                                    david.user = "clusterProfiler@hku.hk")

# annotation = "KEGG_PATHWAY",
#Translate the ENTREZID to SYMBOL
columns(org.Mm.eg.db)
mm_wt_phenotype_david_KEGG_ALL<-setReadable(mm_wt_phenotype_david_KEGG_ALL,OrgDb = org.Mm.eg.db,keyType = "ENTREZID")
mm_wt_phenotype_david_KEGG_DE<-setReadable(mm_wt_phenotype_david_KEGG_DE,OrgDb = org.Mm.eg.db,keyType = "ENTREZID")

#grep(pattern = "wt1",x = mm_wt_phenotype_david_KEGG@result$geneID,ignore.case = T)


pdf("./DAVID/KEGG/Mmusculus/mm_wt_phenotype_DAVID_KEGG_DE_barplot.pdf", height = 20, width = 16)
barplot(mm_wt_phenotype_david_KEGG_DE,drop=T,showCategory = 12,title = "mm DAVID KEGG enrichment wt_phenotype \n genelist Selected by padj<0.05 ")
dev.off()


#All
pdf("./DAVID/KEGG/Mmusculus/mm_wt_phenotype_DAVID_KEGG_ALL_barplot.pdf", height = 20, width = 16)
barplot(mm_wt_phenotype_david_KEGG_ALL,drop=T,showCategory = 12,title = "mm DAVID KEGG enrichment wt_phenotype \n genelist Selected by LFC > 1 ")
dev.off()

#Save the table
write.table(mm_wt_phenotype_david_KEGG_DE@result,"./DAVID/KEGG/Mmusculus/mm_wt_phenotype_DAVID_KEGG_DE.tsv",sep = "\t",row.names = F,col.names = T)
write.table(mm_wt_phenotype_david_KEGG_ALL@result,"./DAVID/KEGG/Mmusculus/mm_wt_phenotype_DAVID_KEGG_ALL.tsv",sep = "\t",row.names = F,col.names = T)







############################
#                          #
#       ZSCORED            #
#                          #
############################

#ALL
#Select the top pathways based on the ZSCORE added to the table.
mm_wt_phenotype_david_KEGG_ALL_topPathWays<-top_paths_DAVID(mm_wt_phenotype_david_KEGG_ALL,anno_mouse_wt_phenotype_stats,20)
head(mm_wt_phenotype_david_KEGG_ALL_topPathWays[order(mm_wt_phenotype_david_KEGG_ALL_topPathWays$zscore,decreasing = T),][,c(1,2,5,6,10)],20)
tail(mm_wt_phenotype_david_KEGG_ALL_topPathWays[order(mm_wt_phenotype_david_KEGG_ALL_topPathWays$zscore,decreasing = T),][,c(1,2,5,6,10)],20)

#Plot in a more visual way the pathways up and down regulated.
#Create the plot of the gsea file using ggplplot

pdf("./DAVID/KEGG/Mmusculus/mm_wt_phenotype_david_KEGG_ALL_topPathWays.pdf",width = 30,height = 25)
ggplot(mm_wt_phenotype_david_KEGG_ALL_topPathWays,aes(x=reorder(Description,zscore),y=zscore))+
  geom_bar(stat = "identity",fill="orange")+
  #geom_text(aes(label= round(p.adjust,3)),hjust=1.2,vjust=-.00001,size=3.5,color="white")+
  #labs(text="padjust")+
  geom_point(aes(label=p.adjust,size=rev(p.adjust)))+
  coord_flip()+
  labs(title = "mm DAVID KEGG ALL genes LFC > 1",subtitle = sampleTable_wt_phenotype$treatment_type,size ="p-adjusted",x="Mouse DAVID:KEGG PATHWAYS",y="Normalized Enrichment Score")
dev.off()


#Save the filtered data table.
fwrite(mm_wt_phenotype_david_KEGG_ALL_topPathWays,"./DAVID/KEGG/Mmusculus/mm_wt_phenotype_david_KEGG_ALL_topPathWays.tsv",sep = "\t",col.names = T,row.names = F)


#DE
#Select the top pathways based on the ZSCORE added to the table.
mm_wt_phenotype_david_KEGG_DE_topPathWays<-top_paths_DAVID(mm_wt_phenotype_david_KEGG_DE,anno_mouse_wt_phenotype_stats)
head(mm_wt_phenotype_david_KEGG_DE_topPathWays[order(mm_wt_phenotype_david_KEGG_DE_topPathWays$zscore,decreasing = T),][,c(1,2,5,6,10)])
tail(mm_wt_phenotype_david_KEGG_DE_topPathWays[order(mm_wt_phenotype_david_KEGG_DE_topPathWays$zscore,decreasing = T),][,c(1,2,5,6,10)])

#Plot in a more visual way the pathways up and down regulated.
#Create the plot of the gsea file using ggplplot
pdf("./DAVID/KEGG/Mmusculus/mm_wt_phenotype_david_KEGG_DE_topPathWays.pdf",width = 30,height = 25)
ggplot(mm_wt_phenotype_david_KEGG_DE_topPathWays,aes(x=reorder(Description,zscore),y=zscore))+
  geom_bar(stat = "identity",fill="orange")+
  #geom_text(aes(label= round(p.adjust,3)),hjust=1.2,vjust=-.00001,size=3.5,color="white")+
  #labs(text="padjust")+
  geom_point(aes(label=p.adjust,size=rev(p.adjust)))+
  coord_flip()+
  labs(title = "mm DAVID KEGG  DE genes selected padj < 0.05",subtitle = sampleTable_wt_phenotype$treatment_type,size ="p.adjusted",x="Mouse DAVID:KEGG PATHWAYS",y="Normalized Enrichment Score")
dev.off()


#Save the filtered data table.
fwrite(mm_wt_phenotype_david_KEGG_DE_topPathWays,"./DAVID/KEGG/Mmusculus/mm_wt_phenotype_david_KEGG_DE_topPathWays.tsv",sep = "\t",col.names = T,row.names = F)










##################################
#                                #
#              GO                #
#                                #
##################################
cat("DAVID GO Enrichments")
print("DAVID GO Enrichments")


##MOUSE
dir.create(recursive = T,"./DAVID/GO/Mmusculus/BP/")
#David enrichment with all genes with LFC > 1
names(mm_wt_phenotype_ranked_genes[abs(mm_wt_phenotype_ranked_genes)> 1])


#David enrichment with DE genes and the significative ones, al lthe padj < 0.05 
head(mm_wt_phenotype_diff_expressed_genes)
head(mm_wt_phenotype_ranked_genes)

###########################
#                         #
#                         #
#     DE padj < 0.05      #
#                         #
###########################
mm_wt_phenotype_david_GO_DE<- enrichDAVID(gene = mm_wt_phenotype_diff_expressed_genes,
                                   idType = "ENTREZ_GENE_ID",
                                   minGSSize = 10,
                                   maxGSSize = 500,
                                   annotation = "GOTERM_BP_ALL",
                                   pvalueCutoff = 0.05,
                                   pAdjustMethod = "BH",
                                   david.user = "clusterProfiler@hku.hk")
mm_wt_phenotype_david_GO_DE<-setReadable(mm_wt_phenotype_david_GO_DE,OrgDb = org.Mm.eg.db,keyType = "ENTREZID")
head(mm_wt_phenotype_david_GO_DE@result)
dim(mm_wt_phenotype_david_GO_DE@result)

# annotation = "GOTERM_BP_ALL","GOTERM_BP_FAT"


###########################
#                         #
#                         #
#     ALL LFC  1          #
#                         #
###########################

mm_wt_phenotype_david_GO_ALL <- enrichDAVID(gene = names(mm_wt_phenotype_ranked_genes[abs(mm_wt_phenotype_ranked_genes)>1]),
                                                   idType = "ENTREZ_GENE_ID",
                                                   minGSSize = 10,
                                                   maxGSSize = 500,
                                                   annotation = "GOTERM_BP_ALL",
                                                   pvalueCutoff = 0.05,
                                                   pAdjustMethod = "BH",
                                                   david.user = "clusterProfiler@hku.hk")


#Translate the ENTREZIDS to SYMBOLS
columns(org.Mm.eg.db)
mm_wt_phenotype_david_GO_DE<-setReadable(mm_wt_phenotype_david_GO_DE,OrgDb = org.Mm.eg.db,keyType = "ENTREZID")
mm_wt_phenotype_david_GO_ALL<-setReadable(mm_wt_phenotype_david_GO_ALL,OrgDb = org.Mm.eg.db,keyType = "ENTREZID")

head(mm_wt_phenotype_david_GO_DE@result,30)
head(mm_wt_phenotype_david_GO_ALL@result,30)


pdf("./DAVID/GO/Mmusculus/BP/mm_wt_phenotype_david_GO_ALL_barplot.pdf", height = 20, width = 16)
barplot(mm_wt_phenotype_david_GO_ALL,drop=T,showCategory = 12,title = "mm DAVID GO enrichment wt_phenotype \n genelist Selected by LFC > 1 ")
dev.off()
pdf("./DAVID/GO/Mmusculus/BP/mm_wt_phenotype_david_GO_DE_barplot.pdf", height = 20, width = 16)
barplot(mm_wt_phenotype_david_GO_DE,drop=T,showCategory = 12,title = "mm DAVID GO enrichment wt_phenotype \n genelist Selected by padj<0.05")
dev.off()



#Dotplot the reuslts.
pdf("./DAVID/GO/Mmusculus/BP/mm_wt_phenotype_DAVID_GO_DE_dotplot.pdf", height = 20, width = 16)
dotplot(mm_wt_phenotype_david_GO_DE, showCategory=12, title="mm DAVID GO enrichment wt_phenotype \n genelist Selected by padj<0.05 ")
dev.off()
pdf("./DAVID/GO/Mmusculus/BP/mm_wt_phenotype_DAVID_GO_ALL_dotplot.pdf", height = 20, width = 16)
dotplot(mm_wt_phenotype_david_GO_ALL, showCategory=12, title="mm DAVID GO enrichment wt_phenotype \n genelist Selected by LFC > 1")
dev.off()

#Save the data
write.table(mm_wt_phenotype_david_GO_ALL@result,"./DAVID/GO/Mmusculus/BP/mm_wt_phenotype_david_GO_ALL.tsv",sep = "\t",row.names = F,col.names = T)
write.table(mm_wt_phenotype_david_GO_DE@result,"./DAVID/GO/Mmusculus/BP/mm_wt_phenotype_david_GO_DE.tsv",sep = "\t",row.names = F,col.names = T)



############################
#                          #
#       ZSCORED            #
#                          #
############################

#ALL
#Select the top pathways based on the ZSCORE added to the table.


mm_wt_phenotype_david_GO_ALL_topPathWays<-top_paths_DAVID(mm_wt_phenotype_david_GO_ALL,anno_mouse_wt_phenotype_stats,20)
head(mm_wt_phenotype_david_GO_ALL_topPathWays[order(mm_wt_phenotype_david_GO_ALL_topPathWays$zscore,decreasing = T),][,c(1,2,5,6,10)],20)
tail(mm_wt_phenotype_david_GO_ALL_topPathWays[order(mm_wt_phenotype_david_GO_ALL_topPathWays$zscore,decreasing = T),][,c(1,2,5,6,10)],20)

#Plot in a more visual way the pathways up and down regulated.
#Create the plot of the gsea file using ggplplot

pdf("./DAVID/GO/Mmusculus/BP/mm_wt_phenotype_david_GO_ALL_topPathWays.pdf",width = 30,height = 25)
ggplot(mm_wt_phenotype_david_GO_ALL_topPathWays,aes(x=reorder(Description,zscore),y=zscore))+
  geom_bar(stat = "identity",fill="darkorchid2")+
  #geom_text(aes(label= round(p.adjust,3)),hjust=1.2,vjust=-.00001,size=3.5,color="white")+
  #labs(text="padjust")+
  geom_point(aes(label=p.adjust,size=rev(p.adjust)))+
  coord_flip()+
  labs(title = "mm DAVID GO ALL genes LFC > 1",subtitle = sampleTable_wt_phenotype$treatment_type,size ="p-adjusted",x="Mouse DAVID:GO Functions",y="Normalized Enrichment Score")
dev.off()


#Save the filtered data table.
fwrite(mm_wt_phenotype_david_GO_ALL_topPathWays,"./DAVID/GO/Mmusculus/BP/mm_wt_phenotype_david_GO_ALL_topPathWays.tsv",sep = "\t",col.names = T,row.names = F)


#DE
#Select the top pathways based on the ZSCORE added to the table.


mm_wt_phenotype_david_GO_DE_topPathWays<-top_paths_DAVID(mm_wt_phenotype_david_GO_DE,anno_mouse_wt_phenotype_stats,20)
head(mm_wt_phenotype_david_GO_DE_topPathWays[order(mm_wt_phenotype_david_GO_DE_topPathWays$zscore,decreasing = T),][,c(1,2,5,6,10)],20)
tail(mm_wt_phenotype_david_GO_DE_topPathWays[order(mm_wt_phenotype_david_GO_DE_topPathWays$zscore,decreasing = T),][,c(1,2,5,6,10)],20)


pdf("./DAVID/GO/Mmusculus/BP/mm_wt_phenotype_david_GO_DE_topPathWays.pdf",width = 30,height = 25)
ggplot(mm_wt_phenotype_david_GO_DE_topPathWays,aes(x=reorder(Description,zscore),y=zscore))+
  geom_bar(stat = "identity",fill="darkorchid2")+
  #geom_text(aes(label= round(p.adjust,3)),hjust=1.2,vjust=-.00001,size=3.5,color="white")+
  #labs(text="padjust")+
  geom_point(aes(label=p.adjust,size=rev(p.adjust)))+
  coord_flip()+
  labs(title = "mm DAVID GO DE genes padj < 0.05",subtitle = sampleTable_wt_phenotype$treatment_type,size ="p.adjusted",x="Mouse DAVID:GO Functions",y="Normalized Enrichment Score")
dev.off()


#Save the filtered data table.
fwrite(mm_wt_phenotype_david_GO_DE_topPathWays,"./DAVID/GO/Mmusculus/BP/mm_wt_phenotype_david_GO_DE_topPathWays.tsv",sep = "\t",col.names = T,row.names = F)



################################################
#                                              #
#                                              #
##          HUMAN GO ALL DAVID                 #
#               padj<0.05                      #  
#               and LFC > 1                    # 
################################################


dir.create(recursive = T,"./DAVID/GO/Hsapiens/BP")

###########################
#                         #
#                         #
#     DE padj < 0.05      #
#                         #
###########################

hs_wt_phenotype_david_GO_DE <- enrichDAVID(gene = hs_wt_phenotype_diff_expressed_genes,
                                                  idType = "ENTREZ_GENE_ID",
                                                  minGSSize = 10,
                                                  maxGSSize = 500,
                                                  annotation = "GOTERM_BP_ALL",
                                                  pvalueCutoff = 0.05,
                                                  pAdjustMethod = "BH",
                                                  david.user = "clusterProfiler@hku.hk")
# annotation = "GOTERM_BP_ALL","GOTERM_BP_FAT"

###########################
#                         #
#                         #
#     ALL LFC  1          #
#                         #
###########################

hs_wt_phenotype_david_GO_ALL <- enrichDAVID(gene = unique(names(hs_wt_phenotype_ranked_genes[abs(hs_wt_phenotype_ranked_genes)>1])),
                                                   idType = "ENTREZ_GENE_ID",
                                                   minGSSize = 10,
                                                   maxGSSize = 500,
                                                   annotation = "GOTERM_BP_ALL",
                                                   pvalueCutoff = 0.05,
                                                   pAdjustMethod = "BH",
                                                   david.user = "clusterProfiler@hku.hk")



#Translate the entezids to symbols.
columns(org.Hs.eg.db)
hs_wt_phenotype_david_GO_DE<-setReadable(hs_wt_phenotype_david_GO_DE,OrgDb = org.Hs.eg.db,keyType = "ENTREZID")
hs_wt_phenotype_david_GO_ALL<-setReadable(hs_wt_phenotype_david_GO_ALL,OrgDb = org.Hs.eg.db,keyType = "ENTREZID")



head(hs_wt_phenotype_david_GO_DE@result,30)
head(hs_wt_phenotype_david_GO_ALL@result,30)


##Search interested genes...
#grep(pattern = "sox",x = hs_wt_phenotype_david_GO_DE@result$geneID,ignore.case = T)
#grep(pattern = "sox",x = hs_wt_phenotype_david_GO_ALL@result$geneID,ignore.case = T)

#paint the plot

#Barplot the results.
pdf("./DAVID/GO/Hsapiens/BP/hs_wt_phenotype_david_GO_DE_barplot.pdf", height = 20, width = 16)
barplot(hs_wt_phenotype_david_GO_DE,drop=T,showCategory = 20,title = "hs DAVID GO enrichment wt_phenotype \n genelist Selected by padj<0.05 ")
dev.off()

pdf("./DAVID/GO/Hsapiens/BP/hs_wt_phenotype_david_GO_ALL_barplot.pdf", height = 20, width = 16)
barplot(hs_wt_phenotype_david_GO_ALL,drop=T,showCategory = 20,title = "hs DAVID GO enrichment wt_phenotype \n genelist Selected by LFC > 1 ")
dev.off()



#Dotplot the reuslts.
pdf("./DAVID/GO/Hsapiens/BP/hs_wt_phenotype_DAVID_GO_DE_dotplot.pdf", height = 20, width = 16)
dotplot(hs_wt_phenotype_david_GO_DE, showCategory=12, title="hs DAVID GO enrichment wt_phenotype \n genelist Selected by padj<0.05 ")
dev.off()

pdf("./DAVID/GO/Hsapiens/BP/hs_wt_phenotype_DAVID_GO_ALL_dotplot.pdf", height = 20, width = 16)
dotplot(hs_wt_phenotype_david_GO_ALL, showCategory=12, title="hs DAVID GO enrichment wt_phenotype \n genelist Selected by LFC > 1 ")
dev.off()

#Save the table
write.table(hs_wt_phenotype_david_GO_DE@result,"./DAVID/GO/Hsapiens/BP/hs_wt_phenotype_DAVID_GO_DE.tsv",sep = "\t",row.names = F,col.names = T)

write.table(hs_wt_phenotype_david_GO_ALL@result,"./DAVID/GO/Hsapiens/BP/hs_wt_phenotype_DAVID_GO_ALL.tsv",sep = "\t",row.names = F,col.names = T)




############################
#                          #
#       ZSCORED            #
#                          #
############################

#ALL
#Select the top pathways based on the ZSCORE added to the table.

hs_wt_phenotype_david_GO_ALL_topPathWays<-top_paths_DAVID(hs_wt_phenotype_david_GO_ALL,anno_human_wt_phenotype_stats,30)
head(hs_wt_phenotype_david_GO_ALL_topPathWays[,c(1,2,5,6,10)]);tail(hs_wt_phenotype_david_GO_ALL_topPathWays[,c(1,2,5,6,10)])


#Plot in a more visual way the pathways up and down regulated.
#Create the plot of the gsea file using ggplplot

pdf("./DAVID/GO/Hsapiens/BP/hs_wt_phenotype_david_GO_ALL_topPathWays.pdf",width = 30,height = 25)
ggplot(hs_wt_phenotype_david_GO_ALL_topPathWays,aes(x=reorder(Description,zscore),y=zscore))+
  geom_bar(stat = "identity",fill="darkorchid2")+
  #geom_text(aes(label= round(p.adjust,3)),hjust=1.2,vjust=-.00001,size=3.5,color="white")+
  #labs(text="padjust")+
  geom_point(aes(label=p.adjust,size=rev(p.adjust)))+
  coord_flip()+
  labs(title = "Hs DAVID GO ALL genes LFC > 1",subtitle = sampleTable_wt_phenotype$treatment_type,size ="p-adjusted",x="Human DAVID:GO Functions",y="Normalized Enrichment Score")
dev.off()


#Save the filtered data table.
fwrite(hs_wt_phenotype_david_GO_ALL_topPathWays,"./DAVID/GO/Hsapiens/BP/hs_wt_phenotype_david_GO_ALL_topPathWays.tsv",sep = "\t",col.names = T,row.names = F)






#DE
#Select the top pathways based on the ZSCORE added to the table.

hs_wt_phenotype_david_GO_DE_topPathWays<-top_paths_DAVID(hs_wt_phenotype_david_GO_DE,anno_human_wt_phenotype_stats,30)
head(hs_wt_phenotype_david_GO_DE_topPathWays[,c(1,2,5,6,10)]);tail(hs_wt_phenotype_david_GO_DE_topPathWays[,c(1,2,5,6,10)])

#Plot in a more visual way the pathways up and down regulated.
#Create the plot of the gsea file using ggplplot

pdf("./DAVID/GO/Hsapiens/BP/hs_wt_phenotype_david_GO_DE_topPathWays.pdf",width = 30,height = 25)
ggplot(hs_wt_phenotype_david_GO_DE_topPathWays,aes(x=reorder(Description,zscore),y=zscore))+
  geom_bar(stat = "identity",fill="darkorchid2")+
  #geom_text(aes(label= round(p.adjust,3)),hjust=1.2,vjust=-.00001,size=3.5,color="white")+
  #labs(text="padjust")+
  geom_point(aes(label=p.adjust,size=rev(p.adjust)))+
  coord_flip()+
  labs(title = "Hs DAVID GO DE genes padj< 0.05" ,subtitle = sampleTable_wt_phenotype$treatment_type,size ="p-adjusted",x="Human DAVID:GO Functions",y="Normalized Enrichment Score")
dev.off()


#Save the filtered data table.
fwrite(hs_wt_phenotype_david_GO_DE_topPathWays,"./DAVID/GO/Hsapiens/BP/hs_wt_phenotype_david_GO_DE_topPathWays.tsv",sep = "\t",col.names = T,row.names = F)



#data("annotationSummary2")
#data("annotationSummary1")
#(annotationSummary2)
#annotationSummary1




##################################
#                                #
#              GO MF             #
#                                #
##################################
#                                #
#           Human                #
#                                #
##################################
dir.create("./DAVID/GO/Hsapiens/MF/")



###########################
#                         #
#                         #
#     DE padj < 0.05      #
#                         #
###########################
hs_wt_phenotype_david_GO_MF_DE <- enrichDAVID(gene = hs_wt_phenotype_diff_expressed_genes,
                                                     idType = "ENTREZ_GENE_ID",
                                                     minGSSize = 10,
                                                     maxGSSize = 500,
                                                     annotation = "GOTERM_MF_ALL",
                                                     pvalueCutoff = 0.05,
                                                     pAdjustMethod = "BH",
                                                     david.user = "clusterProfiler@hku.hk")

###########################
#                         #
#                         #
#     ALL LFC  1          #
#                         #
###########################
hs_wt_phenotype_david_GO_MF_ALL <- enrichDAVID(gene = names(hs_wt_phenotype_ranked_genes[abs(hs_wt_phenotype_ranked_genes)>1]),
                                                      idType = "ENTREZ_GENE_ID",
                                                      minGSSize = 10,
                                                      maxGSSize = 500,
                                                      annotation = "GOTERM_MF_ALL",
                                                      pvalueCutoff = 0.05,
                                                      pAdjustMethod = "BH",
                                                      david.user = "clusterProfiler@hku.hk")

# annotation = "KEGG_PATHWAY",


#Translate the names to entrezIDS.
hs_wt_phenotype_david_GO_MF_DE<-setReadable(hs_wt_phenotype_david_GO_MF_DE,OrgDb = org.Hs.eg.db,keyType = "ENTREZID")
hs_wt_phenotype_david_GO_MF_ALL<-setReadable(hs_wt_phenotype_david_GO_MF_ALL,OrgDb = org.Hs.eg.db,keyType = "ENTREZID")





#Barplot the results.
pdf("./DAVID/GO/Hsapiens/MF/hs_wt_phenotype_david_GO_MF_DE_barplot.pdf", height = 20, width = 16)
barplot(hs_wt_phenotype_david_GO_MF_DE,drop=T,showCategory = 20,title = "hs DAVID  GO MF enrichment wt_phenotype \n genelist Selected by padj<0.05 ")
dev.off()

#ALL
pdf("./DAVID/GO/Hsapiens/MF/hs_wt_phenotype_david_GO_MF_ALL_barplot.pdf", height = 20, width = 16)
barplot(hs_wt_phenotype_david_GO_MF_ALL,drop=T,showCategory = 20,title = "hs DAVID  GO MF enrichment wt_phenotype \n genelist Selected by LFC > 1 ")
dev.off()




#Dotplot the reuslts.
pdf("./DAVID/GO/Hsapiens/MF/hs_wt_phenotype_DAVID_GO_MF_DE_dotplot.pdf", height = 20, width = 16)
dotplot(hs_wt_phenotype_david_GO_MF_DE, showCategory=20, title="hs DAVID GO MF enrichment wt_phenotype \n genelist Selected by padj<0.05 ")
dev.off()


#All
pdf("./DAVID/GO/Hsapiens/MF/hs_wt_phenotype_DAVID_GO_MF_all_dotplot.pdf", height = 20, width = 16)
dotplot(hs_wt_phenotype_david_GO_MF_ALL, showCategory=20, title="hs DAVID GO MF enrichment wt_phenotype \n genelist Selected by LFC>1 ")
dev.off()


#Save the table
write.table(hs_wt_phenotype_david_GO_MF_DE@result,"./DAVID/GO/Hsapiens/MF/hs_wt_phenotype_DAVID_GO_MF_DE.tsv",sep = "\t",row.names = F,col.names = T)
write.table(hs_wt_phenotype_david_GO_MF_ALL@result,"./DAVID/GO/Hsapiens/MF/hs_wt_phenotype_DAVID_GO_MF_ALL.tsv",sep = "\t",row.names = F,col.names = T)



############################
#                          #
#       ZSCORED            #
#                          #
############################

#ALL
#Select the top pathways based on the ZSCORE added to the table.

hs_wt_phenotype_david_GO_MF_ALL_topPathWays<-top_paths_DAVID(hs_wt_phenotype_david_GO_MF_ALL,anno_human_wt_phenotype_stats,20)
head(hs_wt_phenotype_david_GO_MF_ALL_topPathWays[,c(1,2,5,6,10)]);tail(hs_wt_phenotype_david_GO_MF_ALL_topPathWays[,c(1,2,5,6,10)])


#Plot in a more visual way the pathways up and down regulated.
#Create the plot of the gsea file using ggplplot

pdf("./DAVID/GO/Hsapiens/MF/hs_wt_phenotype_david_GO_MF_ALL_topPathWays.pdf",width = 30,height = 25)
ggplot(hs_wt_phenotype_david_GO_MF_ALL_topPathWays,aes(x=reorder(Description,zscore),y=zscore))+
  geom_bar(stat = "identity",fill="darkorchid2")+
  #geom_text(aes(label= round(p.adjust,3)),hjust=1.2,vjust=-.00001,size=3.5,color="white")+
  #labs(text="padjust")+
  geom_point(aes(label=p.adjust,size=rev(p.adjust)))+
  coord_flip()+
  labs(title = "Hs DAVID GO MF ALL genes LFC > 1",subtitle = sampleTable_wt_phenotype$treatment_type,size ="p-adjusted",x="Human DAVID:GO Molecular Functions",y="Normalized Enrichment Score")
dev.off()


#Save the filtered data table.
fwrite(hs_wt_phenotype_david_GO_MF_ALL_topPathWays,"./DAVID/GO/Hsapiens/MF/hs_wt_phenotype_david_GO_MF_ALL_topPathWays.tsv",sep = "\t",col.names = T,row.names = F)

#DE
#Select the top pathways based on the ZSCORE added to the table.

hs_wt_phenotype_david_GO_MF_DE_topPathWays<-top_paths_DAVID(hs_wt_phenotype_david_GO_MF_DE,anno_human_wt_phenotype_stats,20)
head(hs_wt_phenotype_david_GO_MF_DE_topPathWays[,c(1,2,5,6,10)]);tail(hs_wt_phenotype_david_GO_MF_DE_topPathWays[,c(1,2,5,6,10)])


#Plot in a more visual way the pathways up and down regulated.
#Create the plot of the gsea file using ggplplot

pdf("./DAVID/GO/Hsapiens/MF/hs_wt_phenotype_david_GO_MF_DE_topPathWays.pdf",width = 30,height = 25)
ggplot(hs_wt_phenotype_david_GO_MF_DE_topPathWays,aes(x=reorder(Description,zscore),y=zscore))+
  geom_bar(stat = "identity",fill="darkorchid2")+
  #geom_text(aes(label= round(p.adjust,3)),hjust=1.2,vjust=-.00001,size=3.5,color="white")+
  #labs(text="padjust")+
  geom_point(aes(label=p.adjust,size=rev(p.adjust)))+
  coord_flip()+
  labs(title = "Hs DAVID GO MF DE genes padj<0.05",subtitle = sampleTable_wt_phenotype$treatment_type,size ="p.adjusted",x="Human DAVID:GO Molecular Functions",y="Normalized Enrichment Score")
dev.off()


#Save the filtered data table.
fwrite(hs_wt_phenotype_david_GO_MF_DE_topPathWays,"./DAVID/GO/Hsapiens/MF/hs_wt_phenotype_david_GO_MF_DE_topPathWays.tsv",sep = "\t",col.names = T,row.names = F)


##################################
#                                #
#           Mouse                #
#                                #
##################################
dir.create("./DAVID/GO/Mmusculus/MF/")



###########################
#                         #
#                         #
#     DE padj < 0.05      #
#                         #
###########################
mm_wt_phenotype_david_GO_MF_DE <- enrichDAVID(gene = mm_wt_phenotype_diff_expressed_genes,
                                                     idType = "ENTREZ_GENE_ID",
                                                     minGSSize = 10,
                                                     maxGSSize = 500,
                                                     annotation = "GOTERM_MF_ALL",
                                                     pvalueCutoff = 0.05,
                                                     pAdjustMethod = "BH",
                                                     david.user = "clusterProfiler@hku.hk")

###########################
#                         #
#                         #
#     ALL LFC  1          #
#                         #
###########################
mm_wt_phenotype_david_GO_MF_ALL <- enrichDAVID(gene = names(mm_wt_phenotype_ranked_genes[abs(mm_wt_phenotype_ranked_genes)>1]),
                                                      idType = "ENTREZ_GENE_ID",
                                                      minGSSize = 10,
                                                      maxGSSize = 500,
                                                      annotation = "GOTERM_MF_ALL",
                                                      pvalueCutoff = 0.05,
                                                      pAdjustMethod = "BH",
                                                      david.user = "clusterProfiler@hku.hk")

# annotation = "KEGG_PATHWAY",


#Translate the names to entrezIDS.
mm_wt_phenotype_david_GO_MF_DE<-setReadable(mm_wt_phenotype_david_GO_MF_DE,OrgDb = org.Mm.eg.db,keyType = "ENTREZID")
mm_wt_phenotype_david_GO_MF_ALL<-setReadable(mm_wt_phenotype_david_GO_MF_ALL,OrgDb = org.Mm.eg.db,keyType = "ENTREZID")



#Barplot the results.
pdf("./DAVID/GO/Mmusculus/MF/mm_wt_phenotype_david_GO_MF_DE_barplot.pdf", height = 20, width = 16)
barplot(mm_wt_phenotype_david_GO_MF_DE,drop=T,showCategory = 20,title = "mm DAVID  GO MF enrichment wt_phenotype \n genelist Selected by padj<0.05 ")
dev.off()

#ALL
pdf("./DAVID/GO/Mmusculus/MF/mm_wt_phenotype_david_GO_MF_ALL_barplot.pdf", height = 20, width = 16)
barplot(mm_wt_phenotype_david_GO_MF_ALL,drop=T,showCategory = 20,title = "mm DAVID  GO MF enrichment wt_phenotype \n genelist Selected by LFC > 1 ")
dev.off()




#Dotplot the reuslts.
pdf("./DAVID/GO/Mmusculus/MF/mm_wt_phenotype_DAVID_GO_MF_DE_dotplot.pdf", height = 20, width = 16)
dotplot(mm_wt_phenotype_david_GO_MF_DE, showCategory=20, title="mm DAVID GO MF enrichment wt_phenotype \n genelist Selected by padj<0.05 ")
dev.off()


#All
pdf("./DAVID/GO/Mmusculus/MF/mm_wt_phenotype_DAVID_GO_MF_all_dotplot.pdf", height = 20, width = 16)
dotplot(mm_wt_phenotype_david_GO_MF_ALL, showCategory=20, title="mm DAVID GO MF enrichment wt_phenotype \n genelist Selected by LFC>1 ")
dev.off()


#Save the table
write.table(mm_wt_phenotype_david_GO_MF_DE@result,"./DAVID/GO/Mmusculus/MF/mm_wt_phenotype_DAVID_GO_MF_DE.tsv",sep = "\t",row.names = F,col.names = T)
write.table(mm_wt_phenotype_david_GO_MF_ALL@result,"./DAVID/GO/Mmusculus/MF/mm_wt_phenotype_DAVID_GO_MF_ALL.tsv",sep = "\t",row.names = F,col.names = T)



############################
#                          #
#       ZSCORED            #
#                          #
############################

#ALL
#Select the top pathways based on the ZSCORE added to the table.

mm_wt_phenotype_david_GO_MF_ALL_topPathWays<-top_paths_DAVID(mm_wt_phenotype_david_GO_MF_ALL,anno_mouse_wt_phenotype_stats,20)
head(mm_wt_phenotype_david_GO_MF_ALL_topPathWays[,c(1,2,5,6,10)]);tail(mm_wt_phenotype_david_GO_MF_ALL_topPathWays[,c(1,2,5,6,10)])


#Plot in a more visual way the pathways up and down regulated.
#Create the plot of the gsea file using ggplplot

pdf("./DAVID/GO/Mmusculus/MF/mm_wt_phenotype_david_GO_MF_ALL_topPathWays.pdf",width = 30,height = 25)
ggplot(mm_wt_phenotype_david_GO_MF_ALL_topPathWays,aes(x=reorder(Description,zscore),y=zscore))+
  geom_bar(stat = "identity",fill="darkorchid2")+
  #geom_text(aes(label= round(p.adjust,3)),hjust=1.2,vjust=-.00001,size=3.5,color="white")+
  #labs(text="padjust")+
  geom_point(aes(label=p.adjust,size=rev(p.adjust)))+
  coord_flip()+
  labs(title = "Mm DAVID GO MF ALL genes LFC > 1",subtitle = sampleTable_wt_phenotype$treatment_type,size ="p-adjusted",x="Mouse DAVID:GO Molecular Functions",y="Normalized Enrichment Score")
dev.off()


#Save the filtered data table.
fwrite(mm_wt_phenotype_david_GO_MF_ALL_topPathWays,"./DAVID/GO/Mmusculus/MF/mm_wt_phenotype_david_GO_MF_ALL_topPathWays.tsv",sep = "\t",col.names = T,row.names = F)



mm_wt_phenotype_david_GO_CC_DE@result
#DE
#Select the top pathways based on the ZSCORE added to the table.
mm_wt_phenotype_david_GO_MF_DE_topPathWays<-top_paths_DAVID(mm_wt_phenotype_david_GO_MF_DE,anno_mouse_wt_phenotype_stats,30)
head(mm_wt_phenotype_david_GO_MF_DE_topPathWays[,c(1,2,5,6,10)]);tail(mm_wt_phenotype_david_GO_MF_DE_topPathWays[,c(1,2,5,6,10)])


#Plot in a more visual way the pathways up and down regulated.
#Create the plot of the gsea file using ggplplot

pdf("./DAVID/GO/Mmusculus/MF/mm_wt_phenotype_david_GO_MF_DE_topPathWays.pdf",width = 30,height = 25)
ggplot(mm_wt_phenotype_david_GO_MF_DE_topPathWays,aes(x=reorder(Description,zscore),y=zscore))+
  geom_bar(stat = "identity",fill="darkorchid2")+
  #geom_text(aes(label= round(p.adjust,3)),hjust=1.2,vjust=-.00001,size=3.5,color="white")+
  #labs(text="padjust")+
  geom_point(aes(label=p.adjust,size=rev(p.adjust)))+
  coord_flip()+
  labs(title = "Mm DAVID GO MF DE genes padj<0.05",subtitle = sampleTable_wt_phenotype$treatment_type,size ="p.adjusted",x="Mouse DAVID:GO Molecular Functions",y="Normalized Enrichment Score")
dev.off()


#Save the filtered data table.
fwrite(mm_wt_phenotype_david_GO_MF_DE_topPathWays,"./DAVID/GO/Mmusculus/MF/mm_wt_phenotype_david_GO_MF_DE_topPathWays.tsv",sep = "\t",col.names = T,row.names = F)


##################################
#                                #
#              GO CC             #
#                                #
##################################
#                                #
#           Human                #
#                                #
##################################
dir.create("./DAVID/GO/Hsapiens/CC/")

###########################
#                         #
#                         #
#     DE padj < 0.05      #
#                         #
###########################

hs_wt_phenotype_david_GO_CC_DE <- enrichDAVID(gene = hs_wt_phenotype_diff_expressed_genes,
                                                     idType = "ENTREZ_GENE_ID",
                                                     minGSSize = 10,
                                                     maxGSSize = 500,
                                                     annotation = "GOTERM_CC_ALL",
                                                     pvalueCutoff = 0.05,
                                                     pAdjustMethod = "BH",
                                                     david.user = "clusterProfiler@hku.hk")

###########################
#                         #
#                         #
#     ALL LFC  1          #
#                         #
###########################
hs_wt_phenotype_david_GO_CC_ALL <- enrichDAVID(gene = names(hs_wt_phenotype_ranked_genes[abs(hs_wt_phenotype_ranked_genes)>1]),
                                                      idType = "ENTREZ_GENE_ID",
                                                      minGSSize = 10,
                                                      maxGSSize = 500,
                                                      annotation = "GOTERM_CC_ALL",
                                                      pvalueCutoff = 0.05,
                                                      pAdjustMethod = "BH",
                                                      david.user = "clusterProfiler@hku.hk")

#Translate the names to entrezIDS.
hs_wt_phenotype_david_GO_CC_DE<-setReadable(hs_wt_phenotype_david_GO_CC_DE,OrgDb = org.Hs.eg.db,keyType = "ENTREZID")
hs_wt_phenotype_david_GO_CC_ALL<-setReadable(hs_wt_phenotype_david_GO_CC_ALL,OrgDb = org.Hs.eg.db,keyType = "ENTREZID")

#Barplot the results.
pdf("./DAVID/GO/Hsapiens/CC/hs_wt_phenotype_david_GO_CC_DE_barplot.pdf", height = 20, width = 16)
barplot(hs_wt_phenotype_david_GO_CC_DE,drop=T,showCategory = 20,title = "hs DAVID GO CC enrichment wt_phenotype \n genelist Selected by padj<0.05 ")
dev.off()


#All
pdf("./DAVID/GO/Hsapiens/CC/hs_wt_phenotype_david_GO_CC_ALL_barplot.pdf", height = 20, width = 16)
barplot(hs_wt_phenotype_david_GO_CC_ALL,drop=T,showCategory = 20,title = "hs DAVID GO CC enrichment wt_phenotype \n genelist Selected by LFC > 1")
dev.off()


#Dotplot the reuslts.
pdf("./DAVID/GO/Hsapiens/CC/hs_wt_phenotype_DAVID_GO_CC_DE_dotplot.pdf", height = 20, width = 16)
dotplot(hs_wt_phenotype_david_GO_CC_DE, showCategory=20, title="hs DAVID GO CC enrichment wt_phenotype \n genelist Selected by padj<0.05 ")
dev.off()

#All
pdf("./DAVID/GO/Hsapiens/CC/hs_wt_phenotype_david_GO_CC_ALL_dotlot.pdf", height = 20, width = 16)
dotplot(hs_wt_phenotype_david_GO_CC_ALL,showCategory = 20,title = "hs DAVID GO CC enrichment wt_phenotype \n genelist Selected by LFC > 1")
dev.off()



#Save the table
write.table(hs_wt_phenotype_david_GO_CC_DE@result,"./DAVID/GO/Hsapiens/CC/hs_wt_phenotype_DAVID_GO_CC_DE.tsv",sep = "\t",row.names = F,col.names = T)
write.table(hs_wt_phenotype_david_GO_CC_ALL@result,"./DAVID/GO/Hsapiens/CC/hs_wt_phenotype_DAVID_GO_CC_ALL.tsv",sep = "\t",row.names = F,col.names = T)


############################
#                          #
#       ZSCORED            #
#                          #
############################

#ALL
#Select the top pathways based on the ZSCORE added to the table.

hs_wt_phenotype_david_GO_CC_ALL_topPathWays<-top_paths_DAVID(hs_wt_phenotype_david_GO_CC_ALL,anno_human_wt_phenotype_stats,20)
head(hs_wt_phenotype_david_GO_CC_ALL_topPathWays[,c(1,2,5,6,10)]);tail(hs_wt_phenotype_david_GO_CC_ALL_topPathWays[,c(1,2,5,6,10)])


#Plot in a more visual way the pathways up and down regulated.
#Create the plot of the gsea file using ggplplot

pdf("./DAVID/GO/Hsapiens/CC/hs_wt_phenotype_david_GO_CC_ALL_topPathWays.pdf",width = 30,height = 25)
ggplot(hs_wt_phenotype_david_GO_CC_ALL_topPathWays,aes(x=reorder(Description,zscore),y=zscore))+
  geom_bar(stat = "identity",fill="darkorchid2")+
  #geom_text(aes(label= round(p.adjust,3)),hjust=1.2,vjust=-.00001,size=3.5,color="white")+
  #labs(text="padjust")+
  geom_point(aes(label=p.adjust,size=rev(p.adjust)))+
  coord_flip()+
  labs(title = "Hs DAVID GO CC ALL genes LFC > 1",subtitle = sampleTable_wt_phenotype$treatment_type,size ="p-adjusted",x="Human DAVID:GO Cellular Component",y="Normalized Enrichment Score")
dev.off()


#Save the filtered data table.
fwrite(hs_wt_phenotype_david_GO_CC_ALL_topPathWays,"./DAVID/GO/Hsapiens/CC/hs_wt_phenotype_david_GO_CC_ALL_topPathWays.tsv",sep = "\t",col.names = T,row.names = F)

#DE
#Select the top pathways based on the ZSCORE added to the table.

hs_wt_phenotype_david_GO_CC_DE_topPathWays<-top_paths_DAVID(hs_wt_phenotype_david_GO_CC_DE,anno_human_wt_phenotype_stats,20)
head(hs_wt_phenotype_david_GO_CC_DE_topPathWays[,c(1,2,5,6,10)]);tail(hs_wt_phenotype_david_GO_CC_DE_topPathWays[,c(1,2,5,6,10)])


#Plot in a more visual way the pathways up and down regulated.
#Create the plot of the gsea file using ggplplot

pdf("./DAVID/GO/Hsapiens/CC/hs_wt_phenotype_david_GO_CC_DE_topPathWays.pdf",width = 30,height = 25)
ggplot(hs_wt_phenotype_david_GO_CC_DE_topPathWays,aes(x=reorder(Description,zscore),y=zscore))+
  geom_bar(stat = "identity",fill="darkorchid2")+
  #geom_text(aes(label= round(p.adjust,3)),hjust=1.2,vjust=-.00001,size=3.5,color="white")+
  #labs(text="padjust")+
  geom_point(aes(label=p.adjust,size=rev(p.adjust)))+
  coord_flip()+
  labs(title = "Hs DAVID GO CC DE genes padj<0.05",subtitle = sampleTable_wt_phenotype$treatment_type,size ="p.adjusted",x="Human DAVID:GO CC",y="Normalized Enrichment Score")
dev.off()


#Save the filtered data table.
fwrite(hs_wt_phenotype_david_GO_CC_DE_topPathWays,"./DAVID/GO/Hsapiens/CC/hs_wt_phenotype_david_GO_CC_DE_topPathWays.tsv",sep = "\t",col.names = T,row.names = F)

##################################
#                                #
#           Mouse                #
#                                #
##################################
dir.create("./DAVID/GO/Mmusculus/CC/")

###########################
#                         #
#                         #
#     DE padj < 0.05      #
#                         #
###########################

mm_wt_phenotype_david_GO_CC_DE <- enrichDAVID(gene = mm_wt_phenotype_diff_expressed_genes,
                                                     idType = "ENTREZ_GENE_ID",
                                                     minGSSize = 10,
                                                     maxGSSize = 500,
                                                     annotation = "GOTERM_CC_ALL",
                                                     pvalueCutoff = 0.05,
                                                     pAdjustMethod = "BH",
                                                     david.user = "clusterProfiler@hku.hk")

###########################
#                         #
#                         #
#     ALL LFC  1          #
#                         #
###########################
mm_wt_phenotype_david_GO_CC_ALL <- enrichDAVID(gene = names(mm_wt_phenotype_ranked_genes[abs(mm_wt_phenotype_ranked_genes)>1]),
                                                      idType = "ENTREZ_GENE_ID",
                                                      minGSSize = 10,
                                                      maxGSSize = 500,
                                                      annotation = "GOTERM_CC_ALL",
                                                      pvalueCutoff = 0.05,
                                                      pAdjustMethod = "BH",
                                                      david.user = "clusterProfiler@hku.hk")

#Translate the names to entrezIDS.
mm_wt_phenotype_david_GO_CC_DE<-setReadable(mm_wt_phenotype_david_GO_CC_DE,OrgDb = org.Mm.eg.db,keyType = "ENTREZID")
mm_wt_phenotype_david_GO_CC_ALL<-setReadable(mm_wt_phenotype_david_GO_CC_ALL,OrgDb = org.Mm.eg.db,keyType = "ENTREZID")

#Barplot the results.
pdf("./DAVID/GO/Mmusculus/CC/mm_wt_phenotype_david_GO_CC_DE_barplot.pdf", height = 20, width = 16)
barplot(mm_wt_phenotype_david_GO_CC_DE,drop=T,showCategory = 20,title = "mm DAVID GO CC enrichment wt_phenotype \n genelist Selected by padj<0.05 ")
dev.off()


#All
pdf("./DAVID/GO/Mmusculus/CC/mm_wt_phenotype_david_GO_CC_ALL_barplot.pdf", height = 20, width = 16)
barplot(mm_wt_phenotype_david_GO_CC_ALL,drop=T,showCategory = 20,title = "mm DAVID GO CC enrichment wt_phenotype \n genelist Selected by LFC > 1")
dev.off()


#Dotplot the reuslts.
pdf("./DAVID/GO/Mmusculus/CC/mm_wt_phenotype_DAVID_GO_CC_DE_dotplot.pdf", height = 20, width = 16)
dotplot(mm_wt_phenotype_david_GO_CC_DE, showCategory=20, title="mm DAVID GO CC enrichment wt_phenotype \n genelist Selected by padj<0.05 ")
dev.off()

#All
pdf("./DAVID/GO/Mmusculus/CC/mm_wt_phenotype_david_GO_CC_ALL_dotlot.pdf", height = 20, width = 16)
dotplot(mm_wt_phenotype_david_GO_CC_ALL,showCategory = 20,title = "mm DAVID GO CC enrichment wt_phenotype \n genelist Selected by LFC > 1")
dev.off()



#Save the table
write.table(mm_wt_phenotype_david_GO_CC_DE@result,"./DAVID/GO/Mmusculus/CC/mm_wt_phenotype_DAVID_GO_CC_DE.tsv",sep = "\t",row.names = F,col.names = T)
write.table(mm_wt_phenotype_david_GO_CC_ALL@result,"./DAVID/GO/Mmusculus/CC/mm_wt_phenotype_DAVID_GO_CC_ALL.tsv",sep = "\t",row.names = F,col.names = T)


############################
#                          #
#       ZSCORED            #
#                          #
############################

#ALL
#Select the top pathways based on the ZSCORE added to the table.

mm_wt_phenotype_david_GO_CC_ALL_topPathWays<-top_paths_DAVID(mm_wt_phenotype_david_GO_CC_ALL,anno_mouse_wt_phenotype_stats,20)
head(mm_wt_phenotype_david_GO_CC_ALL_topPathWays[,c(1,2,5,6,10)]);tail(mm_wt_phenotype_david_GO_CC_ALL_topPathWays[,c(1,2,5,6,10)])


#Plot in a more visual way the pathways up and down regulated.
#Create the plot of the gsea file using ggplplot

pdf("./DAVID/GO/Mmusculus/CC/mm_wt_phenotype_david_GO_CC_ALL_topPathWays.pdf",width = 30,height = 25)
ggplot(mm_wt_phenotype_david_GO_CC_ALL_topPathWays,aes(x=reorder(Description,zscore),y=zscore))+
  geom_bar(stat = "identity",fill="darkorchid2")+
  #geom_text(aes(label= round(p.adjust,3)),hjust=1.2,vjust=-.00001,size=3.5,color="white")+
  #labs(text="padjust")+
  geom_point(aes(label=p.adjust,size=rev(p.adjust)))+
  coord_flip()+
  labs(title = "mm DAVID GO CC ALL genes LFC > 1",subtitle = sampleTable_wt_phenotype$treatment_type,size ="p-adjusted",x="Mouse DAVID:GO CC",y="Normalized Enrichment Score")
dev.off()


#Save the filtered data table.
fwrite(mm_wt_phenotype_david_GO_CC_ALL_topPathWays,"./DAVID/GO/Mmusculus/CC/mm_wt_phenotype_david_GO_CC_ALL_topPathWays.tsv",sep = "\t",col.names = T,row.names = F)

#DE
#Select the top pathways based on the ZSCORE added to the table.

mm_wt_phenotype_david_GO_CC_DE_topPathWays<-top_paths_DAVID(mm_wt_phenotype_david_GO_CC_DE,anno_mouse_wt_phenotype_stats,20)
head(mm_wt_phenotype_david_GO_CC_DE_topPathWays[,c(1,2,5,6,10)]);tail(mm_wt_phenotype_david_GO_CC_DE_topPathWays[,c(1,2,5,6,10)])


#Plot in a more visual way the pathways up and down regulated.
#Create the plot of the gsea file using ggplplot

pdf("./DAVID/GO/Mmusculus/CC/mm_wt_phenotype_david_GO_CC_DE_topPathWays.pdf",width = 30,height = 25)
ggplot(mm_wt_phenotype_david_GO_CC_DE_topPathWays,aes(x=reorder(Description,zscore),y=zscore))+
  geom_bar(stat = "identity",fill="darkorchid2")+
  #geom_text(aes(label= round(p.adjust,3)),hjust=1.2,vjust=-.00001,size=3.5,color="white")+
  #labs(text="padjust")+
  geom_point(aes(label=p.adjust,size=rev(p.adjust)))+
  coord_flip()+
  labs(title = "mm DAVID GO CC DE genes padj<0.05",subtitle = sampleTable_wt_phenotype$treatment_type,size ="p.adjusted",x="Mouse DAVID:GO CC",y="Normalized Enrichment Score")
dev.off()


#Save the filtered data table.
fwrite(mm_wt_phenotype_david_GO_CC_DE_topPathWays,"./DAVID/GO/Mmusculus/CC/mm_wt_phenotype_david_GO_CC_DE_topPathWays.tsv",sep = "\t",col.names = T,row.names = F)






##################################
#                                #
#       DAVID  REACTOME          #
#                                #
##################################
#                                #
#           HUMAN                #
#                                #
##################################
cat("DAVID REACTOME HUMAN ENRICHMENT")
print("DAVID REACTOME HUMAN ENRICHMENT")


dir.create(recursive = T,"./DAVID/REACTOME/Hsapiens/")
#Differentially expressed ENTREZIDS
head(hs_wt_phenotype_diff_expressed_genes)


#All genes entrez ids, which abs value is over 1.
head(hs_wt_phenotype_ranked_genes)
head(names(hs_wt_phenotype_ranked_genes[abs(hs_wt_phenotype_ranked_genes)>1]))

#DE genes with padj<0.05
###########################
#                         #
#                         #
#     DE padj < 0.05      #
#                         #
###########################

hs_wt_phenotype_david_REACTOME_DE <- enrichDAVID(gene = hs_wt_phenotype_diff_expressed_genes,
                                                        idType = "ENTREZ_GENE_ID",
                                                        minGSSize = 10,
                                                        maxGSSize = 500,
                                                        annotation = "REACTOME_PATHWAY",
                                                        pvalueCutoff = 0.05,
                                                        pAdjustMethod = "BH",
                                                        david.user = "clusterProfiler@hku.hk")



#All LFC > 1
###########################
#                         #
#                         #
#     ALL LFC  1          #
#                         #
###########################
hs_wt_phenotype_david_REACTOME_ALL <- enrichDAVID(gene = names(hs_wt_phenotype_ranked_genes[abs(hs_wt_phenotype_ranked_genes)>1]),
                                                         idType = "ENTREZ_GENE_ID",
                                                         minGSSize = 10,
                                                         maxGSSize = 500,
                                                         annotation = "REACTOME_PATHWAY",
                                                         pvalueCutoff = 0.05,
                                                         pAdjustMethod = "BH",
                                                         david.user = "clusterProfiler@hku.hk")

#Translate the names to entrezIDS.
hs_wt_phenotype_david_REACTOME_ALL<-setReadable(hs_wt_phenotype_david_REACTOME_ALL,OrgDb = org.Hs.eg.db,keyType = "ENTREZID")
hs_wt_phenotype_david_REACTOME_DE<-setReadable(hs_wt_phenotype_david_REACTOME_DE,OrgDb = org.Hs.eg.db,keyType = "ENTREZID")

#View
head(hs_wt_phenotype_david_REACTOME_DE@result,30)
head(hs_wt_phenotype_david_REACTOME_ALL@result,30)



#Barplot the results.
pdf("./DAVID/REACTOME/Hsapiens/hs_wt_phenotype_david_REACTOME_barplot_DE_genes.pdf", height = 20, width = 16)
barplot(hs_wt_phenotype_david_REACTOME_DE,drop=T,showCategory = 12,title = "hs DAVID REACTOME enrichment wt_phenotype \n DE padj < 0.05 genelist")# Selected by padj<0.05 ")
dev.off()
#All
pdf("./DAVID/REACTOME/Hsapiens/hs_wt_phenotype_david_REACTOME_barplot_ALL_genes.pdf", height = 20, width = 16)
barplot(hs_wt_phenotype_david_REACTOME_ALL,drop=T,showCategory = 12,title = "hs DAVID REACTOME enrichment wt_phenotype \n all LFC > 1 genelist")# Selected by padj<0.05 ")
dev.off()


#Dotplot the reuslts.
#All
pdf("./DAVID/REACTOME/Hsapiens/hs_wt_phenotype_david_REACTOMEdotplot_ALL_genes.pdf", height = 20, width = 16)
dotplot(hs_wt_phenotype_david_REACTOME_ALL, showCategory=12, title="hs DAVID REACTOME enrichment wt_phenotype \n all genelist LFC> 1")# Selected by padj<0.05 ")
dev.off()

#DE
pdf("./DAVID/REACTOME/Hsapiens/hs_wt_phenotype_david_REACTOMEdotplot_DE_genes.pdf", height = 20, width = 16)
dotplot(hs_wt_phenotype_david_REACTOME_DE, showCategory=12, title="hs DAVID REACTOME enrichment wt_phenotype \n DE genelist padj<0.05")# Selected by padj<0.05 ")
dev.off()


#Save the table
write.table(hs_wt_phenotype_david_REACTOME_DE@result,"./DAVID/REACTOME/Hsapiens/hs_wt_phenotype_DAVID_REACTOME_DE.tsv",sep = "\t",row.names = F,col.names = T)

write.table(hs_wt_phenotype_david_REACTOME_ALL@result,"./DAVID/REACTOME/Hsapiens/hs_wt_phenotype_DAVID_REACTOME_ALL.tsv",sep = "\t",row.names = F,col.names = T)



############################
#                          #
#       ZSCORED            #
#                          #
############################

#ALL
#Select the top pathways based on the ZSCORE added to the table.

hs_wt_phenotype_david_REACTOME_ALL_topPathWays<-top_paths_DAVID(hs_wt_phenotype_david_REACTOME_ALL,anno_human_wt_phenotype_stats,20)
head(hs_wt_phenotype_david_REACTOME_ALL_topPathWays[,c(1,2,5,6,10)]);tail(hs_wt_phenotype_david_REACTOME_ALL_topPathWays[,c(1,2,5,6,10)])


#Plot in a more visual way the pathways up and down regulated.
#Create the plot of the gsea file using ggplplot

pdf("./DAVID/REACTOME/Hsapiens/hs_wt_phenotype_david_REACTOME_ALL_topPathWays.pdf",width = 30,height = 25)
ggplot(hs_wt_phenotype_david_REACTOME_ALL_topPathWays,aes(x=reorder(Description,zscore),y=zscore))+
  geom_bar(stat = "identity",fill="orange")+
  #geom_text(aes(label= round(p.adjust,3)),hjust=1.2,vjust=-.00001,size=3.5,color="white")+
  #labs(text="padjust")+
  geom_point(aes(label=p.adjust,size=rev(p.adjust)))+
  coord_flip()+
  labs(title = "Hs DAVID REACTOME ALL genes LFC > 1",subtitle = sampleTable_wt_phenotype$treatment_type,size ="p-adjusted",x="Human DAVID:REACTOME",y="Normalized Enrichment Score")
dev.off()


#Save the filtered data table.
fwrite(hs_wt_phenotype_david_REACTOME_ALL_topPathWays,"./DAVID/REACTOME/Hsapiens/hs_wt_phenotype_david_REACTOME_ALL_topPathWays.tsv",sep = "\t",col.names = T,row.names = F)

#ALL
#Select the top pathways based on the ZSCORE added to the table.

hs_wt_phenotype_david_REACTOME_DE_topPathWays<-top_paths_DAVID(hs_wt_phenotype_david_REACTOME_DE,anno_human_wt_phenotype_stats,20)
head(hs_wt_phenotype_david_REACTOME_DE_topPathWays[,c(1,2,5,6,10)]);tail(hs_wt_phenotype_david_REACTOME_DE_topPathWays[,c(1,2,5,6,10)])


#Plot in a more visual way the pathways up and down regulated.
#Create the plot of the gsea file using ggplplot

pdf("./DAVID/REACTOME/Hsapiens/hs_wt_phenotype_david_REACTOME_DE_topPathWays.pdf",width = 30,height = 25)
ggplot(hs_wt_phenotype_david_REACTOME_DE_topPathWays,aes(x=reorder(Description,zscore),y=zscore))+
  geom_bar(stat = "identity",fill="orange")+
  #geom_text(aes(label= round(p.adjust,3)),hjust=1.2,vjust=-.00001,size=3.5,color="white")+
  #labs(text="padjust")+
  geom_point(aes(label=p.adjust,size=rev(p.adjust)))+
  coord_flip()+
  labs(title = "Hs DAVID REACTOME DE genes padj< 0.05",subtitle = sampleTable_wt_phenotype$treatment_type,size ="p.adjusted",x="Human DAVID:REACTOME",y="Normalized Enrichment Score")
dev.off()


#Save the filtered data table.
fwrite(hs_wt_phenotype_david_REACTOME_DE_topPathWays,"./DAVID/REACTOME/Hsapiens/hs_wt_phenotype_david_REACTOME_DE_topPathWays.tsv",sep = "\t",col.names = T,row.names = F)



##################################
#                                #
#           Mouse                #
#                                #
##################################
print("MOUSE DAVID REACTOME")
cat("MOUSE DAVID REACTOME")

dir.create(recursive = T,"./DAVID/REACTOME/Mmusculus/")
#Differentially expressed ENTREZIDS
head(mm_wt_phenotype_diff_expressed_genes)


#All genes entrez ids, which abs value is over 1.
head(mm_wt_phenotype_ranked_genes)
head(names(mm_wt_phenotype_ranked_genes[abs(mm_wt_phenotype_ranked_genes)>1]))

#DE genes with padj<0.05
###########################
#                         #
#                         #
#     DE padj < 0.05      #
#                         #
###########################

mm_wt_phenotype_david_REACTOME_DE <- enrichDAVID(gene = mm_wt_phenotype_diff_expressed_genes,
                                                        idType = "ENTREZ_GENE_ID",
                                                        minGSSize = 10,
                                                        maxGSSize = 500,
                                                        annotation = "REACTOME_PATHWAY",
                                                        pvalueCutoff = 0.05,
                                                        pAdjustMethod = "BH",
                                                        david.user = "clusterProfiler@hku.hk")



#All LFC > 1
###########################
#                         #
#                         #
#     ALL LFC  1          #
#                         #
###########################
mm_wt_phenotype_david_REACTOME_ALL <- enrichDAVID(gene = names(mm_wt_phenotype_ranked_genes[abs(mm_wt_phenotype_ranked_genes)>1]),
                                                         idType = "ENTREZ_GENE_ID",
                                                         minGSSize = 10,
                                                         maxGSSize = 500,
                                                         annotation = "REACTOME_PATHWAY",
                                                         pvalueCutoff = 0.05,
                                                         pAdjustMethod = "BH",
                                                         david.user = "clusterProfiler@hku.hk")

#Translate the names to entrezIDS.
mm_wt_phenotype_david_REACTOME_ALL<-setReadable(mm_wt_phenotype_david_REACTOME_ALL,OrgDb = org.Mm.eg.db,keyType = "ENTREZID")
mm_wt_phenotype_david_REACTOME_DE<-setReadable(mm_wt_phenotype_david_REACTOME_DE,OrgDb = org.Mm.eg.db,keyType = "ENTREZID")

#View
head(mm_wt_phenotype_david_REACTOME_DE@result,30)
head(mm_wt_phenotype_david_REACTOME_ALL@result,30)



#Barplot the results.
pdf("./DAVID/REACTOME/Mmusculus/mm_wt_phenotype_david_REACTOME_barplot_DE_genes.pdf", height = 20, width = 16)
barplot(mm_wt_phenotype_david_REACTOME_DE,drop=T,showCategory = 12,title = "mm DAVID REACTOME enrichment wt_phenotype \n DE padj < 0.05 genelist")# Selected by padj<0.05 ")
dev.off()
#All
pdf("./DAVID/REACTOME/Mmusculus/mm_wt_phenotype_david_REACTOME_barplot_ALL_genes.pdf", height = 20, width = 16)
barplot(mm_wt_phenotype_david_REACTOME_ALL,drop=T,showCategory = 12,title = "mm DAVID REACTOME enrichment wt_phenotype \n all LFC > 1 genelist")# Selected by padj<0.05 ")
dev.off()


#Dotplot the reuslts.
#All
pdf("./DAVID/REACTOME/Mmusculus/mm_wt_phenotype_david_REACTOMEdotplot_ALL_genes.pdf", height = 20, width = 16)
dotplot(mm_wt_phenotype_david_REACTOME_ALL, showCategory=12, title="mm DAVID REACTOME enrichment wt_phenotype \n all genelist LFC> 1")# Selected by padj<0.05 ")
dev.off()

#DE
pdf("./DAVID/REACTOME/Mmusculus/mm_wt_phenotype_david_REACTOMEdotplot_DE_genes.pdf", height = 20, width = 16)
dotplot(mm_wt_phenotype_david_REACTOME_DE, showCategory=12, title="mm DAVID REACTOME enrichment wt_phenotype \n DE genelist padj<0.05")# Selected by padj<0.05 ")
dev.off()


#Save the table
write.table(mm_wt_phenotype_david_REACTOME_DE@result,"./DAVID/REACTOME/Mmusculus/mm_wt_phenotype_DAVID_REACTOME_DE.tsv",sep = "\t",row.names = F,col.names = T)

write.table(mm_wt_phenotype_david_REACTOME_ALL@result,"./DAVID/REACTOME/Mmusculus/mm_wt_phenotype_DAVID_REACTOME_ALL.tsv",sep = "\t",row.names = F,col.names = T)



############################
#                          #
#       ZSCORED            #
#                          #
############################

#ALL
#Select the top pathways based on the ZSCORE added to the table.

mm_wt_phenotype_david_REACTOME_ALL_topPathWays<-top_paths_DAVID(mm_wt_phenotype_david_REACTOME_ALL,anno_mouse_wt_phenotype_stats,20)
head(mm_wt_phenotype_david_REACTOME_ALL_topPathWays[,c(1,2,5,6,10)]);tail(mm_wt_phenotype_david_REACTOME_ALL_topPathWays[,c(1,2,5,6,10)])


#Plot in a more visual way the pathways up and down regulated.
#Create the plot of the gsea file using ggplplot

pdf("./DAVID/REACTOME/Mmusculus/mm_wt_phenotype_david_REACTOME_ALL_topPathWays.pdf",width = 30,height = 25)
ggplot(mm_wt_phenotype_david_REACTOME_ALL_topPathWays,aes(x=reorder(Description,zscore),y=zscore))+
  geom_bar(stat = "identity",fill="orange")+
  #geom_text(aes(label= round(p.adjust,3)),hjust=1.2,vjust=-.00001,size=3.5,color="white")+
  #labs(text="padjust")+
  geom_point(aes(label=p.adjust,size=rev(p.adjust)))+
  coord_flip()+
  labs(title = "Mm DAVID REACTOME ALL genes LFC > 1",subtitle = sampleTable_wt_phenotype$treatment_type,size ="p-adjusted",x="Mouse DAVID:REACTOME",y="Normalized Enrichment Score")
dev.off()


#Save the filtered data table.
fwrite(mm_wt_phenotype_david_REACTOME_ALL_topPathWays,"./DAVID/REACTOME/Mmusculus/mm_wt_phenotype_david_REACTOME_ALL_topPathWays.tsv",sep = "\t",col.names = T,row.names = F)

#ALL
#Select the top pathways based on the ZSCORE added to the table.

mm_wt_phenotype_david_REACTOME_DE_topPathWays<-top_paths_DAVID(mm_wt_phenotype_david_REACTOME_DE,anno_human_wt_phenotype_stats,20)
head(mm_wt_phenotype_david_REACTOME_DE_topPathWays[,c(1,2,5,6,10)]);tail(mm_wt_phenotype_david_REACTOME_DE_topPathWays[,c(1,2,5,6,10)])


#Plot in a more visual way the pathways up and down regulated.
#Create the plot of the gsea file using ggplplot

pdf("./DAVID/REACTOME/Mmusculus/mm_wt_phenotype_david_REACTOME_DE_topPathWays.pdf",width = 30,height = 25)
ggplot(mm_wt_phenotype_david_REACTOME_DE_topPathWays,aes(x=reorder(Description,zscore),y=zscore))+
  geom_bar(stat = "identity",fill="orange")+
  #geom_text(aes(label= round(p.adjust,3)),hjust=1.2,vjust=-.00001,size=3.5,color="white")+
  #labs(text="padjust")+
  geom_point(aes(label=p.adjust,size=rev(p.adjust)))+
  coord_flip()+
  labs(title = "Mm DAVID REACTOME DE genes padj< 0.05",subtitle = sampleTable_wt_phenotype$treatment_type,size ="p.adjusted",x="Mouse DAVID:REACTOME",y="Normalized Enrichment Score")
dev.off()


#Save the filtered data table.
fwrite(mm_wt_phenotype_david_REACTOME_DE_topPathWays,"./DAVID/REACTOME/Mmusculus/mm_wt_phenotype_david_REACTOME_DE_topPathWays.tsv",sep = "\t",col.names = T,row.names = F)



##################################
#                                #
#       DATABASE WITH            #
#                                #
#       enrich REACTOME          #
#                                #
##################################
dir.create(recursive = T,"./REACTOME/Hsapiens")
print("REACTOME ENRICHMENT HUMAN")
cat("REACTOME ENRICHMENT HUMAN")
##################################
#                                #
#             HUMAN              #
#                                #
##################################


#Using  all the expressed genes and their logFC but for the analysis only the names with the pathways enhanced. for example.
head(names(hs_wt_phenotype_ranked_genes[abs(hs_wt_phenotype_ranked_genes)>1]))

#For this use the names of the genes from the list of entrez ids and logfoldchanges.

#Enrich using REACTOME database.
###########################
#                         #
#                         #
#     ALL LFC  1          #
#                         #
###########################
hs_wt_phenotype_enrich_REACTOME_ALL<-enrichPathway(gene = names(hs_wt_phenotype_ranked_genes[abs(hs_wt_phenotype_ranked_genes)>1]),
                                                          #hs_wt_phenotype_diff_expressed_genes,
                                                          pvalueCutoff = 0.05,
                                                          readable = T,
                                                          pAdjustMethod = "BH",
                                                          organism = "human",
                                                          minGSSize = 10,
                                                          maxGSSize = 500)

###########################
#                         #
#                         #
#     DE padj < 0.05      #
#                         #
###########################
hs_wt_phenotype_enrich_REACTOME_DE<-enrichPathway(gene = hs_wt_phenotype_diff_expressed_genes,
                                                         pvalueCutoff = 0.05,
                                                         readable = T,
                                                         pAdjustMethod = "BH",
                                                         organism = "human",
                                                         minGSSize = 10,
                                                         maxGSSize = 500)

#Barplot the results.
pdf("./REACTOME/Hsapiens/hs_wt_phenotype_REACTOME_DE_barplot.pdf", height = 20, width = 16)
barplot(hs_wt_phenotype_enrich_REACTOME_DE,drop=T,showCategory = 20,title = "hs REACTOME enrichment wt_phenotype \n all genelist Selected by padj<0.05 ")
dev.off()
#Barplot the results.
pdf("./REACTOME/Hsapiens/hs_wt_phenotype_REACTOME_ALL_barplot.pdf", height = 20, width = 16)
barplot(hs_wt_phenotype_enrich_REACTOME_ALL,drop=T,showCategory = 20,title = "hs REACTOME enrichment wt_phenotype \n all genelist Selected by LFC >1 ")
dev.off()


#Enrichment map
pdf("./REACTOME/Hsapiens/hs_wt_phenotype_REACTOME_DE_map_plot.pdf")
emapplot(hs_wt_phenotype_enrich_REACTOME_DE,showCategory = 20,color = "p.adjust")
dev.off()

#Enrichment map
pdf("./REACTOME/Hsapiens/hs_wt_phenotype_REACTOME_ALL_map_plot.pdf")
emapplot(hs_wt_phenotype_enrich_REACTOME_ALL,showCategory = 20,color = "p.adjust")
dev.off()



#Cnet plot to see the gene where does it belong to what function,and paint where is it expressed.
pdf("./REACTOME/Hsapiens/hs_wt_phenotype_REACTOME_DE_cnetplot_gene_to_function.pdf")
cnetplot(hs_wt_phenotype_enrich_REACTOME_DE,showCategory = 4,foldChange = hs_wt_phenotype_ranked_genes)
dev.off()


#Cnet plot to see the gene where does it belong to what function,and paint where is it expressed.
pdf("./REACTOME/Hsapiens/hs_wt_phenotype_REACTOME_ALL_cnetplot_gene_to_function.pdf")
cnetplot(hs_wt_phenotype_enrich_REACTOME_ALL,showCategory = 4,foldChange = hs_wt_phenotype_ranked_genes)
dev.off()


############################
#                          #
#       ZSCORED            #
#                          #
############################

#ALL
#Select the top pathways based on the ZSCORE added to the table.
hs_wt_phenotype_enrich_REACTOME_ALL_topPathWays<-top_paths_DAVID(hs_wt_phenotype_enrich_REACTOME_ALL,anno_human_wt_phenotype_stats,15)
head(hs_wt_phenotype_enrich_REACTOME_ALL_topPathWays[,c(1,2,5,6,10)]);tail(hs_wt_phenotype_enrich_REACTOME_ALL_topPathWays[,c(1,2,5,6,10)])


#Plot in a more visual way the pathways up and down regulated.
#Create the plot of the gsea file using ggplplot


pdf("./REACTOME/Hsapiens/hs_wt_phenotype_enrich_REACTOME_ALL_topPathWays.pdf",width = 30,height = 25)
ggplot(hs_wt_phenotype_enrich_REACTOME_ALL_topPathWays,aes(x=reorder(Description,zscore),y=zscore))+
  geom_bar(stat = "identity",fill="orange")+
  #geom_text(aes(label= round(p.adjust,3)),hjust=1.2,vjust=-.00001,size=3.5,color="white")+
  #labs(text="padjust")+
  geom_point(aes(label=p.adjust,size=rev(p.adjust)))+
  coord_flip()+
  labs(title = "Hs REACTOME genes ALL LFC > 1",subtitle = sampleTable_wt_phenotype$treatment_type,size ="p-adjusted",x="Human REACTOME",y="Normalized Enrichment Score")
dev.off()


#Save the filtered data table.
fwrite(hs_wt_phenotype_enrich_REACTOME_ALL_topPathWays,"./REACTOME/Hsapiens/hs_wt_phenotype_enrich_REACTOME_ALL_topPathWays.tsv",sep = "\t",col.names = T,row.names = F)

#DE
#Select the top pathways based on the ZSCORE added to the table.
hs_wt_phenotype_enrich_REACTOME_DE_topPathWays<-top_paths_DAVID(hs_wt_phenotype_enrich_REACTOME_DE,anno_human_wt_phenotype_stats,15)
head(hs_wt_phenotype_enrich_REACTOME_DE_topPathWays[,c(1,2,5,6,10)]);tail(hs_wt_phenotype_enrich_REACTOME_DE_topPathWays[,c(1,2,5,6,10)])


#Plot in a more visual way the pathways up and down regulated.
#Create the plot of the gsea file using ggplplot

pdf("./REACTOME/Hsapiens/hs_wt_phenotype_enrich_REACTOME_DE_topPathWays.pdf",width = 30,height = 25)
ggplot(hs_wt_phenotype_enrich_REACTOME_DE_topPathWays,aes(x=reorder(Description,zscore),y=zscore))+
  geom_bar(stat = "identity",fill="orange")+
  #geom_text(aes(label= round(p.adjust,3)),hjust=1.2,vjust=-.00001,size=3.5,color="white")+
  #labs(text="padjust")+
  geom_point(aes(label=p.adjust,size=rev(p.adjust)))+
  coord_flip()+
  labs(title = "Hs REACTOME genes DE padj<0.05",subtitle = sampleTable_wt_phenotype$treatment_type,size ="p.adjusted",x="Human REACTOME",y="Normalized Enrichment Score")
dev.off()


#Save the filtered data table.
fwrite(hs_wt_phenotype_enrich_REACTOME_DE_topPathWays,"./REACTOME/Hsapiens/hs_wt_phenotype_enrich_REACTOME_DE_topPathWays.tsv",sep = "\t",col.names = T,row.names = F)




##################################
#                                #
#           REACTOME IN          #          
#                                #
#             MOSUE              #
#                                #
##################################
print("REACTOME MOUSE")
cat("REACTOME MOUSE")
dir.create("./REACTOME/Mmusculus")
#Using  all the expressed genes and their logFC but for the analysis only the names with the pathways enhanced. for example.
head(names(mm_wt_phenotype_ranked_genes[abs(mm_wt_phenotype_ranked_genes)>1]))
head(mm_wt_phenotype_diff_expressed_genes)
#For this use the names of the genes from the list of entrez ids and logfoldchanges.

#Enrich using REACTOME database.                     
###########################
#                         #
#                         #
#     ALL LFC  1          #
#                         #
###########################
mm_wt_phenotype_enrich_REACTOME_ALL<-enrichPathway(gene = names(mm_wt_phenotype_ranked_genes[abs(mm_wt_phenotype_ranked_genes)>1]),
                                                          pvalueCutoff = 0.05,
                                                          readable = T,
                                                          pAdjustMethod = "BH",
                                                          organism = "mouse",
                                                          minGSSize = 10,
                                                          maxGSSize = 500)
###########################
#                         #
#                         #
#     DE padj < 0.05      #
#                         #
###########################
mm_wt_phenotype_enrich_REACTOME_DE<-enrichPathway(gene = mm_wt_phenotype_diff_expressed_genes,
                                                         pvalueCutoff = 0.05,
                                                         readable = T,
                                                         pAdjustMethod = "BH",
                                                         organism = "mouse",
                                                         minGSSize = 10,
                                                         maxGSSize = 500)

#Barplot the results.
pdf("./REACTOME/Mmusculus/mm_wt_phenotype_REACTOME_DE_barplot.pdf", height = 20, width = 16)
barplot(mm_wt_phenotype_enrich_REACTOME_DE,drop=T,showCategory = 20,title = "mm REACTOME enrichment wt_phenotype \n all genelist Selected by padj<0.05 ")
dev.off()

pdf("./REACTOME/Mmusculus/mm_wt_phenotype_REACTOME_ALL_barplot.pdf", height = 20, width = 16)
barplot(mm_wt_phenotype_enrich_REACTOME_ALL,drop=T,showCategory = 20,title = "mm REACTOME enrichment wt_phenotype \n all genelist Selected by LFC > 1 ")
dev.off()

#Enrichment map
pdf("./REACTOME/Mmusculus/mm_wt_phenotype_REACTOME_DE_map_plot.pdf")
emapplot(mm_wt_phenotype_enrich_REACTOME_DE,showCategory = 20,color = "p.adjust")
dev.off()

#Enrichment map
pdf("./REACTOME/Mmusculus/mm_wt_phenotype_REACTOME_ALL_map_plot.pdf")
emapplot(mm_wt_phenotype_enrich_REACTOME_ALL,showCategory = 20,color = "p.adjust")
dev.off()



#Cnet plot to see the gene where does it belong to what function,and paint where is it expressed.
pdf("./REACTOME/Mmusculus/mm_wt_phenotype_REACTOME_DE_cnetplot.pdf")
cnetplot(mm_wt_phenotype_enrich_REACTOME_DE,showCategory = 4,foldChange = mm_wt_phenotype_ranked_genes)
dev.off()

pdf("./REACTOME/Mmusculusmm_wt_phenotype_REACTOME_ALL_cnetplot.pdf")
cnetplot(mm_wt_phenotype_enrich_REACTOME_ALL,showCategory = 4,foldChange = mm_wt_phenotype_ranked_genes)
dev.off()

############################
#                          #
#       ZSCORED            #
#                          #
############################

#ALL
#Select the top pathways based on the ZSCORE added to the table.
mm_wt_phenotype_enrich_REACTOME_ALL_topPathWays<-top_paths_DAVID(mm_wt_phenotype_enrich_REACTOME_ALL,anno_mouse_wt_phenotype_stats,15)
head(mm_wt_phenotype_enrich_REACTOME_ALL_topPathWays[,c(1,2,5,6,10)]);tail(mm_wt_phenotype_enrich_REACTOME_ALL_topPathWays[,c(1,2,5,6,10)])


#Plot in a more visual way the pathways up and down regulated.
#Create the plot of the gsea file using ggplplot


pdf("./REACTOME/Mmusculus/mm_wt_phenotype_enrich_REACTOME_ALL_topPathWays.pdf",width = 30,height = 25)
ggplot(mm_wt_phenotype_enrich_REACTOME_ALL_topPathWays,aes(x=reorder(Description,zscore),y=zscore))+
  geom_bar(stat = "identity",fill="orange")+
  #geom_text(aes(label= round(p.adjust,3)),hjust=1.2,vjust=-.00001,size=3.5,color="white")+
  #labs(text="padjust")+
  geom_point(aes(label=p.adjust,size=rev(p.adjust)))+
  coord_flip()+
  labs(title = "mm REACTOME ALL genes LFC > 1",subtitle = sampleTable_wt_phenotype$treatment_type,size ="p-adjusted",x="Mouse REACTOME",y="Normalized Enrichment Score")
dev.off()


#Save the filtered data table.
fwrite(mm_wt_phenotype_enrich_REACTOME_ALL_topPathWays,"./REACTOME/Mmusculus/mm_wt_phenotype_enrich_REACTOME_ALL_topPathWays.tsv",sep = "\t",col.names = T,row.names = F)

#DE
#Select the top pathways based on the ZSCORE added to the table.
mm_wt_phenotype_enrich_REACTOME_DE_topPathWays<-top_paths_DAVID(mm_wt_phenotype_enrich_REACTOME_DE,anno_mouse_wt_phenotype_stats,15)
head(mm_wt_phenotype_enrich_REACTOME_DE_topPathWays[,c(1,2,5,6,10)]);tail(mm_wt_phenotype_enrich_REACTOME_DE_topPathWays[,c(1,2,5,6,10)])


#Plot in a more visual way the pathways up and down regulated.
#Create the plot of the gsea file using ggplplot


pdf("./REACTOME/Mmusculus/mm_wt_phenotype_enrich_REACTOME_DE_topPathWays.pdf",width = 30,height = 25)
ggplot(mm_wt_phenotype_enrich_REACTOME_DE_topPathWays,aes(x=reorder(Description,zscore),y=zscore))+
  geom_bar(stat = "identity",fill="orange")+
  #geom_text(aes(label= round(p.adjust,3)),hjust=1.2,vjust=-.00001,size=3.5,color="white")+
  #labs(text="padjust")+
  geom_point(aes(label=p.adjust,size=rev(p.adjust)))+
  coord_flip()+
  labs(title = "mm REACTOME DE genes padj<0.05",subtitle = sampleTable_wt_phenotype$treatment_type,size ="p-adjusted",x="Mouse REACTOME",y="Normalized Enrichment Score")
dev.off()


#Save the filtered data table.
fwrite(mm_wt_phenotype_enrich_REACTOME_DE_topPathWays,"./REACTOME/Mmusculus/mm_wt_phenotype_enrich_REACTOME_DE_topPathWays.tsv",sep = "\t",col.names = T,row.names = F)




##################################
#                                #
#          GSEA PLOT             #
#                                #
##################################
print("GSEAS plotting")
cat("GSEAS plotting")
dir.create("./gsea_plot")
#Use the preranked list of all genes for human to perform gsea analysis.

?gsePathway
hs_wt_phenotype_enrich_gsea<-gsePathway(geneList = hs_wt_phenotype_ranked_genes,
                                               organism = "human",nPerm = 1000,
                                               pvalueCutoff = 0.11,
                                               pAdjustMethod = "BH")

head(hs_wt_phenotype_enrich_gsea@result)

#Translate the names to entrezIDS.
hs_wt_phenotype_enrich_gsea<-setReadable(hs_wt_phenotype_enrich_gsea,OrgDb = org.Hs.eg.db,keyType = "ENTREZID")
head(hs_wt_phenotype_enrich_gsea@result,30)


#Plot the results from gsea.

#Map plot.
pdf("./gsea_plot/emma_plot_pvlaues.pdf",width = 20,height = 20)
emapplot(hs_wt_phenotype_enrich_gsea,color = "pvalue")
dev.off()



#Net plot of the pathway selected.
pdf("./gsea_plot/gseapathway_enriched_logfoldchange_Mitochondrial_Fatty_Acid_BetaOxidation.pdf")
viewPathway("Mitochondrial Fatty Acid Beta-Oxidation",readable = TRUE,foldChange = hs_wt_phenotype_ranked_genes)
dev.off()




###############################
#                             #
#          GOPLOT             #
#           from              #
#          enrichGO           #
#                             #
###############################
dir.create("./GOplot")
cat("Obtaining the GO PLOT object")
print("Obtaining the GO PLOT object")

table2GOplot<-function(FDR_GO_KEGG_data,anno_specie_used){
  #Select the amount of data out of the table, maximum of 30...
  cat("\n*\n")
  if (length(FDR_GO_KEGG_data$ID)>30){
    FDR_GO_KEGG_data<-FDR_GO_KEGG_data@result[1:30,]
    cat("Data is very large, selecting first 30 rows\n")
  }else {
    FDR_GO_KEGG_data<-FDR_GO_KEGG_data@result[,]
  }#Obtain the objects for circ final object.
  symbols<-list()
  symbols_extracted<-list()
  symbols_extracted_double<-list()
  cat("Creating lists...\n")
  cat("\n**\n")
  for (i in 1:length(FDR_GO_KEGG_data$ID)){
    #Extract from the table all the SYMBOLS involved in the pathways from the gnees.
    target<-paste0("target",i)
    #print(target)
    
    #Clean the genelist of each row and obtain the genes.
    symbols[[target]]<-unique(unlist(strsplit(as.character(FDR_GO_KEGG_data$geneID[i]),split = "\\/")))
    
    #Obtain all the data from the annotated genes in human, as reactome table is done using human genes.
    #Create another list with all the genes that are in the list and in reactome for obtaining all the required information.
    symbols_extracted[[target]]<-anno_specie_used[toupper(FDR_GO_KEGG_data$SYMBOL)%in%toupper(symbols)[[i]],]
    #Select only the SYMBOLS and the logFC for calculating the zscore.
    symbols_extracted_double[[target]]<-symbols_extracted[[i]][,c(6,2)]
    #Now with the genes and the logFC for each row from the selected above, calculate the zscore and add it to the main table as another column.
    FDR_GO_KEGG_data$zscore[i]<-sum(as.numeric(symbols_extracted_double[[i]][,2]))/sqrt(length(symbols[[i]]))}
  cat("Added zscore...\n")
  cat("\n***\n")
  
  FDR_GO_KEGG_data<-FDR_GO_KEGG_data[order(FDR_GO_KEGG_data$zscore,decreasing = T),]
  FDR_GO_KEGG_data_renamed<-FDR_GO_KEGG_data[,c("ID","Description","Description","Count","geneID","p.adjust","zscore")]
  cat("Changing names and formats...\n")
  cat("\n****\n")
  rownames(FDR_GO_KEGG_data_renamed)<-NULL
  colnames(FDR_GO_KEGG_data_renamed)<-c("Category","ID","TERM","Count","GENES","adj_pval","z.score")
  FDR_GO_KEGG_only_symbols<-unique(unlist(strsplit(as.character(FDR_GO_KEGG_data_renamed$GENES),split = "\\/")))
  FDR_GO_KEGG_data_extracted<-anno_specie_used[toupper(anno_specie_used$SYMBOL)%in%toupper(FDR_GO_KEGG_only_symbols),]
  FDR_GO_KEGG_data_extracted<-FDR_GO_KEGG_data_extracted[,c("SYMBOL","log2FoldChange","pvalue","padj")]
  #Change the names, for GOplot functions to understand it.
  colnames(FDR_GO_KEGG_data_extracted)<-c("ID","logFC","P.Value","adj.P.val")
  #Check the names 
  rownames(FDR_GO_KEGG_data_extracted)<-NULL
  #Replace the / for the comm_a to work the function.
  FDR_GO_KEGG_data_renamed$GENES<-gsub(x = FDR_GO_KEGG_data_renamed$GENES,pattern = "\\/",replacement=",")
  cat("Change the format of the data for [,] instead of [/]\n")
  cat("\n*****\n")
  #Create the circ object by mergin both tables.
  cat("Creating the CIRC object for GOPLOT\n")
  cat("\n******\n")
  FDR_GO_KEGG_circ<-circle_dat(FDR_GO_KEGG_data_renamed,FDR_GO_KEGG_data_extracted)
  ######################################
  #Create the processes vector for the names of the funcitons ordered by zscore.
  FDR_GO_KEGG_data_process<-FDR_GO_KEGG_data_renamed$TERM
  FDR_GO_KEGG_only_symbols<-toupper(FDR_GO_KEGG_only_symbols)
  #Create the plot object using the circ generated matrix, the genelis and the processes.
  cat("Creating the CHORD object for GOPLOT\n")
  cat("\n******\n")
  FDR_GO_KEGG_chord_object<-chord_dat(data = FDR_GO_KEGG_circ,genes = FDR_GO_KEGG_data_extracted,process = FDR_GO_KEGG_data_process)
  cat("Completed, next print the object using the GOPlot package...\nDONE!\n")
  cat("\n********\n")
  
  #Plot and save the data obtained in the function
  FDR_GO_KEGG_chord_object<-as.data.frame(FDR_GO_KEGG_chord_object)
  dir.create("./GOplot/")
  pdf("./GOplot/function_legend_used.pdf")
  plot(main="z-score",
       xlab="Functions",
       ylab="Colors",x=rep(1,length(FDR_GO_KEGG_data_renamed$z.score)),col=rainbow(length(FDR_GO_KEGG_data_renamed$z.score)),pch=19,cex=2)
  dev.off()
  
  #Save the tables
  write.table(FDR_GO_KEGG_chord_object,"./GOplot/GO_chord_object.tsv",sep = "\t",col.names = T,row.names = F)
  write.table(FDR_GO_KEGG_data_renamed,"./GOplot/table_annotated_zscore.tsv",sep = "\t",col.names = T,row.names = F)
  write.table(FDR_GO_KEGG_data_extracted[order(FDR_GO_KEGG_data_extracted[2]),],"./GOplot/selected_symbols_from_the_functions.tsv",sep = "\t",col.names = T,row.names = F)
  return(FDR_GO_KEGG_chord_object)
  #return(FDR_GO_KEGG_only_symbols)
  #return(FDR_GO_KEGG_only_process)
  #return(FDR_GO_KEGG_data_renamed)
  #return(FDR_GO_KEGG_data_extracted)
  
}
head(hs_wt_phenotype_BP_genelist_enrichGO@result,2)
head(hs_wt_phenotype_david_GO_ALL@result,2)
test_human<-table2GOplot(hs_wt_phenotype_david_GO_ALL,anno_human_wt_phenotype_stats)
test_human2<-table2GOplot(hs_wt_phenotype_BP_genelist_enrichGO,anno_human_wt_phenotype_stats)


pdf("./GOplot/test_human_david.pdf",width = 30,height = 30)
GOChord(test_human,space = 0.033,gene.order = "logFC",
        lfc.min = min(test_human$logFC),
        lfc.max = max(test_human$logFC),
        gene.size = 7,
        gene.space = 0.2,
        lfc.col =rev(colorRampPalette(brewer.pal(9,"Blues"))(2.6)))

dev.off()




pdf("./GOplot/test_human_GO.pdf",width = 30,height = 30)
GOChord(test_human2,space = 0.033,gene.order = "logFC",
        lfc.min = min(test_human2$logFC),
        lfc.max = max(test_human2$logFC),
        gene.size = 7,
        gene.space = 0.2,
        lfc.col =rev(colorRampPalette(brewer.pal(9,"Blues"))(2.6)))

dev.off()


#####################################
#Save session.
dir.create("./R")
save.image(file = "./R/wt_phenotype_session.RData")

load("/mnt/haus/marius/carol_data/carol_kallisto_deseq2/wt_super_extended/R/wt_phenotype_session.RData")


#Change the permission of the folders if opened with sudo rstudio.
#system("grep marius /etc/passwd")
system("sudo chown -R 1000 $PWD")






#Get the GSEA geneset for OXPHOS.

oxphos_set<-read.csv("/home/marius/Downloads/oxphos_geneset.txt",header = T)
oxphos_set<-oxphos_set[-1,]
oxphos_geneset<-as.data.frame(anno_zebrafish_wt_phenotype_stats[toupper(anno_zebrafish_wt_phenotype_stats$SYMBOL)%in%oxphos_set,])
class(oxphos_geneset)



dev.off()

getwd()
pdf("/mnt/haus/marius/carol_data/carol_kallisto_deseq2/wt_super_extended/Plots/oxphos_annotated_zebrafish.pdf",height = 12,width = 12)
with(anno_zebrafish_wt_phenotype_stats, plot(log2FoldChange,-log10(pvalue), pch=20, 
                                             main="Volcano plot\nWT DD vs SD",
                                             xlab="Log2FoldChange",
                                             ylab="-log10 (pvalue)",
                                             col="72"))
with(subset(anno_zebrafish_wt_phenotype_stats, padj <.05), points(log2FoldChange, -log10(pvalue), pch=20, col="lightblue"))
with(subset(anno_zebrafish_wt_phenotype_stats, padj <=.05 & abs(log2FoldChange) >= 1), points(log2FoldChange, -log10(pvalue), pch=20, col="dodgerblue"))

#with(subset(anno_zebrafish_wt_phenotype_stats, abs(log2FoldChange)<1), points(log2FoldChange, -log10(pvalue), pch=20, col="blue"))

#with(subset(anno_zebrafish_wt_phenotype_stats,padj<.05 & abs(log2FoldChange)>1), points(log2FoldChange, -log10(pvalue), pch=20, col="green3"))


with(anno_zebrafish_wt_phenotype_stats[anno_zebrafish_wt_phenotype_stats$drerio_ensembl_gene_id%in%oxphos_geneset$drerio_ensembl_gene_id,],
     points(log2FoldChange,
            -log10(pvalue),
            pch=1,
            col="red",
            cex = .9))
with(subset(anno_zebrafish_wt_phenotype_stats[anno_zebrafish_wt_phenotype_stats$drerio_ensembl_gene_id%in%oxphos_geneset$drerio_ensembl_gene_id,],
            padj <=.05 & abs(log2FoldChange) > 1),
     points(log2FoldChange,
            -log10(pvalue),
            pch=1,
            col="green",
            cex = 1.3))

with(subset(anno_zebrafish_wt_phenotype_stats[anno_zebrafish_wt_phenotype_stats$drerio_ensembl_gene_id%in%oxphos_geneset$drerio_ensembl_gene_id,],
            abs(log2FoldChange) > 1 & padj <=.05),
     points(log2FoldChange,
            -log10(pvalue),
            pch=1,
            col="green",
            cex = 1.2))
with(subset(anno_zebrafish_wt_phenotype_stats[anno_zebrafish_wt_phenotype_stats$drerio_ensembl_gene_id%in%oxphos_geneset$drerio_ensembl_gene_id,],
            abs(log2FoldChange) > 1 & padj <=.05),
     points(log2FoldChange,
            -log10(pvalue),
            pch=1,
            col="green",
            cex = 0.9))
#grep(x = anno_zebrafish_wt_phenotype_stats[anno_zebrafish_wt_phenotype_stats$drerio_ensembl_gene_id%in%oxphos_geneset$drerio_ensembl_gene_id,]$SYMBOL,pattern = "cox7a2l")[1]
with(anno_zebrafish_wt_phenotype_stats[anno_zebrafish_wt_phenotype_stats$drerio_ensembl_gene_id%in%oxphos_geneset$drerio_ensembl_gene_id,][grep(x = anno_zebrafish_wt_phenotype_stats[anno_zebrafish_wt_phenotype_stats$drerio_ensembl_gene_id%in%oxphos_geneset$drerio_ensembl_gene_id,]$SYMBOL,pattern = "cox7a2l")[1],],
     points(log2FoldChange,
            -log10(pvalue),
            pch=1,
            col="purple",
            cex = 0.9))
#Plot the names
with(anno_zebrafish_wt_phenotype_stats[anno_zebrafish_wt_phenotype_stats$drerio_ensembl_gene_id%in%oxphos_geneset$drerio_ensembl_gene_id,],
     text(subset(anno_zebrafish_wt_phenotype_stats[anno_zebrafish_wt_phenotype_stats$drerio_ensembl_gene_id%in%oxphos_geneset$drerio_ensembl_gene_id,],
                   padj <.05 & abs(log2FoldChange) > 1)$log2FoldChange,
            -log10(subset(anno_zebrafish_wt_phenotype_stats[anno_zebrafish_wt_phenotype_stats$drerio_ensembl_gene_id%in%oxphos_geneset$drerio_ensembl_gene_id,],
                          padj <.05 & abs(log2FoldChange) > 1)$pvalue),
            labels = subset(anno_zebrafish_wt_phenotype_stats[anno_zebrafish_wt_phenotype_stats$drerio_ensembl_gene_id%in%oxphos_geneset$drerio_ensembl_gene_id,],
                          padj <.05 & abs(log2FoldChange) > 1)$SYMBOL,
            cex=.6,
            offset = 0.3,pos = 3))


with(anno_zebrafish_wt_phenotype_stats[anno_zebrafish_wt_phenotype_stats$drerio_ensembl_gene_id%in%oxphos_geneset$drerio_ensembl_gene_id,][46,],
     text(anno_zebrafish_wt_phenotype_stats[anno_zebrafish_wt_phenotype_stats$drerio_ensembl_gene_id%in%oxphos_geneset$drerio_ensembl_gene_id,]$log2FoldChange[46],
     -log10(anno_zebrafish_wt_phenotype_stats[anno_zebrafish_wt_phenotype_stats$drerio_ensembl_gene_id%in%oxphos_geneset$drerio_ensembl_gene_id,]$pvalue[46]),
     labels=anno_zebrafish_wt_phenotype_stats[anno_zebrafish_wt_phenotype_stats$drerio_ensembl_gene_id%in%oxphos_geneset$drerio_ensembl_gene_id,]$SYMBOL[46]),
     cex=.5,
     offset = 0.3,pos=3)


legend("bottomright",
       legend = c("Significative","No Significative","Oxphos Genes"),
       col=c("dodgerblue","grey72","black"),
       bty="o",
       pch=c(19,19,1),
       pt.cex=2,
       cex=1.2,
       text.col=c("dodgerblue","grey72","black"),
       horiz = F,
       inset = c(0.01, 0.05))

dev.off()

#Plot the OXPHOS GENES.



library(ggplot2)

pdf("/mnt/haus/marius/carol_data/carol_kallisto_deseq2/wt_super_extended/Plots/oxphos_annotated_zebrafish_ggplot.pdf",height = 15,width = 15)
ggplot(anno_zebrafish_wt_phenotype_stats,aes(log2FoldChange,-log10(pvalue)))+
  geom_point(alpha=0.6,size=3,aes(color=ifelse(padj <=0.05 & abs(log2FoldChange) >= 1,"Significative","No Significative")))+
  geom_text_repel(data = subset(anno_zebrafish_wt_phenotype_stats[anno_zebrafish_wt_phenotype_stats$drerio_ensembl_gene_id%in%oxphos_geneset$drerio_ensembl_gene_id,],
                   padj <.05 & abs(log2FoldChange) > 1),aes(log2FoldChange,-log10(pvalue),label=SYMBOL),check_overlap = F)+
  scale_color_manual(values=c("Significative"="green","No Significative"="red"))+
  xlim(min(anno_zebrafish_wt_phenotype_stats$log2FoldChange),max(anno_zebrafish_wt_phenotype_stats$log2FoldChange))+
  labs(title = sampleTable_wt_phenotype$treatment_type)+xlab("Log2FoldChange")+ylab("log10-pvalue")
 
#geom_text_repel(data = subset(anno_zebrafish_wt_phenotype_stats[anno_zebrafish_wt_phenotype_stats$drerio_ensembl_gene_id%in%oxphos_geneset$drerio_ensembl_gene_id,],padj >.05 & abs(log2FoldChange) > 1))

dev.off()


pdf("/mnt/haus/marius/carol_data/carol_kallisto_deseq2/wt_super_extended/Plots/oxphos_annotated_zebrafish_enhanced_volcano.pdf",height = 9,width = 9)
EnhancedVolcano(anno_zebrafish_wt_phenotype_stats,x = "log2FoldChange",y="pvalue",lab = anno_zebrafish_wt_phenotype_stats$SYMBOL,
                FCcutoff = 1,
                pCutoff = 0.05,
                title = "Oxphos Dr",
                selectLab = anno_zebrafish_wt_phenotype_stats$SYMBOL[anno_zebrafish_wt_phenotype_stats$drerio_ensembl_gene_id%in%oxphos_geneset$drerio_ensembl_gene_id])

dev.off()



