# BRNAseq
RNA-seq analysis in R.

RNAseq analysis script using Kallisto (pseudoallignment) + tximport(importing h5 files into R) + DESeq2 (analyzing the expression counts of genes and getting statistical data).

Continues with ENRICHMENTS for three species in parallel with Danio rerio homologues genes in Human and Mouse.

  1:Enrichment using statistically relevant data for analysis, always using padj<0.05 for GO (BP,MF,CC) and KEGG.
  
  2:Gen set enrichment analysis (ORA) with all the genes detected and ranked by the Log2FoldChange in GO and KEGG databases.
  
  3:GSEA using Hallmarks and GO:C5 for Mouse and Human, .rData obtained from (http://bioinf.wehi.edu.au/software/MSigDB/)
    This data is /gsea_bundle/ directory called in the Rscript, just save your Hallmarks,GO and use it from there.

Finally concluding with the regular most used and basic plots;
  
  1: Volcano plot.
  2: PCA
  3: Heatmap
  
Save the data filtered and raw data after all, in a txt file.

DONE.

#To be added extra graphs and Reactome pathways.
    
