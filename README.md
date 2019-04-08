# BRNAseq
RNA-seq analysis in R.

RNAseq analysis script using Kallisto (quantification) + tximport(importing h5 files into R) + DESeq2 (analyzing the expression counts of genes and getting statistical data).

Continues with ENRICHMENTS for three species in parallel with Danio rerio homologues genes in Human and Mouse.

  1:Enrichment using statistically relevant data for analysis, always using padj<0.05 for GO (BP,MF,CC) and KEGG.
    1.1: KEGG
    1.2: GO BP
    1.3: GO MF
    1.4: GO CC
      
  2:Gen set enrichment analysis (ORA) with all the genes detected and ranked by the Log2FoldChange in GO and KEGG databases.
    2.1: KEGG
    2.2: GO BP
    2.3: GO MF
    2.4: GO CC
    
  3:GSEA using Hallmarks and GO:C5 for Mouse and Human, .rData obtained from (http://bioinf.wehi.edu.au/software/MSigDB/)
    This data is /gsea_bundle/ directory called in the Rscript, just save your Hallmarks,GO and use it from there.
    3.1: HALLMARKS
    3.2: GO DATABASE

  4:DAVID
    4.1: DAVID KEGG
    4.2: DAVID GO_ALL
    4.3: DAVID GO_MF
    4.4: DAVID GO_CC
    4.5: DAVID REACTOME
    
  5:REACTOME

  6: GSEA PLOTING AND NETWORKS(TODO)
  
Finally concluding with the regular most used and basic plots;
  
  1: Volcano plot.
  2: PCA
  3: Heatmap
  4: GOPloting creating customized table.
  
Save the all the annotated data for each species used and the filtered, in a txt file (tab).

DONE.  
