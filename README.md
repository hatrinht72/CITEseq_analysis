# CITEseq_analysis

## Summary
1. Connect with a HPC for primary analysis
2. Primary analysis : using cellranger
3. Secondary analysis 

### Introduction 

The aim of this git is to analyze CITE-sequencing data from [GSE229018](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE229018)

To get insight into the effect of **aging** on the hematopoiesis, an **inductible lineage tracing systeme** was used to tracing the blood cell differenciation from hematopoietic stem cells (HSC) into more mature progenitor blood cells in old and young mice. The reporter gene Tomato was induced in 18th week for old mice and 3rd week for young mice. **Lin-cKit+Tom+** hematopoetic progenitors were isolated 4 weeks after the induction. Their transcription profile and immunophenotypic profile were analysis with **Cellular Indexing Transcriptomes and Epitopes Sequencing** (CITE-seq). 

Cells are marked by 2 tags : 
  1. **HTO** (Hashtag oligonucleotides) : A pair of antibodies raised against two different targets which are conjugated to an oligo sequence with the same barcode used to run multiple samples on one 10X lane, to distinguish the origins of the reads.
       HTO1, HTO2, HTO3 : old mice
       HTO4, HTO5, HTO6 : young mice
  2. **ADT** (antibody-derived tags) : Antibodies raised agaisnt a specific extracellular target which are conjugated to an oligo sequence with the same barcode used to quantify the expression of that target, to identify the immunophenotypic of each cell 



