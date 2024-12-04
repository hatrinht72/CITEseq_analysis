# Primary Analysis
## Bulk data & Hierarchy 
The experiment performed by 10X Genomics. 
Cells were isolated then tagged by Hashtag and ADT, pooled then divided into 4 sequencing runs. 
Assuming I was already in my citeseq_wk4 folder :

├── fastq
│   ├── MAP17-4wks-ADT
│   │   ├── MAP17-4wks-ADT_S11_I1_001.fastq.gz
│   │   ├── MAP17-4wks-ADT_S11_R1_001.fastq.gz
│   │   ├── MAP17-4wks-ADT_S11_R2_001.fastq.gz
│   │   └── wk4-adt.md5
│   ├── MAP17-4wks-Hashtag
│   │   ├── MAP17-4wks-Hashtag_S12_I1_001.fastq.gz
│   │   ├── MAP17-4wks-Hashtag_S12_R1_001.fastq.gz
│   │   ├── MAP17-4wks-Hashtag_S12_R2_001.fastq.gz
│   │   └── wk4-hashtag.md5
│   ├── MAP17-4wks-cDNA1
│   │   ├── MAP17-4wks-cDNA1_S7_I1_001.fastq.gz
│   │   ├── MAP17-4wks-cDNA1_S7_R1_001.fastq.gz
│   │   ├── MAP17-4wks-cDNA1_S7_R2_001.fastq.gz
│   │   └── wk4-cdna1.md5
│   ├── MAP17-4wks-cDNA2
│   │   ├── MAP17-4wks-cDNA2_S8_I1_001.fastq.gz
│   │   ├── MAP17-4wks-cDNA2_S8_R1_001.fastq.gz
│   │   ├── MAP17-4wks-cDNA2_S8_R2_001.fastq.gz
│   │   └── wk4-cdna2.md5
│   ├── MAP17-4wks-cDNA3
│   │   ├── MAP17-4wks-cDNA3_S9_I1_001.fastq.gz
│   │   ├── MAP17-4wks-cDNA3_S9_R1_001.fastq.gz
│   │   ├── MAP17-4wks-cDNA3_S9_R2_001.fastq.gz
│   │   └── wk4-cdna3.md5
│   └── MAP17-4wks-cDNA4
│       ├── MAP17-4wks-cDNA4_S10_I1_001.fastq.gz
│       ├── MAP17-4wks-cDNA4_S10_R1_001.fastq.gz
│       ├── MAP17-4wks-cDNA4_S10_R2_001.fastq.gz
│       └── wk4-cdna4.md5
├── feature_ref.csv
└── samplesheet.csv

In **fastq**, I have 5 folders, corresponding to 4 sequencing runs, ADT and Hashtag
I also have my feature_ref.csv and my samplesheet.csv where I specify my the list of tags and also which fastq correspond to which type of data (Gene Expression or Antibody Capture) 

**samplesheet.csv**
fastqs                                                               sample             library_type
/home/ha/y1/cahier_labo/241126/citeseq_wk4/fastq/MAP17-4wks-cDNA1/   MAP17-4wks-cDNA1   Gene Expression
/home/ha/y1/cahier_labo/241126/citeseq_wk4/fastq/MAP17-4wks-cDNA2/   MAP17-4wks-cDNA2   Gene Expression
/home/ha/y1/cahier_labo/241126/citeseq_wk4/fastq/MAP17-4wks-cDNA3/   MAP17-4wks-cDNA3   Gene Expression
/home/ha/y1/cahier_labo/241126/citeseq_wk4/fastq/MAP17-4wks-cDNA4/   MAP17-4wks-cDNA4   Gene Expression
/home/ha/y1/cahier_labo/241126/citeseq_wk4/fastq/MAP17-4wks-ADT/     MAP17-4wks-ADT     Antibody Capture
/home/ha/y1/cahier_labo/241126/citeseq_wk4/fastq/MAP17-4wks-Hashtag/ MAP17-4wks-Hashtag Antibody Capture

**feature_ref.csv**
id           name                read pattern sequence        feature_type
B220         B220_TotalA         R2   ^(BC)   CCTACACCTCATAAT Antibody Capture
CD199        CD199_TotalA        R2   ^(BC)   ATTCCTCATTCCTGA Antibody Capture
CD115        CD115_TotalA        R2   ^(BC)   TTCCGTTGTTGTGAG Antibody Capture
CD11b        CD11b_TotalA        R2   ^(BC)   TGAAGGCTCATTTGT Antibody Capture
CD11c        CD11c_TotalA        R2   ^(BC)   GTTATGGACGCTTGC Antibody Capture
CD127        CD127_TotalA        R2   ^(BC)   GTGTGAGGCACTCTT Antibody Capture
CD135        CD135_TotalA        R2   ^(BC)   GTAGCAAGATTCAAG Antibody Capture
CD14         CD14_TotalA         R2   ^(BC)   AACCAACAGTCACGT Antibody Capture
CD150        CD150_TotalA        R2   ^(BC)   CAACGCCTAGAAACC Antibody Capture
CD16_32      CD16_32_TotalA      R2   ^(BC)   TTCGATGCTGGAGCA Antibody Capture
CD3          CD3_TotalA          R2   ^(BC)   GTATGTCCGCTCGAT Antibody Capture
CD34         CD34_TotalA         R2   ^(BC)   GATTCCTTTACGAGC Antibody Capture
CD49d        CD49d_TotalA        R2   ^(BC)   AACAAGACCCTTGAG Antibody Capture
CD41         CD41_TotalA         R2   ^(BC)   ACTTGGATGGACACT Antibody Capture
CD48         CD48_TotalA         R2   ^(BC)   AGAACCGCCGTAGTT Antibody Capture
CD62L        CD62L_TotalA        R2   ^(BC)   TGGGCCTAAGTCATC Antibody Capture
CD71         CD71_TotalA         R2   ^(BC)   ACCGACCAGTAGACA Antibody Capture
CD8          CD8_TotalA          R2   ^(BC)   TACCCGTAATAGCGT Antibody Capture
CD117        CD117_TotalA        R2   ^(BC)   TGCATGTCATCGGTG Antibody Capture
CX3CR1       CX3CR1_TotalA       R2   ^(BC)   CACTCTCAGTCCTAT Antibody Capture
CD201        CD201_TotalA        R2   ^(BC)   TATGATCTGCCCTTG Antibody Capture
Ly6C         Ly6C_TotalA         R2   ^(BC)   AAGTCGTGAGGCATG Antibody Capture
Ly6G         Ly6G_TotalA         R2   ^(BC)   ACATTGACGCAACTA Antibody Capture
NK1.1        NK1.1_TotalA        R2   ^(BC)   GTAACATTACTCGTC Antibody Capture
Sca-1        Sca-1_TotalA        R2   ^(BC)   TTCCTTTCCTACGCA Antibody Capture
CD105        CD105_TotalA        R2   ^(BC)   TATCCCTGCCTTGCA Antibody Capture
isotype_ctrl isotype_ctrl_TotalA R2   ^(BC)   AAGTCAGGTTCGTTT Antibody Capture
HTO1         HTO1_TotalA         R2   ^(BC)   ACCCACCAGTAAGAC Antibody Capture
HTO2         HTO2_TotalA         R2   ^(BC)   GGTCGAGAGCATTCA Antibody Capture
HTO3         HTO3_TotalA         R2   ^(BC)   CTTGCCGCATGTCAT Antibody Capture
HTO4         HTO4_TotalA         R2   ^(BC)   AAAGCATTCTTCACG Antibody Capture
HTO5         HTO5_TotalA         R2   ^(BC)   CTTTGTCTTTGTGAG Antibody Capture
HTO6         HTO6_TotalA         R2   ^(BC)   TATGCTGCCACGGTA Antibody Capture

## Install Cellranger
I download Cellranger version 7 in their website
Make sure that you well installed it by 
'''
$ which cellranger
'''

## Ref data 
Another thing required by cellranger is reference, you can download and extract them according to your model

in my **mm10**:
.
├── fasta
│   ├── genome.fa
│   └── genome.fa.fai
├── genes
│   └── genes.gtf.gz
├── reference.json
└── star
    ├── Genome
    ├── SA
    ├── SAindex
    ├── chrLength.txt
    ├── chrName.txt
    ├── chrNameLength.txt
    ├── chrStart.txt
    ├── exonGeTrInfo.tab
    ├── exonInfo.tab
    ├── geneInfo.tab
    ├── genomeParameters.txt
    ├── sjdbInfo.txt
    ├── sjdbList.fromGTF.out.tab
    ├── sjdbList.out.tab
    └── transcriptInfo.tab



## Launch the code 
Then we can launch our command : 
'''
cellranger count --id "wk4_count" --transcriptome path/to/mm10/ --feature-ref citeseq_wk4/feature_ref.csv --libraries citeseq_wk4/samplesheet.csv
'''

## Output
its will take 5 hours to run this code in plafrim 
the outputs will be in wk4_count/ in the same directory where you launch the code
**wk4_count**

├── SC_RNA_COUNTER_CS
└── outs

the results is in outs/ 
.
├── aggregate_barcodes.csv
├── **analysis**
├── cloupe.cloupe
├── feature_reference.csv
├── **filtered_feature_bc_matrix**
├── filtered_feature_bc_matrix.h5
├── metrics_summary.csv
├── molecule_info.h5
├── possorted_genome_bam.bam
├── possorted_genome_bam.bam.bai
├── **raw_feature_bc_matrix**
├── raw_feature_bc_matrix.h5
└── web_summary.html
 
The folder is in bold 
So depend of your purpose, you can choose the file .h5 to analyze in downtream, whether raw or filtered 



