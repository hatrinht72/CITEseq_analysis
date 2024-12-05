# Primary Analysis
## Raw data & Hierarchy 
The experiment performed by 10X Genomics. 
Cells were isolated then tagged by Hashtag and ADT, pooled then divided into 4 sequencing runs. 
Assuming I was already in my citeseq_wk4 folder :

In **fastq**, I have 5 folders, corresponding to 4 sequencing runs, ADT and Hashtag
I also have my feature_ref.csv and my samplesheet.csv where I specify my the list of tags and also which fastq correspond to which type of data (Gene Expression or Antibody Capture) 
![image](https://github.com/user-attachments/assets/a2dcb03c-fe9f-4551-80fe-329b70544a4a)

## Install Cellranger
I download Cellranger version 7 in their website
Make sure that you well installed it by 
'''
$ which cellranger
'''

## Ref data 
Another thing required by cellranger is reference, you can download and extract them according to your model

in my **mm10**:
![image](https://github.com/user-attachments/assets/5dc677cd-ed16-45f6-a421-b98b618d3241)


## Launch the code 
Then we can launch our command : 
'''
cellranger count --id "wk4_count" --transcriptome path/to/mm10/ --feature-ref citeseq_wk4/feature_ref.csv --libraries citeseq_wk4/samplesheet.csv
'''

## Output
its will take 5 hours to run this code in plafrim 
the outputs will be in wk4_count/ in the same directory where you launch the code
**wk4_count**

![image](https://github.com/user-attachments/assets/e8dcaebe-beb5-475c-9106-524b2e4373e3)

the results is in outs/ 

![image](https://github.com/user-attachments/assets/0d2c3a7e-81ac-471e-8299-809464da5c83)


So depend of your purpose, you can choose the file .h5 to analyze in downtream, whether raw or filtered 



