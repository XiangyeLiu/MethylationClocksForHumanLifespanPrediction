# Methylation Clocks for Human Lifespan Prediction
- Codes for data analyses in the study "Epigenetic estimate of the lifespan limit in modern and archaic humans" and several newly-trained PC-guided clocks are available here.
- For .R or .py files, please refer to code to find out what form of data is needed.
- Codes to train and use new PC-guided clocks are revised from codes provided by "A computational solution for bolstering reliability of epigenetic clocks, implications for clinical trials and longitudinal tracking" (https://doi.org/10.1038/s43587-022-00248-2), which are available at https://github.com/MorganLevineLab/PC-Clocks.
## Software and Operating System Requirements
- R or R Studio for R codes, any version that support R packages used is fine (4.5.0 in the study).
- Python3, any version that support Pandas is fine (3.12.1 in the study).
- Any operating systems in which R and Python3 can run normally.
## General Instruction
- Codes can be downloaded or copied, then run directly, with proper forms of data provided as input. Otherwise, some code modifying may be required. 
- Data series used in the study can be obtained from GEO database at https://www.ncbi.nlm.nih.gov/geo/ or ArrayExpress database at https://www.ebi.ac.uk/biostudies/arrayexpress. Please check the supplementary information for details. Ohter data series availible on the two databases or from other sources can also be used to test codes.
- For reproduction, please use the same data series and process them properly, then run codes step by step in the correct order. Maximum lifespan estimation results of each methylation clock used are expected to be the final output.
