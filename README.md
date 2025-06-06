# RNAseq_Rshiny

This Shiny app allows interactive analysis of RNA-seq data from featureCounts `.rds` files.

## Features (updated as more are added)

- Upload `.rds` count matrix and metadata file in form of `SampleInfo.txt`
- Perform DESeq2 analysis with customizable contrasts
- Visualize PCA, volcano plots, GO analysis
- Clickable plots and real-time feedback

## Usage

1. Upload your files
2. Select groups and contrast
3. Click "Run Analysis"

## Example files

This folder contains 2 example files in the correct format which should be plug-and-play with this app. 

`fc.rds` contains count data from a Rat RNAseq experiment with 2 groups - sham and TBI (from Zhu, Xiaolu, Jin Cheng, Jiangtao Yu, Ruining Liu, Haoli Ma, and Yan Zhao. ‘Nicotinamide Mononucleotides Alleviated Neurological Impairment via Anti-Neuroinflammation in Traumatic Brain Injury’. International Journal of Medical Sciences 20, no. 3 (2023): 307–17. https://doi.org/10.7150/ijms.80942.
).

`SampleInfo.txt` contains metadata about this experiment in the expected format for the app to run properly.
