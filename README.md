# PhyloFunc distance method
Wang and Li et al., PhyloFunc: Phylogeny-informed Functional Distance as a New Ecological Metric for Metaproteomic Data Analysis

## This repository is a series of python codes to:
1) Generate PhyloFunc to incorporate microbiome phylogeny to inform on metaproteomic functional distance.
2) preprocess data,compute PhyloFunc ditance matrix, evaluate the performance by machine learning-based classification algorithms.

## Original dataset and example output
1) The original dataset is given in "./1_datasets_search_results/". Decompress any .zip files.
2) Example outputs are given in "./2_data_processing/".

## System requirements:
1) Python (version >= 3.11.0) and packages pandas (version >= 2.0.3) and numpy(version >= 1.24.3) were required. R (version >= 4.2.0) and Rstudio were required.
2) The codes have been tested on Python (version = 3.11.5) with packages pandas (version = 2.0.3) and numpy(version = 1.24.3) ran by Spyder 5.4.3 on a Windows computer.
3) The time of computing PhyloFunc distance for the human gut microbiome dataset is the most time-consuming, approximately 20-30 min on a Windows System (version = 11) for 48 sample pairs. The calculation of PhyloFunc for the toy dataset and mouse gut microbiome dataset only takes a few minutes because the sample size is relatively small.

## How to run:
1) Download data, example outputs and codes in this repository.
2) Open each .py code fiel: PhyloFunc_toy_dataset.py, PhyloFunc_mouse_gut_dataset.py, or PhyloFunc_human_gut_dataset.py using Spyder.
3) Follow the steps to run through preprocessing raw data, generating branches of tree, creating taxon-function table, and calculating the PhyloFunc distance. The output files will be created automatically during running the code to verify the calculation process of the program.
4) Besides, R codes for computing other three distances (./1_PhyloFunc calculation/3_human gut microbiome and 1_PhyloFunc calculation/2_mouse gut microbiome/four cases/original dataset) and Python codes for evaluating the performance of distances (./2_data_processing/2_Evaluation/) are provided.
