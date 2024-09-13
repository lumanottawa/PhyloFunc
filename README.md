# PhyloFunc distance method
Wang and Li et al., PhyloFunc: Phylogeny-informed Functional Distance as a New Ecological Metric for Metaproteomic Data Analysis

## This repository is a series of python codes to:
1) Generate PhyloFunc to incorporate microbiome phylogeny to inform on metaproteomic functional distance.
2) preprocess data,compute PhyloFunc ditance matrix, evaluate the performance by PCoA and machine learning-based classification algorithms.

## Original dataset and example output
1) The tutorial and its corresponding data are given in "./0_PhyloFunk_package_tutorial".
2) The original dataset is given in "./1_dataset_results/".
3) All output files, including PhyloFunc calculations and distance evaluations, are given in "./2_data_processing/".

## System requirements:
1) Python (version >= 3.11.0) and packages pandas (version >= 2.0.3) and numpy(version >= 1.24.3) were required. R (version >= 4.2.0) and Rstudio were required.
2) The codes have been tested on Python (version = 3.11.5) with packages pandas (version = 2.0.3) and numpy(version = 1.24.3) ran by Spyder 5.4.3 on a Windows computer.
3) Computing the PhyloFunc distance for a sample pair from the human gut microbiome dataset takes nearly 1 minute, while generating the entire distance matrix is more time-consuming, approximately 1 hour for 48 samples on a Windows System (version = 11). In contrast, calculating the PhyloFunc distance matrix for the toy dataset and the mouse gut microbiome dataset takes only a few minutes due to their smaller sample sizes.

## How to run:
1) Download data, example outputs, and codes in this repository.
2) Open ./0_PhyloFunk_package_tutorial/PhyloFunc_Package_Tutorial.ipynb using Jupyter Notebook.
3) Follow the steps to run through calculating PhyloFunc distance, PhyloFunc distance matrix, and generating hierarchical clustering result. The output files will be created automatically during running the code to verify the calculation process of the program.
4) Furthermore, Python codes PhyloFunc_toy_dataset.py, PhyloFunc_mouse_gut_dataset.py, or PhyloFunc_human_gut_dataset.py are provided for calculating PhyloFunc distance matrix for toy dataset, the mouse gut microbiome dataset, and human gut microbiome dataset.
5) Besides, R codes for computing the other three distances (./1_PhyloFunc calculation/3_human gut microbiome and 1_PhyloFunc calculation/2_mouse gut microbiome/four cases/original dataset) and Python codes for evaluating the performance of distances (./2_data_processing/2_Evaluation/) are provided.