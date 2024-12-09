![PhyloFunc_logo2_smallest](https://github.com/user-attachments/assets/90bd7067-c3cb-40dd-a640-8bff9402732c)

# PhyloFunc distance method
Wang and Li et al., PhyloFunc: Phylogeny-informed Functional Distance for Metaproteomics.

## This repository is a series of python codes to:
1. Generate PhyloFunc to incorporate microbiome phylogeny to inform on metaproteomic functional distance.
2. Compute PhyloFunc ditance matrix and visualzie using PCoA.

## System requirements:
1. Python (version >= 3.11.0) and packages pandas (version >= 2.0.3) and numpy(version >= 1.24.3) were required. R (version >= 4.2.0) and Rstudio were required.
2. The codes have been tested on Python (version = 3.11.5) with packages pandas (version = 2.0.3) and numpy(version = 1.24.3) ran by Spyder 5.4.3 on a Windows computer.
3. Computing the PhyloFunc distance for a sample pair from the human gut microbiome dataset takes nearly 1 minute, while generating the entire distance matrix is more time-consuming, approximately 1 hour for 48 samples on a Windows System (version = 11). In contrast, calculating the PhyloFunc distance matrix for the toy dataset and the mouse gut microbiome dataset takes only a few minutes due to their smaller sample sizes.

## How to run:
1. Download all files including PhyloFunc_package_tutorial, datasets, scripts, and results folders in this repository.
2. The PhyloFunc package tutorial and its corresponding data are in the ./1_PhyloFunc_package_tutorial file folder. Run the code /PhyloFunc_Package_Tutorial.ipynb by Jupyter Notebook to calculate PhyloFunc distance and PhyloFunc distance matrix, and to generate hierarchical clustering result.
3. Three datasets are provided: a simulated toy dataset and two real datasets from mouse and human gut microbiomes. The corresponding Taxon-Function tables are located in ./2_datasets. Code for calculating PhyloFunc distances and visualizing PCoA results is available in ./3_scripts. Distance matrices generated by four different methods are stored in ./4_results.
