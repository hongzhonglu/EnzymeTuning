# EnzymeTuning

EnzymeTuning is a computational framework developed to optimize kinetic parameters ($k_{cat}$) in enzyme-constrained genome-scale metabolic models (ecGEMs).

## Table of Contents
* [Installation](#installation)
* [Data Requirements](#data-requirements)
* [Usage Guide](#usage-guide)

## Installation

### 1. Environment Setup
It is recommended to use `conda` to manage the environment:

```bash
# Create and activate the environment
conda create -n enzymetuning python=3.8
conda activate enzymetuning

# Clone the repository
git clone [https://github.com/hongzhonglu/EnzymeTuning.git](https://github.com/hongzhonglu/EnzymeTuning.git)
cd EnzymeTuning

# Install dependencies
pip install -r requirements.txt
```

## Data Requirements


### 1. Metabolic Model: 
COBRA-compatible files (e.g., ecyeast.xml).

### 2. Enzyme Information: 
.mat or .txt files containing kinetic information and molecular weights.

### 3. Proteomics Datasets: 
.csv files with protein abundance data .

### 4. Growth Phenotypes: 
.xlsx or .csv files containing measured exchange rates (e.g., growthdata.xlsx).

## Usage Guide

### 1. Initial Sampling : Step1_data_preprocessing.py
Generates the initial 1,000 $k_{cat}$ samples from Bayesian priors.

### 2. Sensitivity Analysis : sensitivity analysis.py
Identifies the key metabolic reactions and associated enzyme kinetic parameters that significantly influence the model's objective function (e.g., growth rate) or protein abundance predictions. 
This step allows for a more focused optimization.

### 3. Comprehensive Tuning: Step2_gan_training.py or Step4_gan_all.py
Trains the conditional Generative Adversarial Network (cGAN) to learn high-performance parameter distributions.
Step2_gan_training1.py: Focuses on identifying parameter sets that improve the $k_{cat}$ for protein abundance predicitons.
Step4_gan_all.py: Performs an integrated optimization that balances protein abundance prediction, and phenotypic growth prediction ($Root$ $Mean$ $Square$ $Error$, RMSE).

### 4. Turnover Refinement (optional): Step3_turnover.py
Optimizes enzyme turnover rates ($v_{syn}$) in $S.cerevisiae$.

### 5. Cross-Strain Application (optional): Step5_gan_other_yeast_all.py
Applies the framework to other strains such as $K. marxianus$ or $Y. lipolytica$.

#### It is highly recommended to run these scripts on a High-Performance Computing (HPC) cluster.