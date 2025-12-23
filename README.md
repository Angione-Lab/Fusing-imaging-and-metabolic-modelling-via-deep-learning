
# Fusing imaging and metabolic modelling via multimodal deep learning in ovarian cancer

This repository contains the code and data necessary to reproduce the results presented in the paper:

**N. Eftekhari, S. Verma, A. Saha,  G. Zampieri, S. Sawan, A. Occhipinti, C. Angione, "Fusing imaging and metabolic modelling via multimodal deep learning in ovarian cancer, Cell Systems, 2026"**

<img src="Images/Fig1-pipeline.jpg" alt="Pipeline Image" width="720" style="display: block; margin-left: auto; margin-right: auto;"/>

---

## Table of Contents

- [Overview](#overview)
- [Requirements](#requirements)
- [Installation](#installation)
- [Feature Reduction](#feature-reduction)
- [Image Segmentation](#image-segmentation)
- [Metabolic Modelling](#metabolic-modelling)
- [Survival Models](#survival-models)


---

## Overview

This repository includes a Jupyter Notebook, scripts, and data necessary to perform RNA feature extraction and image segmentation for ovarian cancer research. The notebook supports the integration of radiogenomics and metabolic data using deep learning, enabling personalized oncology research.

## Requirements

The repository requires the following Python libraries:

- `numpy`
- `pandas`
- `matplotlib`
- `scikit-learn`
- `nibabel`
- `SimpleITK`
- `tensorflow-gpu==2.1` (or `tensorflow-mkl==2.1` for CPU)

Install these dependencies using the following command:

```bash
pip install -r requirements.txt
```

## Installation
```bash
git clone https://github.com/Angione-Lab/Fusing-imaging-and-metabolic-modelling-via-deep-learning.git
cd Fusing-imaging-and-metabolic-modelling-via-deep-learning
```

## Feature Reduction

1. **Data Preparation**:

Ensure that your RNA data file is in the required format as specified in the notebook.
Place the data file in the appropriate directory, or update the file path in the notebook accordingly.

2. **Running the Notebook**:
 Launch the Jupyter Notebook server and for example open RNA_feature_vector.ipynb
```bash
jupyter notebook RNA_feature_vector.ipynb

```




## Image Segmentation

### Requirements
```bash
conda create -n ovarian_unet_py37 python=3.7
conda activate ovarian_unet_py37

conda install tensorflow-gpu==2.1
```
In case of error:

```bash

conda install tensorflow-estimator=2.1.0

conda install tensorflow-mkl==2.1 (cpu version)
pip install nibabel SimpleItk matplotlib tqdm
```

Training: 
```bash
train/train.py.py
```

Prediction: 
```bash
train/ovarian_unet_predict.py
```


## Metabolic Modelling

This code uses a metabolic model and transcriptomic data to set reaction bounds in accordance with gene expression levels. The code loads a predefined metabolic model (`human1.mat`) and adjusts reaction bounds for specific reactions using gene expression profiles.

Key components:
1. **Loading the Model**: `human1.mat` model is loaded.
2. **Setting Bounds**: Bounds for selected reactions are modified based on user-defined gamma and threshold values.
3. **Transcriptomic Data Integration**: Gene expression data for each patient is mapped to the model.
4. **Flux Analysis**: Reaction fluxes are calculated under different threshold and gamma values.

### Requirements

- **MATLAB** (with COBRA Toolbox installed)
- Metabolic model file `human1.mat`
- Gene expression data files: `gene_exp.mat`, `gene_ids.mat`, and `patient_ids.mat`

Ensure the COBRA Toolbox is installed and properly configured. 

## Survival Models
The notebook contains steps to load and preprocess data, merge multi-omics datasets, and apply several machine learning algorithms for survival analysis. Key components of the analysis include data normalization, survival model fitting, and model evaluation using concordance index.


### Requirements

- Python 3.x
- NumPy
- Pandas
- Seaborn
- Matplotlib
- Scikit-learn
- scikit-survival
- Lifelines











