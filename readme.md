
# DeepCodon: A Deep Learning Codon-optimization Model To Enhance Protein Expression
# Contents
- [DeepCodon: A Deep Learning Codon-optimization Model To Enhance Protein Expression](#deepcodon-a-deep-learning-codon-optimization-model-to-enhance-protein-expression)
- [Contents](#contents)
- [Introduction](#introduction)
- [Directory Structure](#directory-structure)
- [Installation](#installation)
- [Usage](#usage)
  - [Data Collection](#data-collection)
    - [EPCD](#epcd)
    - [RareCodon](#rarecodon)
  - [Deep Learning](#deep-learning)
    - [Training \& Fine-tuning](#training--fine-tuning)
    - [Inference](#inference)
# Introduction
In this study, we leveraged deep learning techniques to analyze 15 billion Enterobacteriaceae gene sequences from natural sources, constructing a translation model for Protein-CDS sequences. The model was further fine-tuned using high-expression gene data. By incorporating a conditional probability-based strategy to preserve conserved rare codon clusters, we developed DeepCodon, a codon optimization tool specifically designed for Escherichia coli. Compared to existing methods, DeepCodon offers three distinct advantages: generating optimized sequences that more closely align with host codon preferences, achieving superior computational evaluation metrics, and preserving rare codons at functionally critical sites. Finally, we experimentally validated the expression of seven poorly-expressed P450 enzymes and thirteen AI-generated G3PDH enzymes. Among the evaluated cases, outperforming traditional codon optimization methods in nine cases, demonstrating the potential of DeepCodon as a practical tool for codon optimization in E. coli.

We also offer web services here: https://deepcodon.biodesign.ac.cn

# Directory Structure
```
.
├── DeepCodon   --> DeepCodon-T&FT
│   ├── config  --> config path
│   ├── data    --> data path
│   ├── __init__.py
│   ├── model   --> model path
│   ├── readme.md
│   └── src --> source code
├── environment.yml
├── EPCD    --> Enterobacteriaceae sequence dataset
│   ├── finalData    --> data path
│   ├── readme.md
│   └── src -->source code
├── RareCodon   -->conserved rare codon clusters dataset
│   ├── finalData    --> data path
│   ├── readme.md
│   ├── slrums --> source code
│   └── src --> source code
└── readme.md
```
# Installation
```bash
git clone https://github.com/JiangLab2020/DeepCodon.git
cd DeepCodon
conda env create -f environment.yml
```
# Usage
## Data Collection
### EPCD
Refer to `EPCD/readme.md`
### RareCodon
Refer to `RareCodon/readme.md`
## Deep Learning
### Training & Fine-tuning
```bash
accelerate config
accelerate launch DeepCodon/src/train/main.py
```
### Inference
```bash
python DeepCodon/src/services/run.py
```
