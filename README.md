# AITeQ (Alzheimer's Identification Tool using RNA-seq) 
![AITeQ_logo](https://github.com/ishtiaque-ahammad/AITeQ/assets/99262870/115e7b7d-0529-4d3e-942a-0326739640ca)


## Introduction
Welcome to AITeQ (Alzheimer's Identification Tool using RNA-seq), a powerful and user-friendly tool that has been developed as part of an integrative machine learning and transcriptomics study. AITeQ is designed for prediction of Alzheimer's disease from RNA-Seq data. This tool takes as input a .csv file containing the output of variance stabilizing transformation (VST) from DESeq2, and it leverages machine learning algorithms to provide accurate predictions of Alzheimer's disease status based on a distinctive 5-gene signature.

## Motivation
Alzheimer's disease is a devastating neurodegenerative disorder that affects millions of individuals worldwide. Accurate diagnosis of the diease is crucial for effective treatment and intervention. AITeQ aims to bridge the gap between RNA-Seq data and Alzheimer's disease diagnosis. By identifying a distinctive 5-gene signature through an integrative approach, we hope to aid researchers and clinicians in the field with accurate detection.

## Features
### Easy Integration: 
AITeQ seamlessly integrates with DESeq2, a widely used tool for differential gene expression analysis, making it easy to transition from data preprocessing to disease prediction.

### Machine Learning: 
AITeQ employs sophisticated machine learning algorithms, guided by a distinctive 5-gene signature, to analyze RNA-Seq data and make predictions about Alzheimer's disease status.

### User-Friendly: 
The tool is designed with user-friendliness in mind, with a simple and intuitive google colab interface that allows users to quickly analyze their data using any operting system.

### Rapid Prediction: 
AITeQ is optimized for rapid predictions, allowing researchers and clinicians to quickly obtain results, which can be crucial for timely diagnosis and intervention.

### Lightweight:
AITeQ is designed to be lightweight and efficient, ensuring that it can easily run on a wide range of computational platforms.

## Usage

### Input
AITeQ takes in .csv file containing the output of variance stabilizing transformation (VST) from DESeq2 (https://bioconductor.org/packages/release/bioc/html/DESeq2.html).

### Output
Once AITeQ is finished running, the Alzheimer's predictions is displayed as a dataframe.

### Step 1: 
Open [AITeQ](https://colab.research.google.com/github/ishtiaque-ahammad/AITeQ/blob/main/AITeQ.ipynb)

### Step 2: 
Upload your input file

### Step 3: 
Make predictions and view the results

## Examples
For understanding how to prepare your input file, please use the example_dataset.csv file in the AITeQ repository.

## Contact
If you have any questions, feedback, or issues, please don't hesitate to contact us at bioinformatics.division.nib.gov.bd@gmail.com

## License
This project is licensed under the GPL-3.0 license.
