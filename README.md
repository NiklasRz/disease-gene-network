# disease-gene-network

## About

This is a script to investigate the networks of genes, diseases and drugs as part of an assigment for the Complex Systems 2 exercise at the Medical University of Vienna.

The assignment is as follows:

### 1
- What are the macroscopic properties of the gene-gene and the disease-disease networks?
- Which two genes are involved in the largest number of disorders?

### 2
- how does the disease-disease network as a projection of the disease-drug network look like?
- how does the drug disease-disease network relate to or differ from the genetic disease-disease  comorbidity network?

## To run:
- create a python virtual environment (optional but recommended)
- install the requirements with "pip3 install -r requirements.txt"
- run with "python3 code/run.py"
- results are in ./results

## Data source
wget -O data/diseases.tsv http://www.complex-systems.meduniwien.ac.at/lectures/DiseaseGeneNet.txt  
wget -O data/drugs.tsv https://www.dgidb.org/data/monthly_tsvs/2022-Feb/interactions.tsv
