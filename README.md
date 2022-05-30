# 3pHLA-score: improved structure-based peptide-HLA binding affinity prediction
 
 ![3pHLA-score and other scoring protocols visualized](https://github.com/KavrakiLab/3pHLA-score/blob/main/FigMethods.png?raw=true)

**Motivation:** Binding affinity prediction of peptide-ligands to Human Leukocyte Antigen (HLA) receptors is an important step in immunotherapy research. Binding of peptides to HLAs is a prerequisite for triggering immune response and therefore crucial for peptide target identification and epitope discovery pipelines, which can be optimized with computational methods. Currently, most computational methods are limited because they rely exclusively on sequence-based data. To overcome this limitation, we propose a novel protocol for using peptide-HLA (pHLA) structures to predict the binding affinity of peptides to HLAs.


**Results:** We apply a non-linear per-peptide-position structure-based training approach to the Rosetta ref2015 scoring function towards the development of a custom score that we call 3pHLA-score. We train 28 per-allele Random Forest Regression models on 77,581 modeled pHLA structures. The 3pHLA-score outperforms widely used scoring functions (AutoDock, Vina, Dope, Vinardo, FoldX, GradDock) in a structural virtual screening task, and shows the ability to generalize well. Finally, the proposed approach outperforms the standard tuning of docking scoring functions on the pHLA system.


**Repo:** In this repo we provide the code used to run the experiments and training as well as trained models which constitute the 3pHLA-score.

## Installation

### Step 1: Install PyRosetta (v2021.07)

Obtain the free academic license for PyRosetta by applying over this link: https://els2.comotion.uw.edu/product/pyrosetta (it's free!). 
You will recieve a username and password for downloading PyRosetta. 
Now you can install PyRosetta using conda (replace the USERNAME and PASSWORD with the ones you received):

```sh
conda install pyrosetta=2021.07 --channel https://USERNAME:PASSWORD@conda.rosettacommons.org
```

Alternalively, if you don't want to go through conda, you can download the python wheel directly and build it following the instructions provided here: https://www.pyrosetta.org/downloads/legacy-pyrosetta3-download 

### Step 2: Install 3pHLA-score

```python
pip install score3pHLA
```

### Step 3: Check out the command line tool

Try out the command line tool by running:

```sh
score3pHLA -h
```

Or running the simple example

```sh
score3pHLA -e
```

### Step 4: Check out the python module

Enter python shell. Type:
```python
from score3pHLA import score
```

## Usage

Try it out on your own file. Replace "your_pHLA_pdb_file.pdb" with the location to your pHLA and "A0201" with your own allele.


```python
from score3pHLA import score

pdb_loc = "your_pHLA_pdb_file.pdb"
allele = "A0201"

score(pdb_loc, allele)
```

Also, check the simple usage in the jupyter notebook: [3pHLA-score:UseCase.ipynb](https://colab.research.google.com/drive/1QbENWaIE-r5AXvUv25IlVEWFPCklMKJ3?usp=sharing). 

## Code for experiments provided
The code used for training the score and running the experiments, along with the datasets can be found in the .\experiments directory.

**.\experiments\Datasets\\** contains all the data used to train the models as well as run the validation experiments. 
Jupyter notebooks used to clean and split the data are provided as well as the source data.

**.\experiments\Featurization\\** contains the code needed to generate ref2015 energies and extract the ref2015 features from PDB files. 
As the PDB files are not provided here, you can simply extract the features from the features.tar.gz file to continue running the experiments.
The pHLA models in PDB file format are available upon request.

**.\experiments\Experiment1 - train ref2015-score, standard-pHLA-score and 3pHLA-score\\** contains the code used to calculate and train the basic scores.
It contains the trained models in the directory final_REGRmodels.

**.\experiments\Results1\\** contains the code used to compare different models trained in the Experiment1 and the resulting figures.

**.\experiments\Experiment2 - Virtual screening\\** contains the code used for running the structural virtual screening validation experiment as well as related results and figures.


**.\experiments\Experiment3 - crystal structures etc** contains the code used for running the second validation experiment with a dataset of crystal structures and Docktope models, as well as related results and figures.



