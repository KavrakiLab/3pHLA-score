# 3pHLA-score

3pHLA-score is a set of Random Forest Regression models trained to predict the binding affinity
of peptides to HLAs given an input structure of the complex. 
Try out our 3pHLA-score to perform structural scoring of pHLA structures!

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

