# 3pHLA-score: improved structure-based peptide-HLA binding affinity prediction

**Motivation:** Binding affinity prediction of peptide-ligands to Human Leukocyte Antigen (HLA) receptors is an important step in immunotherapy research. Binding of peptides to HLAs is a prerequisite for triggering immune response and therefore crucial for peptide target identification and epitope discovery pipelines, which can be optimized with computational methods. Currently, most computational methods are limited because they rely exclusively on sequence-based data. To overcome this limitation, we propose a novel protocol for using peptide-HLA (pHLA) structures to predict the binding affinity of peptides to HLAs.
**Results:** We apply a non-linear per-peptide-position structure-based training approach to the Rosetta ref2015 scoring function towards the development of a custom score that we call 3pHLA-score. We train 28 per-allele Random Forest Regression models on 77,581 modeled pHLA structures. The 3pHLA-score outperforms widely used scoring functions (AutoDock, Vina, Dope, Vinardo, FoldX, GradDock) in a structural virtual screening task, and shows the ability to generalize well. Finally, the proposed approach outperforms the standard tuning of docking scoring functions on the pHLA system.
**Repo:** In this repo we provide the code used to run the experiments and training as well as trained models which constitute the 3pHLA-score.

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

