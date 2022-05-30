from Bio.PDB import PDBParser, PDBIO
import pyrosetta
import numpy as np
import pickle as pk
import pathlib
from score3pHLA.utils import *

"""
3pHLA-score library

Use to score structures of pHLA complexes. 
"""

def score_rosetta(file_name: str,
                  pep_len: int) -> np.array:
    """ Score complex with PyRosetta. Get per-peptide-position energies.

    Call PyRosetta and score the complex with the rosetta scores.

    Args:
        file_name (str): location to the PDB file

    Returns:
        peptide_ene_ar (np.array): numpy array with per-peptide-position energies

    Raises:
        FileNotFoundException: the input file is wrong
        AlleleNotSupportedException: the input allele is not supported.
    """

    # check file
    res = file_pdb_check(file_name)
    if res[0] == False or res[1] == False:
        raise FileNotFoundError("Can not find file {}, \
        or it does not have PDB extension".format(file_name))
    
    # load structure in pyrosetta
    scorefxn=pyrosetta.get_fa_scorefxn()
    pose = pyrosetta.pose_from_pdb(file_name)
    scorefxn(pose)

    # calculate and convert energies
    res_ene = pose.energies().residue_total_energies_array()
    peptide_ene = res_ene[-pep_len:]
    peptide_ene_ar = pppene_to_array(str(peptide_ene))

    return peptide_ene_ar



def score(pdb_loc: str, 
               allele: str) -> float:
    """ Score complex with 3pHLA.

    Call PyRosetta and score the complex to get per-peptide-position features.
    Feed these features into the trained 3pHLA-score.
    Get the predicted score for the given pHLA complex.

    Args:
        file_name (str): location to the PDB file
        allele (str): allele to process

    Returns:
        pred (float): predicted 3pHLA-score (range 0-1 as scaled nM value)

    Raises:
        FileNotFoundException: the input file is wrong
        AlleleNotSupportedException: the input allele is not supported
    """

    # check input allele
    if not allele_check(allele):
        raise AlleleNotSupportedError(allele)
        res = file_pdb_check(file_name)
    
    # check input file
    res = file_pdb_check(pdb_loc)
    if res[0] == False or res[1] == False:
        raise FileNotFoundError("Can not find file {}, \
        or it does not have PDB extension".format(pdb_loc))
    
    # get peptide length
    pdb_str = PDBParser().get_structure("input", pdb_loc)
    # get the length of the smallest chain = peptide length
    pep_len = min([len(list(chain.get_residues())) for chain in pdb_str.get_chains()])

    print("Scoring {} for allele {} \
        and peptide of length {}".format(pdb_loc, allele, pep_len))
    print("-------------------------------------------")
    print("Step 1: extracting 3p energies with pyrosetta")
    print("-------------------------------------------")

    init_rosetta()

    feat = score_rosetta(pdb_loc, pep_len)

    print("-------------------------------------------")
    print("Step 2: scoring the complex with 3pHLA-score")
    print("-------------------------------------------")
    
    # get the root directory
    root_dir = pathlib.Path.joinpath(pathlib.Path(__file__).parent.resolve(), "3pHLAmodels")
    # open the model
    model_9n = pathlib.Path.joinpath(root_dir, allele+"ppp.pkl")
    model_9 = pk.load(open(model_9n, 'rb'))

    # predict the score
    pred = model_9.predict([get_energies(feat)])[0]
    print("-------------------------------------------")
    print("Final score: "+str(pred))
    print("-------------------------------------------")

    # interpret the score
    if pred > 0.5:
        print("predicted as: strong binder (>0.5)")
    print("-------------------------------------------")
    if pred <= 0.1:
        print("predicted as: non-binder (<=0.1)")
    return pred


