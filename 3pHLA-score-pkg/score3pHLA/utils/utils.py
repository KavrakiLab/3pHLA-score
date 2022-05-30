"""
Util functions for running 3pHLA-score.
"""

import numpy as np
import pyrosetta
import argparse
import os

# List of supported alleles
supported_alleles = ['A0101', 'A0201', 'A0301', 'A1101', 'A2402', 'A2902', 'B0702',
       'B0801', 'B1501', 'B2705', 'B3501', 'B4001', 'B4002', 'B4403',
       'B5101', 'B5701', 'C0304', 'C0501', 'C1601']

def allele_check(allele: str) -> bool:
    """ Checks if allele is supported.

    Looks for the given allele in the list of supported alleles.

    Args:
        allele (str): input allele
    
    Returns:
        ret (bool): False if allele is not supported
    """

    ret = True
    if allele not in supported_alleles:
        ret = False
    return ret

def allele_check_p(allele: str,
                    parser: argparse.ArgumentParser) -> str:
    """ Allele check for the command-line parser

    Check allele and notify the parser.

    Args:
        allele (str): input allele
        parser (argparse.ArgumentParser): argument parser

    Returns:
        allele (str): if passes the check, returns allele name
    """

    ret = allele_check(allele)

    if not ret:
        parser.error(\
        "Allele not supported. 3pHLA-score only supports: {}".format(", ".join(supported_alleles))\
        )
        return None
    else:
        return allele


def file_pdb_check(fname: str) -> (bool, bool):
    """ Checks if the file exists and if it has the PDB extension.

    Looks for the file with the PDB extension.

    Args:
        fname (str): location to the PDB file

    Returns:
        ret ((bool, bool)): first bool True if file exists, 
                            second bool True is file has PDB extension
    """
    ret = (True, True)
    ext = os.path.splitext(fname)[1][1:]
    if not os.path.exists(fname):
        ret = (False, ret[1])
    if ext not in ["pdb", "PDB"]:
        ret = (ret[0], False)
    return ret


def file_pdb_check_p(fname, parser):
    """ Check the file format and location and notify parser.

    Args:
        fname (str): input file name
        parser (argparse.ArgumentParser): argument parser

    Returns:
        fname (str): if the checks are passed return file name.
    """

    res = file_pdb_check(fname)

    if not res[0]:
       parser.error("File {} does not exist".format(fname))
       return None
    if not res[1]:
       parser.error("Your input file is not a pdb file")
       return None
    return fname

def init_rosetta():
    """ Initializes PyRosetta.

    Initialize PyRosetta before running analysis with it.
    """

    pyrosetta.init()


def pppene_to_array(str_res: str) -> np.array:
    """ Convert PyRosetta energy string output to a numpy array.

    Take the string as input. Remove all string parenthesis.
    Convert to numpy array. Return the array.

    Args:
        str_res (str): Input string - result of the PyRosetta energy evaluation.

    Returns:
        ret (numpy.array): Converted result array with energy values.

    Raises:
        TypeError: an error occured reformatting the array.
    """

    str_res = str_res.replace("(", "")
    str_res = str_res.replace(")", "")
    str_res = str_res.strip("[]")
    str_res = str_res.replace(" ", "")
    str_res = str_res.replace("\n", ",")
    ret = np.fromstring(str_res, dtype=float, sep=", ")
    if not ret.size == 9*20:
        raise TypeError("Rosetta energy array of wrong dimenstions! Array should have 180 elements.")
    ret = ret.reshape(9,20)
    return ret


def get_energies(X: np.array) -> np.array:
    """ Prepare PyRosetta energies vector for scoring with 3pHLA-score.

    Take the pre-generated vector. 
    Perform the rolling of the vector to extract only nine positions,
    excluding the middle.
    Take only 19 energy terms of the ref2015 energy function.
    Result is the 9x19 vector fit for the 3pHLA-score models.

    Args:
        X (numpy.array): Input array of raw PyRosetta scoring output

    Returns:
        ene (numpy.array): Output array of dimensions  9x19 for scoring with 3pHLA-score

    Raises:
        TypeError: Peptide must have at least 9 amino acids. Input array must have at least 9*19 elements.
    """

    if X.size < 9*19:
        raise TypeError("Peptide must have at least 9 amino acids. Input array must have at least 9*19 elements.")

    ene = np.roll(X, 4, axis = 0)[:9,:19]
    ene = np.roll(ene, -4, axis = 0)
    ene = ene.reshape(9*19)
    return ene

class AlleleNotSupportedError(Exception):
    """ Exception for unsuported allele.
    """
    
    def __init__(self, allele:str):  
        message = "Allele {} not supported.\n \
        Try one of the supported alleles: {}".format(allele, supported_alleles)         
        super().__init__(message)