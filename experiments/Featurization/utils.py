import pyrosetta

# util functions for trating the data
def init_rosetta():
    pyrosetta.init()

def extract_fullc_energies(file_name):
    scorefxn=pyrosetta.get_fa_scorefxn()
    #load rosetta
    pose = pyrosetta.pose_from_pdb(file_name)
    scorefxn(pose)
    #get energies
    res_ene = pose.energies().total_energies_array()
    return res_ene

def extract_ppp_energies(file_name, pep_len):
    scorefxn=pyrosetta.get_fa_scorefxn()
    #load rosetta
    pose = pyrosetta.pose_from_pdb(file_name)
    scorefxn(pose)
    #get energies
    res_ene = pose.energies().residue_total_energies_array()
    peptide_ene = res_ene[-pep_len:]
    return peptide_ene


