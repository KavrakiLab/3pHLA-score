import argparse
import os.path
import pathlib
import sys
from score3pHLA.utils import *
from score3pHLA import score

def main():

   # Instantiate the parser
   parser = argparse.ArgumentParser(description='3pHLA-score for scoring structures of pHLA complexes')
   
   parser.add_argument('-a', '--allele', type=lambda a: allele_check_p(a, parser), 
   required=(not ("--example" in sys.argv or "-e" in sys.argv)), \
   help='3pHLA-score is still per-allele. \
   Please specify which allele you will input. \
   Supported alleles are: A0101, A0201, A0301, A1101, A2402, \
   A2902, B0702, B0801, B1501, B2705, B3501, B4001, B4002, B4403, \
   B5101, B5701. Allele is a required argument.')
   
   parser.add_argument('-i', '--input', type=lambda i: file_pdb_check_p(i, parser), 
   required=(not ("--example" in sys.argv or "-e" in sys.argv)), \
   help='Please specify the location of your input structure. \
   The structure should be in a PDB format.')

   parser.add_argument('-e', '--example', action='store_true', \
   required=False, help='Run the example file.')

   args = parser.parse_args()

   if args.example:
      file_root = pathlib.Path.joinpath(pathlib.Path(__file__).parent.resolve(), "test_pdb_A0201.pdb")
      score(str(file_root), "A0201")
   else:
      score(args.input, args.allele)


if __name__ == "__main__":
   main()
   
