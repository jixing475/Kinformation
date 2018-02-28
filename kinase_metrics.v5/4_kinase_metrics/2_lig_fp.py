#!/usr/bin/python

##########################################################################
#
#
#
#
#
#
#
##########################################################################

import sys,os,re
msg = '''\n  ## Usage: {0}
             [list of knownTemplate MOL] [Compare MOL]  (.smi, .sdf, multi-mol file ok)
             [Rank mol with Fingerprint: str] (ecfp4, dl, maccs, total)
             [No. Processor for Parallel: int]
             [Output Prefix]\n
             e.g.>   x.py a.smi b.sdf.bz2 ecfp4 4 output\n
             return: output.fp.txt (all scores)\n
'''.format(sys.argv[0])
if len(sys.argv) != 6: sys.exit(msg)


from rdkit import Chem
from rdkit import DataStructs
from rdkit_open import *

from rdkit.Chem.Fingerprints import FingerprintMols
from rdkit.Chem import MACCSkeys
from rdkit.Chem import AllChem
from CommonUtility import *
import multiprocessing
import pandas as pd


