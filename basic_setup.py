#! /usr/bin/python3

"""
    Initial setup a structure for Energy evaluation
"""
import argparse
import sys

from Bio.PDB.PDBParser import PDBParser
from Bio.PDB.NACCESS import NACCESS_atomic
from forcefield import VdwParamset
import energies as en

# The specific PATH to naccess binary (in soft) is needed
NACCESS_BIN = '/home/gelpi/DEVEL/BioPhysics/2020-21/soft/NACCESS/naccess'

parser.add_argument(
    '--vdw',
    action='store',
    dest='vdwprm_file',
    default='data/vdwprm',
    help='Vdw parameters'
)

parser.add_argument(
    '--diel',
    action='store',
    dest='diel_type',
    type=int,
    default=2,
    help="Dielectric type: 0: vac, 1: wat, 2: M_S"
    )

parser.add_argument('pdb_file',help='Input PDB', type=open)
parser.add_argument('pdbqt_file',help='Input PDBQT', type=open)

args = parser.parse_args()

print("PDB File:", args.pdb_file.name)
print("PDBQT File:", args.pdbqt_file.name)
print("Parameters File:", args.vdwprm_file)

# loading VdW parameters
ff_params = VdwParamset(args.vdwprm_file)

parser = PDBParser(PERMISSIVE=1)
print('Parsing PDB', args.pdb_file.name)

# load structure from PDB file of PDB ifle handler
st = parser.get_structure('STR', args.pdb_file.name)

# assign data types, and charges from libraries
# We will use the xtra attribute in Bio.PDB.Atom to hold the new data
# Getting Charges and Atom type from PDBQT
print('Parsing PDBQT', args.pdbqt_file.name)
params=[{}]

for line in args.pdbqt_file:
    line = line.rstrip()
    params.append({'charge': line[69:76], 'type': line[77:].replace(' ','')})

total_charge = 0.
for at in st.get_atoms():
    at.xtra['atom_type'] = params[at.serial_number]['type']
    at.xtra['charge'] = float(params[at.serial_number]['charge'])
    at.xtra['vdw'] = ff_params.at_types[at.xtra['atom_type']]
    total_charge += at.xtra['charge']
print('Total Charge:', total_charge)

# Calculating surfaces
# Srf goes to .xtra field directly
srf = NACCESS_atomic(st[0], naccess_binary=NACCESS_BIN)

