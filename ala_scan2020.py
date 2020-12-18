#! /usr/bin/python3

"""
    Initial setup a structure for Energy evaluation
"""
import argparse
import sys
import os

from Bio.PDB.NACCESS import NACCESS_atomic
from forcefield import VdwParamset
from functions import *

#Upper limit for bound atoms
BOND_CUTOFF = 2.0
#Upper limit for int calculations
VDW_CUTOFF = 10.0
INT_CUTOFF = 25.0

# possible names for ala_atoms
ALA_ATOMS = {
    'N', 'H', 'H1','H2','H3',
    'CA', 'HA', 'HA1','HA2',
    'C', 'O', 'OXT',
    'CB', 'HB1', 'HB2', 'HB3'
}

parser = argparse.ArgumentParser(
    prog='structure_setup',
    description='basic structure setup'
)

parser.add_argument(
    '--naccess',
    action='store',
    dest='naccess_bin',
    default=os.path.dirname(os.path.abspath(__file__)) + '/soft/NACCESS/naccess',
    help='Vdw parameters'
)
parser.add_argument(
    '--vdw',
    action='store',
    dest='vdwprm_file',
    default=os.path.dirname(os.path.abspath(__file__)) + '/data/vdwprm',
    help='Vdw parameters'
)

parser.add_argument(
    '--diel_type',
    action='store',
    dest='diel_type',
    default = MHDIEL,
    help="Type of dielectric: 0: Vacuum, 1: water, 2: M.H"
)

parser.add_argument(
    '--int_cutoff',
    action='store',
    type=float,
    default = INT_CUTOFF,
    help="Max distance to calculate interaction energy"
)

parser.add_argument(
    '--vdw_cutoff',
    action='store',
    type=float,
    default = VDW_CUTOFF,
    help="Max distance to calculate VdW interaction energy"
)

parser.add_argument(
    '--debug',
    action='store_true',
    help="Print At-At interactions"
)
parser.add_argument(
    '--aasolv',
    dest='aasolv_file',
    help='Input Solv_UnfoldedAA', 
    type=open, 
    required=False
)

parser.add_argument('pdb_file',help='Input PDB', type=open)
parser.add_argument('pdbqt_file',help='Input PDBQT', type=open)

args = parser.parse_args()

print("PDB File:", args.pdb_file.name)
print("PDBQT File:", args.pdbqt_file.name)
print("Parameters File:", args.vdwprm_file)

# loading VdW parameters
ff_params = VdwParamset(args.vdwprm_file)

#load DDSolv in Unfolded
aasolv={}
if args.aasolv_file:
    for line in args.aasolv_file:
        aa, solv = line.rstrip().split(' ')
        aasolv[aa] = float(solv)
for r in aasolv:
    aasolv[r] = aasolv['ALA'] - aasolv[r]
aasolv['HIE'] = aasolv['HIS']
aasolv['HID'] = aasolv['HIS']
    
# loading PDB and PDBQT files
st = load_protein(args, ff_params)


# Calculating surfaces
# Srf goes to .xtra['EXP_NACCESS'] field
srf = NACCESS_atomic(st[0], naccess_binary=args.naccess_bin)

#TOTAL Energies (evaluating ALA atoms also)
total_energy = 0.
total_elec = 0.
total_vdw = 0.
total_solv = 0.
# Energies per residue as dictionaries, original residue and ALA
elec_res = {}
vdw_res = {}
solv_res = {}
elec_ala = {}
vdw_ala = {}
solv_ala = {}

for res in st.get_residues():
    elec_res[res] = 0.
    vdw_res[res] = 0.
    solv_res[res] = 0.
    elec_ala[res] = 0.
    vdw_ala[res] = 0.
    solv_ala[res] = 0.
    
for at1 in st.get_atoms():
    res1 = at1.get_parent()

    solv = solvation(at1)
   
    solv_res[res1] += solv
    if at1.id in ALA_ATOMS:
        solv_ala[res1] += solv

    total_solv += solv
    
    for at2 in st.get_atoms():
        #Make sure that atoms are from different residues
        res1 = at1.get_parent()
        res2 = at2.get_parent()
        if res1 == res2: 
            continue
        #and they are not bound
        dist = at1 - at2
        if dist < BOND_CUTOFF:
            continue
        #and each pair es counted only once
        if at1.serial_number >= at2.serial_number:
            continue
            
        #Do not calculate to separated atoms
        if dist > INT_CUTOFF:
            continue
            
        elec = electrostatic_int(at1, at2, dist, args.diel_type)
        
        elec_res[res1] += elec
        elec_res[res2] += elec
        if at1.id in ALA_ATOMS:
            elec_ala[res1] += elec
        if at2.id in ALA_ATOMS:
            elec_ala[res2] += elec

        total_energy += elec
        total_elec += elec
        
        if dist > VDW_CUTOFF:
            continue

        vdw = vdw_int(at1, at2, dist)
        
        vdw_res[res1] += vdw
        vdw_res[res2] += vdw
        
        if at1.id in ALA_ATOMS:
            vdw_ala[res1] += vdw
        if at2.id in ALA_ATOMS:
            vdw_ala[res2] += vdw
        if args.debug:
            print("#DEBUG: {:20s} {:20s} {:8.4} {:8.4}".format(
                atom_id(at1),
                atom_id(at2),
                elec, vdw
                )
            )
            
        total_energy += vdw
        total_vdw += vdw

print("Total Elec: {:10.4f}".format(total_elec))
print("Total VdW : {:10.4f}".format(total_vdw))
print("Total Solv: {:10.4f}".format(total_solv))
print("Total     : {:10.4f}".format(total_energy))
print("{:10s}: {:10s} {:10s} {:10s} {:10s} {:10s} {:10s} {:10s} {:10s} {:10s} {:10s} {:10s} ".format(
    'Residue',"Elec","VdW","Solv","Elec_A","VdW_A","Solv_A","DDElec_A","DDVdW_A","DDSolv_A", "DDSolv_U","DDR-A"))
for res in st.get_residues():
    print('{:10s}: {:10.4f} {:10.4f} {:10.4f} {:10.4f} {:10.4f} {:10.4f} {:10.4f} {:10.4f} {:10.4f} {:10.4f} {:10.4f}'.format(
        residue_id(res), 
        elec_res[res],
        vdw_res[res],
        solv_res[res],
        elec_ala[res],
        vdw_ala[res],
        solv_ala[res],
        elec_ala[res] - elec_res[res],
        vdw_ala[res] - vdw_res[res],
        solv_ala[res] - solv_res[res],
        aasolv[res.get_resname()],
        elec_ala[res] + vdw_ala[res] + solv_ala[res] -\
            elec_res[res] - vdw_res[res] - solv_res[res] - aasolv[res.get_resname()]
    )
    )
        
        