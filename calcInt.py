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
DIST_CUTOFF = 20.
COV_CUTOFF = 2.
ala_atoms = {
    'N', 'H', 'CA', 'HA', 'C', 'O', 'CB', 
    'HB', 'HB1', 'HB2', 'HB3', 'HA1', 'HA2', 'HA3'
}

pairs_1_4 = {
    'N' : {'N'},
    'CA': {'N','H','CD','CA'},
    'HA': {'N'},
    'C' : {'N','H','CD','CA','HA','CB','C'},
    'O' : {'N','H','CD','CA'}
}
pairs_1_4_rev = {
    'N' : {'N','CA','HA','C','O'},
    'H' : {'CA','C','O'},
    'CA': {'CA','C','O'},
    'HA': {'C'},
    'CB': {'C'},
    'C' : {'C'},
    'CD': {'CA','C','O'}
}


def is1_4(at1, at2):
    r1 = at1.get_parent()
    r2 = at2.get_parent()
    if abs(r1.id[1] - r2.id[1]) != 1:
        return False
    if at1.serial_number < at2.serial_number:
        return at1.name in pairs_1_4 and at2.name in pairs_1_4[at1.name]
    else:
        if at1.name not in pairs_1_4_rev:
            return False
        return at1.name in pairs_1_4_rev and at2.name in pairs_1_4_rev[at1.name]

def residue_id(res):
    '''Returns readable residue id'''
    return '{} {}{}'.format(res.get_resname(), res.get_parent().id, res.id[1])

def atom_id(at):
    '''Returns readable atom id'''
    return '{}.{}'.format(residue_id(at.get_parent()), at.id)

parser = argparse.ArgumentParser(
    prog='structure_setup',
    description='basic structure setup'
)

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

# evaluating all res - res interaction energies to get total energies. 

en_res = {}
en_res_ala = {}

for r1 in st.get_residues():
    en_res[r1] = {'elec':{'total':0.},'vdw':{'total':0.}, 'solv':0.}
    en_res_ala[r1] = {'elec':{'total':0.},'vdw':{'total':0.}, 'solv':0.}
    for at1 in r1.get_atoms():
        solv = en.calc_solvation(at1)
        en_res[r1]['solv'] += solv
        if at1.name in ala_atoms:
            en_res_ala[r1]['solv'] += solv
        for r2 in st.get_residues():
            if r2 not in en_res[r1]['elec']:
                en_res[r1]['elec'][r2] = 0.
                en_res_ala[r1]['elec'][r2] = 0.
            if r2 not in en_res[r1]['vdw']:
                en_res[r1]['vdw'][r2] = 0.
                en_res_ala[r1]['vdw'][r2] = 0.
            for at2 in r2.get_atoms():
                if r1 == r2:
                    continue
                d = at1 - at2
                if d > DIST_CUTOFF or d < COV_CUTOFF:
                    continue
                # Avoid 1-4 interactions
                if is1_4(at1,at2):
                    continue
                if args.diel_type == 0:
                    elec = en.elec_int(at1, at2, d, 1.)
                elif args.diel_type == 1:
                    elec = en.elec_int(at1, at2, d, 80.)
                else:
                    elec= en.elec_int(at1, at2, d, en.diel_MS(d))
                en_res[r1]['elec'][r2] += elec
                en_res[r1]['elec']['total'] += elec
                if at1.name in ala_atoms:
                    en_res_ala[r1]['elec'][r2] += elec
                    en_res_ala[r1]['elec']['total'] += elec

                vdw = en.vdw_int(at1, at2, d)
                en_res[r1]['vdw'][r2] += vdw
                en_res[r1]['vdw']['total'] += vdw
                if at1.name in ala_atoms:
                    en_res_ala[r1]['vdw'][r2] += vdw
                    en_res_ala[r1]['vdw']['total'] += vdw

for r in st.get_residues():
    print('#ALA {:10}: {:8.3f} {:8.3f} {:8.3f} | {:8.3f} {:8.3f} {:8.3f} | {:8.3f} {:8.3f} {:8.3f} | {:8.3f} '.format(
        residue_id(r), 
        en_res[r]['solv'], 
        en_res_ala[r]['solv'], 
        en_res_ala[r]['solv'] - en_res[r]['solv'],
        en_res[r]['elec']['total'], 
        en_res_ala[r]['elec']['total'], 
        en_res_ala[r]['elec']['total'] - en_res[r]['elec']['total'],
        en_res[r]['vdw']['total'], 
        en_res_ala[r]['vdw']['total'], 
        en_res_ala[r]['vdw']['total'] - en_res[r]['vdw']['total'],
        en_res_ala[r]['solv'] - en_res[r]['solv'] +\
        en_res_ala[r]['elec']['total'] - en_res[r]['elec']['total'] +\
        en_res_ala[r]['vdw']['total'] - en_res[r]['vdw']['total']
        )
    )

for r2 in st.get_residues():
    line += '#    {:8}'.format(residue_id(r2))
print(line)
    
for r1 in st.get_residues():
    line = '#ELEC {:8}'.format(residue_id(r1))
    for r2 in st.get_residues():
        line += ' {:8.3f}'.format(en_res[r1]['elec'][r2])
    print(line)

for r1 in st.get_residues():
    line = '#VDW  {:8}'.format(residue_id(r1))
    for r2 in st.get_residues():
        line += ' {:8.3f}'.format(en_res[r1]['vdw'][r2])
    print(line)
