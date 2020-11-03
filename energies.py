# Basic functions to evaluate energies
# Assume Bio.PDB.Atom objects with charge and vdw 
#

import math



def elec_int(at1, at2, r, diel):
    '''Electrostatic interaction energy between two atoms at r distance'''
    return 332.16 * at1.xtra['charge'] * at2.xtra['charge'] / diel / r

def vdw_int(at1, at2, r):
    '''Vdw interaction energy between two atoms'''
    eps12 = math.sqrt(at1.xtra['vdw'].eps * at2.xtra['vdw'].eps)
    sig12_2 = at1.xtra['vdw'].sig * at2.xtra['vdw'].sig
    return 4 * eps12 * (sig12_2**6/r**12 - sig12_2**3/r**6)

def diel_fix(var):
    return var

def diel_MS(dist):
    '''Mehler-Solmajer dielectric'''
    return 86.9525 / (1 - 7.7839 * math.exp(-0.3153 * dist)) - 8.5525

def calc_solvation(at):
    '''Solvation energy based on ASA'''
    if 'EXP_NACCESS' not in at.xtra:
        return 0
    else:
        return float(at.xtra['EXP_NACCESS'])* at.xtra['vdw'].fsrf
