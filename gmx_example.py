import sys,os
from gmx_job import gmx_pdb
from gmx_top import *
  
def main():
    pdbfile = gmx_pdb('protein.pdb') #reads a pdb file from the pdb
    top = topology() #generates a blank topology object
    #pdbfile.write_pdb('protein_out.pdb') # writes a pdb object 
    residues = pdbfile.assign_residues() # list of residues, each element is a list of atoms of type pdb_atom  
    pdbfile.create_top(pdb_in='protein.pdb',pdb_out='protein_gmx.pdb',forcefield='charmm27') # uses pdb2gmx to generate a gromacs topology
    #pdbfile.update('protein_gmx.pdb')# important! gromacs (pdb2gmx) orders the atoms according to the provided force field! The pdb data must be updated now! see file with _gmx.pdb ending. Also note that the atom type column is missing now.
    top.read_top('topol.top')    # read the topology file that matches the protein_gmx.pdb file from pdb2gmx
    pdbfile.connect_atoms(top)   # creates a list of bonds
    pdbfile.create_simulation_box(pdb_in='protein_gmx.pdb',pdb_out='box.pdb',distance='1.4') # as the name suggest, this creates a simulation box with 1.4 nm distance to the protein
    pdbfile.solvate(pdb_in='box.pdb',pdb_out='water.pdb',watertype='spc216',topology_file='topol.top') # inserts spc water into the simulation box until it is full
    pdbfile.energy_minimize(pdb_in='water.pdb',pdb_out='em.pdb',topology_file='topol.top',nsteps=500)  # runs 500 energy minimization steps on the solvated protein
    os.system('mv em.part0001.pdb em.pdb') 
    pdbfile.update('em.pdb')#this reads in the final, energy minimized solvated simulation box
main()
