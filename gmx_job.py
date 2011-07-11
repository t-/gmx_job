import sys,os
from gmx_top import *

class gmx_pdb:
    def __init__(self,filename):
        self.atoms   = self.read_pdb(filename)
        self.residues= self.assign_residues()
        self.connect = []
        os.system('if [ ! -d "GMX_logfiles" ]; then mkdir GMX_logfiles; fi;')

    class pdb_atom: #creates blank atom 
        def __init__(self,nr=0,name='',res='',seq=0,x=0.0,y=0.0,z=0.0,occup=0.0,tempfac=0.0,element=''):
            self.nr      = ''
            self.name    = ''
            self.res     = ''
            self.seq     = ''
            self.x       = ''
            self.y       = ''
            self.z       = ''
            self.occup   = '' 
            self.tempfac = ''
            self.element = ''

        def s(self,string,length): # prepares string 'string' with length 'length'
            string = str(string)
            while len(string)<length: string=' ' + string;
            if len(string)>length: string=string[0:length];
            return string
        def cprint(self): #creates string in correct PDB format version 3.2 
            return self.s('ATOM  ',6) + self.s(self.nr,5) + ' ' + self.s(self.name,4) + ' ' + self.s(self.res,3) + '  ' + self.s(self.seq,4) + '    ' + self.s(self.x,8) + self.s(self.y,8) + self.s(self.z,8) + self.s(self.occup,6) + self.s(self.tempfac,6) + self.s('',10) + self.s(self.element,2)  
         
    def read_pdb(self,filename):
        print '\n..Reading PDB file',filename
        pdb_atoms=[]        
        for line in open(filename,'r').readlines():
            if line[:6].find('ATOM')>=0:
                atom = self.pdb_atom()
                try:atom.nr      = (int)    (line[6:11].strip(' ' ));
                except Exception:  print 'Error reading PDB file',filename,'atom nr. missing\n-->',line; sys.exit(1); 
                try:atom.name    =          (line[12:16].strip(' '));
                except Exception:  print 'Error reading PDB file',filename,'atom name. missing\n-->',line; sys.exit(1);
                try:atom.res     =          (line[17:20].strip(' '));
                except Exception:  print 'Error reading PDB file',filename,'atom res. missing\n-->',line; sys.exit(1);
                try:atom.seq     = (int)    (line[22:26].strip(' '));
                except Exception: 0 ;
                try:
                    atom.x       = (float)  (line[30:38].strip(' '));
                    atom.y       = (float)  (line[38:46].strip(' '));
                    atom.z       = (float)  (line[46:54].strip(' '));
                except Exception:  print 'Error reading PDB file',filename,'coordinates missing\n-->',line; sys.exit(1);
                try:atom.occup   = (float)  (line[54:60].strip(' '));
                except Exception: 0;
                try:atom.tempfac = (float)  (line[60:66].strip(' '));
                except Exception: 0;
                try:atom.element =          (line[76:78].strip(' '));
                except Exception: 0;
                pdb_atoms.append(atom)
        print '..Identified',len(pdb_atoms),'atoms from PDB file\n'
        return pdb_atoms
    
    def assign_residues(self):
        #creates a residues[] list which contains residue[] lists
        #each residue list contains pdb_atom() objects
        residues = []
        residue  = []
        first_run= 1
        curr_atom= ''
        last_atom= self.atoms[0].name
        curr_res = ''
        last_res = self.atoms[0].res
        for atom in self.atoms:
            curr_res = atom.res
            curr_atom= atom.name
            if curr_res!=last_res or (curr_atom==last_atom and first_run==0):
                residues.append(residue)
                residue=[]
                last_res=curr_res
                last_atom=curr_atom
            if curr_res==last_res :
                residue.append(atom)
            if first_run==1:
                first_run=0
        residues.append(residue)
        #for i in residues:print i[0].res;
        return residues
    
    def write_pdb(self,filename):
        #since each pdb_atom carries a cprint() method, the PDB file can be written very easily
        f=open(filename,'w')
        for atom in self.atoms:
            print >> f, atom.cprint()
            
    def create_top(self,pdb_in,pdb_out,forcefield='charmm27'):
        #this takes a PDB file (pdb_in) and calls gromacs pdb2gmx to create a gromacs topology (topol.top)
        #this file includes the information similar to a namd .psf file
        os.system('pdb2gmx -ignh -f '+pdb_in+' -o '+pdb_out+' -ff '+forcefield+' -water spc> GMX_logfiles/00create_top.log 2>&1 ' )
        success=0
        for line in open('GMX_logfiles/00create_top.log','r'):
            if line.find('successfully generated a topology')>=0:
                success=1
                try:os.system('rm posre.itp \#*');
                except Exception:0;
                print '..Sucessfully generated a topology file for',pdb_in
                print '..Writing PDB file',pdb_out,'in proper Gromacs format\n'
        if success==0:
            print 'Error creating topology for file',pdb_in,'\n-->READ pdb2gmx_error.log'; sys.exit(1);
            
    def create_simulation_box(self,pdb_in,pdb_out,distance):
        #this adds a line to the PDB file with the proper simulation box data
        print '..Creating a simulation box for',pdb_in,'with',distance,'nm spacing between protein and box border'
        print '...writing simulation box to',pdb_out,'\n'
        os.system('editconf -f '+pdb_in+' -o '+pdb_out+' -d '+str(distance)+'> GMX_logfiles/01create_simulation_box.log 2>&1 ')
        #sample create_simulation_box(protein.pdb,protein_in_box.pdb,1.4) # units for distance protein-boxborder is [nm]
         
    def solvate(self,pdb_in,pdb_out,watertype='spc216',topology_file='topol.top'):
        #this calls gromacs to solvate the simulation box using spc216 water
        #tip4p should work just as well
        print '..Solvating',pdb_in
        print '...writing solvated protein to',pdb_out,'\n'
        os.system('genbox -cs '+watertype+' -cp '+pdb_in+' -o '+pdb_out+' -p '+topology_file +'> GMX_logfiles/02solvate.log 2>&1' )
        
    def energy_minimize(self,pdb_in,pdb_out,topology_file,nsteps=500):
        #make config for energy minimization using simple system calls
        #here, we use steepest descent energy minimizer which minimizes to 32bit machine precision for most systems
        os.system('echo "; VARIOUS PREPROCESSING OPTIONS">em.mdp')
        os.system('echo "cpp                      = /lib/cpp">>em.mdp')
        os.system('echo "define                   = -DFLEXIBLE">>em.mdp')
        os.system('echo "; RUN CONTROL PARAMETERS ">>em.mdp')
        os.system('echo "integrator               = steep">>em.mdp')
        os.system('echo "; start time and timestep in ps " >>em.mdp')
        os.system('echo "tinit                    = 0">>em.mdp')
        os.system('echo "dt                       = 0.001">>em.mdp')
        os.system('echo "nsteps                   = '+str(nsteps)+'">>em.mdp')
        os.system('echo "; ENERGY MINIMIZATION OPTIONS " >>em.mdp')
        os.system('echo "emtol                    = 0.00001">>em.mdp')
        os.system('echo "emstep                   = 0.1">>em.mdp')
        os.system('echo "nstcgsteep               = 1000">>em.mdp')
        
        #grompp is used to generate a .tpr file which combines the PDB,topology informations and adds force field parameters as well as simulation parameters, e.g. integrator, duration, etc.
        print '..Preparing',pdb_in,'for energy minimization'
        os.system('grompp -f em.mdp -c '+pdb_in+' -p '+topology_file+' -o em.tpr > GMX_logfiles/03preparing_energy_minimization.log 2>&1 ')
        print '..Running',nsteps,'energy minimization steps for',pdb_in,'and writing output to',pdb_out
        print '..This may take a few minutes for large systems'
        #mdrun is the gromacs integrator. mdrun can be run in parallel on dozens of processors, here it just minimizes the forces on the atoms of the system
        os.system('mdrun -v -s em.tpr -c '+pdb_out+'> GMX_logfiles/04!!IMPORTANT!!energy_minimization_logfile.log 2>&1')           
        print '\n\n_____DO NOT USE THIS SIMULATION BOX WITHOUT READING *ALL* THE LOGFILES IN GMX_logfiles_____\n\n'
        
    def update(self,pdb_in):
        # this method can update an existing gmx_pdb molecule object 
        print '\n ..Updating the PDB object using',pdb_in
        self.atoms=self.read_pdb(pdb_in)
        self.residues= self.assign_residues()
        
    def connect_atoms(self,top):
        #this method creates a list of bonds (don't really need this anymore) as we can create a full psf object
        connect_list=[]
        bonds = []
        for bond in top.bonds_list:
            bonds.append([bond.ai,bond.aj])

        for i in range(1,len(self.atoms)+1):
            connect = [i]
            for atom in bonds:
                if atom[0]==i:
                    connect.append(atom[1])
                if atom[1]==i:
                    connect.append(atom[0])
            connect_list.append(connect) 
