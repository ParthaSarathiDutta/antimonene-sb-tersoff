import os
import subprocess

import numpy as np
from ase.io import read
from math import fabs, sqrt, exp, sin, pi, cos, isnan
import re

idKey2 = re.compile(r'\s+[-+]?\d*\.\d')
magmom_dic = {"Al": 0.6, "Co": 5.0,"W": 0.6}
def shell(command):
    return subprocess.check_output(command,shell=True).strip()


# get dimer distance from poscar(or contcar)
def getDimerDistance(poscar):
    struct = read(poscar, -1, 'vasp')
    pos = struct.get_positions()
    delx = pos[0][0]-pos[1][0]
    dely = pos[0][1]-pos[1][1]
    delz = pos[0][2]-pos[1][2]
    d = sqrt(delx*delx + dely*dely + delz*delz)
    return d

# Total energy from outcar 
def final_energy_from_outcar(filename):
    with open(filename ) as f:
        outcar = f.read()
        energy_re = re.compile( "energy\(sigma->0\) =\s+([-\d\.]+)" )
        energy = float( energy_re.findall( outcar )[-1] )
    return energy



def getPotcar(outcar):
    for line in open(outcar):
        rec = line.strip()
        if rec.startswith('TITEL'):
            potcar = line.split('=')[-1].strip()
    return potcar        

def getSpin(outcar):
    for line in open(outcar):
        rec = line.strip()
        if rec.startswith('ISPIN'):
            spin = int(line.split('=')[-1].strip().split()[-4])
    return spin        
def getEdiff(outcar):
    diff = []
    for line in open(outcar):
        rec = line.strip()
        if rec.startswith('EDIFF'):
            ediff = line.split('=')[-1].strip()
            diff.append(ediff)
    return diff[0].replace('stopping-criterion for ELM', '').strip()        


def checkFile(path):
    if os.path.exists(path) and os.path.getsize(path) > 0: 
        flag = True
    else:
        flag = False   
    return flag

def getNELM(OSZICAR):
    nElm = shell("tail -2  " + OSZICAR + "| head -1 |awk '{print $2}'").decode("utf-8")
    return nElm

def checkSuccess(path):
    search = 'General timing and accounting informations for this job:\n'
    with open(path) as f:
        datafile = f.readlines()
    success = False
    for line in datafile:
        if search in line:
            success = True
            break
    return success

def getTotalEnergy(OUTCAR):
    energy = shell("grep 'free  energy' " +  OUTCAR +  "|awk '{print $5}'|tail -n 1").decode("utf-8")
    return energy


def getPositions(nAtoms, OUTCAR):
    file = open(OUTCAR, 'r')
    lines= file.readlines()
    file.close()
    n=0

def getCella(poscar):
    struct = read(poscar, -1, 'vasp')
    cellVector = struct.cell[:]
    return (cellVector[0][0])

def getPosForces(OUTCAR):
    positions = []
    forces = []
    with open(OUTCAR, "r") as infile:
        lines = infile.readlines()
        n = 0
        pos = -1
        for i, line in enumerate(lines):
            if 'TOTAL-FORCE' in line:
                pos = i
        infile.seek(0)

        for i, line in enumerate(lines):
            if i < pos+2:
                continue
            rawid = idKey2.search(line)
            if rawid is not None:
                x,y,z, fx,fy,fz = tuple([float(x) for x in line.split()])
                positions.append([x,y,z])
                forces.append([fx,fy,fz])
                pass
            else:
                break
    return np.array(positions), np.array(forces)

def getPositions(OUTCAR):
    pos, forces = getPosForces(OUTCAR)
    return pos


def getForces(OUTCAR):
    pos, forces = getPosForces(OUTCAR)
    return forces

def getMaxForceMag(force):
    forceMag = []
    for f in range(len(force)):
        forceMag.append(np.linalg.norm(force[f]))
    return max(forceMag)    


def totalAtoms(poscar):
    with open(poscar, 'r') as poscar:
        poscar_lines = poscar.readlines()
    stoichiometries = poscar_lines[6].split()
    n_atoms = 0
    for stoichiometry in stoichiometries:
        n_atoms += int(stoichiometry)
    return n_atoms    

def getMaxForce(outcar):
    forces = getForces(outcar)
    forceMag = []
    for f in range(len(forces)):
        forceMag.append(np.linalg.norm(forces[f]))
    return max(forceMag)   





magmom_dic = {"Al": 0.6, "Fe":3.0, "Co": 2.0,"W": 0.6}

def getMagmomString(poscar):
    Geometry = read(poscar, -1, 'vasp')
    #Atomic_Symbols = Geometry.get_chemical_symbols()
    #Natoms = Geometry.get_number_of_atoms()
    #Ntypes = len(set(Geometry.get_atomic_numbers()))
    Species={}
    for atmtype in Geometry.get_chemical_symbols():
        if atmtype not in Species:
            Species[atmtype] = 1
        else:
            Species[atmtype] += 1
                                                                                    
    mag = ''
    for k, v in Species.items():
        #print (k, v)
        mag += str(v)+ '*' + str(magmom_dic[k])
        mag +=' '
    print ('MAGMOM=',mag)   
    return (mag)  



def writePotcar(poscar):
    potcarFile = poscar.replace('POSCAR', 'POTCAR')
    struct = read(poscar, -1, 'vasp')
    Species = []
    for x in struct.get_chemical_symbols():
        if x not in Species:
            Species.append(x)

    potcars = []
    potcarSrc = '/global/cfs/cdirs/m1917/smanna/PotcarDirs/'
    
    for sp in Species:
        potcars.append(potcarSrc + sp + '/POTCAR')
    
    with open(potcarFile, 'w') as outfile:
        for potcar in potcars:
            with open(potcar) as infile:
                for line in infile:
                    outfile.write(line)

    print ('Generating potcar with ', Species, 'at ', potcarFile)


def elastic_moduli(path):
    from re import M as multline
    from re import findall
    from numpy import array
    regex = r"\s*TOTAL\s+ELASTIC\s+MODULI\s+\(kBar\)\s*\n"                     \
            r"\s*Direction\s+XX\s*YY\s*ZZ\s*XY\s*YZ\s*ZX\s*\n"                 \
            r"\s*-+\s*\n"                                                      \
            r"\s*XX\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s*\n"      \
            r"\s*YY\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s*\n"      \
            r"\s*ZZ\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s*\n"      \
            r"\s*XY\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s*\n"      \
            r"\s*YZ\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s*\n"      \
            r"\s*ZX\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s*\n"      \
            r"\s*-+\s*\n"
    result = findall(regex,open(path).read(), multline)[0]
    return array(result,dtype='float64').reshape((6,6))

def corrCij(Cij):
    C_ij = np.zeros((6,6))
    C_ij[0,0] = Cij[0,0]*0.1
    C_ij[0,1] = Cij[0,1]*0.1
    C_ij[0,2] = Cij[0,2]*0.1
    C_ij[1,0] = Cij[1,0]*0.1
    C_ij[1,1] = Cij[1,1]*0.1
    C_ij[1,2] = Cij[1,2]*0.1
    C_ij[2,0] = Cij[2,0]*0.1
    C_ij[2,1] = Cij[2,1]*0.1
    C_ij[2,2] = Cij[2,2]*0.1
    C_ij[0,3] = Cij[0,4]*0.1
    C_ij[0,4] = Cij[0,5]*0.1
    C_ij[0,5] = Cij[0,3]*0.1
    C_ij[1,3] = Cij[1,4]*0.1
    C_ij[1,4] = Cij[1,5]*0.1
    C_ij[1,5] = Cij[1,3]*0.1
    C_ij[2,3] = Cij[2,4]*0.1
    C_ij[2,4] = Cij[2,5]*0.1
    C_ij[2,5] = Cij[2,3]*0.1
    C_ij[3,0] = Cij[4,0]*0.1
    C_ij[4,0] = Cij[5,0]*0.1
    C_ij[5,0] = Cij[3,0]*0.1
    C_ij[3,1] = Cij[4,1]*0.1
    C_ij[4,1] = Cij[5,1]*0.1
    C_ij[5,1] = Cij[3,1]*0.1
    C_ij[3,2] = Cij[4,2]*0.1
    C_ij[4,2] = Cij[5,2]*0.1
    C_ij[5,2] = Cij[3,2]*0.1
    C_ij[3,3] = Cij[4,4]*0.1
    C_ij[3,4] = Cij[4,5]*0.1
    C_ij[3,5] = Cij[4,3]*0.1
    C_ij[4,3] = Cij[5,4]*0.1
    C_ij[4,4] = Cij[5,5]*0.1
    C_ij[4,5] = Cij[5,3]*0.1
    C_ij[5,3] = Cij[3,4]*0.1
    C_ij[5,4] = Cij[3,5]*0.1
    C_ij[5,5] = Cij[3,3]*0.1
    
    return C_ij
def final_energy_from_outcar(filename):
    with open( filename ) as f:
        outcar = f.read()
        energy_re = re.compile( "energy\(sigma->0\) =\s+([-\d\.]+)" )
        energy = float( energy_re.findall( outcar )[-1] )
    return energy





if __name__ == "__main__":
    import sys
    filename = sys.argv[1]
#    getPosForces(filename)
    positions, forces = getPosForces(filename)
    print(len(forces))
    for position, force in zip(positions, forces):
        print(position, force)

