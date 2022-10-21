"""
Small example to calculate RESP charges from a SDF file coming from docking
and addition to antechamber.

Inspired from Psikit resp charges calculation example and iwatobipen blog post
https://iwatobipen.wordpress.com/2019/03/11/calculate-atomic-charges-with-psikit-rdkit-psi4/
"""
import pandas as pd
from collections import defaultdict
from psikit import Psikit
pk=Psikit()
from openbabel import pybel
from subprocess import call
import parmed as pmd

def addRESP (inpfile,outfile):
    with open (inpfile, 'r') as fr:
        with open(outfile, 'w')as fw:
            lines= fr.readlines()
            ini=0
            resp_count=0
            init_line=0
            for line in lines:
                ini=ini+1
                if 'UNL' not in line:
                    fw.write(line)
                if '@<TRIPOS>ATOM' in line:
                    init_line=ini
                if init_line > 0 and 'UNL' in line:
                    if '*' not in line:
                        line2=line.split()
                        line2[-1]=molecule_df['RESP'][resp_count]
                        print(line2)
                        for i in range(len(line2)):
                            fw.write(str(line2[i])+' ')
                        fw.write('\n')
                        resp_count=resp_count+1

###Starting from smiles###
metformin='CN(C)C(=N)NC(=N)N'
pk.read_from_smiles(metformin)
pk.optimize('HF/6-31g(d)') #The level of theory can be changed to any methods
#supported by Psi4 like 'B3LIP/def2-QZVPP'  or 'MP2/aug-pVTZ'. Also here you
#control the multiplicity and the maximum number of iterations
respCharges=pk.calc_resp_charges()
res=defaultdict(list)
for i, atom in enumerate(pk.mol.GetAtoms()):
    res['SYMBOL'].append(atom.GetSymbol())
    res['RESP'].append(respCharges[i])
molecule_df = pd.DataFrame(res)

###Starting from a SDF file###
mol2conv=pybel.readfile("sdf", "metformin.sdf")
for mol in mol2conv:
    mol.write("mol", "metformin.mol", overwrite=True)
    mol.write("mol2", "metformin.mol2", overwrite=True)
#SDF or any other format should be convert to MOL. Psikit can read smiles or
#MOL files
pk.read_from_molfile('metformin.mol')
pk.optimize('HF/6-31g(d)')
respCharges=pk.calc_resp_charges()
res=defaultdict(list)
for i, atom in enumerate(pk.mol.GetAtoms()):
    res['SYMBOL'].append(atom.GetSymbol())
    res['RESP'].append(respCharges[i])
molecule_df = pd.DataFrame(res)

####Adding RESP charges###
addRESP('metformin.mol2','metformin-resp.mol2')
"""
To perform MD simulations with gromacs using FFs from the Amber family
an easy option to parameterize organic molecules, not included as standard
residues in the FFs is to do it taking advantage from AmberTools. The
combination of antechamber, parmchk and leap make the work.

As antechamber doesn't recognize mol as input format we use mol2
"""

antechamber="/home/melchor/Software/amber20/bin/antechamber -i metformin-resp.mol2 \
-fi mol2 -o metformin-ach.mol2 -fo mol2 -at amber -pf yes"
call(antechamber.split())

parmchk="/home/melchor/Software/amber20/bin/parmchk2 -i metformin-ach.mol2 -f \
mol2 -o metformin.frcmod"
call(parmchk.split())

#The first lines can be deleted. As there is none protein in the system and
#the smallmolecule is parameteterized from scratch fb15 and amino94 that are used
#to model proteins are not needed. Same for gaff that is used to model known
#small organic molecules and also same for the ions library, as there are not
#ions in the system. However for modelling a protein-ligand complex are need
#so I leave them in the leap template.
LEAP_TEMPLATE="""
source leaprc.protein.fb15
fb15 = loadamberparams frcmod.fb15
source leaprc.gaff
ions1lmtip3p = loadamberparams frcmod.ions1lm_1264_tip3p
ions234lmtip3p= loadamberparams frcmod.ions234lm_1264_tip3p
loadoff phos_amino94.lib

molecule = loadmol2 {}
loadamberparams {}
saveamberparm molecule {} {}
quit
"""

leap = LEAP_TEMPLATE.format("metformin-ach.mol2", "metformin.frcmod", "metformin.prmtop", "metformin.inpcrd")
with open("leap.in", "w") as fd:
    fd.write(leap)
tleap="/home/melchor/Software/amber20/bin/tleap -f leap.in"
call(tleap.split())

##Convert topology to gromacs##
"""
Prmtop and inpcrd can be used to run Amber or NAMD simulations, but these tools
are licenced. NAMD is free for academia. AMBER is not free but the fee
is lower for academia. GROMACS is open source , that's why is interesting to use
it.

Here PARMED is used to perform the topology conversion but other tools like
ACPYPE can also be used. PARMED works as a python library, ACPYPE
has to be called externally.
"""
parm = pmd.load_file("metformin.prmtop", "metformin.inpcrd")
parm.save("metformin.top", format='gromacs')
parm.save("metformin.gro")
