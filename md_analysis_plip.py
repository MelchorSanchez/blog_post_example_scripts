"""
Small script to test PLIP analysis over an MD

Melchor Sanchez-Martinez,2021
"""
import os
from subprocess import call
import shlex
import mdtraj as md
import shutil


#Load the trajectory
traj=md.load('target_ligand_solvent-w-MD.dcd', top='target_ligand_solvent.pdb')
traj=traj.remove_solvent() #To simplify and accelerate analysis
traj_reduced=traj[::int(traj.n_frames/10)] #To simplify and accelerate the
#analysis I only selct 10 frames from the MD
count=0
for frame in traj_reduced:# Saving each frame as a single PDB file
    count+=1
    frame.save('traj_frame'+str(count)+'.pdb')
    #Running PLIP externally. Output (x)xml and (t)txt reports
    plip=('plipcmd.py -f  traj_frame'+str(count)+'.pdb --model 0 -xt')
    call(plip.split())
    #Reaname the PLIP generic report name
    shutil.copy('report.txt', 'traj_frame'+str(count)+'.txt')
    shutil.copy('report.xml', 'traj_frame'+str(count)+'.xml')
