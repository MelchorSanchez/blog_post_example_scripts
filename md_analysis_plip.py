"""
Small script to test PLIP analysis over an MD

Melchor Sanchez-Martinex,2021
"""
import os
import subprocess
import shlex
import pytraj as pt
import shutil

def createPath(path):
    if not os.path.isdir(path):
        os.mkdir(path)

def run_command(command):
    createPath('./logs')
    process = subprocess.Popen(shlex.split(command), stdout=subprocess.PIPE)
    while True:
        output = process.stdout.readline()
        print (output)
        if output == '' and process.poll() is not None:
            break
        if output:
            print (output.strip())
    rc = process.poll()
    return rc

#Load the first 100 frames of the trajectory from 10 to 10
traj=pt.iterload('target_ligand_solvent-w-MD.dcd', top='target_ligand_solvent\
.prmtop', frame_slice=[(0, 100, 10),])
traj=traj.strip(':WAT,:Na+')

count=0
for frame in traj:
    count=count+1
    #Writing each frame as a new PDB file
    pt.write_traj('trajectory_frame'+str(count)+'.pdb', traj=traj, overwrite=True)
    #Writing each frame to the same PDB file
    pt.write_traj('trajectory.pdb', traj=traj, overwrite=True)
    #Running PLIP externally. Output (x)xml and (t)txt reports
    plip=('plipcmd -f  trajectory_frame'+str(count)+'.pdb -xt')
    run_command(plip)
    #Reaname the PLIP generic report name
    shutil.copy('report.txt', 'trajectory_frame'+str(count)+'.txt',)
    shutil.copy('report.xml', 'trajectory_frame'+str(count)+'.xml',)
