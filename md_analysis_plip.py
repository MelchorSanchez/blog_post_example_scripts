"""
Small script to test PLIP analysis over an MD

Melchor Sanchez-Martinex,2021
"""
import subprocess
import shlex
import pytraj as pt

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
traj=pt.iterload('/home/melchor/science/papers/paper_laura/test/target_ligand_\
solvent-w-MD.dcd', top='/home/melchor/science/papers/paper_laura/test/\
target_ligand_solvent.prmtop', frame_slice=[(0, 100, 10),])
traj=traj.strip(':WAT,:Na+')


count=0
for frame in traj:
    count=count+1
    pt.write_traj('trajectory_frame'+str(count)+'.pdb', traj=traj, overwrite=True)
    pt.write_traj('trajectory.pdb', traj=traj, overwrite=True)
    plip=('plipcmd -f ' + trajectory_frame'+str(count)+'.pdb -t')
    run_command(plip)
