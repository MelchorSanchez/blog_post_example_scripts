"""
Small script to test and exemplify PLIP XML output parsing

Melchor Sanchez-Martinez, 2021
"""
import os, sys
from subprocess import call
import shlex
import xml.etree.ElementTree as ET
import pandas as pd

plip=('plipcmd.py -i 1gs5 -xt') #Here 1gs5 is downloaded but with -f can be
#inputed an existing file. With plipcmd.py -h more options can be shown.
call(plip.split()) # The call can be nicer with a better use of subproces

plipoutput='report.xml' # the file name can be replaced by the name of
#your PLIP output. It can be provided as an argyment with sys.argv or with
#argparse.
tree = ET.parse(plipoutput)
root = tree.getroot()

number_detected_ligands=len(root.findall("bindingsite"))
ligands_name=[]
ligtype=[]
longname=[]
for i in range(len(root.findall("./bindingsite/identifiers/"))):
    if root.findall("./bindingsite/identifiers/")[i].tag == 'hetid':
        ligands_name.append(root.findall("./bindingsite/identifiers/")[i].text)
    if root.findall("./bindingsite/identifiers/")[i].tag == 'longname':
        longname.append(root.findall("./bindingsite/identifiers/")[i].text)
    if root.findall("./bindingsite/identifiers/")[i].tag == 'ligtype':
        ligtype.append(root.findall("./bindingsite/identifiers/")[i].text)

desired_ligand='ANP' #this can be provided as an argument with sys.argv or with
#argparse. It can be commented if the information of all the ligands should be
#extracted
interactions={}
for i in range(len(root.findall("./bindingsite/identifiers/longname"))):
    if desired_ligand in root.findall("./bindingsite/identifiers/longname")[i].text: #If all interactions wanted, this line should be commented
        for int in root.findall("./bindingsite/interactions")[i]:
            if int.text == None:
                interactions[int.tag+'-'+root.findall("./bindingsite/identifiers/longname")[i].text]=int.text
            else:
                interactions[int.tag+'-'+root.findall("./bindingsite/identifiers/longname")[i].text]={}
                for j in range(len(root.findall("./bindingsite/interactions/"+int.tag+"/"+int.tag[:-1])) or len(root.findall("./bindingsite/interactions/"+int.tag+"/"+int.tag[:-2]))):
                    if int.tag == 'metal_complexes':
                        if root.findall("./bindingsite/interactions/"+int.tag+"/"+int.tag[:-2]+"/restype")[j].text == desired_ligand:#If all interactions wanted, this line should be commented
                            ident = root.findall("./bindingsite/interactions/"+int.tag+"/metal_complex")[j].attrib['id']
                            metalname=root.findall("./bindingsite/interactions/"+int.tag+"/metal_complex/metal_type")[j].text
                            metalnr=root.findall("./bindingsite/interactions/"+int.tag+"/metal_complex/metal_idx")[j].text
                            metch=root.findall("./bindingsite/interactions/"+int.tag+"/metal_complex/reschain_lig")[j].text
                            interactions[int.tag+'-'+root.findall("./bindingsite/identifiers/longname")[i].text][ident]=metalname+metalnr+'.'+metch
                    elif int.tag == 'water_bridges':
                        if root.findall("./bindingsite/interactions/"+int.tag+"/"+int.tag[:-1]+"/restype_lig")[j].text == desired_ligand: #If all interactions wanted, this line should be commented
                            ident = root.findall("./bindingsite/interactions/"+int.tag+"/"+int.tag[:-1])[j].attrib['id']
                            resname=root.findall("./bindingsite/interactions/"+int.tag+"/"+int.tag[:-1]+"/restype")[j].text
                            resnr=root.findall("./bindingsite/interactions/"+int.tag+"/"+int.tag[:-1]+"/resnr")[j].text
                            reschain=root.findall("./bindingsite/interactions/"+int.tag+"/"+int.tag[:-1]+"/reschain_lig")[j].text
                            watid=root.findall("./bindingsite/interactions/"+int.tag+"/"+int.tag[:-1]+"/water_idx")[j].text
                            interactions[int.tag+'-'+root.findall("./bindingsite/identifiers/longname")[i].text][ident]=resname+resnr+'.'+reschain+'-WAT'+watid
                    else:
                        if root.findall("./bindingsite/interactions/"+int.tag+"/"+int.tag[:-1]+"/restype_lig")[j].text == desired_ligand: #If all interactions wanted, this line should be commented
                            ident = root.findall("./bindingsite/interactions/"+int.tag+"/"+int.tag[:-1])[j].attrib['id']
                            resname=root.findall("./bindingsite/interactions/"+int.tag+"/"+int.tag[:-1]+"/restype")[j].text
                            resnr=root.findall("./bindingsite/interactions/"+int.tag+"/"+int.tag[:-1]+"/resnr")[j].text
                            reschain=root.findall("./bindingsite/interactions/"+int.tag+"/"+int.tag[:-1]+"/reschain_lig")[j].text
                            interactions[int.tag+'-'+root.findall("./bindingsite/identifiers/longname")[i].text][ident]=resname+resnr+'.'+reschain
                    #there are more properties that can be parse as if the
                    #protein is the donor, in the case of HBs, if the bond is
                    #esyablished with the side chain or if there are more than
                    #one desired ligand in different chains
interaction_types={"hydrophobic_interactions":"HI", "hydrogen_bonds":"HB", "halogen_bonds":"HalB"\
,"water_bridges":"WB", "pi_stacks":"PS", "pi_Cation_interactions":"PC", "metal_complexes":"MC", \
"salt_bridges":"SB" }"
with open('interactions_summary.csv', 'w') as fw:
    for key,val in interactions.items():
        for key2 in interaction_types.keys():
            if key2 in key and val != None:
                for ele in val.values():
                    fw.write(interaction_types[key2]+'\t'+ele+os.linesep)
