# -*- coding: utf-8 -*-
"""
Created on Fri Jan 19 10:32:37 2018

@author: Pranav Kodali
"""

import sys, os, time, getpass, errno
from shutil import copyfile
from Bio import pairwise2 as pw2

filename = sys.argv[1]
reference = sys.argv[2]
referenceContents = ' '
pdbCode = []
pdbContents = []
isReference = 0

allXML = ''



with open(filename) as file:
	for line in file.readlines():
		if line.find('>') == -1 and isReference == 0:
			pdbContents[len(pdbContents) - 1] =  pdbContents[len(pdbContents) - 1] + line
		elif line.find('>') == -1 and isReference == 1:
			referenceContents = referenceContents + line
		else:
			if line.find(reference[:4]) > -1:
				isReference = 1
			else:
				isReference = 0
				pdbCode.append(line[1:])
				pdbContents.append(' ')

referenceContents = referenceContents.replace('\n','')
referenceContents = referenceContents.replace('/','-/')

filenames = os.listdir(os.path.dirname(filename))

for filenamed in filenames:
	if filenamed.find('.fasta') > -1 or filenamed.find('.pdb') > -1:
		fileorigin, file_extension = os.path.splitext(filenamed)
		os.rename(os.path.dirname(filename) + '/' + filenamed, os.path.dirname(filename) + '/' + filenamed[:4].lower() + file_extension)

with open('./rosetta_cm.options', 'r') as file:
	Options = file.read()
	newOptions = Options.replace('InsertFastaHere', (reference[:4].lower()) + '.fasta')

with open('./rosetta_cm.xml', 'r') as file:
	XML = file.read()
	XML = XML.replace('3FRAGSFILE', reference[:4].lower() + '_3.frags')
	XML = XML.replace('9FRAGSFILE', reference[:4].lower() + '_9.frags')

with open('./run_array.slurm', 'r') as file:
	RunSlurm = file.read()

with open('./run_relax.slurm', 'r') as file:
	RunRelax = file.read()


with open(os.path.dirname(filename) + '/rosetta_cm.options', 'w') as file:
	file.write(newOptions)
RunAllSlurm = RunSlurm.replace('PATHL', reference[:4] + reference[7:].replace('/',''))
RunAllRelax = RunRelax.replace('PATHL', reference[:4] + reference[7:].replace('/','') + 'relax')
RunAllSlurm = RunAllSlurm.replace('OPTIONSP', os.path.dirname(os.path.abspath(filename)) + '/rosetta_cm.options')
RunAllRelax = RunAllRelax.replace('OPTIONSP', os.path.dirname(os.path.abspath(filename)) + '/relax.options')

if 'all' in reference:
	RunAllSlurm = RunAllSlurm.replace('NUMBEROFTRIALS', '10')
	RunAllSlurm = RunAllSlurm.replace('THREADS', '50')
else:
	RunAllSlurm = RunAllSlurm.replace('NUMBEROFTRIALS', '2')
	RunAllSlurm = RunAllSlurm.replace('THREADS', '10')

if 'all' in reference:
	RunAllRelax = RunAllRelax.replace('THREADS', '100')
else:
	RunAllRelax = RunAllRelax.replace('THREADS', '10')

with open(os.path.dirname(filename) + '/run_array.slurm', 'w') as file:
	file.write(RunAllSlurm)
with open(os.path.dirname(filename) + '/run_relax.slurm', 'w') as file:
	file.write(RunAllRelax)
copyfile('./relax.options', os.path.dirname(filename) + '/relax.options')
copyfile('./sc_parser.bash', os.path.dirname(filename) + '/sc_parser.bash')


for i in range(len(pdbCode) ):
	pdbContents[i] = pdbContents[i].replace('\n','') + '-'
	newFileName = (reference[:4].lower() + "_" + pdbCode[i][:4].lower() + ".grishin")
	f = open(os.path.join(os.path.dirname(filename), newFileName), 'w')
	f.write("## {} {}.pdb\n".format(reference[:4].lower(), pdbCode[i][:4].lower()))
	f.write("#\n")
	f.write("scores_from_program: 0\n")
	f.write("0{}\n".format(referenceContents))
	f.write("0{}\n".format(pdbContents[i]))
	f.close()
	if 'all' not in reference:
		if not os.path.exists(os.path.join(os.path.dirname(filename), (pdbCode[i][:4].lower()) + "_0001.pdb.gz")):
			os.system('( cd ./' + reference + ' && /dors/meilerlab/apps/rosetta/rosetta-3.7/main/source/bin/relax.default.linuxgccrelease -s ' + pdbCode[i][:4].lower() + '.pdb @relax.options -nstruct 1 )')
		while not os.path.exists(os.path.join(os.path.dirname(filename), (pdbCode[i][:4].lower()) + "_0001.pdb.gz")):
			time.sleep(1)
		os.remove(os.path.join(os.path.dirname(filename), (pdbCode[i][:4].lower()) + ".pdb"))
		os.rename(os.path.join(os.path.dirname(filename), (pdbCode[i][:4].lower()) + "_0001.pdb.gz"), os.path.join(os.path.dirname(filename), (pdbCode[i][:4].lower()) + ".pdb.gz"))
		os.system('( cd ./' + reference + ' && /dors/meilerlab/apps/rosetta/rosetta-3.7/main/source/bin/partial_thread.default.linuxgccrelease -in:file:fasta {} -in:file:alignment {} -in:file:template_pdb {} )'.format((reference[:4].lower()) + '.fasta',newFileName,(pdbCode[i][:4].lower()) + ".pdb.gz"))
		while not os.path.exists(os.path.join(os.path.dirname(filename), (pdbCode[i][:4].lower()) + ".pdb.pdb")):
			time.sleep(1)
		copyfile(os.path.join(os.path.dirname(filename), (pdbCode[i][:4].lower()) + ".pdb.pdb"), os.path.join(os.path.dirname(filename)[:-1] + '/all', (pdbCode[i][:4].lower()) + ".pdb.pdb"))
		XMLstring = ' <Template pdb="' + pdbCode[i][:4].lower() + '.pdb.pdb" cst_file="AUTO" weight=   1.000 />'
	else:
		globalalign = pw2.align.globalxx(referenceContents, pdbContents[i])
		score = globalalign[0][2] / len(referenceContents)
		XMLstring = ' <Template pdb="' + pdbCode[i][:4].lower() + '.pdb.pdb" cst_file="AUTO" weight=   ' + str(score) + ' />'
		#XMLstring = ' <Template pdb="' + pdbCode[i][:4].lower() + '.pdb.pdb" cst_file="AUTO" weight=   1.000 />'
	allXML = allXML + XMLstring + '\n             '

with open(os.path.dirname(filename) + '/rosetta_cm.xml', 'w') as file:
	file.write(XML.replace('THREADING', allXML))

