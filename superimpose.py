#!C:\Python27
# sslsa.py
# Structural Superimposition of Local Sequence Alignment
# A program which finds out whether a local sequence
# alignment of two protein sequences also implies structural
# similarity of the aligned parts
import os
import sys
import glob
from Bio.PDB import *
from Bio import pairwise2
from Bio.SubsMat import MatrixInfo as matlist
# Obtain structures directory
str_dir = os.path.dirname(os.path.realpath(__file__)) + "\\structures"
# Create it if it doesn't exist
if not os.path.isdir(str_dir):
	os.makedirs(str_dir)
# Get PDB IDs from the user
if len(sys.argv) > 2:
	pdb_ids = [sys.argv[1], sys.argv[2]]
else:
	sys.exit("Two separate valid PDB IDs must be given")
# Initiate PDB list object
pdb_list = PDBList(server="http://www.rcsb.org/pdb/files")
# Retrieve PDB files from the server to structures directory
pdb_list.retrieve_pdb_file(pdb_ids[0], obsolete=False, pdir=str_dir)
pdb_list.retrieve_pdb_file(pdb_ids[1], obsolete=False, pdir=str_dir)
# Generate PDB file paths
pdb_paths = ["".join(glob.glob(str_dir + "\\*" + pdb_ids[0] + ".ent")), "".join(glob.glob(str_dir + "\\*" + pdb_ids[1] + ".ent"))]
# Initiate PDB parser object
pdb_parser = PDBParser(QUIET=True)
# Generate PDB structures using PDB parser
pdb_strs = [pdb_parser.get_structure(pdb_ids[0], pdb_paths[0]), pdb_parser.get_structure(pdb_ids[1], pdb_paths[1])]
# Initiate an empty list for storing PDB sequences
pdb_seqs = ["", ""]
# Initiate CA polypeptide builder used to get sequences of each protein
ppb = CaPPBuilder()
for i in range(len(pdb_seqs)):
	for pp in ppb.build_peptides(pdb_strs[i]): 
		pdb_seqs[i] += str(pp.get_sequence())
# Get BLOSUM62 matrix
matrix = matlist.blosum62
# Set gap penalties or get from the user
if len(sys.argv) >= 4:
	gap_open = int(sys.argv[3])
else:
	gap_open = -10
if len(sys.argv) == 5:
	gap_extend = int(sys.argv[4])
else:
	gap_extend = -5
# Do the pairwise alignment and get alignments
alns = pairwise2.align.localds(pdb_seqs[0], pdb_seqs[1], matrix, gap_open, gap_extend)
# Obtain the best alignment
best_aln = alns[0]
# Decompose best alignment into its components
aln_first, aln_second, score, begin, end = best_aln
# Print the alignment and alignment length
print aln_first[begin:end] + "\n" + aln_second[begin:end]
print "Alignment length: " + str(end - begin)
# Initiate an empty list to store atom objects
pdb_atms = [[], []]
for i in range(len(pdb_atms)):
	# Get only the first model and use it
	model = pdb_strs[i][0]
	for chain in model:
		for residue in chain:
			# Only if the residue has CA atom
			if "CA" in residue:
				# Append the atom object
				pdb_atms[i].append(residue["CA"])
# Initiate another empty string for mapping the atom objects
pdb_atms_mapped = [[], []]
# i is the index for the two alignments, j is for the first
# atom object list and k is for the other atom object list
i, j, k = 0, 0, 0
while i < len(aln_first[:end]):
	# Check if there is no gap in either part of
	# the alignment because there will be no atom
	# for ones with -
	if aln_first[i] != "-" and aln_second[i] != "-":
		# Check if it's the beginning of the alignment
		# here's where we need to start mapping
		if i >= begin:
			# Append the atom objects accordingly
			pdb_atms_mapped[0].append(pdb_atms[0][j])
			pdb_atms_mapped[1].append(pdb_atms[1][k])
	# Move j to the next amino acid if it wasn't a gap
	# that is we put its atom object in the previous 
	# step. If it's a gap, stay at the same atom object
	if aln_first[i] != "-":
		j += 1
	# Move k to the next amino acid if it wasn't a gap
	# that is we put its atom object in the previous 
	# step. If it's a gap, stay at the same atom object
	if aln_second[i] != "-":
		k += 1
	# Move i to the next amino acid in the alignment
	# because we process it no matter what
	i += 1
# Initiate the superimposer
superimposer = Superimposer()
# Set (translate/rotate) atoms minimizing RMSD
superimposer.set_atoms(pdb_atms_mapped[0], pdb_atms_mapped[1])
# Print RMSD
print "RMSD: " + str(superimposer.rms)