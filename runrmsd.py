#!/bin/env python

import glob, os, sys, fnmatch, getpass, re, subprocess, time, shutil, gzip, pickle, pymol, math

from shutil import copyfile
'''import plotly as py
import plotly.graph_objs as go
from plotly import tools'''
import matplotlib
matplotlib.use("agg")
import matplotlib.pyplot as plt
from pandas import DataFrame
import seaborn as sns
sns.set()
sns.set_context('poster')
sns.set_style("whitegrid")
flatui = ["#9b59b6", "#3498db", "#95a5a6", "#e74c3c", "#34495e", "#2ecc71"]
sns.set_palette(flatui)

from Bio.PDB import *
from Bio import pairwise2
from Bio.PDB.Polypeptide import three_to_one
from Bio.SubsMat import MatrixInfo as matlist

import requests

def read_fasta_file(file_name):
	""" return array of sequences from FASTA formatted file

	file_name -- relative or absolute path to FASTA file
	"""
	seqArray=[]
	seqArrayNames=[]
	seq=""
	for l in open(file_name):
		if l.startswith('>'):
			if not "" == seq:
				seqArray.append(seq)
				seq=""
			seqArrayNames.append(l.rstrip())
		else:
			seq=seq+l.rstrip()
	seqArray.append(seq)
	return seqArray

def write_fasta_file(file_name, data, prefix):
	""" writes fasta file

	file_name -- name of file, will also be header / sequence name
	data -- the single sequence to be written
	prefix -- path prepended to the filename, explicitly have a directory separator if needed
	"""
	with open(prefix+file_name+'.fasta', 'w') as f: f.write('>%s\n' % file_name);  f.write(data); f.write('\n')

def IdentifyCDRs(light_chain, heavy_chain):
	""" Determination of loops in both light and heavy chains

	light_chain -- string representation of light chain
	heavy_chain -- string representation of heavy chain
	"""

	L1_pattern=r'C[A-Z]{1,17}(WYL|WLQ|WFQ|WYQ|WYH|WVQ|WVR|WWQ|WVK|WYR|WLL|WFL|WVF|WIQ|WYR|WNQ|WHL|WHQ|WYM|WYY)'
	L3_pattern=r'C[A-Z]{1,15}(L|F|V|S)G[A-Z](G|Y)'
	H1_pattern=r'C[A-Z]{1,16}(W)(I|V|F|Y|A|M|L|N|G)(R|K|Q|V|N|C|G)(Q|K|H|E|L|R)' # Jeff's mod for ATHM set
   #H1_pattern=r'C[A-Z]{1,16}(W)(I|V|F|Y|A|M|L|N|G)(R|K|Q|V|N|C)(Q|K|H|E|L|R)'
	H3_pattern=r'C[A-Z]{1,33}(W)(G|A|C)[A-Z]{1,2}(Q|S|G|R)'

	''' Identift CDR region and return them as dict with keys: 'FR_H1', 'FR_H2', 'FR_H3', 'FR_H4', 'FR_L1', 'FR_L2', 'FR_L3', 'FR_L4', 'H1', 'H2', 'H3', 'L1', 'L2', 'L3'
	'''
	light_first = light_chain[:65] if len(light_chain) > 130 else light_chain[:60]
	heavy_first = heavy_chain[:70] if len(heavy_chain) > 140 else heavy_chain[:60]

	# LL II GG HH TT

	## L1

	len_FR_L1=0
	res = re.search(L1_pattern,light_first)
	FR_L1 = L1 = False
	if res:
		L1 = res.group()[1:-3]
		L1_start = light_chain.index(L1)
		L1_end = L1_start + len(L1) - 1
		print(("L1 detected: %s (%d residues at positions %d to %d)" % (L1, len(L1), L1_start, L1_end)))
		FR_L1 = light_chain[:L1_start]
		if len(FR_L1) >  24:
			len_FR_L1 = len(FR_L1) - 24
			FR_L1 = light_chain[len_FR_L1:L1_start]
	else:
		print("L1 detected: False")

	light_second = light_chain[L1_end+16+7:L1_end+16+7+80] if len(light_chain) > 130 else light_chain[L1_end+16+7:]

	## L3

	L3 = False
	res = re.search(L3_pattern,light_second)
	if res:
		L3 = res.group()[1:-4]
		L3_start = light_chain.index(L3)
		L3_end = L3_start + len(L3) - 1
		print(("L3 detected: %s ( %d residues at positions %d to %d)" % (L3, len(L3), L3_start, L3_end)))
	else:
		print("L3 detected: False")

	if L1 and L3:
		#L1_start = light_chain.index(L1)
		#L1_end = L1_start + len(L1) - 1

		L2_start = L1_end + 16
		L2_end = L2_start + 7 - 1

		L2 = light_chain[L2_start:L2_start+7]  # L2 is identified here. Current implementation can deal with only 7-resiue L2
		print(("L2 detected: %s (%d residues at positions %d to %d)" % (L2, len(L2), L2_start, L2_end)))

		#FR_L1 = light_chain[:L1_start]
		FR_L2 = light_chain[ L1_end + 1  :  L1_end + 1 + 15					]
		FR_L3 = light_chain[ L2_end + 1  :  L2_end + 1 + L3_start - L2_end - 1 ]
		FR_L4 = light_chain[ L3_end + 1  :  L3_end + 1 + 12					]

		print(("FR_L1: ", FR_L1))
		print(("FR_L2: ", FR_L2))
		print(("FR_L3: ", FR_L3))
		print(("FR_L4: ", FR_L4))
		print(("L segments: ",FR_L1,L1,FR_L2,L2,FR_L3,L3,FR_L4))

		# Light chain sub-type classification. This is useful in the future. But currently this is not used.
		# ... skipped, see Google doc for details
		# FR classification by AHo. This might be useful in the future. But currently this is not used.
		# ... skipped, see Google doc for details


	# HH EE AA VV YY

	## H1
	res = re.search(H1_pattern, heavy_first)
	H1 = False
	len_FR_H1 = 0
	if res:
		H1 = res.group()[4:-4]
		H1_start = heavy_chain.index(H1)
		H1_end = H1_start + len(H1) - 1
		print(("H1 detected: %s (%d residues at positions %d to %d)" % (H1, len(H1), H1_start, H1_end)))
		FR_H1 = heavy_chain[:H1_start]
		if len(FR_H1) >  25:
			len_FR_H1 = len(FR_H1) - 25
			FR_H1 = heavy_chain[len_FR_H1:H1_start]
	else:
		print("H1 detected: False")


	heavy_second = heavy_chain[H1_end+33+15:H1_end+33+15+95+len_FR_H1] if len(heavy_chain) > 140 else heavy_chain[H1_end+33+15:]

	## H3
	H3 = False #H3_and_stem=False
	res = re.search(H3_pattern,heavy_second)
	if res:
		H3 = res.group()[3:-4] #H3_and_stem = res.group()[0:-4]
		H3_start = heavy_chain.index(H3)
		H3_end = H3_start + len(H3) - 1
		print(("H3 detected: %s (%d residues at positions %d to %d)" % (H3, len(H3), H3_start, H3_end)))
	else:
		print("H3 detected: False")


	if H1 and H3:
		#H1_start = heavy_chain.index(H1)
		#H1_end = H1_start + len(H1) - 1
		H2_start = H1_end + 15
		H2_end = H3_start - 33
		H2 = heavy_chain[H2_start:H2_start + H2_end-H2_start+1]
		print(("H2 detected: %s (%d residues at positions %d to %d)" % (H2, len(H2), H2_start, H2_end)))

		#FR_H1 = heavy_chain[:H1_start]

		#if len(FR_H1) >  26:
		#	FR_H1 = light_chain[20:H1_start]

		FR_H2 = heavy_chain[H1_end + 1: H1_end + 1 + H2_start - H1_end - 1]
		FR_H3 = heavy_chain[H2_end + 1: H2_end + 1 + H3_start - H2_end - 1]
		FR_H4 = heavy_chain[H3_end + 1: H3_end + 1 + 12]

		print(("FR_H1: ", FR_H1))
		print(("FR_H2: ", FR_H2))
		print(("FR_H3: ", FR_H3))
		print(("FR_H4: ", FR_H4))
		print(("H segments: ",FR_H1,H1,FR_H2,H2,FR_H3,H3,FR_H4))

	if not (L1 and L3 and H1 and H3):
		if not L1: print(('ERROR: CDR L1 cannot be recognized !!!  L1 pattern: %s' % L1_pattern)) # C[A-Z]{1,17}(WYL|WLQ|WFQ|WYQ|WYH|WVQ|WVR|WWQ|WVK|WYR|WLL|WFL|WVF|WIQ|WYR|WNQ|WHL|WHQ|WYM|WYY)'
		if not L3: print(('ERROR: CDR L3 cannot be recognized !!!  L3 pattern: %s' % L3_pattern)) # C[A-Z]{1,15}(L|F|V|S)G[A-Z](G|Y)
		if not H1: print(('ERROR: CDR H1 cannot be recognized !!!  H1 pattern: %s' % H1_pattern)) # C[A-Z]{1,16}(W)(I|V|F|Y|A|M|L|N|G)(R|K|Q|V|N|C)(Q|K|H|E|L|R)
		if not H3: print(('ERROR: CDR H3 cannot be recognized !!!  H3 pattern: %s' % H3_pattern)) # C[A-Z]{1,33}(W)(G|A|C)[A-Z](Q|S|G|R)
		sys.exit(1)

	res = dict(L1=L1, L2=L2, L3=L3, H1=H1, H2=H2, H3=H3,  FR_L1=FR_L1, FR_L2=FR_L2, FR_L3=FR_L3, FR_L4=FR_L4,  FR_H1=FR_H1, FR_H2=FR_H2, FR_H3=FR_H3, FR_H4=FR_H4)
	#if Options.verbose: print 'L1: %(L1)s\nL2: %(L2)s\nL3: %(L3)s\nH1: %(H1)s\nH2: %(H2)s\nH3: %(H3)s' % res

	return res


def int_(s):
	""" Removes all capital letters from its argument, returns int value of remainder
	"""
	return int( re.sub('[A-Z]', '', s) )  #v = int( re.sub('[A-Z]', '', new_number_FR_L1) )  # $new_number_FR_L1[$i] =~ s/[A-Z]//	 #new_number_FR_L1[i] = string.translate(new_number_FR_L1[i], None, string.ascii_letters)

def Extract_FR_CDR_Sequences(L1='', L2='', L3='', H1='', H2='', H3='', FR_L1='', FR_L2='', FR_L3='', FR_L4='', FR_H1='', FR_H2='', FR_H3='', FR_H4=''):
	""" Find Cothia numbering for loops and conserved regions
	"""

	# LIGHT CHAIN
	# FR_L1	How can we handle missing residue in C/N-terminals?
	if re.search( r'[A-Z][QE][A-Z]{9}[A-Z][A-Z]{4}[LVIMF][A-Z]C', FR_L1): # Change G to [A-Z] (3G04)
		if   len(FR_L1) == 19: new_number_FR_L1="5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23"
		elif len(FR_L1) == 20: new_number_FR_L1="4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23"
		elif len(FR_L1) == 21: new_number_FR_L1="3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23"
		elif len(FR_L1) == 22: new_number_FR_L1="2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23"
		elif len(FR_L1) == 23: new_number_FR_L1="1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23"
		elif len(FR_L1) == 24:
			#new_number_FR_L1="0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23"
			FR_L1 = FR_L1[1:]  # Remove 0th residue 10/24/2012
			new_number_FR_L1="1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23"
		else: print(("ERROR: FR_L1 matches [A-Z][QE][A-Z]{9}[A-Z][A-Z]{4}[LVIMF][A-Z]C but length",len(FR_L1),"is not between 19 and 24"))

	elif re.search( r'[A-Z][QE][A-Z]{8}[A-Z][A-Z]{4}[LVIMF][A-Z]C', FR_L1):  # Change G to [A-Z] (3G04)
		if   len(FR_L1) == 19: new_number_FR_L1="4,5,6,7,8,10,11,12,13,14,15,16,17,18,19,20,21,22,23"
		elif len(FR_L1) == 20: new_number_FR_L1="3,4,5,6,7,8,10,11,12,13,14,15,16,17,18,19,20,21,22,23"
		elif len(FR_L1) == 21: new_number_FR_L1="2,3,4,5,6,7,8,10,11,12,13,14,15,16,17,18,19,20,21,22,23"
		elif len(FR_L1) == 22: new_number_FR_L1="1,2,3,4,5,6,7,8,10,11,12,13,14,15,16,17,18,19,20,21,22,23"
		elif len(FR_L1) == 23:
			FR_L1 = FR_L1[1:]  # Remove 0th residue 10/24/2012
			new_number_FR_L1="1,2,3,4,5,6,7,8,10,11,12,13,14,15,16,17,18,19,20,21,22,23"
		elif len(FR_L1) == 24:
			FR_L1 = FR_L1[2:]  # Remove -1st and 0th residue 10/24/2012
			new_number_FR_L1="1,2,3,4,5,6,7,8,10,11,12,13,14,15,16,17,18,19,20,21,22,23"
		else: print(("ERROR: FR_L1 matches [A-Z][QE][A-Z]{8}[A-Z][A-Z]{4}[LVIMF][A-Z]C but length",len(FR_L1),"is not between 19 and 24"))
	else:
		print('ERROR: Current code could not assign Chothia numbering of FR_L1 in the query sequence!!! Exiting...')

	# L1
	if   len(L1) ==  8: new_number_L1="24,25,26,27,28,29,30,34"
	elif len(L1) ==  9: new_number_L1="24,25,26,27,28,29,30,33,34"
	elif len(L1) == 10: new_number_L1="24,25,26,27,28,29,30,32,33,34"
	elif len(L1) == 11: new_number_L1="24,25,26,27,28,29,30,31,32,33,34"
	elif len(L1) == 12: new_number_L1="24,25,26,27,28,29,30,30A,31,32,33,34"
	elif len(L1) == 13: new_number_L1="24,25,26,27,28,29,30,30A,30B,31,32,33,34"
	elif len(L1) == 14: new_number_L1="24,25,26,27,28,29,30,30A,30B,30C,31,32,33,34"
	elif len(L1) == 15: new_number_L1="24,25,26,27,28,29,30,30A,30B,30C,30D,31,32,33,34"
	elif len(L1) == 16: new_number_L1="24,25,26,27,28,29,30,30A,30B,30C,30D,30E,31,32,33,34"
	elif len(L1) == 17: new_number_L1="24,25,26,27,28,29,30,30A,30B,30C,30D,30E,30F,31,32,33,34"
	else: print(("ERROR: L1 length",len(L1),"is not between 8 and 17"))

	# FR_L2
	if len(FR_L2) == 15: new_number_FR_L2="35,36,37,38,39,40,41,42,43,44,45,46,47,48,49"
	else: print(("ERROR: FR_L2 length",len(FR_L2),"is not 15"))

	# L2
	if   len(L2) ==  7: new_number_L2="50,51,52,53,54,55,56"
	elif len(L2) ==  8: new_number_L2="50,51,52,53,54,54A,55,56"
	elif len(L2) ==  9: new_number_L2="50,51,52,53,54,54A,54B,55,56"
	elif len(L2) == 10: new_number_L2="50,51,52,53,54,54A,54B,54C,55,56"
	elif len(L2) == 11: new_number_L2="50,51,52,53,54,54A,54B,54C,54D,55,56"
	else: print(("ERROR: L2 length",len(L2),"is not between 7 and 11"))

	# FR_L3
	if   len(FR_L3) == 32: new_number_FR_L3="57,58,59,60,61,62,63,64,65,66,67,68,69,70,71,72,73,74,75,76,77,78,79,80,81,82,83,84,85,86,87,88"
	elif len(FR_L3) == 33: new_number_FR_L3="57,58,59,60,61,62,63,64,65,66,66A,67,68,69,70,71,72,73,74,75,76,77,78,79,80,81,82,83,84,85,86,87,88"
	elif len(FR_L3) == 34: new_number_FR_L3="57,58,59,60,61,62,63,64,65,66,66A,66B,67,68,69,70,71,72,73,74,75,76,77,78,79,80,81,82,83,84,85,86,87,88"
	else: print(("ERROR: FR_L3 length",len(FR_L3),"is not between 32 and 34"))

	# L3
	if   len(L3) ==  5: new_number_L3="89,90,91,92,97"
	elif len(L3) ==  6: new_number_L3="89,90,91,92,93,97"
	elif len(L3) ==  7: new_number_L3="89,90,91,92,93,94,97"
	elif len(L3) ==  8: new_number_L3="89,90,91,92,93,94,95,97"
	elif len(L3) ==  9: new_number_L3="89,90,91,92,93,94,95,96,97"
	elif len(L3) == 10: new_number_L3="89,90,91,92,93,94,95,95A,96,97"
	elif len(L3) == 11: new_number_L3="89,90,91,92,93,94,95,95A,95B,96,97"
	elif len(L3) == 12: new_number_L3="89,90,91,92,93,94,95,95A,95B,95C,96,97"
	elif len(L3) == 13: new_number_L3="89,90,91,92,93,94,95,95A,95B,95C,95D,96,97"
	elif len(L3) == 14: new_number_L3="89,90,91,92,93,94,95,95A,95B,95C,95D,95E,96,97"
	elif len(L3) == 15: new_number_L3="89,90,91,92,93,94,95,95A,95B,95C,95D,95E,95F,96,97"
	else: print(("ERROR: L3 length",len(L3),"is not between 5 and 15"))

	# FR_L4
	new_number_FR_L4="98,99,100,101,102,103,104,105,106,107,108,109"

	# HEAVY CHAIN
	# FR_H1	How can we handle missing residue in C/N-terminals?
	if   len(FR_H1) == 16: new_number_FR_H1="10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25"
	elif len(FR_H1) == 17: new_number_FR_H1="9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25"
	elif len(FR_H1) == 18: new_number_FR_H1="8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25"
	elif len(FR_H1) == 19: new_number_FR_H1="7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25"
	elif len(FR_H1) == 20: new_number_FR_H1="6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25"
	elif len(FR_H1) == 21: new_number_FR_H1="5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25"
	elif len(FR_H1) == 22: new_number_FR_H1="4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25"
	elif len(FR_H1) == 23: new_number_FR_H1="3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25"
	elif len(FR_H1) == 24: new_number_FR_H1="2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25"
	elif len(FR_H1) == 25: new_number_FR_H1="1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25"
	elif len(FR_H1) == 26:
		#new_number_FR_H1="0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25"
		new_number_FR_H1="1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25"
		FR_H1 = FR_H1[1:]  # Remove 0th residue 10/24/2012
	else: print(("ERROR: FR_H1 length",len(FR_H1),"is not between 16 and 26"))

	# H1
	if   len(H1) ==  6: new_number_H1="26,27,32,33,34,35"
	elif len(H1) ==  7: new_number_H1="26,27,28,32,33,34,35"
	elif len(H1) ==  8: new_number_H1="26,27,28,29,32,33,34,35"
	elif len(H1) ==  9: new_number_H1="26,27,28,29,30,32,33,34,35"
	elif len(H1) == 10: new_number_H1="26,27,28,29,30,31,32,33,34,35"
	elif len(H1) == 11: new_number_H1="26,27,28,29,30,31,31A,32,33,34,35"
	elif len(H1) == 12: new_number_H1="26,27,28,29,30,31,31A,31B,32,33,34,35"
	elif len(H1) == 13: new_number_H1="26,27,28,29,30,31,31A,31B,31C,32,33,34,35"
	else: print(("ERROR: H1 length",len(H1),"is not between 6 and 13"))

	# FR_H2
	if len(FR_H2) == 14: new_number_FR_H2="36,37,38,39,40,41,42,43,44,45,46,47,48,49"
	else: print(("ERROR: FR_H2 length",len(FR_H2),"is not 14"))

	# H2
	if   len(H2) == 12: new_number_H2="50,51,52,57,58,59,60,61,62,63,64,65"
	elif len(H2) == 13: new_number_H2="50,51,52,56,57,58,59,60,61,62,63,64,65"
	elif len(H2) == 14: new_number_H2="50,51,52,55,56,57,58,59,60,61,62,63,64,65"
	elif len(H2) == 15: new_number_H2="50,51,52,54,55,56,57,58,59,60,61,62,63,64,65"
	elif len(H2) == 16: new_number_H2="50,51,52,53,54,55,56,57,58,59,60,61,62,63,64,65"
	elif len(H2) == 17: new_number_H2="50,51,52,52A,53,54,55,56,57,58,59,60,61,62,63,64,65"
	elif len(H2) == 18: new_number_H2="50,51,52,52A,52B,53,54,55,56,57,58,59,60,61,62,63,64,65"
	elif len(H2) == 19: new_number_H2="50,51,52,52A,52B,52C,53,54,55,56,57,58,59,60,61,62,63,64,65"
	elif len(H2) == 20: new_number_H2="50,51,52,52A,52B,52C,52D,53,54,55,56,57,58,59,60,61,62,63,64,65"
	elif len(H2) == 21: new_number_H2="50,51,52,52A,52B,52C,52D,52E,53,54,55,56,57,58,59,60,61,62,63,64,65"
	elif len(H2) == 22: new_number_H2="50,51,52,52A,52B,52C,52D,52E,52F,53,54,55,56,57,58,59,60,61,62,63,64,65"
	else: print(("ERROR: H2 length",len(H2),"is not between 12 and 22"))

	# FR_H3
	if   len(FR_H3) == 30: new_number_FR_H3="66,67,68,69,70,71,72,73,76,77,78,79,80,81,82,82A,82B,82C,83,84,85,86,87,88,89,90,91,92,93,94"
	elif len(FR_H3) == 31: new_number_FR_H3="66,67,68,69,70,71,72,73,74,76,77,78,79,80,81,82,82A,82B,82C,83,84,85,86,87,88,89,90,91,92,93,94"
	elif len(FR_H3) == 32: new_number_FR_H3="66,67,68,69,70,71,72,73,74,75,76,77,78,79,80,81,82,82A,82B,82C,83,84,85,86,87,88,89,90,91,92,93,94"
	else: print(("ERROR: FR_H3 length",len(FR_H3),"is not between 30 and 32"))

	# H3
	if   len(H3) ==  3: new_number_H3="95,101,102"
	elif len(H3) ==  4: new_number_H3="95,96,101,102"
	elif len(H3) ==  5: new_number_H3="95,96,97,101,102"
	elif len(H3) ==  6: new_number_H3="95,96,97,98,101,102"
	elif len(H3) ==  7: new_number_H3="95,96,97,98,99,101,102"
	elif len(H3) ==  8: new_number_H3="95,96,97,98,99,100,101,102"
	elif len(H3) ==  9: new_number_H3="95,96,97,98,99,100,100A,101,102"
	elif len(H3) == 10: new_number_H3="95,96,97,98,99,100,100A,100B,101,102"
	elif len(H3) == 11: new_number_H3="95,96,97,98,99,100,100A,100B,100C,101,102"
	elif len(H3) == 12: new_number_H3="95,96,97,98,99,100,100A,100B,100C,100D,101,102"
	elif len(H3) == 13: new_number_H3="95,96,97,98,99,100,100A,100B,100C,100D,100E,101,102"
	elif len(H3) == 14: new_number_H3="95,96,97,98,99,100,100A,100B,100C,100D,100E,100F,101,102"
	elif len(H3) == 15: new_number_H3="95,96,97,98,99,100,100A,100B,100C,100D,100E,100F,100G,101,102"
	elif len(H3) == 16: new_number_H3="95,96,97,98,99,100,100A,100B,100C,100D,100E,100F,100G,100H,101,102"
	elif len(H3) == 17: new_number_H3="95,96,97,98,99,100,100A,100B,100C,100D,100E,100F,100G,100H,100I,101,102"
	elif len(H3) == 18: new_number_H3="95,96,97,98,99,100,100A,100B,100C,100D,100E,100F,100G,100H,100I,100J,101,102"
	elif len(H3) == 19: new_number_H3="95,96,97,98,99,100,100A,100B,100C,100D,100E,100F,100G,100H,100I,100J,100K,101,102"
	elif len(H3) == 20: new_number_H3="95,96,97,98,99,100,100A,100B,100C,100D,100E,100F,100G,100H,100I,100J,100K,100L,101,102"
	elif len(H3) == 21: new_number_H3="95,96,97,98,99,100,100A,100B,100C,100D,100E,100F,100G,100H,100I,100J,100K,100L,100M,101,102"
	elif len(H3) == 22: new_number_H3="95,96,97,98,99,100,100A,100B,100C,100D,100E,100F,100G,100H,100I,100J,100K,100L,100M,100N,101,102"
	elif len(H3) == 23: new_number_H3="95,96,97,98,99,100,100A,100B,100C,100D,100E,100F,100G,100H,100I,100J,100K,100L,100M,100N,100O,101,102"
	elif len(H3) == 24: new_number_H3="95,96,97,98,99,100,100A,100B,100C,100D,100E,100F,100G,100H,100I,100J,100K,100L,100M,100N,100O,100P,101,102"
	elif len(H3) == 25: new_number_H3="95,96,97,98,99,100,100A,100B,100C,100D,100E,100F,100G,100H,100I,100J,100K,100L,100M,100N,100O,100P,100Q,101,102"
	elif len(H3) == 26: new_number_H3="95,96,97,98,99,100,100A,100B,100C,100D,100E,100F,100G,100H,100I,100J,100K,100L,100M,100N,100O,100P,100Q,100R,101,102"
	elif len(H3) == 27: new_number_H3="95,96,97,98,99,100,100A,100B,100C,100D,100E,100F,100G,100H,100I,100J,100K,100L,100M,100N,100O,100P,100Q,100R,100S,101,102"
	elif len(H3) == 28: new_number_H3="95,96,97,98,99,100,100A,100B,100C,100D,100E,100F,100G,100H,100I,100J,100K,100L,100M,100N,100O,100P,100Q,100R,100S,100T,101,102"
	elif len(H3) == 29: new_number_H3="95,96,97,98,99,100,100A,100B,100C,100D,100E,100F,100G,100H,100I,100J,100K,100L,100M,100N,100O,100P,100Q,100R,100S,100T,100U,101,102"
	elif len(H3) == 30: new_number_H3="95,96,97,98,99,100,100A,100B,100C,100D,100E,100F,100G,100H,100I,100J,100K,100L,100M,100N,100O,100P,100Q,100R,100S,100T,100U,100V,101,102"
	elif len(H3) == 31: new_number_H3="95,96,97,98,99,100,100A,100B,100C,100D,100E,100F,100G,100H,100I,100J,100K,100L,100M,100N,100O,100P,100Q,100R,100S,100T,100U,100V,100W,101,102"
	elif len(H3) == 32: new_number_H3="95,96,97,98,99,100,100A,100B,100C,100D,100E,100F,100G,100H,100I,100J,100K,100L,100M,100N,100O,100P,100Q,100R,100S,100T,100U,100V,100W,100X,101,102"
	elif len(H3) == 33: new_number_H3="95,96,97,98,99,100,100A,100B,100C,100D,100E,100F,100G,100H,100I,100J,100K,100L,100M,100N,100O,100P,100Q,100R,100S,100T,100U,100V,100W,100X,100Y,101,102"
	elif len(H3) == 34: new_number_H3="95,96,97,98,99,100,100A,100B,100C,100D,100E,100F,100G,100H,100I,100J,100K,100L,100M,100N,100O,100P,100Q,100R,100S,100T,100U,100V,100W,100X,100Y,100Z,101,102"
	else: print(("ERROR: H3 length",len(FR_H3),"is not between 3 and 34"))

	# FR_H4
	new_number_FR_H4="103,104,105,106,107,108,109,110,111,112,113,114"

	try:
		(new_number_L1 and new_number_L2 and new_number_L3 and new_number_H1 and new_number_H2 and new_number_H3
		 and new_number_FR_H1 and new_number_FR_H2 and new_number_FR_H3 and new_number_FR_H4
		 and new_number_FR_L1 and new_number_FR_L2 and new_number_FR_L3 and new_number_FR_L4 )
	except:
		print("Numbering failed.  Exiting.")
		sys.exit(1)

	# Converting all new_number_* vars in to a lists
	new_number_L1=new_number_L1.split(',');  new_number_L2=new_number_L2.split(',');  new_number_L3=new_number_L3.split(',');
	new_number_H1=new_number_H1.split(',');  new_number_H2=new_number_H2.split(',');  new_number_H3=new_number_H3.split(',');
	new_number_FR_L1=new_number_FR_L1.split(',');  new_number_FR_L2=new_number_FR_L2.split(',');  new_number_FR_L3=new_number_FR_L3.split(',');  new_number_FR_L4=new_number_FR_L4.split(',');
	new_number_FR_H1=new_number_FR_H1.split(',');  new_number_FR_H2=new_number_FR_H2.split(',');  new_number_FR_H3=new_number_FR_H3.split(',');  new_number_FR_H4=new_number_FR_H4.split(',');

	print_seq_FR_L1, print_seq_FR_L2, print_seq_FR_L3, print_seq_FR_L4, print_seq_FR_L4_extra = '', '', '', '', ''
	print_seq_FR_H1, print_seq_FR_H2, print_seq_FR_H3, print_seq_FR_H4, print_seq_FR_H4_extra = '', '', '', '', ''
	numbering_L, numbering_H = {}, {}

	# OUTPUT FOR LIGHT CHAIN. This should be save in the file 'numbering_L.txt'.
	for i, s in enumerate(FR_L1): #FR_L1
		numbering_L[ new_number_FR_L1[i] ] = s  #+='%s %s\n' % (s, new_number_FR_L1[i])
		v = int_(new_number_FR_L1[i])
		if v >= 10 and v <= 23: print_seq_FR_L1 += s

	for i, s in enumerate(L1): numbering_L[ new_number_L1[i] ] = s  # +='%s %s\n' % (s, new_number_L1[i])  # L1

	for i, s in enumerate(FR_L2):  #FR_L2
		numbering_L [ new_number_FR_L2[i] ] = s  #+='%s %s\n' % (s, new_number_FR_L2[i])
		v = int_(new_number_FR_L2[i])
		if (v >= 35 and v <= 38) or (v >= 45 and v <= 49): print_seq_FR_L2 += s

	for i, s in enumerate(L2): numbering_L[ new_number_L2[i] ] = s  #+='%s %s\n' % (s, new_number_L2[i])  # L2

	for i, s in enumerate(FR_L3):  #FR_L3
		numbering_L[ new_number_FR_L3[i] ] = s  # +='%s %s\n' % (s, new_number_FR_L3[i])
		v = int_(new_number_FR_L3[i])
		if (v >= 57 and v <= 66) or (v >= 71 and v <= 88): print_seq_FR_L3 += s

	for i, s in enumerate(L3): numbering_L[ new_number_L3[i] ] = s  #+='%s %s\n' % (s, new_number_L3[i])  # L3

	for i, s in enumerate(FR_L4):  #FR_L4
		numbering_L[ new_number_FR_L4[i] ] = s  #+='%s %s\n' % (s, new_number_FR_L4[i])
		v = int_(new_number_FR_L4[i])
		if v >= 98 and v <= 107: print_seq_FR_L4 += s
		if v >= 98 and v <= 101: print_seq_FR_L4_extra += s


	# OUTPUT FOR HEAVY CHAIN. This should be save in the file 'numbering_H.txt'.
	for i, s in enumerate(FR_H1):  #FR_H1
		numbering_H[ new_number_FR_H1[i] ] = s  #+='%s %s\n' % (s, new_number_FR_H1[i])
		v = int_(new_number_FR_H1[i])
		if v >= 10 and v <= 25: print_seq_FR_H1 += s

	for i, s in enumerate(H1): numbering_H[ new_number_H1[i] ] = s  #+='%s %s\n' % (s, new_number_H1[i])  # H1

	for i, s in enumerate(FR_H2):  #FR_H2
		numbering_H[ new_number_FR_H2[i] ] = s  #+='%s %s\n' % (s, new_number_FR_H2[i])
		v = int_(new_number_FR_H2[i])
		if (v >= 36 and v <= 39) or (v >= 46 and v <= 49): print_seq_FR_H2 += s

	for i, s in enumerate(H2): numbering_H[ new_number_H2[i] ] = s  #+='%s %s\n' % (s, new_number_H2[i])  # H2

	for i, s in enumerate(FR_H3):  #FR_H3
		numbering_H[ new_number_FR_H3[i] ] = s  #+='%s %s\n' % (s, new_number_FR_H3[i])
		v = int_(new_number_FR_H3[i])
		if v >= 66 and v <= 94: print_seq_FR_H3 += s

	for i, s in enumerate(H3): numbering_H[ new_number_H3[i] ] = s  #+='%s %s\n' % (s, new_number_H3[i])  # H3

	for i, s in enumerate(FR_H4):  #FR_H4
		numbering_H[ new_number_FR_H4[i] ] = s  #+='%s %s\n' % (s, new_number_FR_H4[i])
		v = int_(new_number_FR_H4[i])
		if v >= 103 and v <= 112: print_seq_FR_H4 += s
		if v >= 103 and v <= 106: print_seq_FR_H4_extra += s

	#if Options.verbose: print 'numbering_L:\n', numbering_L
	#if Options.verbose: print 'numbering_H:\n', numbering_H

	FRL = print_seq_FR_L1 + print_seq_FR_L2 + print_seq_FR_L3 + print_seq_FR_L4
	FRH = print_seq_FR_H1 + print_seq_FR_H2 + print_seq_FR_H3 + print_seq_FR_H4

	if len(FRL) != 58 and len(FRL) != 60 and len(FRL) != 61 and len(FRL) != 63 :
		print("ERROR: Current DB does not cover the length of FRL of your query.")
		print(("ERROR: FRL length of your query:", len(FRL)))
		print("ERROR: DB: 58 or 60")
		sys.exit(1)

	if len(FRH) != 63 and len(FRH) != 65 and len(FRH) != 66 and len(FRH) != 68:
		print("ERORR: Current DB does not cover the length of FRL of your query.")
		print(("ERORR: FRH length of your query:", len(FRH)))
		print("ERORR: DB: 63 or 65")
		sys.exit(1)

	return dict(FRL = print_seq_FR_L1 + print_seq_FR_L2 + print_seq_FR_L3 + print_seq_FR_L4,
				FRH = print_seq_FR_H1 + print_seq_FR_H2 + print_seq_FR_H3 + print_seq_FR_H4,
				light = FR_L1 + L1 + FR_L2 + L2 + FR_L3 + L3 + print_seq_FR_L4,
				heavy = FR_H1 + H1 + FR_H2 + H2 + FR_H3 + H3 + print_seq_FR_H4,
				numbering_L=numbering_L, numbering_H=numbering_H)

def get_RMSD(name1, path1, name2, path2, heavy, light):
	pdb_parser = PDBParser(QUIET=True)
	# Generate PDB structures using PDB parser
	file1 = gzip.open(path1, 'rt')
	file2 = gzip.open(path2, 'rt')
	pdb_strs = [pdb_parser.get_structure(name1, file1), pdb_parser.get_structure(name2, file2)]
	# Initiate an empty list for storing PDB sequences
	# Initiate CA polypeptide builder used to get sequences of each protein
	ppb = CaPPBuilder()
	seqs=['','']
	for i in range(2):
		for pp in ppb.build_peptides(pdb_strs[i]):
			seqs[i] += str(pp.get_sequence())
	pdb_index=[seqs[0].find(heavy[:10]),seqs[1].find(heavy[:10]),seqs[0].find(light[:10]),seqs[1].find(light[:10])]	
	print(seqs)
	# Get BLOSUM62 matrix
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
	pdb_atms_mapped = [[], [],[],[]]
	# i is the index for the two alignments, j is for the first
	# atom object list and k is for the other atom object list
	i, j, k = 0, 0, 0
	while i < len(heavy):
		pdb_atms_mapped[0].append(pdb_atms[0][pdb_index[0]+i])
		pdb_atms_mapped[1].append(pdb_atms[1][pdb_index[1]+i])
		i += 1
	i, j, k = 0, 0, 0
	while i < len(light):
		pdb_atms_mapped[0].append(pdb_atms[0][pdb_index[0]+i])
		pdb_atms_mapped[1].append(pdb_atms[1][pdb_index[1]+i])
		i += 1
	# Do the pairwise alignment and get alignments
	# Initiate the superimposer
	superimposer = Superimposer()
	# Set (translate/rotate) atoms minimizing RMSD
	print(pdb_atms_mapped[0]+pdb_atms_mapped[2])
	print(len(pdb_atms_mapped[0]+pdb_atms_mapped[2]))
	print(pdb_atms_mapped[1]+pdb_atms_mapped[3])
	#pdb_atms_mapped[0].extend(pdb_atms_mapped[2])
	#pdb_atms_mapped[1].extend(pdb_atms_mapped[3])
	superimposer.set_atoms(pdb_atms_mapped[0],pdb_atms_mapped[1])
	return float(superimposer.rms)



folder = sys.argv[1]
ignore_rmsd = int(sys.argv[2])
folderh = sys.argv[1].replace(' ', '_')[0:6].upper()
folderl = sys.argv[1].replace(' ', '_')[0:5].upper() + sys.argv[1][6].upper()
path = os.path.dirname(os.path.realpath(sys.argv[0]))
pdb = sys.argv[1][0:4]

if ignore_rmsd == 0 or ignore_rmsd == 1:
	for the_file in os.listdir(path + '/' + folder):
		file_path = os.path.join(path + '/' + folder, the_file)
		try:
			if os.path.isfile(file_path):
				os.unlink(file_path)
		except Exception as e:
			print(e)



if ignore_rmsd == 0:
	os.system('( cd ' + path + '/' + folder + ' && python3 ../clean_pdb.py ' + folderh.replace('_', ' ') + ' )')
	os.system('( cd ' + path + '/' + folder + ' && python3 ../clean_pdb.py ' + folderl.replace('_', ' ') + ' )')
	heavy_chain = read_fasta_file('./' + folder +'/' + folderh + '.fasta')[0]
	light_chain = read_fasta_file('./' + folder +'/' + folderl + '.fasta')[0]
	
	os.remove(os.path.join(path + '/' + folder +'/', folderh + '.fasta'))
	os.remove(os.path.join(path + '/' + folder +'/', folderl + '.fasta'))
	os.remove(os.path.join(path + '/' + folder +'/', folderh + '.pdb'))
	os.remove(os.path.join(path + '/' + folder +'/', folderl + '.pdb'))

	CDRs = IdentifyCDRs(light_chain, heavy_chain)
	CDRs.update( Extract_FR_CDR_Sequences(**CDRs) )
	#heavy_light = CDRs['heavy'] + '/' + CDRs['light']
	heavy_light = CDRs['light'] + '/' + CDRs['heavy']
	copyfile('./relax.options', path + '/' + folder + '/relax.options')
	write_fasta_file(folder, heavy_light.replace('\n', ''), path + '/' + folder + '/')
	if os.path.isfile('./templates/' + folder[0:4].lower() + '.pdb.gz') is True:
		with gzip.open('./templates/' + folder[0:4].lower() + '.pdb.gz', 'rb') as f_in:
			with open('./' + folder + '/' + folder[0:4].upper() + '.pdb', 'wb') as f_out:
				shutil.copyfileobj(f_in, f_out)
	else:
		os.system('( cd ' + path + '/' + folder + ' && python3 ../clean_pdb.py ' + folder.replace('_', ' ') + ' )')
		os.remove(os.path.join(path + '/' + folder +'/', folderh + '.fasta'))
		os.remove(os.path.join(path + '/' + folder +'/', folderl + '.fasta'))
		with open('./' + folder + '/' + folder[0:4].upper() + '.pdb', 'wb') as f_out:
			files = {'pdb': (folder +'.pdb', open(path + '/' + folder +'/' + folder +'.pdb', 'rb'))}
			values = {'output':' ', 'scheme':'-c', 'dofile':'l'}
			r = requests.post(url = 'http://www.bioinf.org.uk/abs/abnum/abnumpdb.cgi', files=files, data=values)
			f_out.write(r.content)
		
	os.system('python pdb2fasta.py ' + './' + folder +'/' + folder[0:4].upper() + '.pdb > '+ './' + folder +'/' + folder[0:4].upper() + '.fasta')
	os.system('( cd ' + path + '/' + folder + ' && cat *.fasta > alignment.fasta )')
	os.system('clustalo -i ' + path + '/' + folder + '/alignment.fasta -o ' + path + '/' + folder + '/alignment.aln')
	referenceContents = ' '
	pdbCode = []
	pdbContents = []
	isReference = 0
	with open(path + '/' + folder + '/alignment.aln') as file:
		for line in file.readlines():
			if line.find('>') == -1 and isReference == 0:
				pdbContents[len(pdbContents) - 1] =  pdbContents[len(pdbContents) - 1] + line
			elif line.find('>') == -1 and isReference == 1:
				referenceContents = referenceContents + line
			else:
				if line.find(folder) > -1:
					isReference = 1
				else:
					isReference = 0
					pdbCode.append(line[1:].replace('\n', ''))
					pdbContents.append(' ')
	pdbContents[0] = pdbContents[0].replace('\n','')
	pdbContents[0] = pdbContents[0] + '-'
	referenceContents = referenceContents.replace('\n','')
	referenceContents = referenceContents.replace('/','-/')
	pdbCode[0] = 'temp' + pdbCode[0][4:]
	newFileName = (pdb.lower() + "_" + pdbCode[0][:4].lower() + ".grishin")
	f = open(os.path.join(path + '/' + folder, newFileName), 'w')
	f.write("## {} {}.pdb\n".format(pdbCode[0][:4].lower(), pdb.lower()))
	f.write("#\n")
	f.write("scores_from_program: 0\n")
	f.write("0{}\n".format(referenceContents))
	f.write("0{}\n".format(pdbContents[0]))
	f.close()
	os.system('( cd ' + path + '/' + folder + ' && /dors/meilerlab/apps/rosetta/rosetta-3.7/main/source/bin/relax.default.linuxgccrelease -s ' + folder[0:4] + '.pdb @relax.options -nstruct 1 )')
	while not os.path.exists(os.path.join(path + '/' + folder, (pdb + "_0001.pdb.gz"))):
		print(os.path.join(path + '/' + folder, (pdb + "_0001.pdb.gz")))
		time.sleep(1)
	#os.remove(os.path.join(path + '/' + folder, (pdb + ".pdb")))
	#os.rename(os.path.join(path + '/' + folder, (pdb + "_0001.pdb.gz")), os.path.join(path + '/' + folder, (pdbCode[0][:4].lower() + ".pdb.gz")))
	os.rename(os.path.join(path + '/' + folder, (pdb + "_0001.pdb.gz")), os.path.join(path + '/' + folder, (pdb.lower() + ".pdb.gz")))
	os.system('( cd ' + path + '/' + folder + ' && /dors/meilerlab/apps/rosetta/rosetta-3.7/main/source/bin/partial_thread.default.linuxgccrelease -in:file:fasta {} -in:file:alignment {} -in:file:template_pdb {} )'.format((folder.upper()) + '.fasta',newFileName,(pdb.lower()) + ".pdb.gz"))
	while not os.path.exists(os.path.join(path + '/' + folder, pdb.lower() + ".pdb.pdb")):
		time.sleep(1)

	os.rename(os.path.join(path + '/' + folder, (pdb.lower() + ".pdb.pdb")), os.path.join(path + '/' + folder, (pdbCode[0][:4].lower() + ".pdb.pdb")))
	if not os.path.exists(path + '/' + folder + '/Control'):
		os.mkdir(path + '/' + folder + '/Control')
	with open('./rosetta_cm.options', 'r') as file:
		Options = file.read()
		newOptions = Options.replace('InsertFastaHere', pdb + '.fasta')

	with open('./rosetta_cm.xml', 'r') as file:
		XML = file.read()
		XML = XML.replace('3FRAGSFILE', pdb + '_3.frags')
		XML = XML.replace('9FRAGSFILE', pdb + '_9.frags')

	with open('./run_array.slurm', 'r') as file:
		RunSlurm = file.read()

	with open('./run_relax.slurm', 'r') as file:
		RunRelax = file.read()


	with open(path + '/' + folder + '/Control/rosetta_cm.options', 'w') as file:
		file.write(newOptions)
	RunAllSlurm = RunSlurm.replace('PATHL', 'Control')
	RunAllRelax = RunRelax.replace('PATHL', 'ControlRelax')
	RunAllSlurm = RunAllSlurm.replace('OPTIONSP', path + '/' + folder + '/Control' + '/rosetta_cm.options')
	RunAllRelax = RunAllRelax.replace('OPTIONSP', path + '/' + folder + '/Control' + '/relax.options')
	RunAllSlurm = RunAllSlurm.replace('NUMBEROFTRIALS', '2')
	RunAllSlurm = RunAllSlurm.replace('THREADS', '10')
	RunAllRelax = RunAllRelax.replace('THREADS', '10')

	with open(path + '/' + folder + '/Control' + '/run_array.slurm', 'w') as file:
		file.write(RunAllSlurm)
	with open(path + '/' + folder + '/Control' + '/run_relax.slurm', 'w') as file:
		file.write(RunAllRelax)
	copyfile('./relax.options', path + '/' + folder + '/Control' + '/relax.options')
	copyfile('./sc_parser.bash', path + '/' + folder + '/Control' + '/sc_parser.bash')
	copyfile(os.path.join(path + '/' + folder, (pdbCode[0][:4].lower()) + ".pdb.pdb"), os.path.join(path + '/' + folder + '/Control' , (pdbCode[0][:4].lower()) + ".pdb.pdb"))
	copyfile(os.path.join(path + '/' + folder, (folder) + ".fasta"), os.path.join(path + '/' + folder + '/Control' , (pdb[:4].upper()) + ".fasta"))
	XMLstring = ' <Template pdb="temp.pdb.pdb" cst_file="AUTO" weight=   1.000 />'
	allXML = XMLstring + '\n'
	with open(path + '/' + folder  + '/Control/rosetta_cm.xml', 'w') as file:
		file.write(XML.replace('THREADING', allXML))


else:
	if os.path.isfile(folder + '/' + folder[0:4] + '.pdb') is True:
		os.remove(folder + '/' + folder[0:4] + '.pdb')
	os.system('( cd ' + path + '/' + folder + ' && python3 ../clean_pdb.py ' + folderh.replace('_', ' ') + ' )')
	os.system('( cd ' + path + '/' + folder + ' && python3 ../clean_pdb.py ' + folderl.replace('_', ' ') + ' )')

	with open('./' + folder +'/1/'  + folder[0:4].lower() + '.fasta') as reference:
		lines = reference.readlines()
		reference_fasta = lines[1]

	heavy_chain = read_fasta_file('./' + folder +'/' + folderh + '.fasta')[0]
	light_chain = read_fasta_file('./' + folder +'/' + folderl + '.fasta')[0]
	os.remove(os.path.join(path + '/' + folder +'/', folderh + '.fasta'))
	os.remove(os.path.join(path + '/' + folder +'/', folderl + '.fasta'))
	os.remove(os.path.join(path + '/' + folder +'/', folderh + '.pdb'))
	os.remove(os.path.join(path + '/' + folder +'/', folderl + '.pdb'))

	CDRs = IdentifyCDRs(light_chain, heavy_chain)
	CDRs.update( Extract_FR_CDR_Sequences(**CDRs) )
	if ignore_rmsd == 1:
		
		'''if folder[6] > folder[5]:
			heavy_light = CDRs['heavy'] + '/' + CDRs['light']
		else:
			heavy_light = CDRs['light'] + '/' + CDRs['heavy']'''

		#print(reference_fasta)

		if str(reference_fasta[2:10]) in reference_fasta:
			heavy_light = CDRs['heavy'] + '/' + CDRs['light']
		else:
			heavy_light = CDRs['light'] + '/' + CDRs['heavy']

		print(heavy_light)		

		copyfile('./relax.options', path + '/' + folder + '/relax.options')
		os.system('( cd ' + path + '/' + folder + ' && python3 ../clean_pdb.py ' + folder[0:4] + ' ignorechain )')
		write_fasta_file(folder, heavy_light.replace('\n', ''), path + '/' + folder + '/')
		#write_fasta_file(folder, (heavy_chain + '/' + light_chain).replace('\n', ''), path + '/' + folder + '/')

		'''heavy_light = CDRs['heavy'] + '/' + CDRs['light']
		copyfile('./relax.options', path + '/' + folder + '/relax.options')
		write_fasta_file(folder, heavy_light.replace('\n', ''), path + '/' + folder + '/')
		if os.path.isfile('./templates/' + folder[0:4].lower() + '.pdb.gz') is True:
			with gzip.open('./templates/' + folder[0:4].lower() + '.pdb.gz', 'rb') as f_in:
				with open('./' + folder + '/' + folder[0:4].upper() + '.pdb', 'wb') as f_out:
					shutil.copyfileobj(f_in, f_out)
		else:
			os.system('( cd ' + path + '/' + folder + ' && python3 ../clean_pdb.py ' + folder.replace('_', ' ') + ' )')
			os.remove(os.path.join(path + '/' + folder +'/', folderh + '.fasta'))
			os.remove(os.path.join(path + '/' + folder +'/', folderl + '.fasta'))
			with open('./' + folder + '/' + folder[0:4].upper() + '.pdb', 'wb') as f_out:
				files = {'pdb': (folder +'.pdb', open(path + '/' + folder +'/' + folder +'.pdb', 'rb'))}
				values = {'output':' ', 'scheme':'-c', 'dofile':'l'}
				r = requests.post(url = 'http://www.bioinf.org.uk/cgi-bin/abnum/abnumpdb.pl', files=files, data=values)
				f_out.write(r.content)
		
		os.system('( cd ' + path + '/' + folder + ' && python3 ../clean_pdb.py ' + folder + '.pdb  )')'''

		

		pymol.cmd.reinitialize
		pymol.cmd.load(path + '/' + folder + '/' + folder[0:4] + '.pdb', 'temp')
		pymol.cmd.do('run findseq.py')
		pymol.cmd.do('run renumber.py')
		pymol.cmd.do('findseq ' +  CDRs['heavy'] + ',temp and chain ' + folderh[-1] + ', heavy')
		pymol.cmd.do('findseq ' +  CDRs['light'] + ',temp and chain ' + folderl[-1] + ', light')
		pymol.cmd.do("alter light, chain = 'L'")
		pymol.cmd.do("alter heavy, chain = 'H'")
		pymol.cmd.do('create native = (heavy)+(light)')
		pymol.cmd.do('save ' + path + '/' + folder + '/native.pdb, native')

		os.system('( cd ' + path + '/' + folder + ' && /dors/meilerlab/apps/rosetta/rosetta-3.7/main/source/bin/relax.default.linuxgccrelease -s native.pdb @relax.options -nstruct 1 )')
		while not os.path.exists(os.path.join(path + '/' + folder, ("native_0001.pdb.gz"))):
			print(os.path.join(path + '/' + folder, ("native_0001.pdb.gz")))
			time.sleep(1)
		os.remove(os.path.join(path + '/' + folder, ("native.pdb")))
		os.rename(os.path.join(path + '/' + folder, ("native_0001.pdb.gz")), os.path.join(path + '/' + folder, ("native.pdb.gz")))
		#os.system('( cd ' + path + '/' + folder + ' && /dors/meilerlab/apps/rosetta/rosetta-3.7/main/source/bin/partial_thread.default.linuxgccrelease -in:file:fasta {} -in:file:alignment {} -in:file:template_pdb {} )'.format((folder[:4].upper()) + '.fasta',newFileName,(pdbCode[0][:4].lower()) + ".pdb.gz"))
		#while not os.path.exists(os.path.join(path + '/' + folder, (pdbCode[0][:4].lower()) + ".pdb.pdb")):
			#time.sleep(1)
		'''
		q = []
	
		for root, dirnames, filenames in os.walk(path + '/' + folder):
			for filename in fnmatch.filter(filenames, 'run_relax.slurm'):
				#copyfile('./' + folder + '/' + pdbCode[0][:4].lower() + ".pdb.pdb", root + '/' + pdbCode[0][:4].lower() + ".pdb.pdb")
				#os.rename(root + '/' + pdbCode[0][:4].lower() + ".pdb.pdb", root + '/' +  folder[:4].lower() + '_0001.pdb')
				#relaxCommand = '( cd ' + root + ' && /dors/meilerlab/apps/rosetta/rosetta-3.7/main/source/bin/score_jd2.linuxgccrelease -database /dors/meilerlab/apps/rosetta/rosetta-3.7/main/database/ -in:file:l relaxscore.txt -in:file:native ' + folder[:4].lower() + '_0001.pdb' + ' -out:file:scorefile rmsd_scores_relax.sc )'
				#subprocess.Popen(relaxCommand, shell=True)
				#os.system('( cd ' + root + ' && cat *_score.txt > score.txt )')
				#command = '( cd ' + root + ' && /dors/meilerlab/apps/rosetta/rosetta-3.7/main/source/bin/score_jd2.linuxgccrelease -database /dors/meilerlab/apps/rosetta/rosetta-3.7/main/database/ -in:file:l score.txt -in:file:native ' + folder[:4].lower() + '.pdb -out:file:scorefile rmsd_scores_hybridize.sc )'
				command = '( cd ' + root + ' && bash sc_parser.bash relaxscore.fasc && /dors/meilerlab/apps/rosetta/rosetta-3.7/main/source/bin/score_jd2.linuxgccrelease -database /dors/meilerlab/apps/rosetta/rosetta-3.7/main/database/ -in:file:l relaxscore.txt -out:level 300 -out:file:scorefile rmsd_scores_relax.sc )'
				p = subprocess.Popen(command, shell=True)
				q.append(p)

		exitcodes = [p.wait() for p in q]'''

	if os.path.isfile(path + '/' + folder + '/ScorePlots/dumpdict.p') is False:
	
		allresultsrelaxeddict = {}
		for root, dirnames, filenames in os.walk(path + '/' + folder):
			for subdir in dirnames:
				for rootx, dirnamesx, filenamesx in os.walk(root + '/' + subdir):
					for filenamex in fnmatch.filter(filenamesx, 'relaxscore.fasc'):
						with open(rootx + '/relaxscore.fasc') as f:
							f.readline()
							f.readline()
							for line in f:
								try:
									score = float(line[8:20].strip())
									if 'H3_modeling' in rootx.split('/')[-1]:
										pdbFile = line[(line.find('model')-2):-1].strip()
										print(pdbFile)
										pymol.cmd.reinitialize()
										pymol.cmd.load(rootx + '/' + pdbFile + '.pdb', 'Test')
										pymol.cmd.do('run renumber.py')
										pymol.cmd.do('renumber chain H')
										pymol.cmd.do('renumber chain L')
									else:
										pdbFile = line[line.find('relax'):-1].strip()
										pymol.cmd.reinitialize()
										pymol.cmd.load(rootx + '/' + pdbFile + '.pdb.gz', 'Test')
									#print(pdbFile)
									pymol.cmd.load(path + '/' + folder + '/native.pdb.gz', 'Native')
									pymol.cmd.do('run findseq.py')
									#rmsd = get_RMSD('Test', rootx + '/' + pdbFile + '.pdb.gz', folder, path + '/' + folder + '/temp.pdb.gz', CDRs['heavy'], CDRs['light'])
									alignrmsd = pymol.cmd.align('Test', 'Native', cycles=5)[0]
									superrmsd = pymol.cmd.super('Test & backbone', 'Native & backbone', cycles=0)[0]

									pymol.cmd.do('findseq ' +  CDRs['H1'] + ', Native, NativeH1')
									pymol.cmd.do('findseq ' +  CDRs['H1'] + ', Test, TestH1')
									pymol.cmd.do('alter NativeH1, chain="H"')
									pymol.cmd.do('alter TestH1, chain="H"')
									pymol.cmd.do('alter TestH1 & chain H, resi=str(1)')
									pymol.cmd.do('alter NativeH1 & chain H, resi=str(1)')
									H1RMS = pymol.cmd.rms_cur('TestH1 & backbone', 'NativeH1 & backbone')

									pymol.cmd.do('findseq ' +  CDRs['H2'] + ', Native, NativeH2')
									pymol.cmd.do('findseq ' +  CDRs['H2'] + ', Test, TestH2')
									pymol.cmd.do('alter NativeH2, chain="H"')
									pymol.cmd.do('alter TestH2, chain="H"')
									pymol.cmd.do('alter TestH2 & chain H, resi=str(1)')
									pymol.cmd.do('alter NativeH2 & chain H, resi=str(1)')
									H2RMS = pymol.cmd.rms_cur('TestH2 & backbone', 'NativeH2 & backbone')

									pymol.cmd.do('findseq ' +  CDRs['H3'] + ', Native, NativeH3')
									pymol.cmd.do('findseq ' +  CDRs['H3'] + ', Test, TestH3')
									pymol.cmd.do('alter NativeH3, chain="H"')
									pymol.cmd.do('alter TestH3, chain="H"')
									pymol.cmd.do('alter Testh3 & chain H, resi=str(1)')
									pymol.cmd.do('alter Nativeh3 & chain H, resi=str(1)')
									H3RMS = pymol.cmd.rms_cur('Testh3 & backbone', 'Nativeh3 & backbone')

									pymol.cmd.do('findseq ' +  CDRs['L1'] + ', Native, NativeL1')
									pymol.cmd.do('findseq ' +  CDRs['L1'] + ', Test, TestL1')
									pymol.cmd.do('alter NativeL1, chain="L"')
									pymol.cmd.do('alter TestL1, chain="L"')
									pymol.cmd.do('alter TestL1 & chain L, resi=str(1)')
									pymol.cmd.do('alter NativeL1 & chain L, resi=str(1)')
									L1RMS = pymol.cmd.rms_cur('TestL1 & backbone', 'NativeL1 & backbone')

									pymol.cmd.do('findseq ' +  CDRs['L2'] + ', Native, NativeL2')
									pymol.cmd.do('findseq ' +  CDRs['L2'] + ', Test, TestL2')
									pymol.cmd.do('alter NativeL2, chain="L"')
									pymol.cmd.do('alter TestL2, chain="L"')
									pymol.cmd.do('alter TestL2 & chain L, resi=str(1)')
									pymol.cmd.do('alter NativeL2 & chain L, resi=str(1)')
									L2RMS = pymol.cmd.rms_cur('TestL2 & backbone', 'NativeL2 & backbone')

									pymol.cmd.do('findseq ' +  CDRs['L3'] + ', Native, NativeL3')
									pymol.cmd.do('findseq ' +  CDRs['L3'] + ', Test, TestL3')
									pymol.cmd.do('alter NativeL3, chain="L"')
									pymol.cmd.do('alter TestL3, chain="L"')
									pymol.cmd.do('alter TestL3 & chain L, resi=str(1)')
									pymol.cmd.do('alter NativeL3 & chain L, resi=str(1)')
									L3RMS = pymol.cmd.rms_cur('TestL3 & backbone', 'NativeL3 & backbone')


									
									#print(other)
									#rmsd = float(line[33:44].strip())
									#print(rootx.split('/')[-1])
									if 'all' in rootx.split('/')[-1]:
										tname = 'Multi-Template'
									elif len(rootx.split('/')[-1]) < 2:
										tname = 'Rosetta Antibody Inspired'
									elif 'Single' in rootx.split('/')[-1]:
										tname = 'Single Template'
									elif 'H3_modeling' in rootx.split('/')[-1]:
										tname = 'Rosetta Antibody'
									elif 'AbPredict' in rootx.split('/')[-1]:
										tname = 'AbPredict'
										print(superrmsd)
										print(alignrmsd)
										print(H3RMS)	
									else:
										tname = 'Control'
									#print(rootx.split('/')[-1])
									if ('AbPredict' not in tname and 'Rosetta Antibody' not in tname) or ('AbPredict' in tname  and min([H1RMS, H2RMS, H3RMS, L1RMS, L2RMS, L3RMS]) > 0.01) or ('Rosetta Antibody' in tname  and min([H1RMS, H2RMS, H3RMS, L1RMS, L2RMS, L3RMS]) > 0):
										allresultsrelaxeddict.setdefault(tname, [[],[],[],[],[],[],[],[],[],[]])
										allresultsrelaxeddict[tname][0].append(score)
										allresultsrelaxeddict[tname][1].append(superrmsd)
										allresultsrelaxeddict[tname][2].append(alignrmsd)
										allresultsrelaxeddict[tname][3].append(H1RMS)
										allresultsrelaxeddict[tname][4].append(H2RMS)
										allresultsrelaxeddict[tname][5].append(H3RMS)
										allresultsrelaxeddict[tname][6].append(L1RMS)
										allresultsrelaxeddict[tname][7].append(L2RMS)
										allresultsrelaxeddict[tname][8].append(L3RMS)
										allresultsrelaxeddict[tname][9].append(rootx + '/' + pdbFile)
									#print(allresultsrelaxeddict[tname])
								except Exception as e:print(e)
		if not os.path.exists(path + '/' + folder + '/ScorePlots'):
			os.mkdir(path + '/' + folder + '/ScorePlots')
		pickle.dump(allresultsrelaxeddict, open(path + '/' + folder + '/ScorePlots/dumpdict.p', 'wb'))

	else:
		allresultsrelaxeddict = pickle.load(open(path + '/' + folder + '/ScorePlots/dumpdict.p', 'rb'))
	
	
	records = list()                                                                                        
	for k, v in sorted(allresultsrelaxeddict.items(), key = lambda x: x[0], reverse=True):        
		for ii in range(len(v[0])):                                                                                            
			records.append({"Score" : v[0][ii], "Backbone RMSD" : v[1][ii], "RMSDA" : v[2][ii], "CDR3 Backbone RMSD" : v[5][ii], "Experiment" : k})
	dframe = DataFrame(records)
	g = sns.lmplot(data=dframe, x="CDR3 Backbone RMSD", y="Score", hue="Experiment", fit_reg=False, size = 12, aspect = 0.85, scatter_kws={'alpha':0.5, 's':40}, legend_out=False)
	g.set(ylim=(None, 0), xlim=(0, 8))
	xwide = g.axes[0,0].get_xlim()
	legend = plt.legend(loc = 'upper right', frameon=True, framealpha=1)
	frame = legend.get_frame()
	frame.set_facecolor('white')
	frame.set_edgecolor('white')
	for lh in legend.legendHandles:
		lh.set_alpha(1)
		lh.set_sizes([100])
	plt.savefig(path + '/' + folder + '/ScorePlots/combined_H3_RMSD_graph.svg', dpi=600, format='svg', transparent=True)

	plt.clf()
	fig, ax = plt.subplots(figsize=(11, 5))
	histDict = {}
	for index, row in dframe.iterrows():
		histDict.setdefault(row['Experiment'], [])
		histDict[row['Experiment']].append(row['CDR3 Backbone RMSD'])
	i = 0
	for key, value in sorted(histDict.items(), key = lambda x: x[0], reverse=True):
		x = sns.distplot(value, color=sns.color_palette()[i], ax=ax, kde=True, hist=False, norm_hist=True, kde_kws={'clip': xwide})
		x.set_xlim(xwide)
		x.set(xticklabels=[], yticklabels=[])
		x.margins(0.5, 0.5)
		i+=1
	plt.tight_layout()
	plt.savefig(path + '/' + folder + '/ScorePlots/combined_H3_RMSD_graph_Histo.svg', dpi=600, format='svg', transparent=True)

	plt.clf()


	records = list()                                                                                        
	for k, v in sorted(allresultsrelaxeddict.items(), key = lambda x: x[0], reverse=True):        
		for ii in range(len(v[0])):                                                                                            
			records.append({"Score" : v[0][ii], "Backbone RMSD" : v[1][ii], "RMSDA" : v[2][ii], "CDR3 Backbone RMSD" : v[3][ii], "Experiment" : k})
	dframe = DataFrame(records)
	g = sns.lmplot(data=dframe, x="Backbone RMSD", y="Score", hue="Experiment", fit_reg=False, size = 12, aspect = 0.85, scatter_kws={'alpha':0.5, 's':40}, legend_out=False)
	g.set(ylim=(None, 0), xlim=(0, 8))
	xwide = g.axes[0,0].get_xlim()
	legend = plt.legend(loc = 'upper right', frameon=True, framealpha=1)
	frame = legend.get_frame()
	frame.set_facecolor('white')
	frame.set_edgecolor('white')
	for lh in legend.legendHandles:
		lh.set_alpha(1)
		lh.set_sizes([100])
	plt.savefig(path + '/' + folder + '/ScorePlots/combined_Overall_RMSD_graph.svg', dpi=600, format='svg', transparent=True)

	plt.clf()
	fig, ax = plt.subplots(figsize=(11, 5))
	histDict = {}
	for index, row in dframe.iterrows():
		histDict.setdefault(row['Experiment'], [])
		histDict[row['Experiment']].append(row['Backbone RMSD'])
	i = 0
	for key, value in sorted(histDict.items(), key = lambda x: x[0], reverse=True):
		x = sns.distplot(value, color=sns.color_palette()[i], ax=ax, kde=True, hist=False, norm_hist=True, kde_kws={'clip': xwide})
		x.set_xlim(xwide)
		x.set(xticklabels=[], yticklabels=[])
		x.margins(0.5, 0.5)
		i+=1
	plt.tight_layout()
	plt.savefig(path + '/' + folder + '/ScorePlots/combined_Overall_RMSD_graph_Histo.svg', dpi=600, format='svg', transparent=True)
	