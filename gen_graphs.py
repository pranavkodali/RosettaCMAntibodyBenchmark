import os, fnmatch, pickle, operator, math, re, sys, time #, statistics
import random
from Bio import pairwise2 as pw2
import matplotlib
matplotlib.use("agg")
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import pandas as pd
from pandas import DataFrame
import seaborn as sns
import numpy as np
from numpy import nan
from numpy import median
#from numpy import median
sns.set()
sns.set_context('poster')
sns.set_style('whitegrid')
sns.set_palette('Set2')


ScoreResults= {}
ScoreL1 = {}
ScoreL2 = {}
ScoreL3 = {}
ScoreH1 = {}
ScoreH2 = {}
ScoreH3 = {}

H1Sample = {}
H2Sample = {}
H3Sample = {}
L1Sample = {}
L2Sample = {}
L3Sample = {}
BackboneSample = {}

Error = {}

Outliers = {}

bins = ['Short', 'Medium', 'Long', 'Long']


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
		#print(("L1 detected: %s (%d residues at positions %d to %d)" % (L1, len(L1), L1_start, L1_end)))
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
		#print(("L3 detected: %s ( %d residues at positions %d to %d)" % (L3, len(L3), L3_start, L3_end)))
	else:
		print("L3 detected: False")

	if L1 and L3:
		#L1_start = light_chain.index(L1)
		#L1_end = L1_start + len(L1) - 1

		L2_start = L1_end + 16
		L2_end = L2_start + 7 - 1

		L2 = light_chain[L2_start:L2_start+7]  # L2 is identified here. Current implementation can deal with only 7-resiue L2
		#print(("L2 detected: %s (%d residues at positions %d to %d)" % (L2, len(L2), L2_start, L2_end)))

		#FR_L1 = light_chain[:L1_start]
		FR_L2 = light_chain[ L1_end + 1  :  L1_end + 1 + 15					]
		FR_L3 = light_chain[ L2_end + 1  :  L2_end + 1 + L3_start - L2_end - 1 ]
		FR_L4 = light_chain[ L3_end + 1  :  L3_end + 1 + 12					]

		#print(("FR_L1: ", FR_L1))
		#print(("FR_L2: ", FR_L2))
		#print(("FR_L3: ", FR_L3))
		#print(("FR_L4: ", FR_L4))
		#print(("L segments: ",FR_L1,L1,FR_L2,L2,FR_L3,L3,FR_L4))

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
		#print(("H1 detected: %s (%d residues at positions %d to %d)" % (H1, len(H1), H1_start, H1_end)))
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
		#print(("H3 detected: %s (%d residues at positions %d to %d)" % (H3, len(H3), H3_start, H3_end)))
	else:
		print("H3 detected: False")


	if H1 and H3:
		#H1_start = heavy_chain.index(H1)
		#H1_end = H1_start + len(H1) - 1
		H2_start = H1_end + 15
		H2_end = H3_start - 33
		H2 = heavy_chain[H2_start:H2_start + H2_end-H2_start+1]
		#print(("H2 detected: %s (%d residues at positions %d to %d)" % (H2, len(H2), H2_start, H2_end)))

		#FR_H1 = heavy_chain[:H1_start]

		#if len(FR_H1) >  26:
		#	FR_H1 = light_chain[20:H1_start]

		FR_H2 = heavy_chain[H1_end + 1: H1_end + 1 + H2_start - H1_end - 1]
		FR_H3 = heavy_chain[H2_end + 1: H2_end + 1 + H3_start - H2_end - 1]
		FR_H4 = heavy_chain[H3_end + 1: H3_end + 1 + 12]

		#print(("FR_H1: ", FR_H1))
		#print(("FR_H2: ", FR_H2))
		#print(("FR_H3: ", FR_H3))
		#print(("FR_H4: ", FR_H4))
		#print(("H segments: ",FR_H1,H1,FR_H2,H2,FR_H3,H3,FR_H4))

	if not (L1 and L3 and H1 and H3):
		if not L1: print(('ERROR: CDR L1 cannot be recognized !!!  L1 pattern: %s' % L1_pattern)) # C[A-Z]{1,17}(WYL|WLQ|WFQ|WYQ|WYH|WVQ|WVR|WWQ|WVK|WYR|WLL|WFL|WVF|WIQ|WYR|WNQ|WHL|WHQ|WYM|WYY)'
		if not L3: print(('ERROR: CDR L3 cannot be recognized !!!  L3 pattern: %s' % L3_pattern)) # C[A-Z]{1,15}(L|F|V|S)G[A-Z](G|Y)
		if not H1: print(('ERROR: CDR H1 cannot be recognized !!!  H1 pattern: %s' % H1_pattern)) # C[A-Z]{1,16}(W)(I|V|F|Y|A|M|L|N|G)(R|K|Q|V|N|C)(Q|K|H|E|L|R)
		if not H3: print(('ERROR: CDR H3 cannot be recognized !!!  H3 pattern: %s' % H3_pattern)) # C[A-Z]{1,33}(W)(G|A|C)[A-Z](Q|S|G|R)
		sys.exit(1)

	res = dict(L1=L1, L2=L2, L3=L3, H1=H1, H2=H2, H3=H3,  FR_L1=FR_L1, FR_L2=FR_L2, FR_L3=FR_L3, FR_L4=FR_L4,  FR_H1=FR_H1, FR_H2=FR_H2, FR_H3=FR_H3, FR_H4=FR_H4)
	#print 'L1: %(L1)s\nL2: %(L2)s\nL3: %(L3)s\nH1: %(H1)s\nH2: %(H2)s\nH3: %(H3)s' % res

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

def median(lst):
    n = len(lst)
    s = sorted(lst)
    return (sum(s[n//2-1:n//2+1])/2.0, s[n//2])[n % 2] if n else None



if os.path.isfile('./RMSD_Data.p') is False or os.path.isfile('./Sample_Data.p') is False or os.path.isfile('./Error_Data.p') is False:
	for root, dirnames, filenames in os.walk(os.getcwd()):
		for filename in fnmatch.filter(filenames, 'dumpdict.p'):
			allresultsrelaxeddict = pickle.load(open(root + '/dumpdict.p', 'rb'))
			print(root)
			sequences = []
			with open(os.path.dirname(root) + '/Single_Control/alignment.fasta', 'r') as alignment:
				for line in alignment.readlines():
					if '>' in line:
						sequences.append('')
					else:
						sequences[-1] += line.rstrip()
			for sequence in sequences:
				if '/' in sequence:
					chains = sequence.split('/')
					BBlength = len(sequence) + random.random() * .01
					print(BBlength)
					CDRs = IdentifyCDRs(chains[1], chains[0])
					CDRs.update( Extract_FR_CDR_Sequences(**CDRs) )
			H1_IDL = []
			H2_IDL = []
			H3_IDL = []
			L1_IDL = []
			L2_IDL = []
			L3_IDL = []
			for sequence in sequences:
				if '/' not in sequence:
					endres = sequence.find('TVS', 100, 150) + 3
					if endres < 100:
						endres = sequence.find('TLV', 100, 150) + 3
					Templatechains = [sequence[endres + 1:], sequence[:endres]]
					#print(Templatechains)
					TemplateBBlength = len(sequence) + random.random() * .01
					TemplateCDRs = IdentifyCDRs(Templatechains[0], Templatechains[1])
					#TemplateCDRs.update( Extract_FR_CDR_Sequences(**TemplateCDRs) )
			
				global_align = pw2.align.globalxx(sequences[0], sequences[1])
				percent_identity = global_align[0][2]/len(sequences[0])		

				H1Length = len(CDRs['H1'])
				global_align = pw2.align.globalxx(CDRs['H1'], TemplateCDRs['H1'])
				H1_IDL.append(global_align[0][2]/H1Length)

				H2Length = len(CDRs['H2'])
				global_align = pw2.align.globalxx(CDRs['H2'], TemplateCDRs['H2'])
				H2_IDL.append(global_align[0][2]/H2Length)

				H3Length = len(CDRs['H3'])
				global_align = pw2.align.globalxx(CDRs['H3'], TemplateCDRs['H3'])
				H3_IDL.append(global_align[0][2]/H3Length)

				L1Length = len(CDRs['L1'])
				global_align = pw2.align.globalxx(CDRs['L1'], TemplateCDRs['L1'])
				L1_IDL.append(global_align[0][2]/L1Length)
	
				L2Length = len(CDRs['L2'])
				global_align = pw2.align.globalxx(CDRs['L2'], TemplateCDRs['L2'])
				L2_IDL.append(global_align[0][2]/L2Length)

				L3Length = len(CDRs['L3'])
				global_align = pw2.align.globalxx(CDRs['L3'], TemplateCDRs['L3'])
				L3_IDL.append(global_align[0][2]/L3Length)
			H1_ID = max(H1_IDL)
			H2_ID = max(H2_IDL)
			H3_ID = max(H3_IDL)
			L1_ID = max(L1_IDL)
			L2_ID = max(L2_IDL)
			L3_ID = max(L3_IDL)

			for key, value in allresultsrelaxeddict.items():
		
				newlistOverall = []
				newlistH1 = []
				newlistH2 = []
				newlistH3 = []
				newlistL1 = []
				newlistL2 = []
				newlistL3 = []
				for i in range(len(value[0])):
					newlistOverall.append([value[0][i], value[1][i], value[9][i]])
					newlistH1.append([value[0][i], value[3][i], value[9][i]])
					newlistH2.append([value[0][i], value[4][i], value[9][i]])
					newlistH3.append([value[0][i], value[5][i], value[9][i]])
					newlistL1.append([value[0][i], value[6][i], value[9][i]])
					newlistL2.append([value[0][i], value[7][i], value[9][i]])
					newlistL3.append([value[0][i], value[8][i], value[9][i]])
				
				newlistH1.sort(key=operator.itemgetter(0))
				if key == 'Rosetta Antibody': key = 'RosettaAntibody'
				elif key == 'Rosetta Antibody Inspired': key = 'RosettaAntibody Inspired'
				print(key)
				if key == 'Multi-Template' or key == 'RosettaAntibody Inspired' or key == 'RosettaAntibody' or key == 'AbPredict':
					listlength = 10
				else:
					listlength = 2
				if len(newlistH1) < listlength: listlength = len(newlistH1)
				if key == 'Multi-Template':
					ScoreH1.setdefault('Reference', [])
					ScoreH1ReferenceRMSD = []
					for i in range(10):
						ScoreH1ReferenceRMSD.append(newlistH1[i][1])
					ScoreH1['Reference'].extend([[median(ScoreH1ReferenceRMSD), H1_ID, H1Length, newlistH1[0][2], root]] * listlength)
				H1Score = []
				for i in range(listlength):
					if newlistH1[i][1] > 12.5: newlistH1[i][1] = float('NaN')
					H1Score.append([newlistH1[i][1], H1_ID, H1Length, newlistH1[i][2], root])
				H1Error = max([row[0] for row in H1Score]) -  min([row[0] for row in H1Score])							

				H1Sample.setdefault(key, {})
				for i in range(0, len(newlistH1), 5):
					sampleSize = i+1
					H1Sample[key].setdefault(sampleSize, [])
					random.shuffle(newlistH1)
					test = newlistH1[:sampleSize]
					test.sort(key=operator.itemgetter(0))
					H1Sample[key][sampleSize].append(test[0][1])
				
				newlistH2.sort(key=operator.itemgetter(0))
				if key == 'Multi-Template':
					ScoreH2.setdefault('Reference', [])
					ScoreH2ReferenceRMSD = []
					for i in range(10):
						ScoreH2ReferenceRMSD.append(newlistH2[i][1])
					ScoreH2['Reference'].extend([[median(ScoreH2ReferenceRMSD), H2_ID, H2Length, newlistH2[0][2], root]] * listlength)
				H2Score = []
				for i in range(listlength):
					if newlistH2[i][1] > 12.5: newlistH2[i][1] = float('NaN')
					H2Score.append([newlistH2[i][1], H2_ID, H2Length, newlistH2[i][2], root])
				H2Error = max([row[0] for row in H2Score]) -  min([row[0] for row in H2Score])
				
				H2Sample.setdefault(key, {})
				for i in range(0, len(newlistH2), 5):
					sampleSize = i+1
					H2Sample[key].setdefault(sampleSize, [])
					random.shuffle(newlistH2)
					test = newlistH2[:sampleSize]
					test.sort(key=operator.itemgetter(0))
					H2Sample[key][sampleSize].append(test[0][1])
				
				newlistH3.sort(key=operator.itemgetter(0))
				with open(root + '/' + key + 'H3List', 'w') as f:
					for item in newlistH3:
						f.write('%s\n' % item)
				if key == 'Multi-Template':
					ScoreH3.setdefault('Reference', [])
					ScoreH3ReferenceRMSD = []
					for i in range(10):
						ScoreH3ReferenceRMSD.append(newlistH3[i][1])
					ScoreH3['Reference'].extend([[median(ScoreH3ReferenceRMSD), H3_ID, H3Length, newlistH3[0][2], root]] * listlength)
				H3Score = []
				for i in range(listlength):
					if newlistH3[i][1] > 12.5: newlistH3[i][1] = float('NaN')
					H3Score.append([newlistH3[i][1], H3_ID, H3Length, newlistH3[i][2], root])
				H3Error = max([row[0] for row in H3Score]) -  min([row[0] for row in H3Score])
				
				H3Sample.setdefault(key, {})
				for i in range(0, len(newlistH3), 5):
					sampleSize = i+1
					H3Sample[key].setdefault(sampleSize, [])
					random.shuffle(newlistH3)
					test = newlistH3[:sampleSize]
					test.sort(key=operator.itemgetter(0))
					H3Sample[key][sampleSize].append(test[0][1])
				
				newlistL1.sort(key=operator.itemgetter(0))
				if key == 'Multi-Template':
					ScoreL1.setdefault('Reference', [])
					ScoreL1ReferenceRMSD = []
					for i in range(10):
						ScoreL1ReferenceRMSD.append(newlistL1[i][1])
					ScoreL1['Reference'].extend([[median(ScoreL1ReferenceRMSD), L1_ID, L1Length, newlistL1[0][2], root]] * listlength)
				L1Score = []
				for i in range(listlength):
					if newlistL1[i][1] > 12.5: newlistL1[i][1] = float('NaN')
					L1Score.append([newlistL1[i][1], L1_ID, L1Length, newlistL1[i][2], root])
				L1Error = max([row[0] for row in L1Score]) -  min([row[0] for row in L1Score])
				
				L1Sample.setdefault(key, {})
				for i in range(0, len(newlistL1), 5):
					sampleSize = i+1
					L1Sample[key].setdefault(sampleSize, [])
					random.shuffle(newlistL1)
					test = newlistL1[:sampleSize]
					test.sort(key=operator.itemgetter(0))
					L1Sample[key][sampleSize].append(test[0][1])
				
				newlistL2.sort(key=operator.itemgetter(0))
				if key == 'Multi-Template':
					ScoreL2.setdefault('Reference', [])
					ScoreL2ReferenceRMSD = []
					for i in range(10):
						ScoreL2ReferenceRMSD.append(newlistL2[i][1])
					ScoreL2['Reference'].extend([[median(ScoreL2ReferenceRMSD), L2_ID, L2Length, newlistL2[0][2], root]] * listlength)
				L2Score = []
				for i in range(listlength):
					if newlistL2[i][1] > 12.5: newlistL2[i][1] = float('NaN')
					L2Score.append([newlistL2[i][1], L2_ID, L2Length, newlistL2[i][2], root])
				L2Error = max([row[0] for row in L2Score]) -  min([row[0] for row in L2Score])
				
				L2Sample.setdefault(key, {})
				for i in range(0, len(newlistL2), 5):
					sampleSize = i+1
					L2Sample[key].setdefault(sampleSize, [])
					random.shuffle(newlistL2)
					test = newlistL2[:sampleSize]
					test.sort(key=operator.itemgetter(0))
					L2Sample[key][sampleSize].append(test[0][1])
				
				newlistL3.sort(key=operator.itemgetter(0))
				if key == 'Multi-Template':
					ScoreL3.setdefault('Reference', [])
					ScoreL3ReferenceRMSD = []
					for i in range(10):
						ScoreL3ReferenceRMSD.append(newlistL3[i][1])
					ScoreL3['Reference'].extend([[median(ScoreL3ReferenceRMSD), L3_ID, L3Length, newlistL3[0][2], root]] * listlength)
				L3Score = []
				for i in range(listlength):
					if newlistL3[i][1] > 12.5: newlistL3[i][1] = float('NaN')
					L3Score.append([newlistL3[i][1], L3_ID, L3Length, newlistL3[i][2], root])
				L3Error = max([row[0] for row in L3Score]) -  min([row[0] for row in L3Score])
				
				L3Sample.setdefault(key, {})
				for i in range(0, len(newlistL3), 5):
					sampleSize = i+1
					L3Sample[key].setdefault(sampleSize, [])
					random.shuffle(newlistL3)
					test = newlistL3[:sampleSize]
					test.sort(key=operator.itemgetter(0))
					L3Sample[key][sampleSize].append(test[0][1])

				newlistOverall.sort(key=operator.itemgetter(0))
				if key == 'Multi-Template':
        	                        ScoreResults.setdefault('Reference', [])
					ScoreResultsReferenceRMSD = []
					for i in range(10):
						ScoreResultsReferenceRMSD.append(newlistOverall[i][1])
					ScoreResults['Reference'].extend([[median(ScoreResultsReferenceRMSD), percent_identity, BBlength, newlistOverall[0][2], root]] * listlength)
				OverallScore = []
				for i in range(listlength):
					OverallScore.append([newlistOverall[i][1], percent_identity, BBlength, newlistOverall[i][2], root])
				OverallError = max([row[0] for row in OverallScore]) -  min([row[0] for row in OverallScore])


				BackboneSample.setdefault(key, {})
				for i in range(0, len(newlistOverall), 5):
					sampleSize = i+1
					BackboneSample[key].setdefault(sampleSize, [])
					random.shuffle(newlistOverall)
					test = newlistOverall[:sampleSize]
					test.sort(key=operator.itemgetter(0))
					BackboneSample[key][sampleSize].append(test[0][1])



				if key == 'Control' or key == 'Single Template':
					OverallScore = OverallScore * 5
					H1Score = H1Score * 5
					H2Score = H2Score * 5
					H3Score = H3Score * 5
					L1Score = L1Score * 5
					L2Score = L2Score * 5
					L3Score = L3Score * 5
				ScoreResults.setdefault(key, [])
				ScoreResults[key].extend(OverallScore)
				ScoreH1.setdefault(key, [])
				ScoreH1[key].extend(H1Score)
				ScoreH2.setdefault(key, [])
				ScoreH2[key].extend(H2Score)
				ScoreH3.setdefault(key, [])
				ScoreH3[key].extend(H3Score)
				ScoreL1.setdefault(key, [])
				ScoreL1[key].extend(L1Score)
				ScoreL2.setdefault(key, [])
				ScoreL2[key].extend(L2Score)
				ScoreL3.setdefault(key, [])
				ScoreL3[key].extend(L3Score)
				'''
				if key == 'RosettaAntibody':
					for i in range(10):
						if H3Score[(i+1) * -1][0] < 0:
							print(H3Score[(i+1) * -1][0])
							print(H3Score[(i+1) * -1][-2])'''

				Error.setdefault(key, [])
				Error[key].append([root, H1Error, H2Error, H3Error, L1Error, L2Error, L3Error, OverallError])
				#print(key)
				#print(len(ScoreH1[key]))


				#if key == 'AbPredict':
					#print(ScoreL3['AbPredict'][-1])
					#print(CDRs['L3']
			
			ScoreH1.setdefault('RosettaAntibody', [])
			if(len(ScoreH1['RosettaAntibody']) < len(ScoreH1['Multi-Template'])):
				while(len(ScoreH1['RosettaAntibody']) < len(ScoreH1['Multi-Template'])):
					ScoreResults.setdefault('RosettaAntibody', [])
					ScoreResults['RosettaAntibody'].extend([[float('NaN'), float('NaN'), BBlength, float('NaN'), float('NaN')]])
					ScoreH1.setdefault('RosettaAntibody', [])
					ScoreH1['RosettaAntibody'].extend([[float('NaN'), float('NaN'), H1Length, float('NaN'), float('NaN')]])
					ScoreH2.setdefault('RosettaAntibody', [])
					ScoreH2['RosettaAntibody'].extend([[float('NaN'), float('NaN'), H2Length, float('NaN'), float('NaN')]])
					ScoreH3.setdefault('RosettaAntibody', [])
					ScoreH3['RosettaAntibody'].extend([[float('NaN'), float('NaN'), H3Length, float('NaN'), float('NaN')]])
					ScoreL1.setdefault('RosettaAntibody', [])
					ScoreL1['RosettaAntibody'].extend([[float('NaN'), float('NaN'), L1Length, float('NaN'), float('NaN')]])
					ScoreL2.setdefault('RosettaAntibody', [])
					ScoreL2['RosettaAntibody'].extend([[float('NaN'), float('NaN'), L2Length, float('NaN'), float('NaN')]])
					ScoreL3.setdefault('RosettaAntibody', [])
					ScoreL3['RosettaAntibody'].extend([[float('NaN'), float('NaN'), L3Length, float('NaN'), float('NaN')]])
					print('RosettaAntibody Blank Added')
			else:
				H3difference = median([row[0] for row in ScoreH3['RosettaAntibody'][-10:]]) - median([row[0] for row in ScoreH3['Multi-Template'][-10:]])
	
				if H3difference > 1:
					Outliers.setdefault('BadRABH3', [])
					Outliers['BadRABH3'].extend([[root, H3difference, [row[3] for row in ScoreH3['RosettaAntibody'][-10:]], [row[3] for row in ScoreH3['Multi-Template'][-10:]]]])	
	
				if H3difference < -1:
					Outliers.setdefault('BadRCMRABH3', [])
					Outliers['BadRCMRABH3'].extend([[root, H3difference, [row[3] for row in ScoreH3['RosettaAntibody'][-10:]], [row[3] for row in ScoreH3['Multi-Template'][-10:]]]])		

				L3difference = median([row[0] for row in ScoreL3['RosettaAntibody'][-10:]]) - median([row[0] for row in ScoreL3['Multi-Template'][-10:]])		

				if L3difference > 1:
					Outliers.setdefault('BadRABL3', [])
					Outliers['BadRABL3'].extend([[root, L3difference, [row[3] for row in ScoreL3['RosettaAntibody'][-10:]], [row[3] for row in ScoreL3['Multi-Template'][-10:]]]])

				if L3difference < -1:
					Outliers.setdefault('BadRCMRABL3', [])
					Outliers['BadRCMRABL3'].extend([[root, L3difference, [row[3] for row in ScoreL3['RosettaAntibody'][-10:]], [row[3] for row in ScoreL3['Multi-Template'][-10:]]]])
			
				BBdifference = median([row[0] for row in ScoreResults['RosettaAntibody'][-10:]]) - median([row[0] for row in ScoreResults['Multi-Template'][-10:]])

				if BBdifference > 0.5:
					Outliers.setdefault('BadRABBB', [])
					Outliers['BadRABBB'].extend([[root, BBdifference, [row[3] for row in ScoreResults['RosettaAntibody'][-10:]], [row[3] for row in ScoreResults['Multi-Template'][-10:]]]])

				if BBdifference < -0.5:
					Outliers.setdefault('BadRCMRABBB', [])
					Outliers['BadRCMRABBB'].extend([[root, BBdifference, [row[3] for row in ScoreResults['RosettaAntibody'][-10:]], [row[3] for row in ScoreResults['Multi-Template'][-10:]]]])

				if(len(ScoreH1['AbPredict']) == len(ScoreH1['Multi-Template'])):
					H3difference = median([row[0] for row in ScoreH3['AbPredict'][-10:]]) - median([row[0] for row in ScoreH3['RosettaAntibody'][-10:]])
	
					if H3difference > 1:
						Outliers.setdefault('BadABPRABH3', [])
						Outliers['BadABPRABH3'].extend([[root, H3difference, [row[3] for row in ScoreH3['AbPredict'][-10:]], [row[3] for row in ScoreH3['RosettaAntibody'][-10:]]]])	
	
					if H3difference < -1:
						Outliers.setdefault('BadRABABPH3', [])
						Outliers['BadRABABPH3'].extend([[root, H3difference, [row[3] for row in ScoreH3['AbPredict'][-10:]], [row[3] for row in ScoreH3['RosettaAntibody'][-10:]]]])		


			if(len(ScoreH1['AbPredict']) < len(ScoreH1['Multi-Template'])):
				while(len(ScoreH1['AbPredict']) < len(ScoreH1['Multi-Template'])):
					ScoreResults.setdefault('AbPredict', [])
					ScoreResults['AbPredict'].extend([[float('NaN'), float('NaN'), BBlength, float('NaN'), float('NaN')]])
					ScoreH1.setdefault('AbPredict', [])
					ScoreH1['AbPredict'].extend([[float('NaN'), float('NaN'), H1Length, float('NaN'), float('NaN')]])
					ScoreH2.setdefault('AbPredict', [])
					ScoreH2['AbPredict'].extend([[float('NaN'), float('NaN'), H2Length, float('NaN'), float('NaN')]])
					ScoreH3.setdefault('AbPredict', [])
					ScoreH3['AbPredict'].extend([[float('NaN'), float('NaN'), H3Length, float('NaN'), float('NaN')]])
					ScoreL1.setdefault('AbPredict', [])
					ScoreL1['AbPredict'].extend([[float('NaN'), float('NaN'), L1Length, float('NaN'), float('NaN')]])
					ScoreL2.setdefault('AbPredict', [])
					ScoreL2['AbPredict'].extend([[float('NaN'), float('NaN'), L2Length, float('NaN'), float('NaN')]])
					ScoreL3.setdefault('AbPredict', [])
					ScoreL3['AbPredict'].extend([[float('NaN'), float('NaN'), L3Length, float('NaN'), float('NaN')]])
					print('AbPredict Blank Added')
			else:
				H3difference = median([row[0] for row in ScoreH3['AbPredict'][-10:]]) - median([row[0] for row in ScoreH3['Multi-Template'][-10:]])
	
				if H3difference > 1:
					Outliers.setdefault('BadABPH3', [])
					Outliers['BadABPH3'].extend([[root, H3difference, [row[3] for row in ScoreH3['AbPredict'][-10:]], [row[3] for row in ScoreH3['Multi-Template'][-10:]]]])	
	
				if H3difference < -1:
					Outliers.setdefault('BadRCMABPH3', [])
					Outliers['BadRCMABPH3'].extend([[root, H3difference, [row[3] for row in ScoreH3['AbPredict'][-10:]], [row[3] for row in ScoreH3['Multi-Template'][-10:]]]])		

				L3difference = median([row[0] for row in ScoreL3['AbPredict'][-10:]]) - median([row[0] for row in ScoreL3['Multi-Template'][-10:]])		

				if L3difference > 1:
					Outliers.setdefault('BadABPL3', [])
					Outliers['BadABPL3'].extend([[root, L3difference, [row[3] for row in ScoreL3['AbPredict'][-10:]], [row[3] for row in ScoreL3['Multi-Template'][-10:]]]])

				if L3difference < -1:
					Outliers.setdefault('BadRCMABPL3', [])
					Outliers['BadRCMABPL3'].extend([[root, L3difference, [row[3] for row in ScoreL3['AbPredict'][-10:]], [row[3] for row in ScoreL3['Multi-Template'][-10:]]]])
			
				BBdifference = median([row[0] for row in ScoreResults['AbPredict'][-10:]]) - median([row[0] for row in ScoreResults['Multi-Template'][-10:]])

				if BBdifference > 0.5:
					Outliers.setdefault('BadABPBB', [])
					Outliers['BadABPBB'].extend([[root, BBdifference, [row[3] for row in ScoreResults['AbPredict'][-10:]], [row[3] for row in ScoreResults['Multi-Template'][-10:]]]])

				if BBdifference < -0.5:
					Outliers.setdefault('BadRCMABPBB', [])
					Outliers['BadRCMABPBB'].extend([[root, BBdifference, [row[3] for row in ScoreResults['AbPredict'][-10:]], [row[3] for row in ScoreResults['Multi-Template'][-10:]]]])

				
			
			
			




			
			H3difference = median([row[0] for row in ScoreH3['RosettaAntibody Inspired'][-10:]]) - median([row[0] for row in ScoreH3['Multi-Template'][-10:]])

			if H3difference > 1:
				Outliers.setdefault('BadRABLH3', [])
				Outliers['BadRABLH3'].extend([[root, H3difference, [row[3] for row in ScoreH3['RosettaAntibody Inspired'][-10:]], [row[3] for row in ScoreH3['Multi-Template'][-10:]]]])

			if H3difference < -1:
				Outliers.setdefault('BadRCMLH3', [])
				Outliers['BadRCMLH3'].extend([[root, H3difference, [row[3] for row in ScoreH3['RosettaAntibody Inspired'][-10:]], [row[3] for row in ScoreH3['Multi-Template'][-10:]]]])

			L3difference = median([row[0] for row in ScoreL3['RosettaAntibody'][-10:]]) - median([row[0] for row in ScoreL3['Multi-Template'][-10:]])

			if L3difference > 1:
				Outliers.setdefault('BadRABLL3', [])
				Outliers['BadRABLL3'].extend([[root, L3difference, [row[3] for row in ScoreL3['RosettaAntibody Inspired'][-10:]], [row[3] for row in ScoreL3['Multi-Template'][-10:]]]])

			if L3difference < -1:
				Outliers.setdefault('BadRCMLL3', [])
				Outliers['BadRCMLL3'].extend([[root, L3difference, [row[3] for row in ScoreL3['RosettaAntibody Inspired'][-10:]], [row[3] for row in ScoreL3['Multi-Template'][-10:]]]])
			
			BBdifference = median([row[0] for row in ScoreResults['RosettaAntibody'][-10:]]) - median([row[0] for row in ScoreResults['Multi-Template'][-10:]])

			if BBdifference > 0.5:
				Outliers.setdefault('BadRABLBB', [])
				Outliers['BadRABLBB'].extend([[root, BBdifference, [row[3] for row in ScoreResults['RosettaAntibody Inspired'][-10:]], [row[3] for row in ScoreResults['Multi-Template'][-10:]]]])

			if BBdifference < -0.5:
				Outliers.setdefault('BadRCMLBB', [])
				Outliers['BadRCMLBB'].extend([[root, BBdifference, [row[3] for row in ScoreResults['RosettaAntibody Inspired'][-10:]], [row[3] for row in ScoreResults['Multi-Template'][-10:]]]])

			Outliers.setdefault('All', [])
			Outliers['All'].extend([[root, [row[3] for row in ScoreL3['RosettaAntibody Inspired'][-10:]], [row[3] for row in ScoreL3['Multi-Template'][-10:]], [row[3] for row in ScoreL3['RosettaAntibody'][-10:]]]])
			
			print(len(ScoreH1['Reference']))
			print(len(ScoreH1['Multi-Template']))
			print(len(ScoreH1['RosettaAntibody']))
	'''int(round(ScoreH3[j][i][2] * 2) / 2)'''
		
	records = list()
	for i in range(len(ScoreResults['Multi-Template'])):
		for j in sorted(ScoreResults.keys(), key = lambda x: x[0], reverse=True):
			print(i)
			print(j)
			records.append({ "Experiment" : j, "Multi-Template BB RMSD vs Native" : ScoreResults['Reference'][i][0], "Method BB RMSD vs Native" : ScoreResults[j][i][0], "Percent Identity" : ((ScoreResults[j][i][1] * 25) // 5) / 5, "Backbone Length" : round(ScoreResults[j][i][2] * 5) / 5
				, "Multi-Template HCDR1 RMSD vs Native" : ScoreH1['Reference'][i][0], "Method HCDR1 RMSD vs Native" : ScoreH1[j][i][0], "HCDR1 Percent Identity" : ((ScoreH1[j][i][1] * 25) // 5) / 5, "HCDR1 Length" : ScoreH1[j][i][2]
				, "Multi-Template HCDR2 RMSD vs Native" : ScoreH2['Reference'][i][0], "Method HCDR2 RMSD vs Native" : ScoreH2[j][i][0], "HCDR2 Percent Identity" : ((ScoreH2[j][i][1] * 25) // 5) / 5, "HCDR2 Length" : ScoreH2[j][i][2]
				, "Multi-Template HCDR3 RMSD vs Native" : ScoreH3['Reference'][i][0], "Method HCDR3 RMSD vs Native" : ScoreH3[j][i][0], "HCDR3 Percent Identity" : ((ScoreH3[j][i][1] * 25) // 5) / 5, "HCDR3 Length" : bins[(ScoreH3[j][i][2] // 7)], "HCDR3 Length Unbinned" :ScoreH3[j][i][2]
				, "Multi-Template LCDR1 RMSD vs Native" : ScoreL1['Reference'][i][0], "Method LCDR1 RMSD vs Native" : ScoreL1[j][i][0], "LCDR1 Percent Identity" : ((ScoreL1[j][i][1] * 25) // 5) / 5, "LCDR1 Length" : ScoreL1[j][i][2]
				, "Multi-Template LCDR2 RMSD vs Native" : ScoreL2['Reference'][i][0], "Method LCDR2 RMSD vs Native" : ScoreL2[j][i][0], "LCDR2 Percent Identity" : ((ScoreL2[j][i][1] * 25) // 5) / 5, "LCDR2 Length" : ScoreL2[j][i][2]
				, "Multi-Template LCDR3 RMSD vs Native" : ScoreL3['Reference'][i][0], "Method LCDR3 RMSD vs Native" : ScoreL3[j][i][0], "LCDR3 Percent Identity" : ((ScoreL3[j][i][1] * 25) // 5) / 5, "LCDR3 Length" : ScoreL3[j][i][2], "Root" : ScoreH3['Reference'][i][-1], "HCDR3 Percent Identity Unbinned" : ScoreH3[j][i][1], "Percent Identity Unbinned" : ScoreResults[j][i][1], "File" : ScoreResults[j][i][3] })
	dframe = DataFrame(records)
	
	records = list()
	for i in BackboneSample.keys():
		for j in BackboneSample[i].keys():
			records.append({ "Experiment" : i, "Sample Size" : j, "Backbone RMSD" : np.median(BackboneSample[i][j])
				, "HCDR1 RMSD" :  np.median(H1Sample[i][j])
				, "HCDR2 RMSD" :  np.median(H2Sample[i][j])
				, "HCDR3 RMSD" :  np.median(H3Sample[i][j])
				, "LCDR1 RMSD" :  np.median(L1Sample[i][j])
				, "LCDR2 RMSD" :  np.median(L2Sample[i][j])
				, "LCDR3 RMSD" :  np.median(L3Sample[i][j])})
	samplesframe = DataFrame(records)

	errorframe = Error	
	dframe.to_csv('RMSD_data.csv')
	pickle.dump(dframe, open('./RMSD_Data.p', 'wb'))
	pickle.dump(samplesframe, open('./Sample_Data.p', 'wb'))
	pickle.dump(errorframe, open('./Error_Data.p', 'wb'))
	pickle.dump(Outliers, open('./Outliers.p', 'wb'))
	
else:
	dframe = pickle.load(open('./RMSD_Data.p', 'rb'))
	samplesframe = pickle.load(open('./Sample_Data.p', 'rb'))
	errorframe = pickle.load(open('./Error_Data.p', 'rb'))
	Outliers = pickle.load(open('./Outliers.p', 'rb'))

#print(dframe)

#plt.rcParams["font.family"] = "sans-serif"
#plt.rcParams["font.sans-serif"] = "Arial"

BB_lim = 3

g = sns.lmplot(data=dframe, x="Multi-Template BB RMSD vs Native", y="Method BB RMSD vs Native", hue="Experiment", fit_reg=False, x_estimator=np.median, x_ci=None, height = 10, aspect = 1, scatter_kws={'alpha':1, 's':75}, palette=['grey'], legend_out=False)
g.set(ylim=(0, 3), xlim=(0, 3))
legend = plt.legend(loc = 'upper left', frameon=True, framealpha=1, title = 'Method')
frame = legend.get_frame()
frame.set_facecolor('white')
frame.set_edgecolor('black')
for lh in legend.legendHandles:
	lh.set_alpha(1)
	lh.set_sizes([100])
ax = plt.gca()
#ax.set_title('Method Performance of Best Models')
plt.tight_layout()
plt.savefig(os.getcwd() + '/Figures/Backbone_by_Score.png', dpi=600, format='png', transparent=True)
g = sns.lmplot(data=dframe.loc[dframe['Experiment'].isin(['RosettaAntibody Inspired'])], x="Multi-Template BB RMSD vs Native", y="Method BB RMSD vs Native", hue="Experiment", fit_reg=False, height = 10, aspect = 1, x_estimator=np.median, x_ci=None, scatter_kws={'alpha':1, 's':75}, palette=['grey'], legend=False, legend_out=False)
g.set(ylim=(0, 3), xlim=(0, 3))
ax = plt.gca()
X_plot = np.linspace(0, 3, 400)
Y_plot = X_plot
plt.plot(X_plot, Y_plot, color='silver')
plt.plot(X_plot, Y_plot - 0.5, color='silver')
plt.plot(X_plot, Y_plot + 0.5, color='silver')
plt.ylabel('RosettaAntibody Inspired BB RMSD vs Native')
#ax.set_title('Performance of Best RosettaAntibody Inspired Models')
plt.tight_layout()
plt.savefig(os.getcwd() + '/Figures/Backbone_by_Score_Rosetta.png', dpi=600, format='png', transparent=True)
g = sns.lmplot(data=dframe.loc[dframe['Experiment'].isin(['RosettaAntibody'])], x="Multi-Template BB RMSD vs Native", y="Method BB RMSD vs Native", hue="Experiment", fit_reg=False, height = 10, aspect = 1, x_estimator=np.median, x_ci=None, scatter_kws={'alpha':1, 's':75}, palette=['grey'], legend=False, legend_out=False)
g.set(ylim=(0, 3), xlim=(0, 3))
ax = plt.gca()
X_plot = np.linspace(0, 3, 400)
Y_plot = X_plot
plt.plot(X_plot, Y_plot, color='silver')
plt.plot(X_plot, Y_plot - 0.5, color='silver')
plt.plot(X_plot, Y_plot + 0.5, color='silver')
plt.xlabel('Multi-Template Backbone RMSD vs Native', labelpad=20)
plt.ylabel('RosettaAntibody Backbone RMSD vs Native', labelpad=20)
#ax.set_title('Performance of Best RosettaAntibody Models')
plt.tight_layout()
plt.savefig(os.getcwd() + '/Figures/Backbone_by_Score_RosettaAB.png', dpi=600, format='png', transparent=True)
g = sns.lmplot(data=dframe.loc[dframe['Experiment'].isin(['AbPredict'])], x="Multi-Template BB RMSD vs Native", y="Method BB RMSD vs Native", hue="Experiment", fit_reg=False, height = 10, aspect = 1, x_estimator=np.median, x_ci=None, scatter_kws={'alpha':1, 's':75}, palette=['grey'], legend=False, legend_out=False)
g.set(ylim=(0, 3), xlim=(0, 3))
ax = plt.gca()
X_plot = np.linspace(0, 3, 400)
Y_plot = X_plot
plt.plot(X_plot, Y_plot, color='silver')
plt.plot(X_plot, Y_plot - 0.5, color='silver')
plt.plot(X_plot, Y_plot + 0.5, color='silver')
plt.xlabel('Multi-Template Backbone RMSD vs Native', labelpad=20)
plt.ylabel('AbPredict Backbone RMSD vs Native', labelpad=20)
#ax.set_title('Performance of Best AbPredict Models')
plt.tight_layout()
plt.savefig(os.getcwd() + '/Figures/Backbone_by_Score_AbPredict.png', dpi=600, format='png', transparent=True)
g = sns.lmplot(data=dframe.loc[dframe['Experiment'].isin(['Single Template'])], x="Multi-Template BB RMSD vs Native", y="Method BB RMSD vs Native", hue="Experiment", fit_reg=False, height = 10, aspect = 1, x_estimator=np.median, x_ci=None, scatter_kws={'alpha':1, 's':75}, palette=['grey'], legend=False, legend_out=False)
g.set(ylim=(0, 3), xlim=(0, 3))
ax = plt.gca()
X_plot = np.linspace(0, 3, 400)
Y_plot = X_plot
plt.plot(X_plot, Y_plot, color='silver')
plt.plot(X_plot, Y_plot - 0.5, color='silver')
plt.plot(X_plot, Y_plot + 0.5, color='silver')
plt.xlabel('Multi-Template Backbone RMSD vs Native', labelpad=20)
plt.ylabel('Single Template Backbone RMSD vs Native', labelpad=20)
#ax.set_title('Performance of Best Single Template Models')
plt.tight_layout()
plt.savefig(os.getcwd() + '/Figures/Backbone_by_Score_Single.png', dpi=600, format='png', transparent=True)
g = sns.lmplot(data=dframe.loc[dframe['Experiment'].isin(['Control'])], x="Multi-Template BB RMSD vs Native", y="Method BB RMSD vs Native", hue="Experiment", fit_reg=False, height = 10, aspect = 1, x_estimator=np.median, x_ci=None, scatter_kws={'alpha':1, 's':75}, palette=['grey'], legend=False, legend_out=False)
g.set(ylim=(0, 3), xlim=(0, 3))
ax = plt.gca()
X_plot = np.linspace(0, 3, 400)
Y_plot = X_plot
plt.plot(X_plot, Y_plot, color='silver')
plt.plot(X_plot, Y_plot - 0.5, color='silver')
plt.plot(X_plot, Y_plot + 0.5, color='silver')
plt.ylabel('Control Backbone RMSD vs Native', labelpad=20)
plt.xlabel('Multi-Template Backbone RMSD vs Native', labelpad=20)
#ax.set_title('Performance of Best Control Models')
plt.tight_layout()
plt.savefig(os.getcwd() + '/Figures/Backbone_by_Score_Control.png', dpi=600, format='png', transparent=True)


plt.clf()
fig, ax = plt.subplots(figsize=(10, 10))
g = sns.boxplot(data=dframe.loc[dframe['Experiment'].isin(['RosettaAntibody Inspired', 'Single Template', 'Multi-Template', 'Control', 'RosettaAntibody', 'AbPredict'])], x="Percent Identity", y="Method BB RMSD vs Native", hue="Experiment", ax=ax)
g.set(ylim=(0, BB_lim))
legend = plt.legend(loc = 'upper left', frameon=True, framealpha=1, title = 'Method') 
plt.savefig(os.getcwd() + '/Figures/Backbone_by_Percent.png', dpi=600, format='png', transparent=True)

plt.clf()
fig, ax = plt.subplots(figsize=(10, 10))
g = sns.boxplot(data=dframe.loc[dframe['Experiment'].isin(['RosettaAntibody Inspired', 'Single Template', 'Multi-Template', 'Control', 'RosettaAntibody', 'AbPredict'])], x="Percent Identity Unbinned", y="Method BB RMSD vs Native", hue="Experiment", ax=ax)
g.set(ylim=(0, BB_lim))
legend = plt.legend(loc = 'upper left', frameon=True, framealpha=1, title = 'Method') 
plt.savefig(os.getcwd() + '/Figures/Backbone_by_Percent_Unbinned.png', dpi=600, format='png', transparent=True)

plt.clf()
fig, ax = plt.subplots(figsize=(10, 10))
g = sns.boxplot(data=dframe.loc[dframe['Experiment'].isin(['RosettaAntibody Inspired', 'Single Template', 'Multi-Template', 'Control', 'RosettaAntibody', 'AbPredict'])], x="Backbone Length", y="Method BB RMSD vs Native", hue="Experiment", ax=ax)
g.set(ylim=(0, BB_lim))
legend = plt.legend(loc = 'upper left', frameon=True, framealpha=1, title = 'Method')
plt.tight_layout()
plt.savefig(os.getcwd() + '/Figures/Backbone_by_Length.png', dpi=600, format='png', transparent=True)

grouped_dframeBB = dframe.groupby(['Experiment', 'Root'], as_index=False)['Method BB RMSD vs Native'].median()
grouped_dframeBB['Method BB RMSD vs Native'] = (np.floor(grouped_dframeBB['Method BB RMSD vs Native'] / 0.5) * 0.5)
grouped_dframeBB['Method BB RMSD vs Native'].values[grouped_dframeBB['Method BB RMSD vs Native'].values > 5.0] = float('NaN')
grouped_dframeBB['Method BB RMSD vs Native Binned'] = grouped_dframeBB['Method BB RMSD vs Native'].astype(str) + ' -\n' + (grouped_dframeBB['Method BB RMSD vs Native'] + 0.5).astype(str) + '  '
grouped_dframeBB.loc[grouped_dframeBB['Method BB RMSD vs Native Binned'].str.contains('nan'), 'Method BB RMSD vs Native Binned'] = 'Failed to\nModel'
grouped_dframeBB_2 = grouped_dframeBB.groupby(['Experiment', 'Method BB RMSD vs Native Binned']).size().to_frame('Count')
grouped_dframeBB_2 = grouped_dframeBB_2.reset_index()
grouped_dframeBB_2.sort_values(by=['Method BB RMSD vs Native Binned'], inplace=True)
#print(grouped_dframeBB_2)

matplotlib.rcParams["patch.force_edgecolor"] = True

fig, ax = plt.subplots(figsize=(10, 10))
g = sns.barplot(data=grouped_dframeBB_2.loc[grouped_dframeBB_2['Experiment'].isin(['Single Template', 'Multi-Template', 'Control', 'RosettaAntibody', 'AbPredict'])], x="Method BB RMSD vs Native Binned", y="Count", hue="Experiment",hue_order=['Control', 'Single Template', 'Multi-Template', 'RosettaAntibody', 'AbPredict'],linewidth=1.5,edgecolor="black")
g.set_xticklabels(g.get_xticklabels(), rotation=45, horizontalalignment='center')
legend = plt.legend(loc = 'upper right', frameon=True, framealpha=1, title = 'Method')
plt.xlabel('Backbone RMSD vs Native') 
plt.ylabel('Count', labelpad=20)
plt.ylim(0, 90)
plt.tight_layout()
plt.savefig(os.getcwd() + '/Figures/Backbone_by_Bin.png', dpi=600, format='png', transparent=True, bbox_inches = "tight")

grouped_dframeBB = dframe.groupby(['Experiment', 'Root'], as_index=False)['Method BB RMSD vs Native'].median()

fig, ax = plt.subplots(figsize=(10, 10))
g = sns.violinplot(data=grouped_dframeBB.loc[grouped_dframeBB['Experiment'].isin(['Single Template', 'Multi-Template', 'Control', 'RosettaAntibody', 'AbPredict'])], x="Experiment", y="Method BB RMSD vs Native", order=['Control', 'Single Template', 'Multi-Template', 'RosettaAntibody', 'AbPredict'], bw=0.15)
g.set(ylim=(0, BB_lim))
g.set_xticklabels(g.get_xticklabels(), rotation=90, horizontalalignment='center')
plt.xlabel('')
plt.ylabel('Backbone RMSD vs Native', labelpad=20)
#sns.swarmplot(data=grouped_dframeBB.loc[grouped_dframeBB['Experiment'].isin(['Single Template', 'Multi-Template', 'Control', 'RosettaAntibody', 'AbPredict'])], x="Experiment", y="Method BB RMSD vs Native", color="white", order=['Control', 'Single Template', 'Multi-Template', 'RosettaAntibody', 'AbPredict'])
#plt.tight_layout()
plt.savefig(os.getcwd() + '/Figures/Backbone_by_RMSD_Violin.png', dpi=600, format='png', transparent=True, bbox_inches = "tight")

RABABPBBdict = {}
for index, row in grouped_dframeBB.iterrows():
	RABABPBBdict.setdefault(row['Root'], [[],[]])
	if row['Experiment'] == 'RosettaAntibody':
		RABABPBBdict[row['Root']][0] = row['Method BB RMSD vs Native']
	elif row['Experiment'] == 'AbPredict':
		RABABPBBdict[row['Root']][1] = row['Method BB RMSD vs Native']
RABABPBBlist = list()
for key, value in RABABPBBdict.items():
	RABABPBBlist.append({ "Root" : key, "RosettaAntibody BB RMSD": value[0], "AbPredict BB RMSD": value[1]})
RABABPBB = DataFrame(RABABPBBlist)

plt.clf()
g = sns.lmplot(data=RABABPBB, x="RosettaAntibody BB RMSD", y="AbPredict BB RMSD", fit_reg=False, height = 10, aspect = 1, x_ci=None, scatter_kws={'alpha':1, 's':75}, palette=['grey'], legend=False, legend_out=False)
g.set(ylim=(0, 3), xlim=(0, 3))
ax = plt.gca()
X_plot = np.linspace(0, 3, 400)
Y_plot = X_plot
plt.plot(X_plot, Y_plot, color='silver')
plt.plot(X_plot, Y_plot - 0.5, color='silver')
plt.plot(X_plot, Y_plot + 0.5, color='silver')
plt.ylabel('AbPredict Backbone RMSD vs Native', labelpad=20)
plt.xlabel('RosettaAntibody Backbone RMSD vs Native', labelpad=20)
#ax.set_title('Performance of Best Control Models')
plt.tight_layout()
plt.savefig(os.getcwd() + '/Figures/Backbone_by_RAB_ABP.png', dpi=600, format='png', transparent=True)

H1_lim = 4

plt.clf()
fig, ax = plt.subplots(figsize=(10, 10))
g = sns.boxplot(data=dframe.loc[dframe['Experiment'].isin(['RosettaAntibody Inspired', 'Single Template', 'Multi-Template', 'Control', 'RosettaAntibody', 'AbPredict'])], x="HCDR1 Percent Identity", y="Method HCDR1 RMSD vs Native", hue="Experiment", ax=ax)
g.set(ylim=(0, H1_lim)) 
legend = plt.legend(loc = 'upper left', frameon=True, framealpha=1, title = 'Method') 
plt.tight_layout()
plt.savefig(os.getcwd() + '/Figures/H1_by_Percent.png', dpi=600, format='png', transparent=True)

plt.clf()
fig, ax = plt.subplots(figsize=(10, 10))
g = sns.boxplot(data=dframe.loc[dframe['Experiment'].isin(['RosettaAntibody Inspired', 'Single Template', 'Multi-Template', 'Control', 'RosettaAntibody', 'AbPredict'])], x="HCDR1 Length", y="Method HCDR1 RMSD vs Native", hue="Experiment", ax=ax)
g.set(ylim=(0, H1_lim))
legend = plt.legend(loc = 'upper left', frameon=True, framealpha=1, title = 'Method') 
plt.tight_layout()
plt.savefig(os.getcwd() + '/Figures/H1_by_Length.png', dpi=600, format='png', transparent=True)

g = sns.lmplot(data=dframe, x="Multi-Template HCDR1 RMSD vs Native", y="Method HCDR1 RMSD vs Native", hue="Experiment", x_estimator=np.median, x_ci=None, fit_reg=False, height = 10, aspect = 1, scatter_kws={'alpha':1, 's':75}, palette=['grey'], legend_out=False)
g.set(ylim=(0, None), xlim=(0, None))
legend = plt.legend(loc = 'upper left', frameon=True, framealpha=1, title = 'Method')
frame = legend.get_frame()
frame.set_facecolor('white')
frame.set_edgecolor('black')
for lh in legend.legendHandles:
	lh.set_alpha(1)
	lh.set_sizes([100])
ax = plt.gca()
#ax.set_title('Performance of Best Models')
plt.tight_layout()
plt.savefig(os.getcwd() + '/Figures/H1_by_Score.png', dpi=600, format='png', transparent=True)
g = sns.lmplot(data=dframe.loc[dframe['Experiment'].isin(['RosettaAntibody Inspired'])], x="Multi-Template HCDR1 RMSD vs Native", y="Method HCDR1 RMSD vs Native", hue="Experiment", fit_reg=False, height = 10, aspect = 1, x_estimator=np.median, x_ci=None, scatter_kws={'alpha':1, 's':75}, palette=['grey'], legend=False, legend_out=False)
g.set(ylim=(0, H1_lim), xlim=(0, H1_lim))
ax = plt.gca()
X_plot = np.linspace(0, H1_lim, 400)
Y_plot = X_plot
plt.plot(X_plot, Y_plot, color='silver')
plt.plot(X_plot, Y_plot - 0.5, color='silver')
plt.plot(X_plot, Y_plot + 0.5, color='silver')
#ax.set_title('Performance of Best RosettaAntibody Inspired Models')
plt.ylabel('RosettaAntibody Inspired HCDR1 RMSD vs Native', labelpad=20)
plt.xlabel('Multi-Template HCDR1 RMSD vs Native', labelpad=20)
plt.tight_layout()
plt.savefig(os.getcwd() + '/Figures/H1_by_Score_Rosetta.png', dpi=600, format='png', transparent=True)
g = sns.lmplot(data=dframe.loc[dframe['Experiment'].isin(['RosettaAntibody'])], x="Multi-Template HCDR1 RMSD vs Native", y="Method HCDR1 RMSD vs Native", hue="Experiment", fit_reg=False, height = 10, aspect = 1, x_estimator=np.median, x_ci=None, scatter_kws={'alpha':1, 's':75}, palette=['grey'], legend=False, legend_out=False)
g.set(ylim=(0, H1_lim), xlim=(0, H1_lim))
ax = plt.gca()
X_plot = np.linspace(0, H1_lim, 400)
Y_plot = X_plot
plt.plot(X_plot, Y_plot, color='silver')
plt.plot(X_plot, Y_plot - 0.5, color='silver')
plt.plot(X_plot, Y_plot + 0.5, color='silver')
#ax.set_title('Performance of Best RosettaAntibody Models')
plt.ylabel('RosettaAntibody HCDR1 RMSD vs Native', labelpad=20)
plt.xlabel('Multi-Template HCDR1 RMSD vs Native', labelpad=20)
plt.tight_layout()
plt.savefig(os.getcwd() + '/Figures/H1_by_Score_RosettaAB.png', dpi=600, format='png', transparent=True)
g = sns.lmplot(data=dframe.loc[dframe['Experiment'].isin(['AbPredict'])], x="Multi-Template HCDR1 RMSD vs Native", y="Method HCDR1 RMSD vs Native", hue="Experiment", fit_reg=False, height = 10, aspect = 1, x_estimator=np.median, x_ci=None, scatter_kws={'alpha':1, 's':75}, palette=['grey'], legend=False, legend_out=False)
g.set(ylim=(0, H1_lim), xlim=(0, H1_lim))
ax = plt.gca()
X_plot = np.linspace(0, H1_lim, 400)
Y_plot = X_plot
plt.plot(X_plot, Y_plot, color='silver')
plt.plot(X_plot, Y_plot - 0.5, color='silver')
plt.plot(X_plot, Y_plot + 0.5, color='silver')
#ax.set_title('Performance of Best AbPredict Models')
plt.ylabel('AbPredict HCDR1 RMSD vs Native', labelpad=20)
plt.xlabel('Multi-Template HCDR1 RMSD vs Native', labelpad=20)
plt.tight_layout()
plt.savefig(os.getcwd() + '/Figures/H1_by_Score_AbPredict.png', dpi=600, format='png', transparent=True)
g = sns.lmplot(data=dframe.loc[dframe['Experiment'].isin(['Single Template'])], x="Multi-Template HCDR1 RMSD vs Native", y="Method HCDR1 RMSD vs Native", hue="Experiment", fit_reg=False, height = 10, aspect = 1, x_estimator=np.median, x_ci=None, scatter_kws={'alpha':1, 's':75}, palette=['grey'], legend=False, legend_out=False)
g.set(ylim=(0, H1_lim), xlim=(0, H1_lim))
ax = plt.gca()
X_plot = np.linspace(0, H1_lim, 400)
Y_plot = X_plot
plt.plot(X_plot, Y_plot, color='silver')
plt.plot(X_plot, Y_plot - 0.5, color='silver')
plt.plot(X_plot, Y_plot + 0.5, color='silver')
#ax.set_title('Performance of Best Single Template Models')
plt.ylabel('Single Template HCDR1 RMSD vs Native', labelpad=20)
plt.xlabel('Multi-Template HCDR1 RMSD vs Native', labelpad=20)
plt.tight_layout()
plt.savefig(os.getcwd() + '/Figures/H1_by_Score_Single.png', dpi=600, format='png', transparent=True)
g = sns.lmplot(data=dframe.loc[dframe['Experiment'].isin(['Control'])], x="Multi-Template HCDR1 RMSD vs Native", y="Method HCDR1 RMSD vs Native", hue="Experiment", fit_reg=False, height = 10, aspect = 1, x_estimator=np.median, x_ci=None, scatter_kws={'alpha':1, 's':75}, palette=['grey'], legend=False, legend_out=False)
g.set(ylim=(0, H1_lim), xlim=(0, H1_lim))
ax = plt.gca()
X_plot = np.linspace(0, H1_lim, 400)
Y_plot = X_plot
plt.plot(X_plot, Y_plot, color='silver')
plt.plot(X_plot, Y_plot - 0.5, color='silver')
plt.plot(X_plot, Y_plot + 0.5, color='silver')
#ax.set_title('Performance of Best Control Models')
plt.ylabel('Control HCDR1 RMSD vs Native', labelpad=20)
plt.xlabel('Multi-Template HCDR1 RMSD vs Native', labelpad=20)
plt.tight_layout()
plt.savefig(os.getcwd() + '/Figures/H1_by_Score_Control.png', dpi=600, format='png', transparent=True)



H2_lim = 4
Margin = 0.5

plt.clf()
fig, ax = plt.subplots(figsize=(10, 10))
g = sns.boxplot(data=dframe.loc[dframe['Experiment'].isin(['RosettaAntibody Inspired', 'Single Template', 'Multi-Template', 'Control', 'RosettaAntibody', 'AbPredict'])], x="HCDR2 Percent Identity", y="Method HCDR2 RMSD vs Native", hue="Experiment", ax=ax)
g.set(ylim=(0, H2_lim))
legend = plt.legend(loc = 'upper left', frameon=True, framealpha=1, title = 'Method')  
plt.tight_layout()
plt.savefig(os.getcwd() + '/Figures/H2_by_Percent.png', dpi=600, format='png', transparent=True)

plt.clf()
g = sns.boxplot(data=dframe.loc[dframe['Experiment'].isin(['RosettaAntibody Inspired', 'Single Template', 'Multi-Template', 'Control', 'RosettaAntibody', 'AbPredict'])], x="HCDR2 Length", y="Method HCDR2 RMSD vs Native", hue="Experiment", ax=ax)
g.set(ylim=(0, H2_lim))
legend = plt.legend(loc = 'upper left', frameon=True, framealpha=1, title = 'Method') 
plt.tight_layout()
plt.savefig(os.getcwd() + '/Figures/H2_by_Length.png', dpi=600, format='png', transparent=True)

g = sns.lmplot(data=dframe, x="Multi-Template HCDR2 RMSD vs Native", y="Method HCDR2 RMSD vs Native", hue="Experiment", x_estimator=np.median, x_ci=None, fit_reg=False, height = 10, aspect = 1, scatter_kws={'alpha':1, 's':75}, palette=['grey'], legend_out=False)
g.set(ylim=(0, None), xlim=(0, None))
legend = plt.legend(loc = 'upper left', frameon=True, framealpha=1, title = 'Method')
frame = legend.get_frame()
frame.set_facecolor('white')
frame.set_edgecolor('black')
for lh in legend.legendHandles:
	lh.set_alpha(1)
	lh.set_sizes([100])
ax = plt.gca()
#ax.set_title('Performance of Best Models')
plt.tight_layout()
plt.savefig(os.getcwd() + '/Figures/H2_by_Score.png', dpi=600, format='png', transparent=True)
g = sns.lmplot(data=dframe.loc[dframe['Experiment'].isin(['RosettaAntibody Inspired'])], x="Multi-Template HCDR2 RMSD vs Native", y="Method HCDR2 RMSD vs Native", hue="Experiment", fit_reg=False, height = 10, aspect = 1, x_estimator=np.median, x_ci=None, scatter_kws={'alpha':1, 's':75}, palette=['grey'], legend=False, legend_out=False)
g.set(ylim=(0, H2_lim), xlim=(0, H2_lim))
ax = plt.gca()
X_plot = np.linspace(0, H2_lim, 400)
Y_plot = X_plot
plt.plot(X_plot, Y_plot, color='silver')
plt.plot(X_plot, Y_plot - Margin, color='silver')
plt.plot(X_plot, Y_plot + Margin, color='silver')
#ax.set_title('Performance of Best RosettaAntibody Inspired Models')
plt.ylabel('RosettaAntibody Inspired HCDR2 RMSD vs Native', labelpad=20)
plt.xlabel('Multi-Template HCDR2 RMSD vs Native', labelpad=20)
plt.tight_layout()
plt.savefig(os.getcwd() + '/Figures/H2_by_Score_Rosetta.png', dpi=600, format='png', transparent=True)
g = sns.lmplot(data=dframe.loc[dframe['Experiment'].isin(['RosettaAntibody'])], x="Multi-Template HCDR2 RMSD vs Native", y="Method HCDR2 RMSD vs Native", hue="Experiment", fit_reg=False, height = 10, aspect = 1, x_estimator=np.median, x_ci=None, scatter_kws={'alpha':1, 's':75}, palette=['grey'], legend=False, legend_out=False)
g.set(ylim=(0, H2_lim), xlim=(0, H2_lim))
ax = plt.gca()
X_plot = np.linspace(0, H2_lim, 400)
Y_plot = X_plot
plt.plot(X_plot, Y_plot, color='silver')
plt.plot(X_plot, Y_plot - Margin, color='silver')
plt.plot(X_plot, Y_plot + Margin, color='silver')
#ax.set_title('Performance of Best RosettaAntibody Models')
plt.ylabel('RosettaAntibody HCDR2 RMSD vs Native', labelpad=20)
plt.xlabel('Multi-Template HCDR2 RMSD vs Native', labelpad=20)
plt.tight_layout()
plt.savefig(os.getcwd() + '/Figures/H2_by_Score_RosettaAB.png', dpi=600, format='png', transparent=True)
g = sns.lmplot(data=dframe.loc[dframe['Experiment'].isin(['AbPredict'])], x="Multi-Template HCDR2 RMSD vs Native", y="Method HCDR2 RMSD vs Native", hue="Experiment", fit_reg=False, height = 10, aspect = 1, x_estimator=np.median, x_ci=None, scatter_kws={'alpha':1, 's':75}, palette=['grey'], legend=False, legend_out=False)
g.set(ylim=(0, H2_lim), xlim=(0, H2_lim))
ax = plt.gca()
X_plot = np.linspace(0, H2_lim, 400)
Y_plot = X_plot
plt.plot(X_plot, Y_plot, color='silver')
plt.plot(X_plot, Y_plot - Margin, color='silver')
plt.plot(X_plot, Y_plot + Margin, color='silver')
#ax.set_title('Performance of Best AbPredict Models')
plt.ylabel('AbPredict HCDR2 RMSD vs Native', labelpad=20)
plt.xlabel('Multi-Template HCDR2 RMSD vs Native', labelpad=20)
plt.tight_layout()
plt.savefig(os.getcwd() + '/Figures/H2_by_Score_AbPredict.png', dpi=600, format='png', transparent=True)
g = sns.lmplot(data=dframe.loc[dframe['Experiment'].isin(['Single Template'])], x="Multi-Template HCDR2 RMSD vs Native", y="Method HCDR2 RMSD vs Native", hue="Experiment", fit_reg=False, height = 10, aspect = 1, x_estimator=np.median, x_ci=None, scatter_kws={'alpha':1, 's':75}, palette=['grey'], legend=False, legend_out=False)
g.set(ylim=(0, H2_lim), xlim=(0, H2_lim))
ax = plt.gca()
X_plot = np.linspace(0, H2_lim, 400)
Y_plot = X_plot
plt.plot(X_plot, Y_plot, color='silver')
plt.plot(X_plot, Y_plot - Margin, color='silver')
plt.plot(X_plot, Y_plot + Margin, color='silver')
#ax.set_title('Performance of Best Single Template Models')
plt.ylabel('Single Template HCDR2 RMSD vs Native', labelpad=20)
plt.xlabel('Multi-Template HCDR2 RMSD vs Native', labelpad=20)
plt.tight_layout()
plt.savefig(os.getcwd() + '/Figures/H2_by_Score_Single.png', dpi=600, format='png', transparent=True)
g = sns.lmplot(data=dframe.loc[dframe['Experiment'].isin(['Control'])], x="Multi-Template HCDR2 RMSD vs Native", y="Method HCDR2 RMSD vs Native", hue="Experiment", fit_reg=False, height = 10, aspect = 1, x_estimator=np.median, x_ci=None, scatter_kws={'alpha':1, 's':75}, palette=['grey'], legend=False, legend_out=False)
g.set(ylim=(0, H2_lim), xlim=(0, H2_lim))
ax = plt.gca()
X_plot = np.linspace(0, H2_lim, 400)
Y_plot = X_plot
plt.plot(X_plot, Y_plot, color='silver')
plt.plot(X_plot, Y_plot - Margin, color='silver')
plt.plot(X_plot, Y_plot + Margin, color='silver')
#ax.set_title('Performance of Best Control Models')
plt.ylabel('Control HCDR2 RMSD vs Native', labelpad=20)
plt.xlabel('Multi-Template HCDR2 RMSD vs Native', labelpad=20)
plt.tight_layout()
plt.savefig(os.getcwd() + '/Figures/H2_by_Score_Control.png', dpi=600, format='png', transparent=True)




H3_lim = 8


grouped_dframe = dframe.groupby(['Experiment', 'Root','HCDR3 Percent Identity', 'HCDR3 Percent Identity Unbinned','HCDR3 Length'], as_index=False)['Method HCDR3 RMSD vs Native'].median()
#grouped_dframe['Method HCDR3 RMSD vs Native'].values[grouped_dframe['Method HCDR3 RMSD vs Native'].values > 12.5] = float('NaN')
pd.set_option('display.max_rows', 80)
#print(grouped_dframe.loc[grouped_dframe['Experiment'].isin(['AbPredict'])])

plt.clf()
fig, ax = plt.subplots(figsize=(10, 10))
g = sns.boxplot(data=grouped_dframe.loc[grouped_dframe['Experiment'].isin(['RosettaAntibody', 'Single Template', 'Multi-Template', 'Control', 'AbPredict'])], x="HCDR3 Percent Identity", y="Method HCDR3 RMSD vs Native", hue="Experiment", hue_order=['Control', 'Single Template', 'Multi-Template', 'RosettaAntibody', 'AbPredict'], ax=ax)
g.set(ylim=(0, H3_lim + 2))
g.set_xticklabels(['0 - 20%','20 - 40%','40 - 60%','60 - 80%'])
plt.xlabel('HCDR3 Percent Identity', labelpad=20)
legend = plt.legend(loc = 'upper right', frameon=True, framealpha=1, title = 'Method')
#sns.lineplot(data=grouped_dframe.loc[grouped_dframe['Experiment'].isin(['RosettaAntibody', 'Single Template', 'Multi-Template', 'Control', 'AbPredict'])], x="HCDR3 Percent Identity", y="Method HCDR3 RMSD vs Native", hue="Experiment", hue_order=['Control', 'Single Template', 'Multi-Template', 'RosettaAntibody', 'AbPredict'], ax=ax) 
plt.tight_layout()
plt.savefig(os.getcwd() + '/Figures/H3_by_Percent.png', dpi=600, format='png', transparent=True)


plt.clf()
fig, ax = plt.subplots(figsize=(10, 10))
g = sns.lineplot(data=grouped_dframe.loc[grouped_dframe['Experiment'].isin(['RosettaAntibody', 'Single Template', 'Multi-Template', 'Control', 'AbPredict'])], x="HCDR3 Percent Identity", y="Method HCDR3 RMSD vs Native", hue="Experiment", hue_order=['Control', 'Single Template', 'Multi-Template', 'RosettaAntibody', 'AbPredict'], ax=ax)
g.set(ylim=(0, H3_lim + 2))
g.set_xticks([0, 0.2, 0.4, 0.6])
g.set_xticklabels(['0 - 20%','20 - 40%','40 - 60%','60 - 80%'])
plt.xlabel('HCDR3 Percent Identity', labelpad=20)
legend = plt.legend(loc = 'upper right', frameon=True, framealpha=1, title = 'Method') 
plt.tight_layout()
plt.savefig(os.getcwd() + '/Figures/H3_by_Percent_Line.png', dpi=600, format='png', transparent=True)

plt.clf()
fig, ax = plt.subplots(figsize=(10, 10))
g = sns.boxplot(data=dframe.loc[dframe['Experiment'].isin(['RosettaAntibody Inspired', 'Single Template', 'Multi-Template', 'Control', 'RosettaAntibody', 'AbPredict'])], x="HCDR3 Percent Identity Unbinned", y="Method HCDR3 RMSD vs Native", hue="Experiment", ax=ax)
#g.set(ylim=(0, BB_lim))
#legend = plt.legend(loc = 'upper left', frameon=True, framealpha=1, title = 'Method') 
plt.tight_layout()
plt.savefig(os.getcwd() + '/Figures/H3_by_Percent_Unbinned.png', dpi=600, format='png', transparent=True)

plt.clf()
fig, ax = plt.subplots(figsize=(10, 10))
g = sns.boxplot(data=grouped_dframe.loc[grouped_dframe['Experiment'].isin(['Single Template', 'Multi-Template', 'Control', 'RosettaAntibody', 'AbPredict'])], x="HCDR3 Length", y="Method HCDR3 RMSD vs Native", hue="Experiment", hue_order=['Control', 'Single Template', 'Multi-Template', 'RosettaAntibody', 'AbPredict'], order=['Short', 'Medium', 'Long'])
g.set(ylim=(0, H3_lim + 2))
plt.xlabel('HCDR3 Length', labelpad=20)
legend = plt.legend(loc = 'upper left', frameon=True, framealpha=1, title = 'Method') 
plt.tight_layout()
plt.savefig(os.getcwd() + '/Figures/H3_by_Length.png', dpi=600, format='png', transparent=True)

def ax_settings(ax, var_name, x_min, x_max):
	ax.set_xlim(x_min,x_max)
	ax.set_yticks([])
	
	ax.spines['left'].set_visible(False)
	ax.spines['right'].set_visible(False)
	ax.spines['top'].set_visible(False)

	ax.spines['bottom'].set_edgecolor('#444444')
	ax.spines['bottom'].set_linewidth(2)

	ax.text(0.00, 0.5, var_name, fontsize=22, transform = ax.transAxes) 
	return None

plt.clf()
fig = plt.figure(figsize=(10, 10))
gs = gridspec.GridSpec(nrows=5, ncols=1, figure=fig)
ax = [None] * 6
methods = ['Control', 'Single Template', 'Multi-Template', 'RosettaAntibody', 'AbPredict']
for i in range(5):
	ax[i] = fig.add_subplot(gs[i, 0])
	ax_settings(ax[i], methods[i], -2.7, H3_lim + 2.5)
	sns.kdeplot(data=grouped_dframe.loc[(grouped_dframe['Experiment'] == methods[i]) & (grouped_dframe['HCDR3 Length'] == 'Short')]["Method HCDR3 RMSD vs Native"], ax=ax[i], color="blue", bw=0.5, legend=False)
	sns.kdeplot(data=grouped_dframe.loc[(grouped_dframe['Experiment'] == methods[i]) & (grouped_dframe['HCDR3 Length'] == 'Medium')]["Method HCDR3 RMSD vs Native"], ax=ax[i], color="orange", bw=0.5, legend=False)
	sns.kdeplot(data=grouped_dframe.loc[(grouped_dframe['Experiment'] == methods[i]) & (grouped_dframe['HCDR3 Length'] == 'Long')]["Method HCDR3 RMSD vs Native"], ax=ax[i], color="green", bw=0.5, legend=False)
	if i < 4: 
        	ax[i].set_xticks([])
	else: 
        	ax[i].set_xticks([0,2,4,6,8,10])
		ax[i].grid(False)
control_count = grouped_dframe[grouped_dframe['Experiment'] == 'Control']
ax[0].legend(['Short, n={}'.format((control_count['HCDR3 Length'].values == 'Short').sum()), 'Medium, n={}'.format((control_count['HCDR3 Length'].values == 'Medium').sum()), 'Long, n={}'.format((control_count['HCDR3 Length'].values == 'Long').sum())], facecolor = 'w', loc = 'upper right', borderaxespad=0, frameon=True, framealpha=1, title = 'HCDR3 Loop Length')
plt.xlabel('HCDR3 RMSD vs Native', labelpad=20)
plt.tight_layout()
plt.savefig(os.getcwd() + '/Figures/H3_by_Length_New.png', dpi=600, format='png', transparent=True)

plt.clf()
fig = plt.figure(figsize=(10, 10))
gs = gridspec.GridSpec(nrows=5, ncols=1, figure=fig)
ax = [None] * 6
methods = ['Control', 'Single Template', 'Multi-Template', 'RosettaAntibody', 'AbPredict']
for i in range(5):
	ax[i] = fig.add_subplot(gs[i, 0])
	ax_settings(ax[i], methods[i], -2.5, H3_lim + 2.5)
	sns.kdeplot(data=grouped_dframe.loc[(grouped_dframe['Experiment'] == methods[i]) & (grouped_dframe['HCDR3 Percent Identity'] == 0)]["Method HCDR3 RMSD vs Native"], ax=ax[i], color="blue", bw=0.5, legend=False)
	sns.kdeplot(data=grouped_dframe.loc[(grouped_dframe['Experiment'] == methods[i]) & (grouped_dframe['HCDR3 Percent Identity'] == 0.2)]["Method HCDR3 RMSD vs Native"], ax=ax[i], color="orange", bw=0.5, legend=False)
	sns.kdeplot(data=grouped_dframe.loc[(grouped_dframe['Experiment'] == methods[i]) & (grouped_dframe['HCDR3 Percent Identity'] == 0.4)]["Method HCDR3 RMSD vs Native"], ax=ax[i], color="green", bw=0.5, legend=False)
	sns.kdeplot(data=grouped_dframe.loc[(grouped_dframe['Experiment'] == methods[i]) & (grouped_dframe['HCDR3 Percent Identity'] == 0.6)]["Method HCDR3 RMSD vs Native"], ax=ax[i], color="purple", bw=0.5, legend=False)
	if i < 4: 
        	ax[i].set_xticks([])
	else: 
        	ax[i].set_xticks([0,2,4,6,8,10])
		ax[i].grid(False)
ax[0].legend(['  0 - 20%','20 - 40%','40 - 60%','60 - 80%'], facecolor = 'w', loc = 'upper right', borderaxespad=0, ncol=2, columnspacing=0.5, frameon=True, framealpha=1, title = 'HCDR3 Percent Similarity')
plt.xlabel('HCDR3 RMSD vs Native', labelpad=20)
plt.tight_layout()
plt.savefig(os.getcwd() + '/Figures/H3_by_Percent_New.png', dpi=600, format='png', transparent=True)



plt.clf()
g = sns.lmplot(data=dframe, x="Multi-Template HCDR3 RMSD vs Native", y="Method HCDR3 RMSD vs Native", hue="Experiment", x_estimator=np.median, x_ci=None, fit_reg=False, height = 10, aspect = 1, scatter_kws={'alpha':1, 's':75}, palette=['grey'], legend_out=False)
g.set(ylim=(0, None), xlim=(0, None))
legend = plt.legend(loc = 'upper left', frameon=True, framealpha=1, title = 'Method')
frame = legend.get_frame()
frame.set_facecolor('white')
frame.set_edgecolor('black')
for lh in legend.legendHandles:
	lh.set_alpha(1)
	lh.set_sizes([100])
ax = plt.gca()
#ax.set_title('Performance of Best Models')
plt.tight_layout()
plt.savefig(os.getcwd() + '/Figures/H3_by_Score.png', dpi=600, format='png', transparent=True)

H3_lim = 10
g = sns.lmplot(data=dframe.loc[dframe['Experiment'].isin(['RosettaAntibody Inspired'])], x="Multi-Template HCDR3 RMSD vs Native", y="Method HCDR3 RMSD vs Native", hue="Experiment", fit_reg=False, height = 10, aspect = 1, x_estimator=np.median, x_ci=None, scatter_kws={'alpha':1, 's':75}, palette=['grey'], legend=False, legend_out=False)
g.set(ylim=(0, H3_lim), xlim=(0, H3_lim))
ax = plt.gca()
X_plot = np.linspace(0, H3_lim, 400)
Y_plot = X_plot
plt.plot(X_plot, Y_plot, color='silver')
plt.plot(X_plot, Y_plot - 1, color='silver')
plt.plot(X_plot, Y_plot + 1, color='silver')
#ax.set_title('Performance of Best RosettaAntibody Inspired Models')
plt.ylabel('RosettaAntibody Inspired HCDR3 RMSD vs Native', labelpad=20)
plt.tight_layout()
plt.savefig(os.getcwd() + '/Figures/H3_by_Score_Rosetta.png', dpi=600, format='png', transparent=True)
g = sns.lmplot(data=dframe.loc[dframe['Experiment'].isin(['RosettaAntibody'])], x="Multi-Template HCDR3 RMSD vs Native", y="Method HCDR3 RMSD vs Native", hue="Experiment", fit_reg=False, height = 10, aspect = 1, x_estimator=np.median, x_ci=None, scatter_kws={'alpha':1, 's':75}, palette=['grey'], legend=False, legend_out=False)
g.set(ylim=(0, H3_lim), xlim=(0, H3_lim))
ax = plt.gca()
X_plot = np.linspace(0, H3_lim, 400)
Y_plot = X_plot
plt.plot(X_plot, Y_plot, color='silver')
plt.plot(X_plot, Y_plot - 1, color='silver')
plt.plot(X_plot, Y_plot + 1, color='silver')
#ax.set_title('Performance of Best RosettaAntibody Models')
plt.xlabel('Multi-Template HCDR3 RMSD vs Native', labelpad=20)
plt.ylabel('RosettaAntibody HCDR3 RMSD vs Native', labelpad=20)
plt.tight_layout()
plt.savefig(os.getcwd() + '/Figures/H3_by_Score_RosettaAB.png', dpi=600, format='png', transparent=True)
g = sns.lmplot(data=dframe.loc[dframe['Experiment'].isin(['AbPredict'])], x="Multi-Template HCDR3 RMSD vs Native", y="Method HCDR3 RMSD vs Native", hue="Experiment", fit_reg=False, height = 10, aspect = 1, x_estimator=np.median, x_ci=None, scatter_kws={'alpha':1, 's':75}, palette=['grey'], legend=False, legend_out=False)
g.set(ylim=(0, H3_lim), xlim=(0, H3_lim))
ax = plt.gca()
X_plot = np.linspace(0, H3_lim, 400)
Y_plot = X_plot
plt.plot(X_plot, Y_plot, color='silver')
plt.plot(X_plot, Y_plot - 1, color='silver')
plt.plot(X_plot, Y_plot + 1, color='silver')
#ax.set_title('Performance of Best AbPredict Models')
plt.xlabel('Multi-Template HCDR3 RMSD vs Native', labelpad=20)
plt.ylabel('AbPredict HCDR3 RMSD vs Native', labelpad=20)
plt.tight_layout()
plt.savefig(os.getcwd() + '/Figures/H3_by_Score_AbPredict.png', dpi=600, format='png', transparent=True)
g = sns.lmplot(data=dframe.loc[dframe['Experiment'].isin(['Single Template'])], x="Multi-Template HCDR3 RMSD vs Native", y="Method HCDR3 RMSD vs Native", hue="Experiment", fit_reg=False, height = 10, aspect = 1, x_estimator=np.median, x_ci=None, scatter_kws={'alpha':1, 's':75}, palette=['grey'], legend=False, legend_out=False)
g.set(ylim=(0, H3_lim), xlim=(0, H3_lim))
ax = plt.gca()
X_plot = np.linspace(0, H3_lim, 400)
Y_plot = X_plot
plt.plot(X_plot, Y_plot, color='silver')
plt.plot(X_plot, Y_plot - 1, color='silver')
plt.plot(X_plot, Y_plot + 1, color='silver')
#ax.set_title('Performance of Best Single Template Models')
plt.xlabel('Multi-Template HCDR3 RMSD vs Native', labelpad=17)
plt.ylabel('Single Template HCDR3 RMSD vs Native', labelpad=16)
plt.tight_layout()
plt.savefig(os.getcwd() + '/Figures/H3_by_Score_Single.png', dpi=600, format='png', transparent=True)
g = sns.lmplot(data=dframe.loc[dframe['Experiment'].isin(['Control'])], x="Multi-Template HCDR3 RMSD vs Native", y="Method HCDR3 RMSD vs Native", hue="Experiment", fit_reg=False, height = 10, aspect = 1, x_estimator=np.median, x_ci=None, scatter_kws={'alpha':1, 's':75}, palette=['grey'], legend=False, legend_out=False)
g.set(ylim=(0, H3_lim), xlim=(0, H3_lim))
ax = plt.gca()
X_plot = np.linspace(0, H3_lim, 400)
Y_plot = X_plot
plt.plot(X_plot, Y_plot, color='silver')
plt.plot(X_plot, Y_plot - 1, color='silver')
plt.plot(X_plot, Y_plot + 1, color='silver')
#ax.set_title('Performance of Best Control Models')
plt.xlabel('Multi-Template HCDR3 RMSD vs Native', labelpad=20)
plt.ylabel('Control HCDR3 RMSD vs Native', labelpad=20)
plt.tight_layout()
plt.savefig(os.getcwd() + '/Figures/H3_by_Score_Control.png', dpi=600, format='png', transparent=True)

plt.clf()

grouped_dframeH3 = dframe.groupby(['Experiment', 'Root'], as_index=False)['Method HCDR3 RMSD vs Native'].median()
grouped_dframeH3['Method HCDR3 RMSD vs Native'] = (grouped_dframeH3['Method HCDR3 RMSD vs Native'] // 2.5) * 2.5
grouped_dframeH3['Method HCDR3 RMSD vs Native'].values[grouped_dframeH3['Method HCDR3 RMSD vs Native'].values > 12.5] = float('NaN')
grouped_dframeH3['Method HCDR3 RMSD vs Native Binned'] = grouped_dframeH3['Method HCDR3 RMSD vs Native'].astype(str) + ' -\n' + (grouped_dframeH3['Method HCDR3 RMSD vs Native'] + 2.5).astype(str)
grouped_dframeH3.loc[grouped_dframeH3['Method HCDR3 RMSD vs Native Binned'].str.contains('nan'), 'Method HCDR3 RMSD vs Native Binned'] = 'Failed to\nModel'
grouped_dframeH3_2 = grouped_dframeH3.groupby(['Experiment', 'Method HCDR3 RMSD vs Native Binned']).size().to_frame('Count')
grouped_dframeH3_2 = grouped_dframeH3_2.reset_index()
grouped_dframeH3_2.sort_values(by=['Method HCDR3 RMSD vs Native Binned'], inplace=True)
#print(grouped_dframeH3_2)

matplotlib.rcParams["patch.force_edgecolor"] = True

fig, ax = plt.subplots(figsize=(10, 10))
g = sns.barplot(data=grouped_dframeH3_2.loc[grouped_dframeH3_2['Experiment'].isin(['RosettaAntibody', 'Single Template', 'Multi-Template', 'Control', 'RosettaAntibody', 'AbPredict'])], x="Method HCDR3 RMSD vs Native Binned", y="Count", hue="Experiment", order=['0.0 -\n2.5','2.5 -\n5.0','5.0 -\n7.5','7.5 -\n10.0','10.0 -\n12.5','Failed to\nModel'], hue_order=['Control', 'Single Template', 'Multi-Template', 'RosettaAntibody', 'AbPredict'],linewidth=1.5,edgecolor="black")
g.set_xticklabels(g.get_xticklabels(), rotation=45, horizontalalignment='center')
plt.xlabel('HCDR3 RMSD vs Native')
plt.ylabel('Count', labelpad=20)
legend = plt.legend(loc = 'upper right', frameon=True, framealpha=1, title = 'Method') 
plt.ylim(0, 90)
plt.tight_layout()
plt.savefig(os.getcwd() + '/Figures/H3_by_Bin.png', dpi=600, format='png', transparent=True, bbox_inches = "tight")

grouped_dframeH3 = dframe.groupby(['Experiment', 'Root'], as_index=False)['Method HCDR3 RMSD vs Native'].median()
fig, ax = plt.subplots(figsize=(10, 10))
for a in ['Single Template', 'Multi-Template', 'Control', 'RosettaAntibody', 'AbPredict']:
	sns.distplot(grouped_dframeH3.loc[grouped_dframeH3['Experiment'] == a]['Method HCDR3 RMSD vs Native'], bins=range(0, 16, 1), ax=ax, norm_hist=False, kde=False, hist_kws={"histtype": "bar"})
#plt.ylim([0,0.14])
plt.xlim([0,15])
plt.xlabel('HCDR3 RMSD vs Native', labelpad=20)
plt.ylabel('Count', labelpad=20)
#plt.title('', fontweight='bold', pad=20)
plt.legend(['Single Template', 'Multi-Template', 'Control', 'RosettaAntibody', 'AbPredict'], loc = 'upper left')

sns.despine()
plt.tight_layout()
plt.savefig(os.getcwd() + '/Figures/H3_by_Bin_Histo.png', dpi=600, format='png', transparent=True, bbox_inches = "tight")

fig, ax = plt.subplots(figsize=(10, 10))
g = sns.violinplot(data=grouped_dframeH3.loc[grouped_dframeH3['Experiment'].isin(['Single Template', 'Multi-Template', 'Control', 'RosettaAntibody', 'AbPredict'])], x="Experiment", y="Method HCDR3 RMSD vs Native", order=['Control', 'Single Template', 'Multi-Template', 'RosettaAntibody', 'AbPredict'], bw=0.15)
g.set(ylim=(0, H3_lim))
g.set_xticklabels(g.get_xticklabels(), rotation=90, horizontalalignment='center')
plt.xlabel('')
plt.ylabel('HCDR3 RMSD vs Native', labelpad=20)
#plt.tight_layout()
plt.savefig(os.getcwd() + '/Figures/H3_by_RMSD_Violin.png', dpi=600, format='png', transparent=True, bbox_inches = "tight") 


RABABPH3dict = {}
for index, row in grouped_dframeH3.iterrows():
	RABABPH3dict.setdefault(row['Root'], [[],[]])
	if row['Experiment'] == 'RosettaAntibody':
		RABABPH3dict[row['Root']][0] = row['Method HCDR3 RMSD vs Native']
	elif row['Experiment'] == 'AbPredict':
		RABABPH3dict[row['Root']][1] = row['Method HCDR3 RMSD vs Native']
RABABPH3list = list()
for key, value in RABABPH3dict.items():
	RABABPH3list.append({ "Root" : key, "RosettaAntibody H3 RMSD": value[0], "AbPredict H3 RMSD": value[1]})
RABABPH3 = DataFrame(RABABPH3list)

plt.clf()
g = sns.lmplot(data=RABABPH3, x="RosettaAntibody H3 RMSD", y="AbPredict H3 RMSD", fit_reg=False, height = 10, aspect = 1, x_ci=None, scatter_kws={'alpha':1, 's':75}, palette=['grey'], legend=False, legend_out=False)
g.set(ylim=(0, H3_lim), xlim=(0, H3_lim))
ax = plt.gca()
X_plot = np.linspace(0, H3_lim, 400)
Y_plot = X_plot
plt.plot(X_plot, Y_plot, color='silver')
plt.plot(X_plot, Y_plot - 1, color='silver')
plt.plot(X_plot, Y_plot + 1, color='silver')
plt.ylabel('AbPredict HCDR3 RMSD vs Native', labelpad=20)
plt.xlabel('RosettaAntibody HCDR3 RMSD vs Native', labelpad=20)
plt.tight_layout()
plt.savefig(os.getcwd() + '/Figures/H3_by_RAB_ABP.png', dpi=600, format='png', transparent=True)

L1_lim = 4

plt.clf()
fig, ax = plt.subplots(figsize=(10, 10))
g = sns.boxplot(data=dframe.loc[dframe['Experiment'].isin(['RosettaAntibody Inspired', 'Single Template', 'Multi-Template', 'Control', 'RosettaAntibody', 'AbPredict'])], x="LCDR1 Percent Identity", y="Method LCDR1 RMSD vs Native", hue="Experiment", ax=ax)
g.set(ylim=(0, L1_lim))
legend = plt.legend(loc = 'upper left', frameon=True, framealpha=1, title = 'Method')  
plt.tight_layout()
plt.savefig(os.getcwd() + '/Figures/L1_by_Percent.png', dpi=600, format='png', transparent=True)

plt.clf()
fig, ax = plt.subplots(figsize=(10, 10))
g = sns.boxplot(data=dframe.loc[dframe['Experiment'].isin(['RosettaAntibody Inspired', 'Single Template', 'Multi-Template', 'Control', 'RosettaAntibody', 'AbPredict'])], x="LCDR1 Length", y="Method LCDR1 RMSD vs Native", hue="Experiment", ax=ax)
g.set(ylim=(0, L1_lim))
legend = plt.legend(loc = 'upper left', frameon=True, framealpha=1, title = 'Method') 
plt.tight_layout()
plt.savefig(os.getcwd() + '/Figures/L1_by_Length.png', dpi=600, format='png', transparent=True)

g = sns.lmplot(data=dframe, x="Multi-Template LCDR1 RMSD vs Native", y="Method LCDR1 RMSD vs Native", hue="Experiment", x_estimator=np.median, x_ci=None, fit_reg=False, height = 10, aspect = 1, scatter_kws={'alpha':1, 's':75}, palette=['grey'], legend_out=False)
g.set(ylim=(0, L1_lim), xlim=(0, L1_lim))
legend = plt.legend(loc = 'upper left', frameon=True, framealpha=1, title = 'Method')
frame = legend.get_frame()
frame.set_facecolor('white')
frame.set_edgecolor('black')
for lh in legend.legendHandles:
	lh.set_alpha(1)
	lh.set_sizes([100])
ax = plt.gca()
#ax.set_title('Performance of Best Models')
plt.tight_layout()
plt.savefig(os.getcwd() + '/Figures/L1_by_Score.png', dpi=600, format='png', transparent=True)
g = sns.lmplot(data=dframe.loc[dframe['Experiment'].isin(['RosettaAntibody Inspired'])], x="Multi-Template LCDR1 RMSD vs Native", y="Method LCDR1 RMSD vs Native", hue="Experiment", fit_reg=False, height = 10, aspect = 1, x_estimator=np.median, x_ci=None, scatter_kws={'alpha':1, 's':75}, palette=['grey'], legend=False, legend_out=False)
g.set(ylim=(0, L1_lim), xlim=(0, L1_lim))
ax = plt.gca()
X_plot = np.linspace(0, L1_lim, 400)
Y_plot = X_plot
plt.plot(X_plot, Y_plot, color='silver')
plt.plot(X_plot, Y_plot - Margin, color='silver')
plt.plot(X_plot, Y_plot + Margin, color='silver')
#ax.set_title('Performance of Best RosettaAntibody Inspired Models')
plt.ylabel('RosettaAntibody Inspired LCDR1 RMSD vs Native', labelpad=20)
plt.xlabel('Multi-Template LCDR1 RMSD vs Native', labelpad=20)
plt.tight_layout()
plt.savefig(os.getcwd() + '/Figures/L1_by_Score_Rosetta.png', dpi=600, format='png', transparent=True)
g = sns.lmplot(data=dframe.loc[dframe['Experiment'].isin(['RosettaAntibody'])], x="Multi-Template LCDR1 RMSD vs Native", y="Method LCDR1 RMSD vs Native", hue="Experiment", fit_reg=False, height = 10, aspect = 1, x_estimator=np.median, x_ci=None, scatter_kws={'alpha':1, 's':75}, palette=['grey'], legend=False, legend_out=False)
g.set(ylim=(0, L1_lim), xlim=(0, L1_lim))
ax = plt.gca()
X_plot = np.linspace(0, L1_lim, 400)
Y_plot = X_plot
plt.plot(X_plot, Y_plot, color='silver')
plt.plot(X_plot, Y_plot - Margin, color='silver')
plt.plot(X_plot, Y_plot + Margin, color='silver')
#ax.set_title('Performance of Best RosettaAntibody Models')
plt.ylabel('RosettaAntibody LCDR1 RMSD vs Native', labelpad=20)
plt.xlabel('Multi-Template LCDR1 RMSD vs Native', labelpad=20)
plt.tight_layout()
plt.savefig(os.getcwd() + '/Figures/L1_by_Score_RosettaAB.png', dpi=600, format='png', transparent=True)
g = sns.lmplot(data=dframe.loc[dframe['Experiment'].isin(['AbPredict'])], x="Multi-Template LCDR1 RMSD vs Native", y="Method LCDR1 RMSD vs Native", hue="Experiment", fit_reg=False, height = 10, aspect = 1, x_estimator=np.median, x_ci=None, scatter_kws={'alpha':1, 's':75}, palette=['grey'], legend=False, legend_out=False)
g.set(ylim=(0, L1_lim), xlim=(0, L1_lim))
ax = plt.gca()
X_plot = np.linspace(0, L1_lim, 400)
Y_plot = X_plot
plt.plot(X_plot, Y_plot, color='silver')
plt.plot(X_plot, Y_plot - Margin, color='silver')
plt.plot(X_plot, Y_plot + Margin, color='silver')
#ax.set_title('Performance of Best AbPredict Models')
plt.ylabel('AbPredict LCDR1 RMSD vs Native', labelpad=20)
plt.xlabel('Multi-Template LCDR1 RMSD vs Native', labelpad=20)
plt.tight_layout()
plt.savefig(os.getcwd() + '/Figures/L1_by_Score_AbPredict.png', dpi=600, format='png', transparent=True)
g = sns.lmplot(data=dframe.loc[dframe['Experiment'].isin(['Single Template'])], x="Multi-Template LCDR1 RMSD vs Native", y="Method LCDR1 RMSD vs Native", hue="Experiment", fit_reg=False, height = 10, aspect = 1, x_estimator=np.median, x_ci=None, scatter_kws={'alpha':1, 's':75}, palette=['grey'], legend=False, legend_out=False)
g.set(ylim=(0, L1_lim), xlim=(0, L1_lim))
ax = plt.gca()
X_plot = np.linspace(0, L1_lim, 400)
Y_plot = X_plot
plt.plot(X_plot, Y_plot, color='silver')
plt.plot(X_plot, Y_plot - Margin, color='silver')
plt.plot(X_plot, Y_plot + Margin, color='silver')
#ax.set_title('Performance of Best Single Template Models')
plt.ylabel('Single Template LCDR1 RMSD vs Native', labelpad=20)
plt.xlabel('Multi-Template LCDR1 RMSD vs Native', labelpad=20)
plt.tight_layout()
plt.savefig(os.getcwd() + '/Figures/L1_by_Score_Single.png', dpi=600, format='png', transparent=True)
g = sns.lmplot(data=dframe.loc[dframe['Experiment'].isin(['Control'])], x="Multi-Template LCDR1 RMSD vs Native", y="Method LCDR1 RMSD vs Native", hue="Experiment", fit_reg=False, height = 10, aspect = 1, x_estimator=np.median, x_ci=None, scatter_kws={'alpha':1, 's':75}, palette=['grey'], legend=False, legend_out=False)
g.set(ylim=(0, L1_lim), xlim=(0, L1_lim))
ax = plt.gca()
X_plot = np.linspace(0, L1_lim, 400)
Y_plot = X_plot
plt.plot(X_plot, Y_plot, color='silver')
plt.plot(X_plot, Y_plot - Margin, color='silver')
plt.plot(X_plot, Y_plot + Margin, color='silver')
#ax.set_title('Performance of Best Control Models')
plt.ylabel('Control LCDR1 RMSD vs Native', labelpad=20)
plt.xlabel('Multi-Template LCDR1 RMSD vs Native', labelpad=20)
plt.tight_layout()
plt.savefig(os.getcwd() + '/Figures/L1_by_Score_Control.png', dpi=600, format='png', transparent=True)





L2_lim = 4

plt.clf()
fig, ax = plt.subplots(figsize=(10, 10))
g = sns.boxplot(data=dframe.loc[dframe['Experiment'].isin(['RosettaAntibody Inspired', 'Single Template', 'Multi-Template', 'Control', 'RosettaAntibody', 'AbPredict'])], x="LCDR2 Percent Identity", y="Method LCDR2 RMSD vs Native", hue="Experiment", ax=ax)
g.set(ylim=(0, L2_lim))
legend = plt.legend(loc = 'upper left', frameon=True, framealpha=1, title = 'Method')  
plt.tight_layout()
plt.savefig(os.getcwd() + '/Figures/L2_by_Percent.png', dpi=600, format='png', transparent=True)

plt.clf()
fig, ax = plt.subplots(figsize=(10, 10))
g = sns.boxplot(data=dframe.loc[dframe['Experiment'].isin(['RosettaAntibody Inspired', 'Single Template', 'Multi-Template', 'Control', 'RosettaAntibody', 'AbPredict'])], x="LCDR2 Length", y="Method LCDR2 RMSD vs Native", hue="Experiment", ax=ax)
g.set(ylim=(0, L2_lim))
legend = plt.legend(loc = 'upper left', frameon=True, framealpha=1, title = 'Method') 
plt.tight_layout()
plt.savefig(os.getcwd() + '/Figures/L2_by_Length.png', dpi=600, format='png', transparent=True)

g = sns.lmplot(data=dframe, x="Multi-Template LCDR2 RMSD vs Native", y="Method LCDR2 RMSD vs Native", hue="Experiment", x_estimator=np.median, x_ci=None, fit_reg=False, height = 10, aspect = 1, scatter_kws={'alpha':1, 's':75}, palette=['grey'], legend_out=False)
g.set(ylim=(0, L2_lim), xlim=(0, L2_lim))
legend = plt.legend(loc = 'upper left', frameon=True, framealpha=1, title = 'Method')
frame = legend.get_frame()
frame.set_facecolor('white')
frame.set_edgecolor('black')
for lh in legend.legendHandles:
	lh.set_alpha(1)
	lh.set_sizes([100])
ax = plt.gca()
#ax.set_title('Performance of Best Models')
plt.tight_layout()
plt.savefig(os.getcwd() + '/Figures/L2_by_Score.png', dpi=600, format='png', transparent=True)
g = sns.lmplot(data=dframe.loc[dframe['Experiment'].isin(['RosettaAntibody Inspired'])], x="Multi-Template LCDR2 RMSD vs Native", y="Method LCDR2 RMSD vs Native", hue="Experiment", fit_reg=False, height = 10, aspect = 1, x_estimator=np.median, x_ci=None, scatter_kws={'alpha':1, 's':75}, palette=['grey'], legend=False, legend_out=False)
g.set(ylim=(0, L2_lim), xlim=(0, L2_lim))
ax = plt.gca()
X_plot = np.linspace(0, L2_lim, 400)
Y_plot = X_plot
plt.plot(X_plot, Y_plot, color='silver')
plt.plot(X_plot, Y_plot - Margin, color='silver')
plt.plot(X_plot, Y_plot + Margin, color='silver')
#ax.set_title('Performance of Best RosettaAntibody Inspired Models')
plt.ylabel('RosettaAntibody Inspired LCDR2 RMSD vs Native', labelpad=20)
plt.xlabel('Multi-Template LCDR2 RMSD vs Native', labelpad=20)
plt.tight_layout()
plt.savefig(os.getcwd() + '/Figures/L2_by_Score_Rosetta.png', dpi=600, format='png', transparent=True)
g = sns.lmplot(data=dframe.loc[dframe['Experiment'].isin(['RosettaAntibody'])], x="Multi-Template LCDR2 RMSD vs Native", y="Method LCDR2 RMSD vs Native", hue="Experiment", fit_reg=False, height = 10, aspect = 1, x_estimator=np.median, x_ci=None, scatter_kws={'alpha':1, 's':75}, palette=['grey'], legend=False, legend_out=False)
g.set(ylim=(0, L2_lim), xlim=(0, L2_lim))
ax = plt.gca()
X_plot = np.linspace(0, L2_lim, 400)
Y_plot = X_plot
plt.plot(X_plot, Y_plot, color='silver')
plt.plot(X_plot, Y_plot - Margin, color='silver')
plt.plot(X_plot, Y_plot + Margin, color='silver')
#ax.set_title('Performance of Best RosettaAntibody Models')
plt.ylabel('RosettaAntibody LCDR2 RMSD vs Native', labelpad=20)
plt.xlabel('Multi-Template LCDR2 RMSD vs Native', labelpad=20)
plt.tight_layout()
plt.savefig(os.getcwd() + '/Figures/L2_by_Score_RosettaAB.png', dpi=600, format='png', transparent=True)
g = sns.lmplot(data=dframe.loc[dframe['Experiment'].isin(['AbPredict'])], x="Multi-Template LCDR2 RMSD vs Native", y="Method LCDR2 RMSD vs Native", hue="Experiment", fit_reg=False, height = 10, aspect = 1, x_estimator=np.median, x_ci=None, scatter_kws={'alpha':1, 's':75}, palette=['grey'], legend=False, legend_out=False)
g.set(ylim=(0, L2_lim), xlim=(0, L2_lim))
ax = plt.gca()
X_plot = np.linspace(0, L2_lim, 400)
Y_plot = X_plot
plt.plot(X_plot, Y_plot, color='silver')
plt.plot(X_plot, Y_plot - Margin, color='silver')
plt.plot(X_plot, Y_plot + Margin, color='silver')
#ax.set_title('Performance of Best AbPredict Models')
plt.ylabel('AbPredict LCDR2 RMSD vs Native', labelpad=20)
plt.xlabel('Multi-Template LCDR2 RMSD vs Native', labelpad=20)
plt.tight_layout()
plt.savefig(os.getcwd() + '/Figures/L2_by_Score_AbPredict.png', dpi=600, format='png', transparent=True)
g = sns.lmplot(data=dframe.loc[dframe['Experiment'].isin(['Single Template'])], x="Multi-Template LCDR2 RMSD vs Native", y="Method LCDR2 RMSD vs Native", hue="Experiment", fit_reg=False, height = 10, aspect = 1, x_estimator=np.median, x_ci=None, scatter_kws={'alpha':1, 's':75}, palette=['grey'], legend=False, legend_out=False)
g.set(ylim=(0, L2_lim), xlim=(0, L2_lim))
ax = plt.gca()
X_plot = np.linspace(0, L2_lim, 400)
Y_plot = X_plot
plt.plot(X_plot, Y_plot, color='silver')
plt.plot(X_plot, Y_plot - Margin, color='silver')
plt.plot(X_plot, Y_plot + Margin, color='silver')
#ax.set_title('Performance of Best Single Template Models')
plt.ylabel('Single Template LCDR2 RMSD vs Native', labelpad=20)
plt.xlabel('Multi-Template LCDR2 RMSD vs Native', labelpad=20)
plt.tight_layout()
plt.savefig(os.getcwd() + '/Figures/L2_by_Score_Single.png', dpi=600, format='png', transparent=True)
g = sns.lmplot(data=dframe.loc[dframe['Experiment'].isin(['Control'])], x="Multi-Template LCDR2 RMSD vs Native", y="Method LCDR2 RMSD vs Native", hue="Experiment", fit_reg=False, height = 10, aspect = 1, x_estimator=np.median, x_ci=None, scatter_kws={'alpha':1, 's':75}, palette=['grey'], legend=False, legend_out=False)
g.set(ylim=(0, L2_lim), xlim=(0, L2_lim))
ax = plt.gca()
X_plot = np.linspace(0, L2_lim, 400)
Y_plot = X_plot
plt.plot(X_plot, Y_plot, color='silver')
plt.plot(X_plot, Y_plot - Margin, color='silver')
plt.plot(X_plot, Y_plot + Margin, color='silver')
#ax.set_title('Performance of Best Control Models')
plt.ylabel('Control LCDR2 RMSD vs Native', labelpad=20)
plt.xlabel('Multi-Template LCDR2 RMSD vs Native', labelpad=20)
plt.tight_layout()
plt.savefig(os.getcwd() + '/Figures/L2_by_Score_Control.png', dpi=600, format='png', transparent=True)



L3_lim = 4

plt.clf()
fig, ax = plt.subplots(figsize=(10, 10))
g = sns.boxplot(data=dframe.loc[dframe['Experiment'].isin(['RosettaAntibody Inspired', 'Single Template', 'Multi-Template', 'Control', 'RosettaAntibody', 'AbPredict'])], x="LCDR3 Percent Identity", y="Method LCDR3 RMSD vs Native", hue="Experiment", ax=ax)
g.set(ylim=(0, L3_lim))
legend = plt.legend(loc = 'upper left', frameon=True, framealpha=1, title = 'Method')  
plt.tight_layout()
plt.savefig(os.getcwd() + '/Figures/L3_by_Percent.png', dpi=600, format='png', transparent=True)

plt.clf()
fig, ax = plt.subplots(figsize=(10, 10))
g = sns.boxplot(data=dframe.loc[dframe['Experiment'].isin(['RosettaAntibody Inspired', 'Single Template', 'Multi-Template', 'Control', 'RosettaAntibody', 'AbPredict'])], x="LCDR3 Length", y="Method LCDR3 RMSD vs Native", hue="Experiment", ax=ax)
g.set(ylim=(0, L3_lim))
legend = plt.legend(loc = 'upper left', frameon=True, framealpha=1, title = 'Method') 
plt.tight_layout()
plt.savefig(os.getcwd() + '/Figures/L3_by_Length.png', dpi=600, format='png', transparent=True)

g = sns.lmplot(data=dframe, x="Multi-Template LCDR3 RMSD vs Native", y="Method LCDR3 RMSD vs Native", hue="Experiment", x_estimator=np.median, x_ci=None, fit_reg=False, height = 10, aspect = 1, scatter_kws={'alpha':1, 's':75}, palette=['grey'], legend_out=False)
g.set(ylim=(0, L3_lim), xlim=(0, L3_lim))
legend = plt.legend(loc = 'upper left', frameon=True, framealpha=1, title = 'Method')
frame = legend.get_frame()
frame.set_facecolor('white')
frame.set_edgecolor('black')
for lh in legend.legendHandles:
	lh.set_alpha(1)
	lh.set_sizes([100])
ax = plt.gca()
#ax.set_title('Performance of Best Models')
plt.tight_layout()
plt.savefig(os.getcwd() + '/Figures/L3_by_Score.png', dpi=600, format='png', transparent=True)
g = sns.lmplot(data=dframe.loc[dframe['Experiment'].isin(['RosettaAntibody Inspired'])], x="Multi-Template LCDR3 RMSD vs Native", y="Method LCDR3 RMSD vs Native", hue="Experiment", fit_reg=False, height = 10, aspect = 1, x_estimator=np.median, x_ci=None, scatter_kws={'alpha':1, 's':75}, palette=['grey'], legend=False, legend_out=False)
g.set(ylim=(0, L3_lim), xlim=(0, L3_lim))
ax = plt.gca()
X_plot = np.linspace(0, L3_lim, 400)
Y_plot = X_plot
plt.plot(X_plot, Y_plot, color='silver')
plt.plot(X_plot, Y_plot - Margin, color='silver')
plt.plot(X_plot, Y_plot + Margin, color='silver')
#ax.set_title('Performance of Best RosettaAntibody Inspired Models')
plt.ylabel('RosettaAntibody Inspired LCDR3 RMSD vs Native', labelpad=20)
plt.xlabel('Multi-Template LCDR3 RMSD vs Native', labelpad=20)
plt.tight_layout()
plt.savefig(os.getcwd() + '/Figures/L3_by_Score_Rosetta.png', dpi=600, format='png', transparent=True)
g = sns.lmplot(data=dframe.loc[dframe['Experiment'].isin(['RosettaAntibody'])], x="Multi-Template LCDR3 RMSD vs Native", y="Method LCDR3 RMSD vs Native", hue="Experiment", fit_reg=False, height = 10, aspect = 1, x_estimator=np.median, x_ci=None, scatter_kws={'alpha':1, 's':75}, palette=['grey'], legend=False, legend_out=False)
g.set(ylim=(0, L3_lim), xlim=(0, L3_lim))
ax = plt.gca()
X_plot = np.linspace(0, L3_lim, 400)
Y_plot = X_plot
plt.plot(X_plot, Y_plot, color='silver')
plt.plot(X_plot, Y_plot - Margin, color='silver')
plt.plot(X_plot, Y_plot + Margin, color='silver')
#ax.set_title('Performance of Best RosettaAntibody Models')
plt.ylabel('RosettaAntibody LCDR3 RMSD vs Native', labelpad=20)
plt.xlabel('Multi-Template LCDR3 RMSD vs Native', labelpad=20)
plt.tight_layout()
plt.savefig(os.getcwd() + '/Figures/L3_by_Score_RosettaAB.png', dpi=600, format='png', transparent=True)
g = sns.lmplot(data=dframe.loc[dframe['Experiment'].isin(['AbPredict'])], x="Multi-Template LCDR3 RMSD vs Native", y="Method LCDR3 RMSD vs Native", hue="Experiment", fit_reg=False, height = 10, aspect = 1, x_estimator=np.median, x_ci=None, scatter_kws={'alpha':1, 's':75}, palette=['grey'], legend=False, legend_out=False)
g.set(ylim=(0, L3_lim), xlim=(0, L3_lim))
ax = plt.gca()
X_plot = np.linspace(0, L3_lim, 400)
Y_plot = X_plot
plt.plot(X_plot, Y_plot, color='silver')
plt.plot(X_plot, Y_plot - Margin, color='silver')
plt.plot(X_plot, Y_plot + Margin, color='silver')
#ax.set_title('Performance of Best AbPredict Models')
plt.ylabel('AbPredict LCDR3 RMSD vs Native', labelpad=20)
plt.xlabel('Multi-Template LCDR3 RMSD vs Native', labelpad=20)
plt.tight_layout()
plt.savefig(os.getcwd() + '/Figures/L3_by_Score_AbPredict.png', dpi=600, format='png', transparent=True)
g = sns.lmplot(data=dframe.loc[dframe['Experiment'].isin(['Single Template'])], x="Multi-Template LCDR3 RMSD vs Native", y="Method LCDR3 RMSD vs Native", hue="Experiment", fit_reg=False, height = 10, aspect = 1, x_estimator=np.median, x_ci=None, scatter_kws={'alpha':1, 's':75}, palette=['grey'], legend=False, legend_out=False)
g.set(ylim=(0, L3_lim), xlim=(0, L3_lim))
ax = plt.gca()
X_plot = np.linspace(0, L3_lim, 400)
Y_plot = X_plot
plt.plot(X_plot, Y_plot, color='silver')
plt.plot(X_plot, Y_plot - Margin, color='silver')
plt.plot(X_plot, Y_plot + Margin, color='silver')
#ax.set_title('Performance of Best Single Template Models')
plt.ylabel('Single Template LCDR3 RMSD vs Native', labelpad=20)
plt.xlabel('Multi-Template LCDR3 RMSD vs Native', labelpad=20)
plt.tight_layout()
plt.savefig(os.getcwd() + '/Figures/L3_by_Score_Single.png', dpi=600, format='png', transparent=True)
g = sns.lmplot(data=dframe.loc[dframe['Experiment'].isin(['Control'])], x="Multi-Template LCDR3 RMSD vs Native", y="Method LCDR3 RMSD vs Native", hue="Experiment", fit_reg=False, height = 10, aspect = 1, x_estimator=np.median, x_ci=None, scatter_kws={'alpha':1, 's':75}, palette=['grey'], legend=False, legend_out=False)
g.set(ylim=(0, L3_lim), xlim=(0, L3_lim))
ax = plt.gca()
X_plot = np.linspace(0, L3_lim, 400)
Y_plot = X_plot
plt.plot(X_plot, Y_plot, color='silver')
plt.plot(X_plot, Y_plot - Margin, color='silver')
plt.plot(X_plot, Y_plot + Margin, color='silver')
#ax.set_title('Performance of Best Control Models')
plt.ylabel('Control LCDR3 RMSD vs Native', labelpad=20)
plt.xlabel('Multi-Template LCDR3 RMSD vs Native', labelpad=20)
plt.tight_layout()
plt.savefig(os.getcwd() + '/Figures/L3_by_Score_Control.png', dpi=600, format='png', transparent=True)



g = sns.lmplot(data=samplesframe, x="Sample Size", y="HCDR3 RMSD", hue="Experiment", fit_reg=False, height = 10, aspect = 1, scatter_kws={'alpha':1, 's':10}, legend_out=False)
legend = plt.legend(loc = 'lower right', frameon=True, framealpha=1, title = 'Method')
frame = legend.get_frame()
frame.set_facecolor('white')
frame.set_edgecolor('black')
for lh in legend.legendHandles:
	lh.set_alpha(1)
	lh.set_sizes([100])
ax = plt.gca()
#ax.set_title('Method Performance of Best Models by Sample Size')
plt.tight_layout()
plt.savefig(os.getcwd() + '/Figures/H3_By_Sample_Size.png', dpi=600, format='png', transparent=True)

g = sns.lmplot(data=samplesframe, x="Sample Size", y="Backbone RMSD", hue="Experiment", fit_reg=False, height = 10, aspect = 1, scatter_kws={'alpha':1, 's':10}, legend_out=False)
legend = plt.legend(loc = 'lower right', frameon=True, framealpha=1, title = 'Method')
frame = legend.get_frame()
frame.set_facecolor('white')
frame.set_edgecolor('black')
for lh in legend.legendHandles:
        lh.set_alpha(1)
        lh.set_sizes([100])
ax = plt.gca()
#ax.set_title('Method Performance of Best Models by Sample Size')
plt.tight_layout()
plt.savefig(os.getcwd() + '/Figures/Backbone_By_Sample_Size.png', dpi=600, format='png', transparent=True)
