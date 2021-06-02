import os, fnmatch, pickle, operator, math, re
import random
from Bio import pairwise2 as pw2
import matplotlib
matplotlib.use("agg")
import matplotlib.pyplot as plt
import pandas as pd
from pandas import DataFrame
import seaborn as sns
import numpy as np

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
				CDRs = IdentifyCDRs(chains[1], chains[0])
				CDRs.update( Extract_FR_CDR_Sequences(**CDRs) )
