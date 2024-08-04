#!/usr/bin/python3

__author__ = "Giorgia Del Missier"

import sys
import numpy as np


def parsePDB(pdbfile):

	fpdb = open(pdbfile)
	chains_dict, hemo_dict = {}, {}

	for line in fpdb:
		if line[0:4] == 'ATOM': 
			x, y, z = float(line[30:38]), float(line[38:46]), float(line[46:54])
			at = line[12:16].strip()
			ch = line[21]
			res = line[17:20]
			seq_n = line[22:26].strip()
			coord = [x,y,z]

			chains_dict[ch] = chains_dict.get(ch, {})
			chains_dict[ch][(res, seq_n)] = chains_dict[ch].get((res, seq_n), {})
			chains_dict[ch][(res, seq_n)][at] = coord

		elif line[0:6] == 'HETATM':
			res = line[17:20]
			if res == 'HEM' or res == 'OXY':
				x, y, z = float(line[30:38]), float(line[38:46]), float(line[46:54])
				at = line[12:16].strip()
				ch = line[21]
				seq_n = line[22:26].strip()
				coord = [x,y,z]

				hemo_dict[ch] = hemo_dict.get(ch, {})
				hemo_dict[ch][(res, seq_n)] = hemo_dict[ch].get((res, seq_n), {})
				hemo_dict[ch][(res, seq_n)][at] = coord

	return chains_dict, hemo_dict


def computeDistance(at1, at2):
	distance = np.sqrt((at1[0]-at2[0])**2 + (at1[1]-at2[1])**2 + (at1[2]-at2[2])**2)
	return distance


def storeDistances(hemo, chain, minimum_dist = 3.5):
	all_dists = {}

	for hres in hemo:
		for hat in hemo[hres]:
			for cres in chain:
				for cat in chain[cres]:
					dist = computeDistance(hemo[hres][hat], chain[cres][cat])
					if dist <= minimum_dist:
						all_dists[hres] = all_dists.get(hres, {})
						all_dists[hres][cres] = all_dists[hres].get(cres, {})
						all_dists[hres][cres][hat] = (dist, cat)

	return all_dists


def printAllInteractions(distancedict, chain):								#to be formatted
	print("CHAIN\tRESIDUE\tHETERO\tATOMS(<=3.5A)\tDISTANCE")
	for hetero in distancedict:
		for residue in distancedict[hetero]:
			for hat in distancedict[hetero][residue]:
				cat = distancedict[hetero][residue][hat][1]
				dist = round(distancedict[hetero][residue][hat][0], 6)
				print(chain + "\t" + residue[0]+residue[1] + "\t" + hetero[0]+hetero[1] + "\t" + cat + "-" + hat + "  \t" + str(dist))


def printBestInteractions(distancedict, chain):
	print("CHAIN\tRESIDUE\tHETERO\tBEST DIST\tATOMS")
	for hetero in distancedict:
		for residue in distancedict[hetero]:
			atoms, distances = '', []
			for hat in distancedict[hetero][residue]:
				cat = distancedict[hetero][residue][hat][1]
				atoms += cat + "-" + hat + ", "
				distances.append(round(distancedict[hetero][residue][hat][0],6))
			print(chain + "\t" + residue[0]+residue[1] + "\t" + hetero[0]+hetero[1] + "\t" + str(min(distances)) + "  \t" + atoms)



if __name__ == '__main__':

	if len(sys.argv) == 2:
		pdbfile = sys.argv[1]
		chains_dict, hemo_dict = parsePDB(pdbfile)

		for chain in 'ABCD':
			distances = storeDistances(hemo_dict[chain], chains_dict[chain])
			printAllInteractions(distances, chain)
			print('\n')
			printBestInteractions(distances, chain)
			print('\n')
	
	else:	
		print('parsepdb.py pdbfile')

