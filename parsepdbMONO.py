#!/usr/bin/python3

__author__ = "Giorgia Del Missier"

import sys
import numpy as np


def parsePDB(pdbfile):

	fpdb = open(pdbfile)
	chains_dict = {}

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

	return chains_dict


def computeDistance(at1, at2):
	distance = np.sqrt((at1[0]-at2[0])**2 + (at1[1]-at2[1])**2 + (at1[2]-at2[2])**2)
	return distance


def storeDistances(chain1, chain2, minimum_dist):
	all_dists = {}

	for res1 in chain1:
		for at1 in chain1[res1]:
			for res2 in chain2:
				for at2 in chain2[res2]:
					dist = computeDistance(chain1[res1][at1], chain2[res2][at2])
					if dist <= minimum_dist:
						all_dists[res1] = all_dists.get(res1, {})
						all_dists[res1][res2] = all_dists[res1].get(res2, {})
						all_dists[res1][res2][at1] = (dist, at2)

	return all_dists


def printAllInteractions(distancedict, chain1, chain2):
	print("CHAIN1\tRES1\tCHAIN2\tRES2\tATOMS(<=3.5A)\tDISTANCE")
	for res1 in distancedict:
		for res2 in distancedict[res1]:
			for at1 in distancedict[res1][res2]:
				at2 = distancedict[res1][res2][at1][1]
				dist = round(distancedict[res1][res2][at1][0],6)
				print(chain1 + "\t" + res1[0]+res1[1] + "\t" + chain2 + "\t" + res2[0]+res2[1] + "\t" + at1 + "-" + at2 + "\t" + str(dist))


def printBestInteractions(distancedict, chain1, chain2):
	print("CHAIN1\tRES1\tCHAIN2\tRES2\tBEST DIST\tATOMS")
	for res1 in distancedict:
		for res2 in distancedict[res1]:
			atoms, distances = '', []
			for at1 in distancedict[res1][res2]:
				at2 = distancedict[res1][res2][at1][1]
				atoms += at1 + "-" + at2 + ", "
				distances.append(round(distancedict[res1][res2][at1][0],6))
			print(chain1 + "\t" + res1[0]+res1[1] + "\t" + chain2 + "\t" + res2[0]+res2[1] + "\t" + str(min(distances)) + "  \t" + atoms)



if __name__ == '__main__':

	if len(sys.argv) == 3:
		pdbfile = sys.argv[1]
		distance = sys.argv[2]
		chains_dict = parsePDB(pdbfile)

		for chain1 in 'ABCD':
			for chain2 in 'ABCD':
				if chain1 != chain2:
					distances = storeDistances(chains_dict[chain1], chains_dict[chain2], float(distance))
					printAllInteractions(distances, chain1, chain2)
					print('\n')
					printBestInteractions(distances, chain1, chain2)
					print('\n')

	else:	
		print('parsepdb2.py pdbfile distance')

		