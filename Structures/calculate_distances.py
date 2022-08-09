from pymol import cmd
from itertools import combinations
import os

def calculate_distances(residues):
	res = residues.split('+')
	y = combinations(list(res), 2)
	i = 0
	for combo in y:
		x = cmd.get_distance("".join(["i. ", combo[0], " and n. CA"]), "".join(["i. ", combo[1], " and n. CA"]))
		#print(combo, x)
		print(",".join([str(combo[0]), str(combo[1]), str(x)]), file=open("output.txt", "a"))

cmd.extend('calculate_distances', calculate_distances)