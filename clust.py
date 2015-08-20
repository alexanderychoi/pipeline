# -*-coding:Utf-8 -*

def to_coord(string,nucleotide_integer):
	coord = "".join(nucleotide_integer.get(base, base) for base in string)
	ints=[]
	ints.extend(coord)
	ints=map(int,ints)
	return ints

bc_file = open('./shortbc.txt','r')
nucleotide_integer = {'A': '1', 'T': '2', 'C': '3', 'G': '4'}
array=[]

while True:
	line = bc_file.readline()
	if not line:
		break
	else:
		coords = to_coord(line,complement)
		array.append(coords)