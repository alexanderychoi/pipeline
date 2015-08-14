# -*-coding:Utf-8 -*

samfile = open('./Day-0_S1_L001_R1_001_noTA.sam', 'r')
qual = open('./qual.txt','w+')
while True:
	line = samfile.readline()
	if not line:
		break
	else:
		columns = line.split("\t")
		gene = columns[2]
		if gene != '*':
			AS_score = int(columns[11][5:])
			qual.write(str(AS_score))
			qual.write('\n')
