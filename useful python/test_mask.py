# -*-coding:Utf-8 -*
from collections import defaultdict
import random
from pympler import asizeof

seq_length = 50

complement = {'A': '00', 'C': '01', 'T': '10', 'G': '11'}

seq1 = "TCTTCCTGCCCAGTAGCGGATGATAATCGTTGTTGCCAGCCGGCGTGGAA"
umi1 = "TAGAC"
bit_seq1 = "".join(complement.get(base, base) for base in seq1)
print bit_seq1

seq1copie = "TCTTCCTGCCCAGTAGCGGATGATAATCGTTGTTGCCAGCCGGCGTGGAA"
umi1copie = "TAGAC"
bit_seq1copie = "".join(complement.get(base, base) for base in seq1copie)
print bit_seq1copie

seq1false = "TCTTCCTGCCCAGTAGCGGATGATAATCGTTGTTGCCAGCCGGCGTGGAA"
umi1false = "TAAAC"
bit_seq1false = "".join(complement.get(base, base) for base in seq1false)
print bit_seq1false

seq2 = "TCTTCCTGCCCAGTAGTGGATGATAATCGTTGTTGCCAGCCGGCGTGGCA"
umi2 = "TAGAC"
bit_seq2 = "".join(complement.get(base, base) for base in seq2)
print bit_seq2

number_of_masks = 30
random_bits_per_mask = 25
random.seed()
matrix = [[random.randint(0,(seq_length*2)-1) for x in range(random_bits_per_mask)] for x in range(number_of_masks)]

bit_seq_list = [bit_seq1, bit_seq1copie, bit_seq1false, bit_seq2]
umi_list = [umi1, umi1copie, umi1false, umi2]

for i in range(number_of_masks):
	print matrix[i]

seq_dictionary = defaultdict(list)
for bit_sequence,umi in zip(bit_seq_list,umi_list):
	bit_code = ''
	for i in range(number_of_masks):
		bit_list = []
		for j in range(random_bits_per_mask):
			bit_list.append(bit_sequence[matrix[i][j]])
		print bit_list
		sum_bits = 0
		for bit in bit_list:
			sum_bits=sum_bits^int(bit)
		print sum_bits
		bit_code=bit_code+str(sum_bits)
	print bit_code
	int_code = int(bit_code,2)
	print int_code
	if int_code in seq_dictionary:
		if umi not in seq_dictionary[int_code]:
			seq_dictionary[int_code].append(umi)
	else:
		seq_dictionary[int_code].append(umi)
print seq_dictionary
print asizeof.asizeof(seq_dictionary)
seq_dictionary.update(dict(seq_dictionary))
print asizeof.asizeof(seq_dictionary)
