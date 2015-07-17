# -*-coding:Utf-8 -*
from collections import defaultdict
import random

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
matrix = [[random.randint(0,seq_length*2) for x in range(random_bits_per_mask)] for x in range(number_of_masks)]

bit_seq_list = [bit_seq1, bit_seq1copie, bit_seq1false, bit_seq2]

for i in range(number_of_masks):
	print matrix[i]

for j in matrix[0]:
	print matrix[0][j]

'''
for bit_sequence in bit_seq_list
	for i in range(number_of_masks):
		for j in
'''






seq_dictionary = defaultdict(list)