import re

string = 'ATTTTTTTTTTTAAAAAAAAAAAAAAAAATTTTTTTTTTTTGHTGMJTTTTMKGKAAAAAAAAGMTKKTTTTAAAAAAAAAAAAAAAAAAAAAATTTTTTTTTTTTTTTTTTTTTTLLLLVGG'
error_str_t = 'TTTTTTTTTT'
error_str_a = 'AAAAAAAAAA'

index_t =  [m.start() for m in re.finditer('(?='+error_str_t+')', string)]
reversed_index_t = list(reversed(index_t))
print reversed_index_t
for x in reversed_index_t:
	while string[x] == 'T':
		string = string[:x] + string[x+1:]
	
print string

index_a =  [m.start() for m in re.finditer('(?='+error_str_a+')', string)]
reversed_index_a = list(reversed(index_a))
print reversed_index_a
for x in reversed_index_a:
	while string[x] == 'A':
		string = string[:x] + string[x+1:]
	
print string