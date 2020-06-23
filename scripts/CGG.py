## count how many arginine - CGG in P/A sites

import sys

f = open(sys.argv[1], 'r')

eCGG = 0
pCGG = 0
epCGG = 0
for line in f:
	seq = line.strip()
	e = seq[12:15]
	p = seq[15:18]
	if e[0:3] == 'cgg' and p[0:3] == 'cgg':
		epCGG += 1
	if e[0:3] == 'cgg':
		eCGG += 1
	if p[0:3] == 'cgg':
		pCGG += 1

f.close()

print "e-site CGG:"
print eCGG
print "p-site CGG:"
print pCGG
print "e- and p-site CGG:"
print epCGG

