## count how many prolines in P/A sites

import sys

f = open(sys.argv[1], 'r')

#ePro = 0
pPro = 0
aPro = 0
#epPro = 0
apPro = 0
for line in f:
	seq = line.strip()
	#e = seq[12:15]
	p = seq[15:18]
	a = seq[18:21]
	#if e[0:2] == 'cc' and p[0:2] == 'cc':
	#	epPro += 1
	if a[0:2] == 'cc' and p[0:2] == 'cc':
		apPro += 1
	#if e[0:2] == 'cc':
	#	ePro += 1
	if p[0:2] == 'cc':
		pPro += 1
	if a[0:2] == 'cc':
		aPro += 1

f.close()

#print "e-site Prolines:"
#print ePro
print "p-site Prolines:"
print pPro
print "a-site Prolines:"
print aPro
#print "e- and p-site Prolines:"
#print epPro
print "a- and p-site Prolines:"
print apPro

aORp = aPro + pPro

print "---------------------"
print "a- or p-site Prolines"
print aORp
print "---------------------"


