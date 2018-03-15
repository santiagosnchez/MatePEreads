import sys
import gzip
from itertools import izip
import os

basefile = sys.argv[1]
pair1 = basefile+"_R1.fastq.gz"
pair2 = basefile+"_R2.fastq.gz"
size1 = os.stat(pair1).st_size
size2 = os.stat(pair2).st_size
ln = 0
c = 0
lookup1 = {}
lookup2 = {}
with gzip.open(basefile+"_paired_R1"+".fastq.gz","w") as p1_fqgz, gzip.open(basefile+"_paired_R2"+".fastq.gz","w") as p2_fqgz:
	with gzip.open(pair1,"r") as p1, gzip.open(pair2,"r") as p2:
		for line_p1,line_p2 in izip(p1,p2):
			ln += 1
			if ln % 4 == 1:
				h1 = line_p1.split(" ")
				h2 = line_p2.split(" ")
				lookup1[h1[0]] = ['','',h1[1]]
				lookup2[h2[0]] = ['','',h2[1]]
			elif ln % 4 == 2:
				lookup1[h1[0]][0] = line_p1
				lookup2[h2[0]][0] = line_p2
			elif ln % 4 == 0:
				lookup1[h1[0]][1] = line_p1
				lookup2[h2[0]][1] = line_p2
				if h1[0] != h2[0]:
					if lookup1.get(h1[0]) and lookup2.get(h1[0]):
						p1_fqgz.write(' '.join(h1))
						p1_fqgz.write(lookup1[h1[0]][0]+"+\n"+lookup1[h1[0]][1])
						p2_fqgz.write(' '.join([h1[0],lookup2[h1[0]][2]]))
						p2_fqgz.write(lookup2[h1[0]][0]+"+\n"+lookup2[h1[0]][1])
						del lookup1[h1[0]]
						del lookup2[h1[0]]
					elif lookup1.get(h2[0]) and lookup2.get(h2[0]):
						p2_fqgz.write(' '.join(h2))
						p2_fqgz.write(lookup2[h2[0]][0]+"+\n"+lookup2[h2[0]][1])
						p1_fqgz.write(' '.join([h2[0],lookup1[h2[0]][2]]))
						p1_fqgz.write(lookup1[h2[0]][0]+"+\n"+lookup1[h2[0]][1])
						del lookup1[h2[0]]
						del lookup2[h2[0]]
				else:
					p1_fqgz.write(' '.join(h1))
					p1_fqgz.write(lookup1[h1[0]][0]+"+\n"+lookup1[h1[0]][1])
					p2_fqgz.write(' '.join(h2))
					p2_fqgz.write(lookup2[h2[0]][0]+"+\n"+lookup2[h2[0]][1])
					del lookup1[h1[0]]
					del lookup2[h2[0]]
			sys.stdout.write(" {:<18s}done for {} lines\r".format("[MatePEreads.py]",ln)),
                        sys.stdout.flush()
	if size1 > size2:
		incomplete,pair = pair1,1
	else:
		incomplete,pair = pair2,2
	with os.popen('gzcat {} | tail -n +{}'.format(incomplete,ln+1)) as leftover:
		for line in leftover:
			ln += 1
			if ln % 4 == 1:
				h = line.split(" ")
				if pair == 1:
					lookup1[h[0]] = ['','',h[1]]
				elif pair == 2:
					lookup2[h[0]] = ['','',h[1]]
			elif ln % 4 == 2:
				if pair == 1:
					lookup1[h[0]][0] = line
				elif pair == 2:
					lookup2[h[0]][0] = line
			elif ln % 4 == 0:
				if pair == 1:
					lookup1[h[0]][1] = line
					if lookup2.get(h[0]):
						p1_fqgz.write(' '.join(h))
						p1_fqgz.write(lookup1[h[0]][0]+"+\n"+lookup1[h[0]][1])
						p2_fqgz.write(' '.join([h[0],lookup2[h[0]][2]]))
						p2_fqgz.write(lookup2[h[0]][0]+"+\n"+lookup2[h[0]][1])
						del lookup1[h[0]]
						del lookup2[h[0]]
				elif pair == 2:
					lookup2[h[0]][1] = line
					if lookup1.get(h[0]):
						p1_fqgz.write(' '.join([h[0],lookup1[h[0]][2]]))
						p1_fqgz.write(lookup1[h[0]][0]+"+\n"+lookup1[h[0]][1])
						p2_fqgz.write(' '.join(h))
						p2_fqgz.write(lookup2[h[0]][0]+"+\n"+lookup2[h[0]][1])
						del lookup1[h[0]]
						del lookup2[h[0]]
				sys.stdout.write(" {:<18s}done for {} lines\r".format("[MatePEreads.py]",ln)),
				sys.stdout.flush()
with gzip.open(basefile+"_unpaired"+".fastq.gz","w") as unp:
	sys.stdout.write("\n {:<18s}writing unpaired reads".format("[MatePEreads.py]"))
	for key in lookup1.keys():
		unp.write(' '.join([key,lookup1[key][2]]))
		unp.write(lookup1[key][0]+"+\n"+lookup1[key][1])
                c += 1
	for key in lookup2.keys():
		unp.write(' '.join([key,lookup2[key][2]]))
		unp.write(lookup2[key][0]+"+\n"+lookup2[key][1])
                c += 1
        sys.stdout.write("\n {:<18}{} reads were unpaired\n".format("[MatePEreads.py]",c))

