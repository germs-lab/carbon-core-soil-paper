import sys

l = []
for line in open(sys.argv[1]):
    l.append(line.rstrip())

for n, line in enumerate(open(sys.argv[2], 'rU')):
    if n > 2:
        dat = line.rstrip().split('\t')
        if dat[0] in l:
            print line, 
    else:
        print line,
