#!/usr/bin/python
import sys

filename = sys.argv[1]
f = open(filename, 'r')
line_num = 100
step = 0.2
for line in f.readlines():
	print line_num, line,
	line_num += step
