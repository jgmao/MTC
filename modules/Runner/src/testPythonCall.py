#!/usr/bin/env python
# -*- coding: utf-8 -*-
import subprocess

def step1():
    global p
    print 'before open'
    p = subprocess.Popen("/home/guoxin/Projects/MTC/Debug/bin/testPythonCall",stdin=subprocess.PIPE, stdout=subprocess.PIPE)
    print p    

def step2(): 
    global p
    print p
    outline = p.stdout.readline()
    print outline
    return outline


p = None

step1()
while True:
#for x in xrange(11):
    #outline = p.stdout.readline()
    #print outline
    outline = step2()
    print outline=='waiting input:\n'
    if outline=='waiting input:\n':
        y = raw_input("input to python: ")
    #print y
        line = '%d\n' % int(y)
    #line = line + '\n'
    #line = '%d\n' % x
        try:
            p.stdin.write(line)
        except IOError as e:
            if e.errno == errno.EPIPE or e.errno == errno.EINVAL:
            # Stop loop on "Invalid pipe" or "Invalid argument".
            # No sense in continuing with broken pipe.
              break
            else:
            # Raise any other error.
              raise
    elif outline=='done!\n':
        break
    else:
        print 'error'
        continue
p.stdin.close()
p.wait()
output = p.stdout.read()
print 'after open'
print output

