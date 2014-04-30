#!/usr/bin/env python
# -*- coding: utf-8 -*-
import subprocess

print 'before open'
p = subprocess.Popen("/home/guoxin/Projects/MTC/Debug/bin/testPythonCall",stdin=subprocess.PIPE, stdout=subprocess.PIPE);
while True:
#for x in xrange(11):
    outline = p.stdout.readline()
    print outline
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
    else:
        print 'error'
        continue
p.stdin.close()
p.wait()
output = p.stdout.read()
print 'after open'
print output
