# Create your views here.
from django.shortcuts import render
from django.conf import settings
from models import SubCodingResult
from django.http import Http404,HttpResponseNotFound
import re, glob, os
#from django.conf.settings import PROJECT_ROOT
import subprocess
from dajax.core import Dajax
import threading

#threadlocal = threading.local()
p= None
def startsession(request):
  global p
  cfg_file = open(os.path.join(settings.ROOT_DIR    ,'../../encoder_configure.cfg'))
  print('done========\n')
  print cfg_file
  #parse file
  for line in cfg_file:
    if line[0]=='#':
        continue
    if 'ImageName' in line:
        regex = re.compile('/([A-Za-z0-9-_]+)\.(.+)$') #match filename.ext
        r = regex.search(line)
        temp =  r.groups()
        #print r.group(0)
        #print r.group(1)
        print str(temp[0])+'.'+str(temp[1])
        #print r.group(3)
        request.session['imagename']=str(temp[0])#+'.'+str(temp[1]);#use regex to get file name xxx.xxx$
    if 'InitBlockSize' in line:
        print line
        regex = re.compile('\((\d+),(\d+),(\d+)\)')
        r = regex.search(line)
        temp = r.groups()
        request.session['initsize']=temp[0]
    #get block size
  #request.session['dbname']=str(database)+str(bsize)
  request.session['count']=0
  request.session['accepted']=-1
  #request.session['total']=30
#  threadlocal.p =[]
#  threadlocal.org=[]
#  threadlocal.cands=[]
  return render(request,'SubCoding/test_begin1.html')

def codingprocess(request):
    global p
    if 'count' not in request.session: #invalid entry
      return startsession(request)
    elif request.session['count']==0: #first step run the coder
      p = subprocess.Popen(os.path.join(settings.ROOT_DIR, '../../Debug/bin/testPythonCall'),stdin=subprocess.PIPE, stdout=subprocess.PIPE)
    #wait until get subjective test ready respone from mtc coder
    while True:
      outline = p.stdout.readline()
      print outline
      if outline=='waiting input:\n':
        #display images , make selection and send response to mtc
        #let mtc send the file names and just use it in to the template
        request.session['org'] = 'everything/org_64_(256_384).png'
        request.session['cands'] = ['everything/cand_plc_64_(256_384)_(86_247)_(3)_(2).png',
                 'everything/cand_plc_64_(256_384)_(125_247)_(4)_(2).png',
                 'everything/cand_plc_64_(256_384)_(137_247)_(5)_(2).png',
                 'everything/cand_plc_64_(256_384)_(148_247)_(6)_(2).png']
        localview = 'everything/seam_rst_0335.png'
        #test input
        #line = '%d\n' % int(3)
        #p.stdin.write(line)
        #outline = p.stdout.readline()
        #print outline
        print p
        return render(request, 'SubCoding/Coding.html',{'org':request.session['org'] , 'cands':request.session['cands'] , 'localview':localview})
#              #y = raw_input("input to python:")
#              #line = '%d\n' % int(y)
#              line = '%d\n' % int(request.session['count'])
#              request.session['count']= request.session['count']+1
#              try:
#                  p.stdin.write(line)
#              except IOError as e:
#                  if e.errno == errno.EPIPE or e.errno == errno.EINVAL:
#                  # Stop loop on "Invalid pipe" or "Invalid argument".
#                  # No sense in continuing with broken pipe.
#                    break
#                  else:
#                  # Raise any other error.
#                    raise
#          elif outline=='done!\n':
#               break
      else:
               #print 'error'
        continue
      p.stdin.close()
      p.wait()
      output = p.stdout.read()
      print 'after open'
      print output
      #render finish page
      return render(request, 'SubCoding/finish_session.html',{'output':output})

def cand_selected(request):
    #find out which candidate is selected
    global p
    localcount=0;
    print request.session['cands']
    for cand in request.session['cands'] :
        if cand in request.POST:
           request.session['accepted']  = localcount;
           print 'FOUND '+ str(localcount)
           break
        else:
           localcount=localcount+1
    line = '%d\n' % int(localcount)
    print p
    print "THE SELECT IS " +line
    try:
      p.stdin.write(line)
      outline = p.stdout.readline()
      print outline
    except IOError as e:
      if not e.errno == errno.EPIPE or e.errno == errno.EINVAL:
        raise
    return codingprocess(request)

