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
      print 'starting coder...'
      p = subprocess.Popen([os.path.join(settings.ROOT_DIR, '../../temprun/Release/bin/mtcmain'),os.path.join(settings.ROOT_DIR, '../../temprun/encoder_configure.cfg')],stdin=subprocess.PIPE, stdout=subprocess.PIPE)
      #p = subprocess.Popen([os.path.join(settings.ROOT_DIR, '../../Debug/bin/testPythonCall'),os.path.join(settings.ROOT_DIR, '../../data/baboon.pgm')],stdin=subprocess.PIPE, stdout=subprocess.PIPE)
      print 'coder started.'
    #wait until get subjective test ready respone from mtc coder
    while True:
      outline = p.stdout.readline()
      print 'now printing: '+outline
      #raw_input("Press Enter to continue...")

      #start from sentence 'Subject Select:'
      #the org is printed in the first line
      #the localview is in the second line
      #and each candidate name are printed in each following lines
      #end with 'End of Candidate.'
      if outline=='Subject Select:\n':
        flag = 0
        request.session['org']=[]
        request.session['cands']=[]
        request.session['localview']=[]
        request.session['localviewplus']=[]
        request.session['localorg']=[]
        while outline:
          outline = p.stdout.readline()
          #print outline
          #print 'flag: '+str(flag)
          #raw_input("Press Enter to continue...")
          if outline != 'End of Candidate.\n':
            regex = re.compile('(/everything/.*.png)')
            r = regex.search(outline)
            temp = r.groups()
            name = temp[0]
            if flag is 0:
              request.session['org'] = name
            elif flag is 1:
              request.session['localview']= name
              request.session['localorg']=name[:11]+'/localorg'+name[17:]
            else:
              request.session['cands'].append(name)
              request.session['localviewplus'].append(name[:11]+'/localview'+name[20:])
            flag=flag+1
          else:
            print "Reading list of candidates end"
            break
        #test input
        #line = '%d\n' % int(3)
        #p.stdin.write(line)
        #outline = p.stdout.readline()
        #print outline
        print request.session['org']
        print request.session['localview']
        print request.session['localviewplus']
        print request.session['cands']
        #raw_input("Press Enter to continue...")
        noneitem = '/data/none'+str(request.session['initsize'])+'.png'

        return render(request, 'SubCoding/Coding.html',{'org':request.session['org'] , 'cands':request.session['cands'] , 'localorg':request.session['localorg'],'localview':request.session['localview'],'localviewplus':request.session['localviewplus'], 'noneitem':noneitem})
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
    x_forwarded_for = request.META.get('HTTP_X_FORWARDED_FOR')
    if x_forwarded_for:
        ip = x_forwarded_for.split(',')[-1].strip()
    else:
        ip = request.META.get('REMOTE_ADDR')
    localcount=0;
    print request.session['cands']
    request.session['accepted']=-1
    for cand in request.session['cands'] :
        thistest = SubCodingResult()
        thistest.orgName = request.session['org']
        thistest.ip = ip
        thistest.candName = cand
        if cand in request.POST:
           request.session['accepted']  = localcount;
           print 'FOUND '+ str(localcount)
           thistest.score = 10
        else:
           localcount=localcount+1
           thistest.score = 0
        thistest.save()
    print p
    print "THE SELECT IS " + str(request.session['accepted'])
    line = '%d\n' % request.session['accepted']
    #raw_input("Press Enter to continue...")
    try:
      p.stdin.write(line)
      request.session['count']=request.session['count']+1
      outline = p.stdout.readline()
      print outline
    except IOError as e:
      if not e.errno == errno.EPIPE or e.errno == errno.EINVAL:
        raise
    #raw_input("Press Enter to continue...")
    return codingprocess(request)

