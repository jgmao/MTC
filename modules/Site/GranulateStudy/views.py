# Create your views here.
from django.shortcuts import render
from django.conf import settings
from GranulateStudy.models import Gran
from GranulateStudy.forms import GranTestForm
import re
import glob
import os
from random import randint

orgs=[]
ROOT_DIR = os.path.dirname(__file__)
texture_dir = ROOT_DIR+'/../../../data/totest'
print 'texture dir is '+texture_dir
def initTestData(request):
  global orgs
  if len(orgs)==0: 
    fcount=0
    path = texture_dir+'/gran'
    print 'gran path is '+path
    for file in os.listdir(path):
      if file.endswith("org.png"):
        print file
        orgs.append(file)
        fcount=fcount+1
  request.session['bsize']=3 #always start from 32x32
  request.session['org_idx'] = randint(0,len(orgs)-1)
  request.session['total']=25
  request.session['count']=0
  request.session['pre_select']='NA'
  print "Num of images: "+str(len(orgs))+'or '+str(fcount)
  return render(request,'GranulateStudy/test_begin2.html')

def show_image(request):
  global orgs
  if request.session['count']==request.session['total']:
    request.session['org_idx'] = randint(0,len(orgs)-1)
    request.session['count']=0
    request.session['bsize']=3 #always start from 32x32
    return render(request, 'GranulateStudy/test_finish.html')
  if request.method == 'POST':
    cur_select = 'NA'
    if 'button_YES' in request.POST: #if select YES it is a texture
      request.session['bsize'] -=1 #decrease size
      cur_select='YES'
    elif 'button_NO' in request.POST:#if select it is not a texture
      request.session['bsize'] +=1 #increase size
      cur_select = 'NO'
    else:
      print 'error neither yes or no selected'

    thistest = Gran()
    thistest.orgName = orgs[request.session['org_idx']]
    nextflag = True
    #determine status
    #if previous select yes and now select know.
    if (request.session['pre_select'] is 'YES' and cur_select is 'NO'):
        thistest.gran = request.session['bsize']+1
        thistest.save()
    elif (request.session['pre_select'] is 'NO' and cur_select is 'YES'):
        thistest.gran = request.session['bsize']-1
        thistest.save()
    else:
        request.session['pre_select']=cur_select
        nextflag = False

    if nextflag is True:
       request.session['org_idx']=randint(0,len(orgs)-1)
       request.session['bsize'] = 3
       request.session['count'] += 1

  else:
    form = GranTestForm()
    #get the image path

  prefix = texture_dir+'/gran/'
  curname = orgs[request.session['org_idx']]
  temppos = curname.rfind('_org.png')
  curname = curname[0:temppos]
  cursize = Gran.SIZE[request.session['bsize']][1]
  outfiles = [prefix+orgs[request.session['org_idx']], prefix+curname+'_g'+ str(cursize)+'.png']
  #outfiles = [prefix+orgs[request.session['org_idx']],prefix+orgs[request.session['org_idx']]+Gran.SIZE[request.session['bsize'][1]]
  print outfiles

  return render(request, 'GranulateStudy/show_images.html',{'form':form,'images':outfiles})
