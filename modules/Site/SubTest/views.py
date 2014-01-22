# Create your views here.
from django.shortcuts import render
from django.conf import settings
from forms import SubTestForm
from models import SubTestResult
from django.core.paginator import Paginator, EmptyPage, PageNotAnInteger
from django.http import Http404,HttpResponseNotFound
import re
import glob
import os
from random import randint
dict={}
orgs=[]
org128=[]
org64=[]
orgd2=[]
orgi64=[]
orgi128=[]
cans=[]
log=[]
imagetype = 0
dist_type  = ('Micro Shift','Micro Rotate','Blur','Blocking','Light Change Const', 'Light Change Linear')
dist_suffix = ('shift','rotate','blur','block','light_const','light_linear')
ROOT_DIR = os.path.dirname(__file__)
texture_dir = ROOT_DIR+'/../../../data/totest'
print 'texture dir is '
print texture_dir
dbpath = {'coding64_inter':texture_dir+'/dist_coding_64_inter',
          'coding64_intra':texture_dir+'/dist_coding_64_intra',
          'corbis128_intra':texture_dir+'/dist_corbis_128_intra',
          'corbis64_intra':texture_dir+'/dist_corbis_64_intra',
          'corbis128_inter':texture_dir+'/dist_corbis_128_inter',
          'corbis64_inter':texture_dir+'/dist_corbis_64_inter',}
dbdata = {'coding64_inter':[],'coding64_intra':[],'corbis128_intra':[],'corbis64_intra':[],'corbis128_inter':[],'corbis64_inter':[],}
dbtype = {'coding64_inter':3,'coding64_intra':4,'corbis128_intra':0,'corbis64_intra':1,'corbis128_inter':5,'corbis64_inter':6}
savename = ('corbis128_intra','corbis64_intra','corbis_d2_intra','coding64_inter','coding64_intra','corbis128_inter','corbis64_inter')
#def show_inter(request,count):
#    global orgs
#    global org_idx_old
#    global dist_idx_old
#    global total
#    global imagetype
#def show_inter(request,count):
#    global orgclasses
#    global org_idx_old_A
#    global org_idx_old_B
#    global cls_idx_old
#    global dist_idx_old
#    global total
#    global imagetype
#    if int(count)==total: #finish
#        org_idx_oldA = randint(0,len(orgs)-1)
#        org_idx_oldB = randint(0,len(
#        dist_idx_old = randint(0,len(dist_type)-1)
#        return render(request,'test_finish.html',{'type':imagetype,})


def show_image(request,count):
    global orgs
    global org_idx_old
    global org_idx_old_A
    global org_idx_old_B
    global cls_idx_old
    global dist_idx_old
    global total
    global imagetype
    global dist_type
    total=30
    if request.session['count']==total: #finish
        org_idx_old = randint(0,len(orgs)-1)
        dist_idx_old = genDistIdx() #randint(0,len(dist_type)-1)
        return render(request,'SubTest/test_finish.html',{'type':imagetype,})
    if request.method == 'POST':
        for i in range(0,11):
            if  'button '+str(i) in request.POST:       
                thistest = SubTestResult()
                thistest.score =  i
                thistest.distortion = dist_type[dist_idx_old]
                if imagetype<3:
                    thistest.orgName = orgs[org_idx_old]+'_org.png'
                    thistest.candName = orgs[org_idx_old]+'_replace_'+dist_suffix[dist_idx_old]+'.png'
                else:
                    thistest.orgName = orgs[cls_idx_old][0]+"_"+orgs[cls_idx_old][1][org_idx_old_B]+'_org.png'
                    thistest.candName = orgs[cls_idx_old][0]+"_"+orgs[cls_idx_old][1][org_idx_old_A]+'_'+orgs[cls_idx_old][1][org_idx_old_B]+'_replace_'+dist_suffix[dist_idx_old]+'.png'
                thistest.imagetype=imagetype
                x_forwarded_for = request.META.get('HTTP_X_FORWARDED_FOR')
                if x_forwarded_for:
                    ip = x_forwarded_for.split(',')[-1].strip()
                else:
                    ip = request.META.get('REMOTE_ADDR')
                thistest.ip = ip
                thistest.dbname=request.session['dbname']
                thistest.save()
                request.session['count']+=1
                if imagetype<3:
                    org_idx = randint(0,len(orgs)-1)
                    org_idx_old = org_idx
                else:
                    cls_idx = randint(0,len(orgs)-1)
                    print orgs[cls_idx]
                    org_idx_A = randint(0,len(orgs[cls_idx][1])-1)
                    org_idx_B = randint(0,len(orgs[cls_idx][1])-1)
                    while org_idx_B == org_idx_A:
                        org_idx_B = randint(0,len(orgs[cls_idx][1])-1)
                    cls_idx_old = cls_idx
                    org_idx_old_A = org_idx_A
                    org_idx_old_B = org_idx_B
                dist_idx = genDistIdx() #randint(0,len(dist_type)-1)
                dist_idx_old = dist_idx
                #return redirct('subtest/'+str(count))
                break
        else:
            print 'no data input'
            if imagetype<3:
                org_idx = org_idx_old
            else:
                org_idx_A = org_idx_old_A
                org_idx_B = org_idx_old_B
                cls_idx = cls_idx_old_B
            dist_idx = dist_idx_old
        form = SubTestForm();
    else:
        if imagetype<3:
            org_idx = org_idx_old
        else:
            org_idx_A = org_idx_old_A
            org_idx_B = org_idx_old_B
            cls_idx = cls_idx_old
        dist_idx = dist_idx_old #randint(0,len(orgs))
        form = SubTestForm()
    for elem in dbtype:
        if dbtype[elem]==imagetype:
            prefix = dbpath[elem]
            break
    temppos = prefix.rfind('/')
    prefix = prefix[temppos+1:len(prefix)]+'/'
    #prefix = 'dist';
    #if imagetype == 1:
    #    prefix = prefix+'64/'
    #elif imagetype == 2:
    #    prefix = prefix+'_d2/'
    #elif imagetype == 3: 
    #    prefix = prefix+'_inter/'
    #else:
    #    prefix = prefix+'/'
    position = randint(0,1)#swap position randomly
    if imagetype >2:
        outfiles = [prefix+orgs[cls_idx][0]+'_'+str(org_idx_B)+'_org.png',prefix+orgs[cls_idx][0]+'_'+str(org_idx_A)+'_'+str(org_idx_B)+'_replace_'+dist_suffix[dist_idx]+'.png'] 
    else:
        outfiles = [prefix+orgs[org_idx]+'_org.png',prefix+orgs[org_idx]+'_replace_'+dist_suffix[dist_idx]+'.png'] 
    #there is a 10% possiblity do not use replacement but show the original
    chance = randint(0,100)
    if chance <10:
        outfiles[1]=outfiles[0]

    if position==1:
        outfiles[0],outfiles[1] = outfiles[1],outfiles[0]
    return render(request,'display_images.html',{'images':outfiles,'form':form,'count':request.session['count'],'total':total})

def startcompare(request,database,bsize,pair):
    kword = database+bsize+'_'+pair
    if (kword in dbpath):
        initTestData(database,bsize,pair)
        request.session['dbname']=str(database)+str(bsize)
        request.session['count']=0
        request.session['total']=30
        return render(request,'SubTest/test_begin1.html')
    else:
        #print 'print select one of valid database:'
        #print dbpath.keys()
        return HttpResponseNotFound('<h1> database not found </h1>')
    #<p> please select one of the valide database below </p> <p>{{dpath.keys()}}</p>')
    ##print path

def genDistIdx():
    temp = randint(0,len(dist_type)-1)
    if temp==2:
        temp = randint(0,len(dist_type)-1)
    while(temp==4): #if the result is LCC or blur, try again and to make the possiblity much smaller
        temp = randint(0,len(dist_type)-1)

    return temp

def initTestData(database,bsize,pair):
    global orgs
    global dbdata
    global dbpath
    global dbtype
    global cls_idx_old 
    global org_idx_old
    global org_idx_old_A
    global org_idx_old_B
    global dist_idx_old
    global total
    global count
    global size
    global class_name
    global imagetype
    kword = database+bsize+'_'+pair
    tempdata = dbdata[kword]
    path = dbpath[kword]
    total = 30
    count=1
    imagetype=dbtype[kword]

    if pair=='inter':
        if len(tempdata)==0:
            tempdict={}
            for file in os.listdir(path):
                if file.endswith("org.png"):          
                    print file
                    if database=='coding':
                        subnames = re.findall(r'(.*)_(\d+)_(\d+)',file)
                        if tempdict.get(subnames[0][0]+'_'+subnames[0][1])==None:
                            templist=[]
                        else:
                            templist=tempdict.get(subnames[0][0]+'_'+subnames[0][1])
                        templist.append(subnames[0][2])
                        tempdict[subnames[0][0]+'_'+subnames[0][1]] = templist
                    else:
                        subnames = re.findall(r'(.*)_(\d+)_org',file)
                        if tempdict.get(subnames[0][0])==None:
                            templist=[]
                        else:
                            templist=tempdict.get(subnames[0][0])
                        templist.append(subnames[0][1])
                        tempdict[subnames[0][0]]=templist
            tempdata = tempdict.items();
            cls_idx_old = randint(0,len(tempdata)-1)
            org_idx_old_A = randint(0,len(tempdata[cls_idx_old][1])-1)
            org_idx_old_B = randint(0,len(tempdata[cls_idx_old][1])-1)
            while org_idx_old_B == org_idx_old_A:
                org_idx_old_B = randint(0,len(tempdata[cls_idx_old][1])-1)
            dist_idx_old = genDistIdx() #randint(0,len(dist_type)-1)
    else:
        for file in os.listdir(path):
            print file
            if file.endswith("org.png"):          
                subname = file[0:file.find("_org.png")]
                org64.append(subname)
                #dict[count] = subname
                #count=count+1
                #now get the modified name for distorted image 
                #because in generating distortion, the file name format is changed
                org_idx_old = randint(0,len(org64)-1)
                dist_idx_old = genDistIdx()  #randint(0,len(dist_type)-1)
    orgs = tempdata
    dbdata[kword]=tempdata




#def startinter(request):
#    global orgs;
#    global orgi64;
#    global cls_idx_old 
#    global org_idx_old_A
#    global org_idx_old_B
#    global dist_idx_old
#    global total
#    global count
#    global size
#    global class_name
#    global imagetype
#    imagetype = 3
#    total=30
#    count=1     
#    path = 'H:/texture/data64/dist_inter/'
#    if len(orgi64)==0:
#        tempdict={}
#        for file in os.listdir(path):
#            if file.endswith("org.png"):          
#                print file
#                subnames = re.findall(r'(.*)_(\d+)_(\d+)',file)
#                if tempdict.get(subnames[0][0]+'_'+subnames[0][1])==None:
#                    templist=[]
#                else:
#                    templist=tempdict.get(subnames[0][0]+'_'+subnames[0][1])
#                templist.append(subnames[0][2])
#                tempdict[subnames[0][0]+'_'+subnames[0][1]] = templist
#        orgi64 = tempdict.items();
#        cls_idx_old = randint(0,len(orgi64)-1)
#        org_idx_old_A = randint(0,len(orgi64[cls_idx_old][1])-1)
#        org_idx_old_B = randint(0,len(orgi64[cls_idx_old][1])-1)
#        while org_idx_old_B == org_idx_old_A:
#             org_idx_old_B = randint(0,len(orgi64[cls_idx_old][1])-1)
#        dist_idx_old = randint(0,len(dist_type)-1)
#    orgs = orgi64
#    return render(request,'test_begin1.html')
#def startpage(request):
#    global orgs
#    global org128
#    global org_idx_old
#    global dist_idx_old
#    global total
#    global count
#    global size
#    global imagetype
#    imagetype = 0
#    total=30
#    count=1     
#    path = 'H:/texture/totest/dist/'
#    if len(org128)==0:
#        for file in os.listdir(path):
#            print file
#            if file.endswith("org.png"):          
#                subname = file[0:file.find("_org.png")]
#                org128.append(subname)
#                #dict[count] = subname
#                #count=count+1
#                #now get the modified name for distorted image 
#                #because in generating distortion, the file name format is changed
#                org_idx_old = randint(0,len(org128)-1)
#                dist_idx_old = randint(0,len(dist_type)-1)
#    orgs = org128
#    return render(request,'test_begin1.html')

#def startpage64(request):
#    global orgs
#    global org64
#    global org_idx_old
#    global dist_idx_old
#    global total
#    global count
#    global imagetype
#    imagetype = 1
#    total=30
#    count=1     
#    path = 'H:/texture/totest/dist64/'
#    if len(org64)==0:
#        for file in os.listdir(path):
#            print file
#            if file.endswith("org.png"):          
#                subname = file[0:file.find("_org.png")]
#                org64.append(subname)
#                #dict[count] = subname
#                #count=count+1
#                #now get the modified name for distorted image 
#                #because in generating distortion, the file name format is changed
#                org_idx_old = randint(0,len(org64)-1)
#                dist_idx_old = randint(0,len(dist_type)-1)
#    orgs = org64
#    return render(request,'test_begin1.html')

#def startpaged2(request):
#    global orgs
#    global orgd2
#    global org_idx_old
#    global dist_idx_old
#    global total
#    global count
#    global imagetype
#    imagetype=2
#    total=30
#    count=1     
#    path = 'H:/texture/totest/dist_d2/'
#    if len(orgd2)==0:
#        for file in os.listdir(path):
#            print file
#            if file.endswith("org.png"):          
#                subname = file[0:file.find("_org.png")]
#                orgd2.append(subname)
#                #dict[count] = subname
#                #count=count+1
#                #now get the modified name for distorted image 
#                #because in generating distortion, the file name format is changed
#                org_idx_old = randint(0,len(orgd2)-1)
#                dist_idx_old = randint(0,len(dist_type)-1)
#    orgs = orgd2
#    return render(request,'test_begin1.html')

testc=0
def threadtest(request):
    global testc
    if request.method == 'POST':
        testc=testc+1
        print testc
    return render(request,'threadtest.html',{'count':testc,})

def output(request):
    global savename
    out = SubTestResult.objects.all()
    prefix = 'subtestoutput_'
    f=[file(prefix+savename[0]+".txt",'w'),file(prefix+savename[1]+".txt",'w'),file(prefix+savename[2]+".txt",'w'),file(prefix+savename[3]+".txt",'w'),file(prefix+savename[4]+".txt",'w'),file(prefix+savename[5]+".txt",'w'),file(prefix+savename[6]+".txt",'w')]
    
    #write to file
    for item in out:
        mystr =  item.orgName+", "+item.candName+", "+str(item.score)+", "+item.distortion+str(item.create_date)+";\n"
        f[item.imagetype].write(mystr)
    for fp in f:
        fp.close()
    return render(request,'SubTest/output.html',{'data':out,})
