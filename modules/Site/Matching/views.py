# Create your views here.
from django.shortcuts import render
from django.conf import settings
from Matching.forms import MatchForm
from django.core.paginator import Paginator, EmptyPage, PageNotAnInteger
import re
import glob
import os
class Rec(object):
    def __init__(self,x,y,size):
        self.x=x
        self.y=y
        self.size=size
    def __str__(self):
        return '('+str(self.x)+','+str(self.y)+')_'+str(self.size)
def read_log(filename):
    f = open(filename,'r')
    for line in f:
        print line

def show_image(request):
    dict={}
    orgs16=[]
    orgs32=[]
    cans=[]
    log=[]
    for file in os.listdir(settings.ROOT_DIR+"../Tester/everything"):
        if file.endswith(".txt"):
                log.append(file)
        else:
            regex = re.compile('_(\d+)_') #find size
            r = regex.search(file)
            size = r.group(1)
            g = re.findall(r'((?<=(?:[(]))[^(].*?(?=[)]))',file) 
            mykey = g[0]+'_'+str(size)
            if mykey in dict:
                if file.startswith("org"):
                    temp =dict[mykey]
                    temp.insert(0,file) 
                elif file.startswith("cand"):
                    temp = dict[mykey]
                    temp.append(file)
            else:
                dict[mykey]=[file]
    
    return render(request,'show_image.html',{'item':dict.items[0],})

def matching(request,x,y,size):
    global files
    global positions
    global scores
    global plcs
    score_in_match={}
    file_in_match={}
    pos_in_match={}
    plc_in_match={}
    light_in_match={}
    global out_files
    global out_scores
    global out_plcs
    global out_pos
    global out_light
    jumpto = "TOP"
    w=100
    h=100
    item_per_page = 5
    os.chdir(settings.ROOT_DIR+"../Tester/everything")
    log = glob.glob("*.txt")
    f = open(log[0])
    if request.method=='POST': #search particular block, display it, may display patches around it
        form = MatchForm(request.POST)
        if form.is_valid(): 
            size = form.cleaned_data['size']
            x = str((form.cleaned_data['x']/size)*size)
            y = str((form.cleaned_data['y']/size)*size)
            size=str(size) 
            #pattern = "_"+str(form.cleaned_data['size'])+"_("+str(form.cleaned_data['x'])+"_"+str(form.cleaned_data['y'])+")"
            count=0  
            for pos in positions:
                if pos['org'][0] == x and pos['org'][1]==y:
                    #check size
                    orgname = files[count]['org']
                    g = re.findall(r'org_(\d+)',orgname)
                    if g[0] == size:
                         pno = count/item_per_page +1
                         out_files = out_files.paginator.page(pno)
                         out_scores = out_scores.paginator.page(pno)
                         out_plcs = out_plcs.paginator.page(pno)
                         out_pos = out_pos.paginator.page(pno)
                         jumpto = str(count%item_per_page)
                         break
                count=count+1
    else:    
        files=[]
        plcs = []
        scores=[]
        light=[]
        positions=[]
        out_files=[]
        out_scores=[]
        out_plcs=[]
        out_pos=[]
        out_light=[]
        form = MatchForm()
        size=0
        for line in f:
            isOrg = False
            #find target location
            g = re.findall(r'(?<=\()(\d+)(?=,)',line) 
            x=g[0]
            g = re.findall(r'(?<=,)(\d+)(?=\))',line)
            y=g[0]
            #find org and size
            regex = re.compile('org\s+size\s+(\d+)')
            r = regex.search(line)
            if r is not None: # a org line
                size = r.group(1)
                isOrg = True
                org_x = x
                org_y = y
            regex = re.compile('candid\s+(\d+)')
            r = regex.search(line)
            if r is not None: # a cand line
                idx = r.group(1)

            if isOrg:
                if len(score_in_match) !=0: 
                    scores.append(score_in_match)
                if len(pos_in_match) !=0:
                    positions.append(pos_in_match)
                if len(file_in_match) !=0:
                    files.append(file_in_match)
                if len(plc_in_match) !=0:
                    plcs.append(plc_in_match)
                if len(light_in_match)!=0:
                    light.append(light_in_match)
                score_in_match={}
                file_in_match={}
                pos_in_match={}
                plc_in_match={}
                orgpattern = "org_"+str(size)+"_("+str(org_x)+"_"+str(org_y)+").png"
                file_in_match["org"] ="everything/"+orgpattern
                pos_in_match['org'] = (x,y)
                score_in_match['org']=1 #dummy to make the order in the dict same as pos_in_match
                light_in_match['org']=1
                plc_in_match['org']="everything/"+orgpattern #dummy same as above
            else: #store the scores
                regex = re.compile('score:\s+(.*)')
                r = regex.search(line)
                score_in_match[idx]=r.group(1)
                regex = re.compile('light:\s+(.[0-9.]*)')
                r = regex.search(line)
                light_in_match[idx]=r.group(1)
                pos_in_match[idx]=(x,y)
                candpattern = "cand_"+str(size)+"_("+str(org_x)+"_"+str(org_y)+")_("+str(x)+"_"+str(y)+")_("+str(idx)+").png"
                plcpattern="cand_plc_"+str(size)+"_("+str(org_x)+"_"+str(org_y)+")_("+str(x)+"_"+str(y)+")_("+str(idx)+").png"
                file_in_match[idx] = "everything/"+candpattern
                plc_in_match[idx] = "everything/"+plcpattern
        if len(score_in_match) !=0: #the last one need to do manually
            scores.append(score_in_match)
            if len(pos_in_match) !=0:
                positions.append(pos_in_match)
            if len(file_in_match) !=0:
                files.append(file_in_match)
            if len(plc_in_match) !=0:
                plcs.append(plc_in_match)
        #get file page
        out_files = getpage(request,files,item_per_page)
        #get score page
        out_scores = getpage(request,scores,item_per_page)
        #get pos page
        out_pos = getpage(request,positions,item_per_page)
        #get plc page
        out_plcs = getpage(request,plcs,item_per_page)
        out_light = getpage(request,light,item_per_page)
        #paginator = Paginator(files,5)
        #page = request.GET.get('page')
        #try: 
        #    out_files =  paginator.page(page)
        #except PageNotAnInteger:
        #    out_files = paginator.page(1)
        #except EmptyPage:
        #    out_files = paginator.page(paginator.num_pages)
    return render(request,'matching.html',{'form':form,'files':out_files, 'plcs':out_plcs, 'scores':out_scores,'light':out_light,'positions':out_pos,'jump':jumpto})

def getpage(request, list, perpage):
    paginator = Paginator(list,perpage)
    page = request.GET.get('page')
    try: 
        out_files =  paginator.page(page)
    except PageNotAnInteger:
        out_files = paginator.page(1)
    except EmptyPage:
        out_files = paginator.page(paginator.num_pages)
    return out_files

def gotopage(request,list,pno):
    paginator=Paginator(list,perpage)
    page = request.GET.get('page')
    try: 
        out_files =  paginator.page(pno)
    except PageNotAnInteger:
        print "not a valid page number",pno 
        out_files = paginator.page(1)
    except EmptyPage:
        out_files = paginator.page(paginator.num_pages)
        print "empty"
    return out_files

def show_result(request):
    item_per_page =1 
    os.chdir(settings.ROOT_DIR+"../Tester/everything")
    g= glob.glob("*rst*.bmp")
    results = [None] * (len(g))
    print len(results)
    for r in g:
       print re.findall(r'^rst_(\d+)',r)[0]
       results[ int(re.findall(r'^rst_(\d+)',r)[0])] = "everything/"+ r
    file = getpage(request,results,item_per_page)
    return render(request,'result.html',{'file':file})