#include "Granulate.h"
#include <string>
#include <dirent.h>
#include <sys/stat.h>
namespace metric {

  const vector<Size3> Granulate::blkSizes({Size3(4,4,1),Size3(8,8,1),Size3(16,16,1),Size3(32,32,1),Size3(64,64,1),Size3(128,128,1)});
  const vector<int> Granulate::sz_idx({4,8,16,32,64,128});

  bool Granulate::readFiles(string path, string searchExt)
  {


    //check if path end with /
    char endwith = *path.rbegin();
    if (endwith!='/')
      path=path+"/";


#ifdef WIN32

    string fullSearchPath = path + "*."+ext;
    WIN32_FIND_DATA FindData;
    HANDLE hFind;
    hFind = FindFirstFile( fullSearchPath.c_str(), &FindData );

    if( hFind == INVALID_HANDLE_VALUE )
    {
      cout << "Error searching directory\n";
      return false;
    }
#else

    DIR *dir;
    class dirent *ent;
    dir = opendir(path.c_str());
    cout<<"path is: "<<path<<endl;
    cout<<"search ext is"<<searchExt<<endl;
#endif
    do
    {
#     ifdef WIN32
      string filename = FindData.cFileName;
      filenames.push_back(FindData.cFileName);
#     else
      if (!strcmp(ent->d_name,"." )) continue;
      if (!strcmp(ent->d_name,"..")) continue;
      string file_name = ent->d_name;
      unsigned long found = file_name.rfind("."+searchExt);
    //  cout<<found<<endl;
     // cout<<std::string::npos;
      if (found!=std::string::npos)
      {
          string filename = path + "/" + file_name;
          filenames.push_back(filename);
          cout<<"Got: "<<file_name<<endl;
      }
#     endif

    }
#ifdef WIN32
    while( FindNextFile(hFind, &FindData) > 0 );
    FindClose(path);
#else
    while ((ent = readdir(dir)) != NULL);
    closedir(dir);
#endif

    return true;
  }

  bool Granulate::generateGrid()
  {
    string outdir = "/home/guoxin/Projects/MTC/data/totest/gran";
    mkdir(outdir.c_str(),0755);
    for (string& f : this->filenames)
    {
        cout<<f<<endl;
        Tensor<uchar,1> img(f);
        //Tensor<uchar,1> imgc=img.Clone();
        //img.Display();
        cout<<f<<endl;
        int pos1 = f.find_last_of('.');
        int pos2 = f.find_last_of("/");
        //cout<<pos1<<","<<pos2<<endl;
        //cout<<f.substr(pos2,pos1)<<endl;
        // cout<<outdir+"/"+f.substr(pos2+1,pos1-pos2-1)+outstring<<endl;
        //write original
        img.SaveBlock(outdir+"/"+f.substr(pos2+1,pos1-pos2-1)+"_org.png");
        for (int i=(int)sz_idx.size()-1; i>=0;i--)
        {

            int stepx = blkSizes[i].height;
            int stepy = blkSizes[i].width;
            cout<<stepx<<","<<stepy<<endl;
            for (int x=0; x< img.size().height; x++)
              for (int y=0; y< img.size().width; y++)
              {
                //draw vertical lines
                  if ((y%stepy)==0)
                    img(x,y,0)=255;
                //draw horizontal lines
                  if ((x%stepx)==0)
                    img(x,y,0)=255;
              }
          //   img.Display();
            string outstring="_g"+std::to_string(sz_idx[i])+".png";
            img.SaveBlock(outdir+"/"+f.substr(pos2+1,pos1-pos2-1)+outstring);
            //img=imgc.Clone();
        }

    }
    return true;
  }

}
