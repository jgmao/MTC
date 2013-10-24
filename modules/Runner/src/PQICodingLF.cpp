using namespace std;
#include "MTC.h"
int main(int argc, char* argv[])
{
    int bsize = atoi(argv[1]);
    cout<<bsize<<endl;
    string filename = argv[2];
    string outname = argv[3];
    mtc::MTC coder(filename,1,Size3(2*bsize,2*bsize,1),Size3(1,1,1),4,8,Size3_<double>(0.25,0.25,0));
    coder.SetPQILFBlockSize(bsize);
    Tensor<double,1> out = coder.PQICodingLFComponent();
    out.SaveBlock(outname);
    return 0;

}

