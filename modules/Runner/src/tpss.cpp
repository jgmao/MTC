
using namespace std;
#include <TensorLite.h>
#include <Lighting.h>
int main(int argc, char* argv[])
{
    string filename = argv[1];
    string outname = argv[2];
    Tensor<double,1> im(filename);
    Tensor<double,1> out = lighting::ComputeTPSS(im,0.01);
    out.SaveBlock(outname);
    return 0;
}
