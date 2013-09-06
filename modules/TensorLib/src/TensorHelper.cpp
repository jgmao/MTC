#include <TensorHelper.h>
namespace tensor{
void CombineImage(const vector<string>& infilenames, string& outfilename, int step)
{
	int count = infilenames.size();
	vector<Tensor<uchar,3>> t(count);
	int maxheight = 0;
	int totalwidth = 0;
	for (int i=0; i< count; i++)
	{
		t[i].Load(infilenames[i]);
		if (t[i].size().height > maxheight)
			maxheight = t[i].size().height;
		totalwidth+= t[i].size().width;
	}
	Tensor<uchar,3> r(Size3(maxheight,totalwidth + (count-1)*step,1));
	for (int i=0; i< count; i++)
		r.SetBlock(Point3i(0,i*step+i*t[i].size().width,0),t[i]);
	r.SaveBlock(outfilename);
	r.Display(3000,1);

}
}
