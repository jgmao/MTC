#include <iostream>
#include <TensorLite.h>
using namespace tensor;
using namespace std;
int main(int argc, char *argv[])
{

  string temp;
  //read file
  temp = string(argv[1]);
  cout<<"reading: "<<temp<<endl;
  Tensor<double,1> test(temp);
  //test.Print();
  for (int i=0; i<10; i++)
    {
  cout<<"waiting input:"<<endl;
  cin >>temp;
  //int key=std::stoi(temp);
  //int key  =1;
  cout<<"the input is: "<<temp<<endl;
    }
  cout<<"done!"<<endl;
  return 0;
}
