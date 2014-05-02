#include <iostream>
using namespace std;
int main(int argc, char *argv[])
{

  string temp;
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
