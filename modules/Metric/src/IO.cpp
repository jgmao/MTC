#include "IO.h"
#include <fstream>
namespace metric {
  bool readMahalanobis(vector<Mat> &covars, string filename)
  {
        covars.clear();
        std::ifstream file(filename);
        vector<double> temp;
        int i=-1;
        for (string line; getline(file,line);)
        {
          if (string::npos==line.find('<'))
          {
                temp.push_back(std::stod(line));
          }
          else
          {
              if (i<0)
              {
                i++;
                continue;
              }
              int sz = temp.size();
              covars.push_back(Mat::zeros(sz,sz,CV_64F));
              for (int j=0; j<sz;j++)
              {
                  covars[i].at<double>(j,j) = temp[j];
              }
              temp.clear();
              i++;
          }
        }
        int sz = temp.size();
        covars.push_back(Mat::zeros(sz,sz,CV_64F));
        for (int j=0; j<sz;j++)
        {
            covars[i].at<double>(j,j) = temp[j];
        }
        temp.clear();
        return true;
  }

  bool readMahalanobis2(vector<Mat> &covars, string filename)
  {
        covars.clear();
        std::ifstream file(filename);
        vector<vector<double> > temp;
        vector<double> temp2;
        int i=-1;
        int m=0;
        for (string line; getline(file,line);)
        {
          if (string::npos==line.find('<'))
          {
              temp.push_back(vector<double>());
              std::istringstream iss(line);
              double value;
              while (iss>>value)
              {
                temp[m].push_back(value);
              }
              m++;
          }
          else
          {

              if (i<0)
              {
                i++;
                m=0;
                continue;
              }
              int sz = temp.size();
              covars.push_back(Mat::zeros(sz,sz,CV_64F));

              for (int j=0; j<sz;j++)
                for (int k=0; k<sz; k++)
                {
                  covars[i].at<double>(j,k) = temp[j][k];
                }
              m=0;
              temp.clear();
              i++;
          }

        }
        int sz = temp.size();
        covars.push_back(Mat::zeros(sz,sz,CV_64F));
        for (int j=0; j<sz;j++)
          for (int k=0; k<sz; k++)
          {
            covars[i].at<double>(j,k) = temp[j][k];
          }
        temp.clear();
        for (unsigned int i=0; i<covars.size(); i++)
          cout<<covars[i]<<endl;
        return true;
  }

}
