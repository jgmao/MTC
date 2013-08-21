#include "LRI.h"
namespace metric{
LRI::LRI()
{
  a = new int[8]();
  gr = new double[16];
}
LRI::~LRI()
{
  delete [] a;
  delete [] gr;
}
double LRI::var(const Mat& im)
{
    int n = im.size().height;
	int m = im.size().width;
    double mean = 0;
    double variance = 0;
    for(int i = 0; i < n; i++)
    {
        for (int j = 0; j < m; j++)
        {
            mean = mean + im.at<double>(i,j) / (double)n / (double)m;
            variance = variance + im.at<double>(i,j) * im.at<double>(i,j) / (double)n / (double)m;
        }
    }
    variance = variance - mean * mean;
    return variance;
}


int* LRI::lri8k(double image[][9], double th)
    {
	int k = 4;
        int flag[8] = {1,1,1,1,1,1,1,1};
	//int *a = new int [8]();// = {0,0,0,0,0,0,0,0};
        for (int ai=0; ai<8;ai++)
          a[ai]=0;

        double c = image[k][k];

        if(image[k][k+1] > (c+th))
            a[0] = 1;
        else if(image[k][k+1] < (c-th))
            a[0] = -1;
        else
            flag[0] = 0;
        if(image[k-1][k+1] > (c+th))
            a[1] = 1;
        else if(image[k-1][k+1] < (c-th))
            a[1] = -1;
        else
            flag[1] = 0;
        if(image[k-1][k] > (c+th))
            a[2] = 1;
        else if(image[k-1][k] < (c-th))
            a[2] = -1;
        else
            flag[2] = 0;
        if(image[k-1][k-1] > (c+th))
            a[3] = 1;
        else if(image[k-1][k-1] < (c-th))
            a[3] = -1;
        else
            flag[3] = 0;
        if(image[k][k-1] > (c+th))
            a[4] = 1;
        else if(image[k][k-1] < (c-th))
            a[4] = -1;
        else
            flag[4] = 0;
        if(image[k+1][k-1] > (c+th))
            a[5] = 1;
        else if(image[k+1][k-1] < (c-th))
            a[5] = -1;
        else
            flag[5] = 0;
        if(image[k+1][k] > (c+th))
            a[6] = 1;
        else if(image[k+1][k] < (c-th))
            a[6] = -1;
        else
            flag[6] = 0;
        if(image[k+1][k+1] > (c+th))
            a[7] = 1;
        else if(image[k+1][k+1] < (c-th))
            a[7] = -1;
        else
            flag[7] = 0;

        for(int i = 2;i <= k; i++)
        {
            if((flag[0] == 1) && (image[k][k+i] > (c+th)) && (a[0] > 0))
                a[0] = a[0] + 1;
            else if((flag[0] == 1) && (image[k][k+i] < (c-th)) && (a[0] < 0))
                a[0] = a[0] - 1;
            else
                flag[0] = 0;
            if((flag[1] == 1) && (image[k-i][k+i] > (c+th)) && (a[1] > 0))
                a[1] = a[1] + 1;
            else if((flag[1] == 1) && (image[k-i][k+i] < (c-th)) && (a[1] < 0))
                a[1] = a[1] - 1;
            else
                flag[1] = 0;  
            if((flag[2] == 1) && (image[k-i][k] > (c+th)) && (a[2] > 0))
                a[2] = a[2] + 1;
            else if((flag[2] == 1) && (image[k-i][k] < (c-th)) && (a[2] < 0))
                a[2] = a[2] - 1;
            else
                flag[2] = 0;
            if((flag[3] == 1) && (image[k-i][k-i] > (c+th)) && (a[3] > 0))
                a[3] = a[3] + 1;
            else if((flag[3] == 1) && (image[k-i][k-i] < (c-th)) && (a[3] < 0))
                a[3] = a[3] - 1;
            else
                flag[3] = 0;
            if((flag[4] == 1) && (image[k][k-i] > (c+th)) && (a[4] > 0))
                a[4] = a[4] + 1;
            else if((flag[4] == 1) && (image[k][k-i] < (c-th)) && (a[4] < 0))
                a[4] = a[4] - 1;
            else
                flag[4] = 0;
            if((flag[5] == 1) && (image[k+i][k-i] > (c+th)) && (a[5] > 0))
                a[5] = a[5] + 1;
            else if((flag[5] == 1) && (image[k+i][k-i] < (c-th)) && (a[5] < 0))
                a[5] = a[5] - 1;
            else
                flag[5] = 0; 
            if((flag[6] == 1) && (image[k+i][k] > (c+th)) && (a[6] > 0))
                a[6] = a[6] + 1;
            else if((flag[6] == 1) && (image[k+i][k] < (c-th)) && (a[6] < 0))
                a[6] = a[6] - 1;
            else
                flag[6] = 0;
            if((flag[7] == 1) && (image[k+i][k+i] > (c+th)) && (a[7] > 0))
                a[7] = a[7] + 1;
            else if((flag[7] == 1) && (image[k+i][k+i] < (c-th)) && (a[7] < 0))
                a[7] = a[7] - 1;
            else
                flag[7] = 0;
        }
        return a;
}

double LRI::computeLRI(const Mat& im1, const Mat& im2)
{
	int K = 4;
	double T = 2;
	double tempim[9][9];
    int* temp1;
    int* temp2;

	for(int i = 0; i < 9; i++)
		for(int j = 0; j < 9; j++)
			tempim[i][j] = 0;
	//for(int i = 0; i < 8; i++)
	//{
	//	temp1[i] = 0;
	//	temp2[i] = 0;
	//}

    double sum = 0;
    double th1;
    double th2;
    th1 = sqrt(var(im1)) / T;
	th2 = sqrt(var(im2)) / T;
	double h1[8][9];
	double h2[8][9];

	for(int i = 0; i < 8; i++)
	{
		for(int j = 0; j < 9; j++)
		{
			h1[i][j] = 0;
			h2[i][j] = 0;
		}
	}

	double score = 0;

	for(int i = K; i < im1.size().height-K; i++)
	{
		for(int j = K; j < im1.size().width-K; j++)
		{
            for (int u = -K; u <= K; u++)
                for (int v = -K; v <= K; v++)
                    tempim[u + K][v + K] = im1.at<double>(i+u,j+v);
            temp1 = lri8k(tempim, th1);
            for (int k = 0; k < 8; k++)
                h1[k][temp1[k] + K] = h1[k][temp1[k] + K] + 1;
        }
    }

            for (int i = 0; i < 8; i++)
            {
                for (int j = 0; j < 9; j++)
                {
                    if (h1[i][j] == 0)
                        h1[i][j] = 1;
                }
            }

            for (int i = 0; i < 8; i++)
            {
                sum = 0;
                for (int j = 0; j < 9; j++)
                    sum = sum + h1[i][j];
                for (int j = 0; j < 9; j++)
                    h1[i][j] = h1[i][j] / sum;
            }


	for(int i = K; i < im2.size().height-K; i++)
	{
		for(int j = K; j < im2.size().width-K; j++)
		{
            for (int u = -K; u <= K; u++)
                for (int v = -K; v <= K; v++)
                    tempim[u + K][v + K] = im2.at<double>(i+u,j+v);
            temp2 = lri8k(tempim, th2);
            for (int k = 0; k < 8; k++)
                h2[k][temp2[k] + K] = h2[k][temp2[k] + K] + 1;
        }
    }

            for (int i = 0; i < 8; i++)
            {
                for (int j = 0; j < 9; j++)
                {
                    if (h2[i][j] == 0)
                        h2[i][j] = 1;
                }
            }

            for (int i = 0; i < 8; i++)
            {
                sum = 0;
                for (int j = 0; j < 9; j++)
                    sum = sum + h2[i][j];
                for (int j = 0; j < 9; j++)
                    h2[i][j] = h2[i][j] / sum;
            }
            

    double value[8];
	for(int i = 0; i < 8; i++)
		value[i] = 0;
    for(int i = 0; i < 8; i++)
    {
		for (int j = 0; j < 9; j++)
		{
            value[i] = value[i] + h1[i][j] * log(h1[i][j] / h2[i][j]);
		}
        score = score + value[i] / 8.0;
    }
    //delete [] temp1;
    //delete [] temp2;
	return score;
}


int LRI::lbp81(double image[][3])
{
	int index = 0;
	double a[8] = {0 ,0, 0, 0, 0, 0, 0, 0};
	a[0] = image[1][2];
	a[2] = image[0][1];
	a[4] = image[1][0];
	a[6] = image[2][1];
	a[1] = image[0][2];
	a[3] = image[0][0];
	a[5] = image[2][0];
	a[7] = image[2][2];

	for(int i = 0; i < 8; i++)
	{
		if(a[i] >= image[1][1])
			a[i] = 1;
		else
			a[i] = 0;
	}
	int t = 0;
	for(int i = 0; i < 8; i++)
	{
		if(i < 7)
			t = t + abs(a[i] - a[i+1]);
		else
			t = t + abs(a[7] - a[0]);
	}

	if(t > 2)
		index = 9;
	else
		index = a[0] + a[1] + a[2] + a[3] + a[4] + a[5] + a[6] + a[7];
	
	return index;
}


double LRI::computeLBP(const Mat& im1, const Mat& im2)
    {
        double tempim[3][3];
	for(int i = 0; i < 3; i++)
		for(int j = 0; j < 3; j++)
			tempim[i][j] = 0;

        int temp1 = 0;
        int temp2 = 0;
        double sum = 0;
        double value = 0;

        double h1[10];
	for(int i = 0; i < 10; i++)
		h1[i] = 0;

        for (int i = 1; i < im1.size().height-1; i++)
        {
            for (int j = 1; j < im1.size().width-1; j++)
            {
                for (int u = -1; u <= 1; u++)
                    for (int v = -1; v <= 1; v++)
                        tempim[u + 1][v + 1] = im1.at<double>(i + u, j + v);
                temp1 = lbp81(tempim);
                h1[temp1] = h1[temp1] + 1;
            }
        }

        for (int i = 0; i < 10; i++)
        {
            if (h1[i] == 0)
                h1[i] = 1;
        }

        sum = 0;
        for (int j = 0; j < 10; j++)
            sum = sum + h1[j];
        for (int j = 0; j < 10; j++)
            h1[j] = h1[j] / sum;

	double h2[10];
	for(int i = 0; i < 10; i++)
		h2[i] = 0;

        for (int i = 1; i < im2.size().height-1; i++)
        {
            for (int j = 1; j < im2.size().width-1; j++)
            {
                for (int u = -1; u <= 1; u++)
                    for (int v = -1; v <= 1; v++)
                        tempim[u + 1][v + 1] = im2.at<double>(i + u, j + v);
                temp2 = lbp81(tempim);
                h2[temp2] = h2[temp2] + 1;
            }
        }

        for (int i = 0; i < 10; i++)
        {
            if (h2[i] == 0)
                h2[i] = 1;
        }

        sum = 0;
        for (int j = 0; j < 10; j++)
            sum = sum + h2[j];
        for (int j = 0; j < 10; j++)
            h2[j] = h2[j] / sum;
               
        value = 0;
        sum = 0;
        for(int j = 0; j < 10; j++)
            value = value + h1[j] * log(h1[j] / h2[j]);
        return value;
    }


double* LRI::grad(const Mat& im)
{
	//double* gr = new double[16];// = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
  for (int gi=0; gi<16; gi++)
    gr[gi]=0;
	Mat g1(im.size().height-4,im.size().width,CV_64F);
	Mat g2(im.size().height,im.size().width-4,CV_64F);
	Mat g3(im.size().height-4,im.size().width-4,CV_64F);
	Mat g4(im.size().height-4,im.size().width-4,CV_64F);

	for(int i = 1;i <= 4 ; i++)
	{
		for(int u = 0; u < im.size().height-4; u++)
			for(int v = 0; v < im.size().width; v++)
				g1.at<double>(u,v) = im.at<double>(4+u,v) - im.at<double>(4+u-i,v);
		gr[4*(i-1)] = var(g1);

		for(int u = 0; u < im.size().height; u++)
			for(int v = 0; v < im.size().width-4; v++)
				g2.at<double>(u,v) = im.at<double>(u,4+v) - im.at<double>(u,4+v-i);
		gr[4*(i-1)+1] = var(g2);

		for(int u = 0; u < im.size().height-4; u++)
			for(int v = 0; v < im.size().width-4; v++)
				g3.at<double>(u,v) = im.at<double>(4+u,4+v) - im.at<double>(4+u-i,4+v-i);
		gr[4*(i-1)+2] = var(g3);

		for(int u = 0; u < im.size().height-4; u++)
			for(int v = 0; v < im.size().width-4; v++)
				g4.at<double>(u,v) = im.at<double>(4+u,v) - im.at<double>(4+u-i,v+i);
		gr[4*(i-1)+3] = var(g4);
	}
	return gr;
}


double LRI::sedest(const Mat& im1, const Mat& im2)
{
	double* gr;
	double grad1[16];
	double grad2[16];
	for(int i = 0; i < 16; i++)
	{
		grad1[i] = 0;
		grad2[i] = 0;
	}

	gr = grad(im1);
	for(int i = 0; i < 16; i++)
		grad1[i] = gr[i];

	gr = grad(im2);
	for(int i = 0; i < 16; i++)
		grad2[i] = gr[i];

	double C = 10;
	double score = 1;
	for(int i = 0; i < 16; i++)
		score = score * (2*sqrt((grad1[i]+C)*(grad2[i]+C))) / (grad1[i]+grad2[i]+2*C);
 // delete [] gr;
	return score;
}


double LRI::lp(const Mat& im1, const Mat& im2)
{
	double mvalue = 0;;
	double m1=0,m2=0;
	for(int i=0;i < im1.size().height;i++)
	{
		for(int j=0;j < im1.size().width;j++)
		{
			m1 = m1 + im1.at<double>(i,j);
			m2 = m2 + im2.at<double>(i,j);
		}
	}
	m1 = m1 / im1.size().height / im1.size().width;
	m2 = m2 / im2.size().height / im2.size().width;

	mvalue = pow(max(11.0,abs(m1 - m2))/256 , 1.3);

	return mvalue;
}


double LRI::computeNewMetric(const Mat& im1, const Mat& im2)
{
	double lrivalue = 0;
	lrivalue = computeLRI(im1, im2);

	double lbpvalue = 0;
	lbpvalue = computeLBP(im1, im2);

	double sedvalue = 0;
	sedvalue = sedest(im1, im2);

	double lpvalue = 0;
	lpvalue = lp(im1, im2);

	double metricvalue = 0;
	metricvalue = lbpvalue * lrivalue * lpvalue * pow((1 - sedvalue),0.5);
	metricvalue = -log10(metricvalue);//output non negative 

	return metricvalue;
}
}
