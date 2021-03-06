#include <iostream>                   // Console I/O
#include <sstream>                    // String to number conversion
#if (CV_MINOR_VERSION > 5)
#include <opencv2/core.hpp>      // Basic OpenCV structures
#include <opencv2/core/utility.hpp>
#include <opencv2/imgproc.hpp>// Image processing methods for the CPU
#include <opencv2/highgui.hpp>// Read images
// GPU structures and methods
#include <opencv2/gpuarithm.hpp>
#include <opencv2/gpufilters.hpp>

#else
#include <opencv2/opencv.hpp>
#include <opencv2/gpu/gpu.hpp>
#endif

using namespace std;
using namespace cv;

double getPSNR(const Mat& I1, const Mat& I2);      // CPU versions
Scalar getMSSIM( const Mat& I1, const Mat& I2);

double getPSNR_GPU(const Mat& I1, const Mat& I2);  // Basic GPU versions
Scalar getMSSIM_GPU( const Mat& I1, const Mat& I2);

struct BufferPSNR                                     // Optimized GPU versions
{   // Data allocations are very expensive on GPU. Use a buffer to solve: allocate once reuse later.
    gpu::GpuMat gI1, gI2, gs, t1,t2;

    gpu::GpuMat buf;
};
double getPSNR_GPU_optimized(const Mat& I1, const Mat& I2, BufferPSNR& b);

struct BufferMSSIM                                     // Optimized GPU versions
{   // Data allocations are very expensive on GPU. Use a buffer to solve: allocate once reuse later.
    gpu::GpuMat gI1, gI2, gs, t1,t2;

    gpu::GpuMat I1_2, I2_2, I1_I2;
    vector<gpu::GpuMat> vI1, vI2;

    gpu::GpuMat mu1, mu2;
    gpu::GpuMat mu1_2, mu2_2, mu1_mu2;

    gpu::GpuMat sigma1_2, sigma2_2, sigma12;
    gpu::GpuMat t3;

    gpu::GpuMat ssim_map;

    gpu::GpuMat buf;
};
Scalar getMSSIM_GPU_optimized( const Mat& i1, const Mat& i2, BufferMSSIM& b);

static void help()
{
    cout
        << "\n--------------------------------------------------------------------------" << endl
        << "This program shows how to port your CPU code to GPU or write that from scratch." << endl
        << "You can see the performance improvement for the similarity check methods (PSNR and SSIM)."  << endl
        << "Usage:"                                                               << endl
        << "./gpu-basics-similarity referenceImage comparedImage numberOfTimesToRunTest(like 10)." << endl
        << "--------------------------------------------------------------------------"   << endl
        << endl;
}

int main(int, char *argv[])
{
    help();
    Mat I1 = imread(argv[1]);           // Read the two images
    Mat I2 = imread(argv[2]);
    namedWindow("test1",WINDOW_NORMAL );
    imshow("test1",I1);
    waitKey();
    imshow("test1",I2);
    waitKey();
    if (!I1.data || !I2.data)           // Check for success
    {
        cout << "Couldn't read the image";
        return 0;
    }

    BufferPSNR bufferPSNR;
    BufferMSSIM bufferMSSIM;

    int TIMES;
    stringstream sstr(argv[3]);
    sstr >> TIMES;
    double time, result = 0;

    //------------------------------- PSNR CPU ----------------------------------------------------
    time = (double)getTickCount();

    for (int i = 0; i < TIMES; ++i)
        result = getPSNR(I1,I2);

    time = 1000*((double)getTickCount() - time)/getTickFrequency();
    time /= TIMES;

    cout << "Time of PSNR CPU (averaged for " << TIMES << " runs): " << time << " milliseconds."
        << " With result of: " << result << endl;

    //------------------------------- PSNR GPU ----------------------------------------------------
    time = (double)getTickCount();

    for (int i = 0; i < TIMES; ++i)
        result = getPSNR_GPU(I1,I2);

    time = 1000*((double)getTickCount() - time)/getTickFrequency();
    time /= TIMES;

    cout << "Time of PSNR GPU (averaged for " << TIMES << " runs): " << time << " milliseconds."
        << " With result of: " <<  result << endl;

    //------------------------------- PSNR GPU Optimized--------------------------------------------
    time = (double)getTickCount();                                  // Initial call
    result = getPSNR_GPU_optimized(I1, I2, bufferPSNR);
    time = 1000*((double)getTickCount() - time)/getTickFrequency();
    cout << "Initial call GPU optimized:              " << time  <<" milliseconds."
        << " With result of: " << result << endl;

    time = (double)getTickCount();
    for (int i = 0; i < TIMES; ++i)
        result = getPSNR_GPU_optimized(I1, I2, bufferPSNR);

    time = 1000*((double)getTickCount() - time)/getTickFrequency();
    time /= TIMES;

    cout << "Time of PSNR GPU OPTIMIZED ( / " << TIMES << " runs): " << time
        << " milliseconds." << " With result of: " <<  result << endl << endl;


    //------------------------------- SSIM CPU -----------------------------------------------------
    Scalar x;
    time = (double)getTickCount();

    for (int i = 0; i < TIMES; ++i)
        x = getMSSIM(I1,I2);

    time = 1000*((double)getTickCount() - time)/getTickFrequency();
    time /= TIMES;

    cout << "Time of MSSIM CPU (averaged for " << TIMES << " runs): " << time << " milliseconds."
        << " With result of B" << x.val[0] << " G" << x.val[1] << " R" << x.val[2] << endl;

    //------------------------------- SSIM GPU -----------------------------------------------------
    time = (double)getTickCount();

    for (int i = 0; i < TIMES; ++i)
        x = getMSSIM_GPU(I1,I2);

    time = 1000*((double)getTickCount() - time)/getTickFrequency();
    time /= TIMES;

    cout << "Time of MSSIM GPU (averaged for " << TIMES << " runs): " << time << " milliseconds."
        << " With result of B" << x.val[0] << " G" << x.val[1] << " R" << x.val[2] << endl;

    //------------------------------- SSIM GPU Optimized--------------------------------------------
    time = (double)getTickCount();
    x = getMSSIM_GPU_optimized(I1,I2, bufferMSSIM);
    time = 1000*((double)getTickCount() - time)/getTickFrequency();
    cout << "Time of MSSIM GPU Initial Call            " << time << " milliseconds."
        << " With result of B" << x.val[0] << " G" << x.val[1] << " R" << x.val[2] << endl;

    time = (double)getTickCount();

    for (int i = 0; i < TIMES; ++i)
        x = getMSSIM_GPU_optimized(I1,I2, bufferMSSIM);

    time = 1000*((double)getTickCount() - time)/getTickFrequency();
    time /= TIMES;

    cout << "Time of MSSIM GPU OPTIMIZED ( / " << TIMES << " runs): " << time << " milliseconds."
        << " With result of B" << x.val[0] << " G" << x.val[1] << " R" << x.val[2] << endl << endl;
    return 0;
}


double getPSNR(const Mat& I1, const Mat& I2)
{
    Mat s1;
    absdiff(I1, I2, s1);       // |I1 - I2|
    s1.convertTo(s1, CV_32F);  // cannot make a square on 8 bits
    s1 = s1.mul(s1);           // |I1 - I2|^2

    Scalar s = sum(s1);         // sum elements per channel

    double sse = s.val[0] + s.val[1] + s.val[2]; // sum channels

    if( sse <= 1e-10) // for small values return zero
        return 0;
    else
    {
        double  mse =sse /(double)(I1.channels() * I1.total());
        double psnr = 10.0*log10((255*255)/mse);
        return psnr;
    }
}



double getPSNR_GPU_optimized(const Mat& I1, const Mat& I2, BufferPSNR& b)
{
    b.gI1.upload(I1);
    b.gI2.upload(I2);

    b.gI1.convertTo(b.t1, CV_32F);
    b.gI2.convertTo(b.t2, CV_32F);

    gpu::absdiff(b.t1.reshape(1), b.t2.reshape(1), b.gs);
    gpu::multiply(b.gs, b.gs, b.gs);

    double sse = gpu::sum(b.gs, b.buf)[0];

    if( sse <= 1e-10) // for small values return zero
        return 0;
    else
    {
        double mse = sse /(double)(I1.channels() * I1.total());
        double psnr = 10.0*log10((255*255)/mse);
        return psnr;
    }
}

double getPSNR_GPU(const Mat& I1, const Mat& I2)
{
    gpu::GpuMat gI1, gI2, gs, t1,t2;

    gI1.upload(I1);
    gI2.upload(I2);

    gI1.convertTo(t1, CV_32F);
    gI2.convertTo(t2, CV_32F);

    gpu::absdiff(t1.reshape(1), t2.reshape(1), gs);
    gpu::multiply(gs, gs, gs);

    Scalar s = gpu::sum(gs);
    double sse = s.val[0] + s.val[1] + s.val[2];

    if( sse <= 1e-10) // for small values return zero
        return 0;
    else
    {
        double  mse =sse /(double)(gI1.channels() * I1.total());
        double psnr = 10.0*log10((255*255)/mse);
        return psnr;
    }
}

Scalar getMSSIM( const Mat& i1, const Mat& i2)
{
    const double C1 = 6.5025, C2 = 58.5225;
    /***************************** INITS **********************************/
    int d     = CV_32F;

    Mat I1, I2;
    i1.convertTo(I1, d);           // cannot calculate on one byte large values
    i2.convertTo(I2, d);

    Mat I2_2   = I2.mul(I2);        // I2^2
    Mat I1_2   = I1.mul(I1);        // I1^2
    Mat I1_I2  = I1.mul(I2);        // I1 * I2

    /*************************** END INITS **********************************/

    Mat mu1, mu2;   // PRELIMINARY COMPUTING
    GaussianBlur(I1, mu1, Size(11, 11), 1.5);
    GaussianBlur(I2, mu2, Size(11, 11), 1.5);

    Mat mu1_2   =   mu1.mul(mu1);
    Mat mu2_2   =   mu2.mul(mu2);
    Mat mu1_mu2 =   mu1.mul(mu2);

    Mat sigma1_2, sigma2_2, sigma12;

    GaussianBlur(I1_2, sigma1_2, Size(11, 11), 1.5);
    sigma1_2 -= mu1_2;

    GaussianBlur(I2_2, sigma2_2, Size(11, 11), 1.5);
    sigma2_2 -= mu2_2;

    GaussianBlur(I1_I2, sigma12, Size(11, 11), 1.5);
    sigma12 -= mu1_mu2;

    ///////////////////////////////// FORMULA ////////////////////////////////
    Mat t1, t2, t3;

    t1 = 2 * mu1_mu2 + C1;
    t2 = 2 * sigma12 + C2;
    t3 = t1.mul(t2);              // t3 = ((2*mu1_mu2 + C1).*(2*sigma12 + C2))

    t1 = mu1_2 + mu2_2 + C1;
    t2 = sigma1_2 + sigma2_2 + C2;
    t1 = t1.mul(t2);               // t1 =((mu1_2 + mu2_2 + C1).*(sigma1_2 + sigma2_2 + C2))

    Mat ssim_map;
    divide(t3, t1, ssim_map);      // ssim_map =  t3./t1;

    Scalar mssim = mean( ssim_map ); // mssim = average of ssim map
    return mssim;
}

Scalar getMSSIM_GPU( const Mat& i1, const Mat& i2)
{
    const float C1 = 6.5025f, C2 = 58.5225f;
    /***************************** INITS **********************************/
    gpu::GpuMat gI1, gI2, gs1, tmp1,tmp2;

    gI1.upload(i1);
    gI2.upload(i2);

    gI1.convertTo(tmp1, CV_MAKE_TYPE(CV_32F, gI1.channels()));
    gI2.convertTo(tmp2, CV_MAKE_TYPE(CV_32F, gI2.channels()));

    vector<gpu::GpuMat> vI1, vI2;
    gpu::split(tmp1, vI1);
    gpu::split(tmp2, vI2);
    Scalar mssim;
#if (CV_MINOR_VERSION >5)
    Ptr<gpu::Filter> gauss = gpu::createGaussianFilter(vI2[0].type(), -1, Size(11, 11), 1.5);
#else
    cv::Ptr<gpu::FilterEngine_GPU> gauss=gpu::createGaussianFilter_GPU(vI2[0].type(), Size(11, 11), 1.5);
#endif
    for( int i = 0; i < gI1.channels(); ++i )
    {
        gpu::GpuMat I2_2, I1_2, I1_I2;

        gpu::multiply(vI2[i], vI2[i], I2_2);        // I2^2
        gpu::multiply(vI1[i], vI1[i], I1_2);        // I1^2
        gpu::multiply(vI1[i], vI2[i], I1_I2);       // I1 * I2

        /*************************** END INITS **********************************/
        gpu::GpuMat mu1, mu2;   // PRELIMINARY COMPUTING
        gauss->apply(vI1[i], mu1);
        gauss->apply(vI2[i], mu2);

        gpu::GpuMat mu1_2, mu2_2, mu1_mu2;
        gpu::multiply(mu1, mu1, mu1_2);
        gpu::multiply(mu2, mu2, mu2_2);
        gpu::multiply(mu1, mu2, mu1_mu2);

        gpu::GpuMat sigma1_2, sigma2_2, sigma12;

        gauss->apply(I1_2, sigma1_2);
        gpu::subtract(sigma1_2, mu1_2, sigma1_2); // sigma1_2 -= mu1_2;

        gauss->apply(I2_2, sigma2_2);
        gpu::subtract(sigma2_2, mu2_2, sigma2_2); // sigma2_2 -= mu2_2;

        gauss->apply(I1_I2, sigma12);
        gpu::subtract(sigma12, mu1_mu2, sigma12); // sigma12 -= mu1_mu2;

        ///////////////////////////////// FORMULA ////////////////////////////////
        gpu::GpuMat t1, t2, t3;

        mu1_mu2.convertTo(t1, -1, 2, C1); // t1 = 2 * mu1_mu2 + C1;
        sigma12.convertTo(t2, -1, 2, C2); // t2 = 2 * sigma12 + C2;
        gpu::multiply(t1, t2, t3);        // t3 = ((2*mu1_mu2 + C1).*(2*sigma12 + C2))

        gpu::addWeighted(mu1_2, 1.0, mu2_2, 1.0, C1, t1);       // t1 = mu1_2 + mu2_2 + C1;
        gpu::addWeighted(sigma1_2, 1.0, sigma2_2, 1.0, C2, t2); // t2 = sigma1_2 + sigma2_2 + C2;
        gpu::multiply(t1, t2, t1);                              // t1 =((mu1_2 + mu2_2 + C1).*(sigma1_2 + sigma2_2 + C2))

        gpu::GpuMat ssim_map;
        gpu::divide(t3, t1, ssim_map);      // ssim_map =  t3./t1;

        Scalar s = gpu::sum(ssim_map);
        mssim.val[i] = s.val[0] / (ssim_map.rows * ssim_map.cols);

    }
    return mssim;
}

Scalar getMSSIM_GPU_optimized( const Mat& i1, const Mat& i2, BufferMSSIM& b)
{
    const float C1 = 6.5025f, C2 = 58.5225f;
    /***************************** INITS **********************************/

    b.gI1.upload(i1);
    b.gI2.upload(i2);

    gpu::Stream stream;

    b.gI1.convertTo(b.t1, CV_32F, stream);
    b.gI2.convertTo(b.t2, CV_32F, stream);

    gpu::split(b.t1, b.vI1, stream);
    gpu::split(b.t2, b.vI2, stream);
    Scalar mssim;
#if (CV_MINOR_VERSION >5)
    Ptr<gpu::Filter> gauss = gpu::createGaussianFilter(b.vI1[0].type(), -1, Size(11, 11), 1.5);
#else
    cv::Ptr<gpu::FilterEngine_GPU> gauss= gpu::createGaussianFilter_GPU(b.vI1[0].type(), Size(11, 11), 1.5);
#endif
#if (CV_MINOR_VERSION >5)
    for( int i = 0; i < b.gI1.channels(); ++i )
    {
        gpu::multiply(b.vI2[i], b.vI2[i], b.I2_2, 1, -1, stream);        // I2^2
        gpu::multiply(b.vI1[i], b.vI1[i], b.I1_2, 1, -1, stream);        // I1^2
        gpu::multiply(b.vI1[i], b.vI2[i], b.I1_I2, 1, -1, stream);       // I1 * I2

        gauss->apply(b.vI1[i], b.mu1, stream);
        gauss->apply(b.vI2[i], b.mu2, stream);

        gpu::multiply(b.mu1, b.mu1, b.mu1_2, 1, -1, stream);
        gpu::multiply(b.mu2, b.mu2, b.mu2_2, 1, -1, stream);
        gpu::multiply(b.mu1, b.mu2, b.mu1_mu2, 1, -1, stream);

        gauss->apply(b.I1_2, b.sigma1_2, stream);
        gpu::subtract(b.sigma1_2, b.mu1_2, b.sigma1_2, gpu::GpuMat(), -1, stream);
        //b.sigma1_2 -= b.mu1_2;  - This would result in an extra data transfer operation

        gauss->apply(b.I2_2, b.sigma2_2, stream);
        gpu::subtract(b.sigma2_2, b.mu2_2, b.sigma2_2, gpu::GpuMat(), -1, stream);
        //b.sigma2_2 -= b.mu2_2;
s
        gauss->apply(b.I1_I2, b.sigma12, stream);
        gpu::subtract(b.sigma12, b.mu1_mu2, b.sigma12, gpu::GpuMat(), -1, stream);
        //b.sigma12 -= b.mu1_mu2;

        //here too it would be an extra data transfer due to call of operator*(Scalar, Mat)
        gpu::multiply(b.mu1_mu2, 2, b.t1, 1, -1, stream); //b.t1 = 2 * b.mu1_mu2 + C1;
        gpu::add(b.t1, C1, b.t1, gpu::GpuMat(), -1, stream);
        gpu::multiply(b.sigma12, 2, b.t2, 1, -1, stream); //b.t2 = 2 * b.sigma12 + C2;
        gpu::add(b.t2, C2, b.t2, gpu::GpuMat(), -12, stream);

        gpu::multiply(b.t1, b.t2, b.t3, 1, -1, stream);     // t3 = ((2*mu1_mu2 + C1).*(2*sigma12 + C2))

        gpu::add(b.mu1_2, b.mu2_2, b.t1, gpu::GpuMat(), -1, stream);
        gpu::add(b.t1, C1, b.t1, gpu::GpuMat(), -1, stream);

        gpu::add(b.sigma1_2, b.sigma2_2, b.t2, gpu::GpuMat(), -1, stream);
        gpu::add(b.t2, C2, b.t2, gpu::GpuMat(), -1, stream);


        gpu::multiply(b.t1, b.t2, b.t1, 1, -1, stream);     // t1 =((mu1_2 + mu2_2 + C1).*(sigma1_2 + sigma2_2 + C2))
        gpu::divide(b.t3, b.t1, b.ssim_map, 1, -1, stream);      // ssim_map =  t3./t1;

        stream.waitForCompletion();

        Scalar s = gpu::sum(b.ssim_map, b.buf);
        mssim.val[i] = s.val[0] / (b.ssim_map.rows * b.ssim_map.cols);

    }
#endif
    return mssim;
}
