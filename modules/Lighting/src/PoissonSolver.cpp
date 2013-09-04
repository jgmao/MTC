#include "PoissonSolver.h"
namespace lighting{
PoissonSolver::PoissonSolver()
{
}
 
void PoissonSolver::copyTo(PoissonSolver &b) const
{
  b.setImages(_cand, _tar, _mask);
}
 
inline PoissonSolver PoissonSolver::clone() const
{
  PoissonSolver b;
  copyTo(b);
  return b; 
}
 
// constructor
PoissonSolver::PoissonSolver(const cv::Mat &cand, const cv::Mat &target, const cv::Mat &mask,const cv::Mat& ref,  const double &foot,const double& upfoot, const double& leftfoot)
  : _cand(cand), _tar(target), _mask(mask),foot(foot),upfoot(upfoot),leftfoot(leftfoot)
{
  _ref = ref;
  CV_Assert(_mask.channels()==1);
  CV_Assert(_cand.cols==_mask.cols && _cand.rows==_mask.rows);
  
}
PoissonSolver::PoissonSolver(const cv::Mat &cand, const cv::Mat &target, const cv::Mat &mask,const cv::Mat& ref,  const vector<FootType>& feet)
: _cand(cand), _tar(target), _mask(mask),feet(feet)
{
  _ref = ref;
  if(feet.size()>0)
  {
  upfoot = feet[0].value;
  leftfoot= feet[(feet.size()-1)/2].value;
  foot = (feet.end()-1)->value;

  }
  CV_Assert(_mask.channels()==1);
  CV_Assert(_cand.cols==_mask.cols && _cand.rows==_mask.rows);
  
}
PoissonSolver::PoissonSolver(const cv::Mat &cand, const cv::Mat &target, const cv::Mat &mask)
  : _cand(cand), _tar(target), _mask(mask)
{
  CV_Assert(_mask.channels()==1);
  CV_Assert(_cand.cols==_mask.cols && _cand.rows==_mask.rows);
}
 
// set source, tareget and destination images
bool PoissonSolver::setImages(const cv::Mat &cand, const cv::Mat &target, const cv::Mat &mask=cv::Mat())
{
  _cand = cand;
  _tar = target;
  _mask = mask;
   
  CV_Assert(_mask.channels()==1);
  CV_Assert(_cand.cols==_mask.cols && _cand.rows==_mask.rows);
   
  return true;
}
 
// solver sparse linear system
template <typename T> 
bool PoissonSolver::solve(const Eigen::SparseMatrix<T> &A, const Eigen::Matrix<T, Eigen::Dynamic, 1> &b, 
                           Eigen::Matrix<T, Eigen::Dynamic, 1> &u)
{
 // Eigen::SparseLU<Eigen::SparseMatrix<T>> lu_of_A(A);
  Eigen::SimplicialLDLT<Eigen::SparseMatrix<T>> lu_of_A(A);
  //if(!lu_of_A.succeeded()) {
  //  std::cerr<< "decomposition failed" << std::endl;
  //  return false;
  //}
  /*if(!lu_of_A.solve(b,&u)) {
    std::cerr<< "solving failed" << std::endl;
    return false;
  }*/
  u=lu_of_A.solve(b); 
  return true;
}
 
// build matrix as linear system
template <typename T>
bool PoissonSolver::buildLCMatrix(Eigen::SparseMatrix<T> &A, Eigen::Matrix<T, Eigen::Dynamic,1> &b, 
                                 Eigen::Matrix<T, Eigen::Dynamic, 1> &u)
{
  int w = mask_roi1.width;
  int h = mask_roi1.height;
  int nz=0;
  for(int y=0; y<h; ++y) {
    double *p = mask1.ptr<double>(y);
    for(int x=0; x<w; ++x, ++p) {
      if(*p==0) continue;
 
      int id = y*(w*ch)+(x*ch);
      mp[id] = nz++;   // r
     // mp[++id] = nz++; // g
     // mp[++id] = nz++; // b
    }
  }
  //mylib::DisplayMat(ref1,"ref");
 /* if (feet.size()==0)
  {
  feet.push_back(FootType(foot,h-1,w-1));
  feet.push_back(FootType(leftfoot,h-1, 0));
  feet.push_back(FootType(upfoot,0,w-1));
  }*/
  A = Eigen::SparseMatrix<double>(nz, nz);
  b = Eigen::VectorXd(nz);
  u = Eigen::VectorXd(nz);
  int rowA = 0;
  //int bh = h-2;
  //int bw = w-2;
  A.reserve(5*nz);
  //mylib::DisplayMat(mask1,"mask1");
  //mylib::DisplayMat(laplacian_cand,"lap");
  //mylib::DisplayMat(tar1,"tar1");
  // Tensor<double,1>(tar1).Print("tar1",true);
   vector<double> interpRight,interpButtom;
  if(feet.size()>0)
  {
   interpRight=interpFeet(feet, Point(0,w-1),Point(h-1,w-1));
   interpButtom=interpFeet(feet, Point(h-1,0),Point(h-1,w-1));
  }
  for(int y=1; y<=h-1; ++y) //01062013 change y<h-1 to y<=h-1
  {
    double *p = mask1.ptr<double>(y)+1;
    //cv::Vec3d *drv = laplacian.ptr<cv::Vec3d>(y)+1;
    double* drv = laplacian_cand.ptr<double>(y)+1;

  //mylib::DisplayMat(mask1);
  //mylib::DisplayMat(_tar);
  //feet.push_back(FootType(50,32,49));
  //feet.push_back(FootType(30,49,24));

	for(int x=1; x<=w-1; ++x, ++p, ++drv)//01062013 change x<w-1 to x<=w-1
  {
    if(*p==0) continue;
    int id = y*(w*ch)+(x*ch);
    int tidx=id-ch*w, lidx=id-ch, ridx=id+ch, bidx=id+ch*w; //index of up, left, right, down in sparse matrix A
    //double temp=0;
    uchar tlrb = 15; // 0b1111
    //disabled//20130528, solve lighting itself, and add it to candid
    //disabled//*drv = 0; 
    if(mask1.at<double>(y-1,x)==0)//up = 8
    {
     /* int tempX,tempY;
      int tempS = _tar.size().height/8;
      tempX = x -  tempS/2+tempS-1;
      tempY = y -1;
      int tempW;*/
	  	//if (tar1.at<double>(y-1,x)==0)
		  //if (x > w-2-tempS*2)
      //if (
      //{

			  //if (tar1.at<double>(y,x)>0)
        //{
          //_tar(Rect(
          //use ave of that part
				  //*drv -= tar1.at<double>(y,x);
        //}
			  //else
        //{
          //double d ;
          //d = ref1.at<double>(y,x);
				  //*drv-=upfoot;//gj 112612 cheating use org
       // cout<<y<<", "<<x<<":"<<d<<endl;
        //double ttt;
      //  mylib::DisplayMat(tar1,"tar1");
        //cin>>ttt;
        //}
		  //}
		  //else
      //{
      //  if (x+tempS/2 > w-2-tempS*2)
      //    tempW = w-2-tempS*2-x;
      //  else
      //    tempW = tempS/2;
      //  tempW += tempS/2+1;
      //Mat tempT = _tar(Rect(tempX,tempY,tempW,tempS+1));
        //mylib::DisplayMat(_tar);
       // mylib::DisplayMat(tempT);
        //*drv -= cv::mean(tempT)[0];
			*drv -= (tar1.at<double>(y-1,x)-cand1.at<double>(y-1,x));//20130528 try Sadeghi's method
      //*drv -= tar1.at<double>(y-1,x); //from rst
      //cout<<"up:"<<tar1.at<double>(y-1,x)<<endl;
      //}
      tlrb &= 7; //0b0111
    }

    if(mask1.at<double>(y,x-1)==0) //left = 4
    {
      //int tempX,tempY;
      //int tempS = _tar.size().height/8;
      //tempY = y -  tempS/2+tempS-1;
      //tempX = x -1;
      //int tempH;
		  //if (tar1.at<double>(y,x-1)==0)
		  //if (y > h-2-tempS*2)
      //{
			  //if (tar1.at<double>(y,x)>0)
				//  *drv -= tar1.at<double>(y,x);
			  //else
        //{
          //double d;
          //d = ref1.at<double>(y,x); 
  			//	*drv -=  leftfoot; //gj 112612 cheating use org
         //cout<<y<<", "<<x<<":"<<tar1.at<double>(y,x)<<endl;
        //double ttt;
       // mylib::DisplayMat(tar1,"tar1");

        //cin>>ttt;
        //}
		  //}
	  	//else
      //{
      //  if (y+tempS/2 > h-2-tempS*2)
      //    tempH = h-2-tempS*2-y;
      //  else
      //    tempH = tempS/2;
      //  tempH += tempS/2+1;
      //  Mat tempT = _tar(Rect(tempX,tempY,tempS+1,tempH));
        //mylib::DisplayMat(_tar);
        //mylib::DisplayMat(tempT);
        //*drv -= cv::mean(tempT)[0];
			  *drv -= (tar1.at<double>(y,x-1)-cand1.at<double>(y,x-1));//20130528 try Sadeghi's method
       // *drv -= tar1.at<double>(y,x-1); //from the result edge
        // cout<<"left:"<<tar1.at<double>(y,x-1)<<endl;
      //}
      tlrb &= 11; //0b1011
    }
	  //this part is tricky, let it be default for test only
	  //if (mask1.at<uchar>(y,x+1)==0 && mask1.at<uchar>(y+1,x)==0)
	  
    if(mask1.at<double>(y,x+1)==0) //right = 2
    {
		  //double alpha = double(h-2-y)/double(h-3);
		  //double beta = double(h-2*y)/double(h-2);
    //    //*drv -= tar1.at<UINT8>(y,x+1);
		  //double temp;
		  //double delta= abs(upfoot-foot);
		//if (upfoot<0)
		//	if (y<h/2)
		//		temp = foot;
		//	else
		//		temp = tar1.at<double>(0,w-2) *alpha + (1-alpha)*foot;
		//else
    // comment out for testing , Nov 14 2012
	/*	  if (delta/foot >0.3)
		  {
			  if (y>=h/2)
				  temp = foot;
			  else
			  {
				  temp = upfoot *beta + (1-beta)*foot;
			  }
		  }
		  else
			  temp = upfoot*alpha + (1-alpha)*foot;*/
      //change to generic interpolation function May 2 2013
		  //old////*drv -= temp;
      /*disable this, use Neumann condition 20130528*/
      if (feet.size()>0)
      {
      if (interpRight[y]>0)
        *drv -= interpRight[y];
      else
        if (tar1.at<double>(y,x+1)>0)
          *drv -= tar1.at<double>(y,x+1);
        else
          *drv -= tar1.at<double>(y,x);
      *drv += cand1.at<double>(y,x); //20130529 add difference back
      }
      //*drv -=  -(cand1.at<double>(y,x+1)-cand1.at<double>(y,x-1));//Neumann//20130528 try Sadeghi's method
      //*drv -= -cand1.at<double>(y,x-1);
      //cout<<"right:"<<tar1.at<double>(y,x+1)<<endl;
      //cout<<"cur:"<<tar1.at<double>(y,x)<<endl;
      //cout<<ref1.at<double>(y,x)<<endl;
      //*drv -= ref1.at<double>(y,x);
      tlrb &= 13; //0b1101
    }
    if(mask1.at<double>(y+1,x)==0) //down =1
    {
  //      //*drv -= tar1.at<UINT8>(y+1,x);
		//double alpha = double(w-2-x)/double(w-3);
		//double beta = double(w-2*x)/double(w-2);
		//double temp;
		//double delta = abs(leftfoot-foot);
		////if (leftfoot<0)
		////	if (leftfoot - foot)
		////	//if (x<w/2)
		////	//	temp = foot;
		////	//else
		////		temp = tar1.at<double>(h-2,0)*alpha + (1-alpha)*foot;
  //  // comment out for testing, Nov 12 2012
		//if (delta/foot >0.3)
		//{
		//	if (x>=w/2)
		//		temp = foot;
		//	else
		//	{
		//		temp = leftfoot *beta + (1-beta)*foot;
		//	}
		//}	
		//else 
		//	temp = leftfoot*alpha + (1-alpha)*foot;
    //change to generic interpolation function May 2 2013
		//old////*drv -= temp;
    if (feet.size()>0)
    {
    /*disable this, use Neumann condition 20130528*/
    if (interpButtom[x]>0)
      *drv -= interpButtom[x];
    else
      if (tar1.at<double>(y+1,x)>0)
        *drv -= tar1.at<double>(y+1,x);
      else
        *drv -= tar1.at<double>(y,x);
     *drv += cand1.at<double>(y,x); //20130529 add difference
    }
    //*drv -=  (cand1.at<double>(y+1,x)-cand1.at<double>(y-1,x));//Neumann//20130528 try Sadeghi's method
    //*drv -= -cand1.at<double>(y-1,x);
    //cout<<"down:"<<tar1.at<double>(y+1,x)<<endl;
    //cout<<"cur:"<<tar1.at<double>(y,x)<<endl;
    //*drv -= ref1.at<double>(y,x);
		tlrb &= 14; //0b1110
  }
    //for Neumann BC the value of equation will change
  for(int k=0; k<ch; ++k) {
    A.startVec(rowA+k);
      if(tlrb&8)
			  A.insertBack(mp[tidx+k], rowA+k) = 1.0; //not top
      if(tlrb&4)
			  A.insertBack(mp[lidx+k], rowA+k) = 1.0; // not left
      A.insertBack(mp[id  +k], rowA+k) = -4.0;// center
      if(tlrb&2)
			  A.insertBack(mp[ridx+k], rowA+k) = 1.0; // not right
      if(tlrb&1)
		 A.insertBack(mp[bidx+k], rowA+k) = 1.0; // not bottom
  //ofstream of;
  //of.open("Aout.txt",ios::out);
  //of<<A.innerVector(rowA);
    //////////////////////////////////////////////////////
    if(!(tlrb&1))// if bottom boundary
    {
      A.coeffRef(mp[tidx+k],rowA+k)=2.0;//be careful , if the element isn't exist, this will insert a value and make it slow
    }
    ///////////////////////////////////
    if (!(tlrb&2)) // if right boundary
    {
		    A.coeffRef(mp[lidx+k], rowA+k) = 2.0;
    }
    //if(!(tlrb&4))// if left boundary
    //{
    //  A.coeffRef(mp[ridx+k],rowA+k)=2.0;//be careful , if the element isn't exist, this will insert a value and make it slow
    //}
    ///////////////////////////////////////
    //if (!(tlrb&8)) // if top boundary
    //{
		  //A.coeffRef(mp[bidx+k], rowA+k) = 2.0;
    //}
    //of<<"---------------------------------------------\n";
    //of<<A.innerVector(rowA);
    //of.close();
  }
	b(rowA) = *drv;
     // b(rowA+0) = cv::saturate_cast<double>((*drv)[0]);
    //  b(rowA+1) = cv::saturate_cast<double>((*drv)[1]);
     // b(rowA+2) = cv::saturate_cast<double>((*drv)[2]);
      rowA+=ch;
    }
  }
  A.finalize();
  CV_Assert(nz==rowA);
 
  return true;
}

template <typename T>
bool PoissonSolver::buildPostLCMatrix(Eigen::SparseMatrix<T> &A, Eigen::Matrix<T, Eigen::Dynamic,1> &b, 
                                 Eigen::Matrix<T, Eigen::Dynamic, 1> &u)
{
  int w = mask_roi1.width;
  int h = mask_roi1.height;
  int nz=0;
  for(int y=0; y<h; ++y) {
    double *p = mask1.ptr<double>(y);
    for(int x=0; x<w; ++x, ++p) {
      if(*p==0) continue;
 
      int id = y*(w*ch)+(x*ch);
      mp[id] = nz++;   // r
     // mp[++id] = nz++; // g
     // mp[++id] = nz++; // b
    }
  }
  //mylib::DisplayMat(ref1,"ref");

  A = Eigen::SparseMatrix<double>(nz, nz);
  b = Eigen::VectorXd(nz);
  u = Eigen::VectorXd(nz);
  int rowA = 0;
  //int bh = h-2;
  //int bw = w-2;
  A.reserve(5*nz);
  //mylib::DisplayMat(mask1,"mask1");
  //mylib::DisplayMat(laplacian_cand,"lap");
  //mylib::DisplayMat(tar1,"tar1");
  for(int y=1; y<=h-1; ++y) {//01062013 change y<h-1 to y<=h-1
    double *p = mask1.ptr<double>(y)+1;
    //cv::Vec3d *drv = laplacian.ptr<cv::Vec3d>(y)+1;
    double* drv = laplacian_cand.ptr<double>(y)+1;//+1????
  //mylib::DisplayMat(mask1);
  //mylib::DisplayMat(_tar);
	for(int x=1; x<=w-1; ++x, ++p, ++drv)//01062013 change x<w-1 to x<=w-1
  {
    if(*p==0) continue;
    int id = y*(w*ch)+(x*ch);
    int tidx=id-ch*w, lidx=id-ch, ridx=id+ch, bidx=id+ch*w; //index of up, left, right, down in sparse matrix A
    //double temp=0;
    uchar tlrb = 15; // 0b1111
    if(mask1.at<double>(y-1,x)==0)//up
    { 
			*drv -= tar1.at<double>(y-1,x);
      *drv += cand1.at<double>(y-1,x);
     // cout<<"up "<<tar1.at<double>(y-1,x)<<endl;
      tlrb &= 7; //0b0111
    }

    if(mask1.at<double>(y,x-1)==0) //left
    {
			*drv -= tar1.at<double>(y,x-1);
      *drv += cand1.at<double>(y-1,x);
     // cout<<"left "<<tar1.at<double>(y,x-1)<<endl;
      tlrb &= 11; //0b1011
    }
	  
    if(mask1.at<double>(y,x+1)==0) //right
    {
      *drv -= tar1.at<double>(y,x+1);
      *drv += cand1.at<double>(y+1,x);
     // cout<<"right "<<tar1.at<double>(y,x+1)<<endl;
      tlrb &= 13; //0b1101
    }
    if(mask1.at<double>(y+1,x)==0) //down
    {
      *drv -= tar1.at<double>(y+1,x);
      *drv += cand1.at<double>(y+1,x);
      //cout<<"down "<<tar1.at<double>(y+1,x)<<endl;
		  tlrb &= 14; //0b1110
    }
    for(int k=0; k<ch; ++k) {
      A.startVec(rowA+k);
      if(tlrb&8)
			  A.insertBack(mp[tidx+k], rowA+k) = 1.0; // top
      if(tlrb&4)
			  A.insertBack(mp[lidx+k], rowA+k) = 1.0; // left
      A.insertBack(mp[id  +k], rowA+k) = -4.0;// center
      if(tlrb&2)
			  A.insertBack(mp[ridx+k], rowA+k) = 1.0; // right
      if(tlrb&1)
			  A.insertBack(mp[bidx+k], rowA+k) = 1.0; // bottom
    
    //if(!(tlrb&1))// if bottom boundary
    //{
    //  A.coeffRef(mp[tidx+k],rowA+k)=2.0;//be careful , if the element isn't exist, this will insert a value and make it slow
    //}
    /////////////////////////////////////
    //if (!(tlrb&2)) // if right boundary
    //{
		  //  A.coeffRef(mp[lidx+k], rowA+k) = 2.0;
    //}
    // if(!(tlrb&4))// if left boundary
    //{
    //  A.coeffRef(mp[ridx+k],rowA+k)=2.0;//be careful , if the element isn't exist, this will insert a value and make it slow
    //}
    /////////////////////////////////////
    //if (!(tlrb&8)) // if top boundary
    //{
		  //  A.coeffRef(mp[bidx+k], rowA+k) = 2.0;
    //}


    }

	  b(rowA) = *drv;
     // b(rowA+0) = cv::saturate_cast<double>((*drv)[0]);
    //  b(rowA+1) = cv::saturate_cast<double>((*drv)[1]);
     // b(rowA+2) = cv::saturate_cast<double>((*drv)[2]);
      rowA+=ch;
    }
  }
  A.finalize();
  CV_Assert(nz==rowA);
 
  return true;
}

 template <typename T>
bool PoissonSolver::buildStichingMatrix(Eigen::SparseMatrix<T> &A, Eigen::Matrix<T, Eigen::Dynamic,1> &b, 
                                 Eigen::Matrix<T, Eigen::Dynamic, 1> &u)
{
  int w = mask_roi1.width;
  int h = mask_roi1.height;
  int nz=0;
 // mylib::DisplayMat(mask1,"mask1");
  for(int y=0; y<h; ++y) {
    double *p = mask1.ptr<double>(y);
    for(int x=0; x<w; ++x, ++p) {
      if(*p==0) continue;
 
      int id = y*(w*ch)+(x*ch);
      mp[id] = nz++;   // r
     // mp[++id] = nz++; // g
     // mp[++id] = nz++; // b
    }
  }
 
  A = Eigen::SparseMatrix<double>(nz, nz);
  b = Eigen::VectorXd(nz);
  u = Eigen::VectorXd(nz);
  int rowA = 0;
  //int bh = h-2;
  //int bw = w-2;
  A.reserve(5*nz);
  for(int y=1; y<h-1; ++y) {
    double *p = mask1.ptr<double>(y)+1;
    //cv::Vec3d *drv = laplacian.ptr<cv::Vec3d>(y)+1;
    double* drv = NI.ptr<double>(y)+1;
	for(int x=1; x<w-1; ++x, ++p, ++drv) {
      if(*p==0) continue;
 
      int id = y*(w*ch)+(x*ch);
      int tidx=id-ch*w, lidx=id-ch, ridx=id+ch, bidx=id+ch*w; //index of up, left, right, down in sparse matrix A
          //double temp=0;
//	  cout<<tar1;
      uchar tlrb = 15; // 0b1111
      if(mask1.at<double>(y-1,x)==0) {
		*drv -= tar1.at<double>(y-1,x);
        tlrb &= 7; //0b0111
      }

      if(mask1.at<double>(y,x-1)==0) {
		*drv -= tar1.at<double>(y,x-1);
        tlrb &= 11; //0b1011
      }
	  
      if(mask1.at<double>(y,x+1)==0) {
		*drv -= dst1.at<double>(y,x+1);
        tlrb &= 13; //0b1101
      }
      if(mask1.at<double>(y+1,x)==0) {
        *drv -= dst1.at<double>(y+1,x);
		tlrb &= 14; //0b1110
      }
      for(int k=0; k<ch; ++k) {
        A.startVec(rowA+k);
        if(tlrb&8)
			A.insertBack(mp[tidx+k], rowA+k) = 1.0; // top
        if(tlrb&4)
			A.insertBack(mp[lidx+k], rowA+k) = 1.0; // left
        A.insertBack(mp[id  +k], rowA+k) = -4.0;// center
        if(tlrb&2)
			A.insertBack(mp[ridx+k], rowA+k) = 1.0; // right
        if(tlrb&1)
			A.insertBack(mp[bidx+k], rowA+k) = 1.0; // bottom
      }
	  b(rowA) = *drv;
     // b(rowA+0) = cv::saturate_cast<double>((*drv)[0]);
    //  b(rowA+1) = cv::saturate_cast<double>((*drv)[1]);
     // b(rowA+2) = cv::saturate_cast<double>((*drv)[2]);
      rowA+=ch;
    }
  }
  A.finalize();
  CV_Assert(nz==rowA);
 
  return true;
}
template <typename T>
bool PoissonSolver::copyResult(Eigen::Matrix<T, Eigen::Dynamic, 1> &u)
{
  int w = mask_roi1.width;
  int h = mask_roi1.height;
  for(int y=1; y<h-1; ++y) {
    double *pd = dst1.ptr<double>(y)+1;
    double *pm = mask1.ptr<double>(y)+1;
    for(int x=1; x<w-1; ++x, ++pm) {
      if(*pm==0) {
        pd ++;
      } else {
        int idx = mp[y*(w*ch)+(x*ch)];
        //*pd++ = cv::saturate_cast<double>(u[idx]);
        *pd = u[idx]+cand1.at<double>(y,x); //20130528, solve lighting itself, and add it to candid
		pd++;
      }
    }
  }
 
  return true;
}
void PoissonSolver::computeLaplacian(cv::Mat& _dst, int offx, int offy)
{
  ch = _tar.channels();
  cv::Point offset(offx, offy);
  cv::Point tl(_mask.size()), br(-1,-1);
  // calc bounding box
  for(int y=0; y<_mask.rows; ++y) {
    double *p = _mask.ptr<double>(y);
    for(int x=0; x<_mask.cols; ++x,++p) {
	  //cout<<*p;
      if(*p==0) continue;
      if(tl.x>x) tl.x=x;
      if(tl.y>y) tl.y=y;
      if(br.x<x) br.x=x;
      if(br.y<y) br.y=y;
    }
  }
  br.x += 1;
  br.y += 1;
 //
  cv::Rect mask_roi(tl, br); //roi of cand and tar
  mask_roi1 = cv::Rect(tl-cv::Point(1,1), br+cv::Point(1,1));  //roi of mask box extend by 1
  cv::Mat _candExt2, _tarExt2, _maskExt1, _dstExt2, _refExt2;
  cv::copyMakeBorder(_cand, _candExt2, 2,2,2,2, cv::BORDER_REPLICATE); //get extend by 2 cand
  cv::copyMakeBorder(_tar, _tarExt2, 2,2,2,2, cv::BORDER_REPLICATE); //get extend by 2 tag
  if (_ref.data!=NULL)//if there is a reference data
    cv::copyMakeBorder(_ref, _refExt2, 2,2,2,2, cv::BORDER_REPLICATE);
 //mylib::DisplayMat(_cand,"cand");
 //mylib::DisplayMat(_tar,"tar");
 //mylib::DisplayMat(_ref,"ref");
 //mylib::DisplayMat(_candExt2,"candext");
 //mylib::DisplayMat(_tarExt2,"tarext");
 cv::Rect mask_roi_ext(tl,br+cv::Point(4,4)); //roi use for locate cand and tar in ext version

 // // allocate destination image
 //_dstExt2 = _tarExt2.clone(); //copy extend by 2 target as dstUp
 //change to 
 _dstExt2 = _candExt2.clone();
 _dst = cv::Mat(_dstExt2, cv::Rect(2,2,_dstExt2.cols-4, _dstExt2.rows-4)); //dst is center by 2 of dstUp
 //
 mask1 = cv::Mat(_mask, mask_roi1);
 tar1 = cv::Mat(_tarExt2, mask_roi1+cv::Point(2,2)); //extend 1 of mask box + offset of target
 cand1 = cv::Mat(_candExt2,mask_roi1+cv::Point(2,2)); 
//mylib::DisplayMat(tar1,"tar1");
//mylib::DisplayMat(cand1,"cand1");
 if (_ref.data!=NULL)
  ref1 =  cv::Mat(_refExt2, mask_roi1+cv::Point(2,2));
 dst1 = cv::Mat(_dstExt2, mask_roi1+cv::Point(2,2)); //extend 1of mask box + offset of dst
// mylib::DisplayMat(dst1,"dst1");
 cv::Mat cand(_candExt2, mask_roi_ext);
 cv::Mat tar(_tarExt2, mask_roi_ext);
 // // calc differential image
 //mylib::DisplayMat(cand,"cand");
 //mylib::DisplayMat(tar,"tar");
 cv::Mat cand64, target64;
 int pw = mask_roi_ext.width-1, ph = mask_roi_ext.height-1;
 cand64 = Mat::zeros(cand.size(),cand.type());
 //cv::copyMakeBorder(cand(Rect(2,2,pw-3,ph-3)).clone(),cand64,2,2,2,2,BORDER_CONSTANT,0);
 //cand64 = cand;
 //apply Dirchlet CV on up and left
 
 //Mat((tar.row(0)-cand.row(0))).copyTo(cand64.row(0)); 
 //Mat((tar.col(0)-cand.col(0))).copyTo(cand64.col(0)); 
 //apply Neumann BC on the right
 //cand64.col(cand64.size().width-2).copyTo(cand64.col(cand64.size().width-1));
 //apply Nueman BC on the buttom
 //cand64.row(cand64.size().height-2).copyTo(cand64.row(cand64.size().height-1));
 //mylib::DisplayMat(cand64);
 target64 = tar;
 //cand.convertTo(cand64, CV_64F);
 //tar.convertTo(target64, CV_64F);
 cv::Rect roi00(0,0,pw,ph), roi10(1,0,pw,ph), roi01(0,1,pw,ph);
 cv::Mat _cand64_00(cand64, roi00), _tar64_00(target64, roi00);
 cv::Mat cand_dx = cv::Mat(cand64, roi10) - _cand64_00;
 cv::Mat cand_dy = cv::Mat(cand64, roi01) - _cand64_00;
 //mylib::DisplayMat(cand64,"cand64",true);
// mylib::DisplayMat(cand_dy,"cand_dy");
 //mylib::DisplayMat(cand_dy,"cand_dy");
 cv::Mat tar_dx = cv::Mat(target64,roi10) - _tar64_00;
 cv::Mat tar_dy = cv::Mat(target64,roi01) - _tar64_00;
 //
 // // laplacian
 int w = pw-1, h = ph-1;
// mylib::DisplayMat(cand1,"cand1");
 //mylib::DisplayMat(tar1,"tar1");
 //laplacian_cand = cv::Mat(cand_dx,cv::Rect(1,0,w,h)) - cv::Mat(cand_dx,cv::Rect(0,0,w,h))
 //  + cv::Mat(cand_dy,cv::Rect(0,1,w,h)) - cv::Mat(cand_dy,cv::Rect(0,0,w,h));
 
 Mat ker = (Mat_<double>(3,3) <<0,1,0,1,-4,1,0,1,0);
 cv::filter2D(cand64,laplacian_cand,-1,ker);
 laplacian_cand = laplacian_cand(Rect(2,2,w,h));
 //mylib::DisplayMat(laplacian_cand,"lap_cand");
 laplacian_tar = cv::Mat(tar_dx,cv::Rect(1,0,w,h)) - cv::Mat(tar_dx,cv::Rect(0,0,w,h))
 
	 + cv::Mat(tar_dy,cv::Rect(0,1,w,h))- cv::Mat(tar_dy,cv::Rect(0,0,w,h));
 if (_mask.size().height/2==DEBUG_SIZE&&offx == DEBUG_X && offy == DEBUG_Y)
 {
	 mylib::DisplayMat(laplacian_cand,"candlap",true);
	 mylib::DisplayMat(laplacian_tar,"tarlap",true);
 }

}

void PoissonSolver::gradientStiching(cv::Mat& _dst, int offx, int offy)
{
	computeLaplacian(_dst,offx,offy);
	int T=1000;
	cv::Mat N_old,W,temp;
	NI = cv::Mat::ones(laplacian_cand.size(),laplacian_cand.type());
	W=cv::Mat::zeros(laplacian_cand.size(),laplacian_cand.type());
//	double dev=DBL_MAX;
	//cv::Mat dstCln = dst1.clone();
	//mylib::DisplayMat(mask1);
	for (int iter=1; iter< T; iter++)
	{
	    N_old = NI.clone();
		for (int i=2; i< laplacian_cand.size().height-1; i++)
		  for( int j=2; j< laplacian_cand.size().width-1;j++)
		  {
            if (mask1.at<double>(i,j)>0)
			{
			  W.at<double>(i,j)= (laplacian_cand.at<double>(i,j)*laplacian_cand.at<double>(i,j)-NI.at<double>(i,j)*laplacian_cand.at<double>(i,j)+NI.at<double>(i,j)*laplacian_tar.at<double>(i,j))/(laplacian_tar.at<double>(i,j)*laplacian_tar.at<double>(i,j)+laplacian_cand.at<double>(i,j)*laplacian_cand.at<double>(i,j));
			}
		  }
		for (int i=2; i< laplacian_cand.size().height-1; i++)
		  for( int j=2; j< laplacian_cand.size().width-1;j++)
		  {
            if (mask1.at<double>(i,j)>0)
			{
			  NI.at<double>(i,j)=0.5*((1-W.at<double>(i,j))*laplacian_cand.at<double>(i,j)+W.at<double>(i,j)*laplacian_tar.at<double>(i,j));
			}
		  }
		//mylib::DisplayMat(NI);
		//mylib::DisplayMat(N_old);
		double dev = cv::norm(N_old,NI,NORM_L2);
		if (dev <10)
			break;
	}
	//solve possion equation with full boundary condition
	//mylib::DisplayMat(dst);
	//mylib::DisplayMat(NI,"NI");
	//mylib::DisplayMat(W,"W");
	 Eigen::SparseMatrix<double> A;
	 Eigen::VectorXd b;
	 Eigen::VectorXd u;
		
	buildStichingMatrix<double>(A, b, u);
  //for (int i=0; i< 256;i++)
  //{
	 // if (i%16==0)
		//  cout<<";\n";
	 // cout<<u[i]<<", ";
  //}

  // solve sparse linear system
	solve<double>(A, b, u);
	copyResult<double>(u);
	// dst1(Rect() = dstCln(mask1);
}
bool PoissonSolver::PostPLC(cv::Mat &_dst)
{
  computeLaplacian(_dst);
 // // solve an poisson's equation
  Eigen::SparseMatrix<double> A;
  Eigen::VectorXd b;
  Eigen::VectorXd u;
// mylib::DisplayMat(mask1,"mask1");
  // build right-hand and left-hand matrix

  buildPostLCMatrix<double>(A, b, u);
  //for (int i=0; i< 256;i++)
  //{
	 // if (i%16==0)
		//  cout<<";\n";
	 // cout<<u[i]<<", ";
  //}

  // solve sparse linear system
  solve<double>(A, b, u);

  //for (int i=0; i< 256;i++)
  //{
	 // if (i%16==0)
		//  cout<<";\n";
	 // cout<<u[i]<<", ";
  //}
  // copy computed result to destination image
  copyResult<double>(u);
 
  return true;
}
bool
PoissonSolver::poissonLightCorrection(cv::Mat &_dst, const int offx, const int offy)
{
  // 
//  ch = _tar.channels();
//  cv::Point offset(offx, offy);
//  cv::Point tl(_mask.size()), br(-1,-1);
//  //this ptr++ operation only move by 8bit ...
//  // calc bounding box
//  for(int y=0; y<_mask.rows; ++y) {
//    double *p = _mask.ptr<double>(y);
//    for(int x=0; x<_mask.cols; ++x,++p) {
//	  //cout<<*p;
//      if(*p==0) continue;
//      if(tl.x>x) tl.x=x;
//      if(tl.y>y) tl.y=y;
//      if(br.x<x) br.x=x;
//      if(br.y<y) br.y=y;
//    }
//  }
//  br.x += 1;
//  br.y += 1;
// //
// cv::Rect mask_roi(tl, br); //roi of cand and tar
// mask_roi1 = cv::Rect(tl-cv::Point(1,1), br+cv::Point(1,1));  //roi of mask box extend by 1
// cv::Mat _candExt2, _tarExt2, _maskExt1, _dstExt2;
// cv::copyMakeBorder(_cand, _candExt2, 2,2,2,2, cv::BORDER_REPLICATE); //get extend by 2 cand
// cv::copyMakeBorder(_tar, _tarExt2, 2,2,2,2, cv::BORDER_REPLICATE); //get extend by 2 tag
// //mylib::DisplayMat(_cand,"cand");
// //mylib::DisplayMat(_tar,"tar");
// //mylib::DisplayMat(_candExt2,"candext");
// //mylib::DisplayMat(_tarExt2,"tarext");
// cv::Rect mask_roi_ext(tl,br+cv::Point(4,4)); //roi use for locate cand and tar in ext version
//
// // // allocate destination image
// //_dstExt2 = _tarExt2.clone(); //copy extend by 2 target as dstUp
// //change to 
// _dstExt2 = _candExt2.clone();
// _dst = cv::Mat(_dstExt2, cv::Rect(2,2,_dstExt2.cols-4, _dstExt2.rows-4)); //dst is center by 2 of dstUp
// //
// mask1 = cv::Mat(_mask, mask_roi1);
// tar1 = cv::Mat(_tarExt2, mask_roi1+cv::Point(2,2)); //extend 1 of mask box + offset of target
// //mylib::DisplayMat(tar1,"tar1");
// dst1 = cv::Mat(_dstExt2, mask_roi1+cv::Point(2,2)); //extend 1of mask box + offset of dst
//// mylib::DisplayMat(dst1,"dst1");
// cv::Mat cand(_candExt2, mask_roi_ext);
// cv::Mat tar(_tarExt2, mask_roi_ext);
// // // calc differential image
////mylib::DisplayMat(cand,"cand");
// //mylib::DisplayMat(tar,"tar");
// cv::Mat cand64, target64;
// int pw = mask_roi_ext.width-1, ph = mask_roi_ext.height-1;
// cand64 = cand;
// target64 = tar;
// //cand.convertTo(cand64, CV_64F);
// //tar.convertTo(target64, CV_64F);
// cv::Rect roi00(0,0,pw,ph), roi10(1,0,pw,ph), roi01(0,1,pw,ph);
// cv::Mat _cand64_00(cand64, roi00), _tar64_00(target64, roi00);
// cv::Mat cand_dx = cv::Mat(cand64, roi10) - _cand64_00;
// cv::Mat cand_dy = cv::Mat(cand64, roi01) - _cand64_00;
// //mylib::DisplayMat(cand_dx,"cand_dx");
// //mylib::DisplayMat(cand_dy,"cand_dy");
// cv::Mat tar_dx = cv::Mat(target64,roi10) - _tar64_00;
// cv::Mat tar_dy = cv::Mat(target64,roi01) - _tar64_00;
// //
// // // laplacian
// int w = pw-1, h = ph-1;
// laplacian_cand = cv::Mat(cand_dx,cv::Rect(1,0,w,h)) - cv::Mat(cand_dx,cv::Rect(0,0,w,h))
//   + cv::Mat(cand_dy,cv::Rect(0,1,w,h)) - cv::Mat(cand_dy,cv::Rect(0,0,w,h));
// laplacian_tar = cv::Mat(tar_dx,cv::Rect(1,0,w,h)) - cv::Mat(tar_dx,cv::Rect(0,0,w,h))
//	 + cv::Mat(tar_dy,cv::Rect(0,1,w,h))- cv::Mat(tar_dy,cv::Rect(0,0,w,h));
// if (offx == DEBUG_X && offy == DEBUG_Y)
// {
//	 mylib::DisplayMat(laplacian_cand,"candlap");
//	 mylib::DisplayMat(tar1,"tar1");
// }

	// compute laplacian and prepare mask
  computeLaplacian(_dst,offx,offy);
 // // solve an poisson's equation
  Eigen::SparseMatrix<double> A;
  Eigen::VectorXd b;
  Eigen::VectorXd u;
// mylib::DisplayMat(mask1,"mask1");
  // build right-hand and left-hand matrix

  buildLCMatrix<double>(A, b, u);
  //for (int i=0; i< 256;i++)
  //{
	 // if (i%16==0)
		//  cout<<";\n";
	 // cout<<u[i]<<", ";
  //}

  // solve sparse linear system
  solve<double>(A, b, u);

  //for (int i=0; i< 256;i++)
  //{
	 // if (i%16==0)
		//  cout<<";\n";
	 // cout<<u[i]<<", ";
  //}
  // copy computed result to destination image
  copyResult<double>(u);
 
  return true;
}

void PoissonSolver::test(void)
{
	Eigen::SparseMatrix<double> mat(100,100);
	Eigen::SparseVector<double> vec(1000);

	for (int k=0; k<mat.outerSize(); ++k)
		for (Eigen::SparseMatrix<double>::InnerIterator it(mat,k); it; ++it)
		{
			it.value();
			it.row();
			it.col();
			it.index();
		}
		//Eigen::SuperLU<Eigen::SparseMatrix<double>> slu;
}

void PoissonSolver::poissonLightChange(cv::Mat &_dst)
{
  computeLaplacian(_dst,0,0);
 // // solve an poisson's equation
  Eigen::SparseMatrix<double> A;
  Eigen::VectorXd b;
  Eigen::VectorXd u;
  buildPostLCMatrix<double>(A, b, u);
  // solve sparse linear system
  solve<double>(A, b, u);
  copyResult<double>(u);
}

vector<double> PoissonSolver::interpFeet(vector<FootType>& feet, Point spos, Point epos)//spos start point, epos end point, the interpolate line is between two points
{
  
  int len = (int)ceil(sqrt(double(spos.x - epos.x)* double(spos.x - epos.x)+ double(spos.y - epos.y)* double(spos.y - epos.y)))+1;
  double a = double(epos.x-spos.x);
  double c = double(epos.y*spos.x-spos.y*epos.x); 
  double b= double(spos.y-epos.y);
  //double d = sqrt(double(spos.x)*double(spos.x)+double(spos.y)*double(spos.y));
  //double r = sqrt(double(epos.x)*double(epos.x)+double(epos.y)*double(epos.y));
  //double theta = atan2(double(spos.x),double(spos.y));
 // double phi = atan2(double(epos.x),double(epos.y));
  vector<double> rst(len);
  for (auto it = rst.begin();it!=rst.end(); it++)
    *it = -255;
  vector<int> idx;
  for (FootType& f : feet)
  {
    //cout<<f.x<<","<<f.y<<endl;
    if (abs(double(f.y)*a+b*double(f.x)+c) < 0.1) //  a foot in line
    {
      int index = (int)ceil(sqrt(double(f.x - spos.x)* double(f.x - spos.x)+ double(f.y - spos.y)* double(f.y - spos.y)));
      idx.push_back(index);
      rst[index]=f.value;
    }
  }
  //sort
  //cout<<"sorting...."<<idx.size()<<"\n";
  std::sort(idx.begin(),idx.end());
  for (unsigned int i=0; i<idx.size()-1; i++)
  {
   // cout<<rst[idx[i+1]]<<endl;
    //cout<<rst[idx[i]]<<endl;
    double alpha = (rst[idx[i+1]]-rst[idx[i]])/double(idx[i+1]-idx[i]);
    for (int j=idx[i]+1; j<idx[i+1]; j++)
    {
       rst[j] = alpha*(j-idx[i])+rst[idx[i]];
    }
  }
  return rst;
}

void PoissonSolver::addFoot(const FootType& f)
{
  //the coordinate x, y must match the coordinate of h and w in the buildLCMatrix
  //that means it is the content size + 2 mask size
  CV_Assert(f.x >= mask_roi1.y && f.x < mask_roi1.y+mask_roi1.height);
  CV_Assert(f.y >=mask_roi1.x && f.y < mask_roi1.x+mask_roi1.width);
  feet.push_back(f);
}

}
