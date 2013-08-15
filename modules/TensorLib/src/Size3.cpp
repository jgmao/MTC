#include "Size3.h"


template <typename Tp> Size2_< Tp >::Size2_():cv::Size_<Tp>() 
{

}

template <typename Tp> 
Size2_<Tp>::Size2_( Tp height, Tp width ):cv::Size_<Tp>(width, height)
{

}

template <typename Tp> inline
	Size2_<Tp>::Size2_(const Size2_<Tp>& sz)
{
	this->height = sz.height;
	this->width = sz.width;
}

template <typename Tp>
Size2_<Tp>::Size2_(const Size_<Tp>& sz) : cv::Size_<Tp>(sz)
{
}

template <typename Tp> inline
	Size2_<Tp>::Size2_(const cv::Point_<Tp>& pt)
{
	this->height = pt.x;
	this->width = pt.y;
}

template <typename Tp> inline
	Size2_<Tp>& Size2_<Tp>::operator= (const Size2_<Tp>& sz)
{
	this->height = sz.height;
	this->width = sz.width;
	return *this;
}

template <typename Tp> inline
	Size2_<Tp>& Size2_<Tp>::operator= (const Size_<Tp> & sz)
{
	this->height = sz.height;
	this->width = sz.width;
	return *this;
}
//

template <typename Tp> Tp Size2_<Tp>::area(void) const
{
	return cv::Size_<Tp>::area();
}

template <typename Tp> inline Size2_<Tp> Size2_<Tp>::operator* (Tp b)
{
	return Size2_<Tp>(this->height*b, this->width * b); 
}

template <typename Tp> inline Size2_<Tp> Size2_<Tp>::operator* (const Size2_<Tp>& b)
{
	return Size2_<Tp>(this->height*b.height,this->width*b.width);
}

template <typename Tp> inline Size2_<Tp> Size2_<Tp>::operator/ (Tp b)
{
	return Size2_<Tp>(this->height/b, this->width/b);
}

template <typename Tp> inline Size2_<Tp> Size2_<Tp>::operator/ (const Size2_<Tp>& b)
{
	return Size2_<Tp>(this->height/b.height, this->width/b.width);
}

template <typename Tp> inline Size2_<Tp> Size2_<Tp>::operator+ (const Size2_<Tp>& b)
{
	return Size2_<Tp>(this->height+b.height, this->width+b.width);
}

template <typename Tp> inline Size2_<Tp> Size2_<Tp>::operator- (const Size2_<Tp>& b)
{
	return Size2_<Tp>(this->height-b.height, this->width-b.width);
}

template <typename Tp> inline bool Size2_<Tp>::operator== (const Size2_<Tp>& b)
{
	return this->height==b.height && this->width == b.width;
}

template <typename Tp> inline Size2_<Tp> Size2_<Tp>::Floor() const
{
	return Size2_<Tp>(
	saturate_cast<Tp>(floor(saturate_cast<double>(this->height))),
	saturate_cast<Tp>(floor(saturate_cast<double>(this->width))));
}
template <typename Tp> inline Size2_<Tp> Size2_<Tp>::Ceil() const
{
	return Size2_<Tp>(
		saturate_cast<Tp>(ceil(saturate_cast<double>(this->height))),
		saturate_cast<Tp>(ceil(saturate_cast<double>(this->width))));
}
////////////////////////////
template <typename Tp> Size3_<Tp>::Size3_():Size2_<Tp>()
{
	this->depth = 0;
}

template <typename Tp> Size3_<Tp>::Size3_(Tp height, Tp width, Tp depth):Size2_<Tp>(height, width)
{
	this->depth =this->depth;
}

template <typename Tp> inline Size3_<Tp>::Size3_(const Size3_<Tp>&sz)
{
	this->height = sz.height;
	this->width = sz.width;
	this->depth = sz.depth;
}

template <typename Tp> inline Size3_<Tp>::Size3_(const Point3_<Tp>& pt)
{
	this->height = pt.x;
	this->width = pt.y;
	this->depth = pt.z;
}

template <typename Tp> inline Size3_<Tp>& Size3_<Tp>::operator= (const Size3_<Tp>& sz)
{
	this->height =sz.height;
	this->width = sz.width;
	this->depth = sz.depth;
	return *this;
}

template <typename Tp> Size3_<Tp>::Size3_(const Size_<Tp>& cvsz)
{
	this->height = cvsz.height;
	this->width = cvsz.width;
	this->depth = 1;
}

template <typename Tp> Size3_<Tp>::Size3_(const Size2_<Tp>& sz2): Size2_<Tp>(sz2)
{
	this->depth = 1;
}


template<typename Tp> Tp Size3_<Tp>::volumn() const
{
	return this->area()*this->depth;
}

template <typename Tp> inline Size3_<Tp> Size3_<Tp>::operator* (Tp b)
{
	return Size3_<Tp>(this->height*b, this->width * b, this->depth*b); 
}

template <typename Tp> inline Size3_<Tp> Size3_<Tp>::operator* (const Size3_<Tp>& b)
{
	return Size3_<Tp>(this->height*b.height, this->width * b.width, this->depth*b.depth); 
}

template <typename Tp> inline Size3_<Tp> Size3_<Tp>::operator/ (Tp b)
{
	return Size3_<Tp>(this->height/b, this->width / b, this->depth/b); 
}

template <typename Tp> inline Size3_<Tp> Size3_<Tp>::operator/ (const Size3_<Tp>& b)
{
	return Size3_<Tp>(this->height/b.height, this->width/b.width, this->depth/b.depth); 
}

template <typename Tp> inline Size3_<Tp> Size3_<Tp>::operator+ (const Size3_<Tp>& b)
{
	return Size3_<Tp>(this->height+b.height, this->width + b.width, this->depth+b.depth); 
}

template <typename Tp> inline Size3_<Tp> Size3_<Tp>::operator - (const Size3_<Tp>& b) const
{
	return Size3_<Tp>(this->height-b.height, this->width - b.width, this->depth-b.depth); 
}


template <typename Tp> inline bool Size3_<Tp>::operator== (const Size3_<Tp>& b) 
{
	return this->height==b.height && this->width == b.width && this->depth == b.depth;
}

//template <typename Tp> inline Point3_<Tp> Size3_<Tp>::operator-( const Point3_<Tp>& p)
//{
//	return Point3_<Tp>(this->height - p.x, this->width - p.y, this->depth - p.z);
//}

template <typename Tp> inline Point3_<Tp> Size3_<Tp>::Point3(void) const
{
	return Point3_<Tp>(this->height, this->width, this->depth);
}

template <typename Tp> inline void Size3_<Tp>::Print(void) const
{
	cout<<"("<<this->height<<","<<this->width<<","<<this->depth<<")"<<endl;
}

template <typename Tp> inline Size3_<Tp> Size3_<Tp>::Floor(void) const
{
	return Size3_<Tp>(
		saturate_cast<Tp>(floor(saturate_cast<double>(this->height))),
		saturate_cast<Tp>(floor(saturate_cast<double>(this->width))),
		saturate_cast<Tp>(floor(saturate_cast<double>(this->depth))));
}

template <typename Tp> inline Size3_<Tp> Size3_<Tp>::Ceil(void) const
{
	return Size3_<Tp>(
		saturate_cast<Tp>(ceil(saturate_cast<double>(this->height))),
		saturate_cast<Tp>(ceil(saturate_cast<double>(this->width))),
		saturate_cast<Tp>(ceil(saturate_cast<double>(this->depth))));
}
//void operator << (ostream & fd, const Point3i& pos)
//{
//	fd << "(" << pos.x<<","<<pos.y<<","<<pos.z<<")\n";
//}


//explicity instantiation
template class Size2_<int>;
template class Size3_<int>;
template class Size2_<float>;
template class Size3_<float>;
template class Size3_<double>;
template class Size2_<double>;
