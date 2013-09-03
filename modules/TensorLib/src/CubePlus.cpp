#include "CubePlus.h"
namespace tensor{
template<typename T, size_t cn> CubePlus_<T,cn>::CubePlus_(void):Cube()
{

}

template<typename T, size_t cn> CubePlus_<T,cn>::CubePlus_(const Cube& c):Cube(c)
{
	content = Tensor<T,cn>(c.size());	
}

template<typename T, size_t cn> CubePlus_<T,cn>::CubePlus_(const Cube& c, const Tensor<T,cn>& src):Cube(c)
{
	content = src.Crop(c.offset(),c.size());
}

template<typename T, size_t cn> CubePlus_<T,cn>::CubePlus_(const Tensor<T,cn>& t):Cube(t.offset(),t.size())
{
	content = t;
}

template<typename T, size_t cn> CubePlus_<T,cn>::CubePlus_(const Tensor<T,cn>& t, const Tensor<T,cn>& extra):Cube(t.offset(),t.size())
{
	content = t;
  extraContent = extra;
}
template<typename T, size_t cn>
void CubePlus_<T,cn>::SetContent(const Tensor<T,cn>& t)
{
	content = t;
}
template<typename T, size_t cn>
void CubePlus_<T,cn>::SetExtraContent(const Tensor<T,cn>& t)
{
	extraContent = t;
}
template<typename T, size_t cn>
const Tensor<T,cn>& CubePlus_<T,cn>::GetContent(void) const
{
	return content;
}
template<typename T, size_t cn>
const Tensor<T,cn>& CubePlus_<T,cn>::GetExtraContent(void) const
{
	return extraContent;
}

template<typename T, size_t cn>
CubePlus_<T,cn>& CubePlus_<T,cn>::operator=(const CubePlus_<T,cn>& c)
{
	this->height=c.height;
	this->width=c.width;
	this->depth=c.depth;
	this->x = c.x;
	this->y = c.y;
	this->z = c.z;
	content = c.content;
  extraContent = c.extraContent;
	return *this;
}

template<typename T, size_t cn>
CubePlus_<T,cn>& CubePlus_<T,cn>::operator=(const Cube& c)
{
	this->height=c.height;
	this->width=c.width;
	this->depth=c.depth;
	this->x = c.x;
	this->y = c.y;
	this->z = c.z;
	return *this;
}

template class CubePlus_<double,1>;
template class CubePlus_<double,2>;
template class CubePlus_<double,3>;
template class CubePlus_<uchar,1>;
template class CubePlus_<uchar,2>;
template class CubePlus_<uchar,3>;
template class CubePlus_<float,1>;
template class CubePlus_<float,2>;
template class CubePlus_<float,3>;
}
