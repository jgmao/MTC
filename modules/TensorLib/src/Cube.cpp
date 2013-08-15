#include "Cube.h"

template <typename Tp>
Cube_<Tp>::Cube_(void):Point3_<Tp>(),Size3_<Tp>()
{
}

template <typename Tp>
Cube_<Tp>::~Cube_(void)
{
}

template <typename Tp>
Cube_<Tp>::Cube_(Tp x, Tp y, Tp z, Tp height, Tp width, Tp depth):Point3_<Tp>(x,y,z),Size3_<Tp>(height,width,depth)
{

}
template <typename Tp>
Cube_<Tp>::Cube_(const Point3_<Tp>& pos, const Size3_<Tp>& sz):Point3_<Tp>(pos),Size3_<Tp>(sz)
{

}

template <typename Tp>
Cube_<Tp>::Cube_(const Rect_<Tp>& r): Point3_<Tp>(r.y,r.x,0), Size3_<Tp>(r.height,r.width,1)
{
}

template <typename Tp>
Cube_<Tp>::Cube_(const Cube_<Tp>& c): Point3_<Tp>(c.x,c.y,c.z),Size3_<Tp>(c.height,c.width,c.depth)
{

}

template <typename Tp>
Cube_<Tp>& Cube_<Tp>::operator = (const Cube_<Tp>& c)
{
	this->x = c.x;
	this->y = c.y;
	this->z = c.z;
	this->height = c.height;
	this->width = c.width;
	this->depth = c.depth;
	return *this;
}

template <typename Tp>
Point3_<Tp> Cube_<Tp>::offset(void) const
{
	return Point3_<Tp>(this->x,this->y,this->z);
}

template <typename Tp>
Size3_<Tp> Cube_<Tp>::size(void) const
{
	return Size3_<Tp>(this->height,this->width,this->depth);
}

template <typename Tp> inline
	void Cube_<Tp>::Print(void) const
{
	cout << this->offset();
	this->size().Print();
}

template class Cube_<int>;
