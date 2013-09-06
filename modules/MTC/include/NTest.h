#include "ConsoleColor.h"
#include <exception>

template< class T> class NTest
{
public:
	//EQ, NEAR, NEQ
	class NTEexception: public exception
	{
		virtual const char* what() const throw()
		{
			return "NTest exception happened";
		}
	} nte;

	NTest& ExpectEQ(const T& target, const T& ref, const string& comment="") 
	{
		cout<<yellow<<endl<<"TESTING: "<<comment<<endl;
		if (target!=ref)
			cout<< red<< "[---FAIL!---]: "<<"EXPECT " << ref << ", but result is "<< target;
		else
			cout<<green<<"[---PASS!---]: "<<"EXPECT " << ref;
		cout<< white<<endl;
		return *this;
	}
	NTest& operator<<(const string comment)//, std::ostream& (*fp)(std::ostream &)) const
	{
		cout<<blue<<comment;
		cout<< white;
		return *this;
	}
	NTest& operator<<(const T& comment)//, std::ostream& (*fp)(std::ostream &)) const
	{
		cout<<blue<<comment;
		cout<< white;
		return *this;
	}

	NTest& AssertEQ(const T& target, const T& ref, const string& comment="")
	{
		cout<<yellow<<endl<<"TESTING: "<<comment<<endl;
		if (target!=ref)
		{
			cout<< red<< "[---FAIL!---]: "<<"EXPECT " << ref << " but result is "<< target;
			cout<< white;
			CV_DbgAssert(ref==target);
		}
		else
			cout<<green<<"[---PASS!---]: "<<"EXPECT " << ref;
		cout<< white<<endl;
		return *this;
	} 

	NTest& ExpectNR(const T& target, const T& ref, const string& comment="", const T& eps = 0.001)
	{
		//only work on float/double numbers; 
		cout<<yellow<<endl<<"TESTING: "<<comment<<endl;
		float ft = 0;
		double dt = 0;
		CV_DbgAssert(typeid(target) == typeid(ft) || typeid(target)==typeid(dt));
		if ( ref-target > eps || target-ref > eps)
		{
			cout<< red<< "[---FAIL!---]: "<<"EXPECT NEAR " << ref << ", but diff about"<< target-ref;
		}
		else
		{
			cout<<green<<"[---PASS!---]: "<<"EXPECT NEAR " << ref << ", and diff about"<< target-ref;
		}
		cout<< white<<endl;
		return *this;
	}
	
	//NTest& ASSERT_NR(const T& target, const T& ref, const string& comment="", const T& eps = 0.001)
	//{
	//	//only work on float/double numbers; 
	//	assert(typeof(target)==typeof(float) || typeof(target)==typeof(double));
	//	cout<<yellow<<"TESTING: "<<comment<<endl;
	//	if ( ref-target > eps || target-ref > eps)
	//	{
	//		cout<< red<< "[---FAIL!---]: "<<"EXPECT NEAR" << ref << "but result is different about"<< target-ref <<endl;
	//		assert(ref-target > eps || target-ref > eps);
	//	}
	//	else
	//	{
	//		cout<<green<<"[---PASS!---]: "<<"EXPECT NEAR" << ref <<endl;
	//	}
	//}




};