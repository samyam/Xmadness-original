
//
//
//
#ifndef __MADNESS_MRA_FUSET_FUSETCONTAINER__INCLUDED__
#define __MADNESS_MRA_FUSET_FUSETCONTAINER__INCLUDED__

#include <madness/world/MADworld.h>
#include <madness/world/buffer_archive.h>
#include <madness/tensor/gentensor.h>

using namespace madness;

//
//
//
template<typename T>
class FuseTContainer;

template<typename T>
class Base;

template<typename T>
class FuseT_Type;

template<typename T>
class FuseT_CoeffT;

template<typename T>
class FuseT_VParameter;

template<typename T>
class FuseT_VCoeffT;

// Map types to integers
enum class WHAT_AM_I : int {FuseT_VCoeffT, FuseT_CoeffT, FuseT_VParameter, FuseT_Type, EMPTY};

//
//
//
template <typename T>
struct WhatAmI {};

template<typename T>
struct WhatAmI<FuseT_Type<T>> {static const WHAT_AM_I t=WHAT_AM_I::FuseT_Type;};

template<typename T>
struct WhatAmI<FuseT_CoeffT<T>> {static const WHAT_AM_I t=WHAT_AM_I::FuseT_CoeffT;};

template<typename T>
struct WhatAmI<FuseT_VCoeffT<T>> {static const WHAT_AM_I t=WHAT_AM_I::FuseT_VCoeffT;};

template<typename T>
struct WhatAmI<FuseT_VParameter<T>> {static const WHAT_AM_I t=WHAT_AM_I::FuseT_VParameter;};

// Simple base class
// Do I have to add get() methods for a variety of return types?
template <typename T>
class Base
{
public:
	virtual WHAT_AM_I what() const = 0;
	virtual ~Base() {};

	template<typename Archive>
	void serialize(Archive& ar) {  }
};

//
template<typename T>
class FuseT_Type : public Base<T>
{
public:
	T value;

	FuseT_Type() { }
	FuseT_Type(const FuseT_Type& other) : value (other.value) { }
	FuseT_Type(T &v) { value = v; }
	~FuseT_Type() { }	

	WHAT_AM_I what() const { return WhatAmI<FuseT_Type>::t; };

	template<typename Archive>
	void serialize(Archive& ar) { ar & value; }
};

//
template<typename T>
class FuseT_CoeffT : public Base<T>
{
	typedef GenTensor<T> coeffT;
public:
	coeffT value;
	
	FuseT_CoeffT() { }
	FuseT_CoeffT(const FuseT_CoeffT& other) : value (other.value) { }
	FuseT_CoeffT(coeffT &v) { value = v; }
	~FuseT_CoeffT() { }
	
	WHAT_AM_I what() const { return WhatAmI<FuseT_CoeffT>::t; };

	template<typename Archive>
	void serialize(Archive& ar) { ar & value; }
};

//
template<typename T>
class FuseT_VCoeffT : public Base<T>
{
	typedef GenTensor<T> coeffT;
public:
	std::vector<coeffT> value;

	FuseT_VCoeffT() { value = std::vector<coeffT>(); } 
	FuseT_VCoeffT(int size) { value = std::vector<coeffT>(size); }
	FuseT_VCoeffT(const FuseT_VCoeffT& other) : value (other.value) { }
	FuseT_VCoeffT(std::vector<coeffT> &v) { value = v; }
	~FuseT_VCoeffT() { }
	
	WHAT_AM_I what() const { return WhatAmI<FuseT_VCoeffT>::t; };	

	template<typename Archive>
	void serialize(Archive& ar) { ar & value; }
};

//
template<typename T>
class FuseT_VParameter : public Base<T>
{
public:
	std::vector<FuseTContainer<T>> value;

	FuseT_VParameter() { value = std::vector<FuseTContainer<T>>(); }
	FuseT_VParameter(int size) { value = std::vector<FuseTContainer<T>>(size); }
	FuseT_VParameter(const FuseT_VParameter& other) : value (other.value) { }
	FuseT_VParameter(std::vector<FuseTContainer<T> > &v) { value = v; }
	~FuseT_VParameter() { }
	FuseTContainer<T> operator[](int i){return value[i];}
	WHAT_AM_I what() const { return WhatAmI<FuseT_VParameter>::t; };	

	template<typename Archive>
	void serialize(Archive& ar) { ar & value; }
};



//
//	FuseTContainer
//
template <typename T>
class FuseTContainer
{
	void allocate(WHAT_AM_I t)
	{
		if (t == WHAT_AM_I::FuseT_CoeffT)
			data = static_cast<Base<T>*>(new FuseT_CoeffT<T>);
		else if (t == WHAT_AM_I::FuseT_VCoeffT)
			data = static_cast<Base<T>*>(new FuseT_VCoeffT<T>);
		else if (t == WHAT_AM_I::FuseT_VParameter)
			data = static_cast<Base<T>*>(new FuseT_VParameter<T>);
		else if (t == WHAT_AM_I::FuseT_Type)
			data = static_cast<Base<T>*>(new FuseT_Type<T>);
		else
			data = 0;
	}

public:
	Base<T>* data;

	// Default constructor makes an empty wrapper
	FuseTContainer() : data(0) {}
	
	// Need to handle "data" pointer
	FuseTContainer<T>& operator=(FuseTContainer<T> other)
	{
		FuseT_Type<T>*			copiedFuseTType;
		FuseT_CoeffT<T>*		copiedFuseTCoeffT;
		FuseT_VCoeffT<T>*		copiedFuseTVCoeffT;
		FuseT_VParameter<T>*	copiedFuseTVParameter;

		switch(other.what())
		{	
			case WHAT_AM_I::FuseT_Type:
				//std::cout<<"=Type"<<std::endl;
				copiedFuseTType			= new FuseT_Type<T>();
				copiedFuseTType->value	= ((FuseT_Type<T>*)(other.data))->value;
				this->data				= static_cast<Base<T>*>(copiedFuseTType);
				break;
			case WHAT_AM_I::FuseT_CoeffT:
				//std::cout<<"=CoeffT"<<std::endl;
				copiedFuseTCoeffT			= new FuseT_CoeffT<T>();
				copiedFuseTCoeffT->value	= ((FuseT_CoeffT<T>*)(other.data))->value;
				this->data					= static_cast<Base<T>*>(copiedFuseTCoeffT);
				break;
			case WHAT_AM_I::FuseT_VCoeffT:
				//std::cout<<"=VCoeffT"<<std::endl;
				copiedFuseTVCoeffT			= new FuseT_VCoeffT<T>();
				copiedFuseTVCoeffT->value	= ((FuseT_VCoeffT<T>*)(other.data))->value;								
				this->data					= static_cast<Base<T>*>(copiedFuseTVCoeffT);
				break;
			case WHAT_AM_I::FuseT_VParameter:
				//std::cout<<"=VParameter"<<std::endl;
				copiedFuseTVParameter			= new FuseT_VParameter<T>();
			    copiedFuseTVParameter->value	= ((FuseT_VParameter<T>*)(other.data))->value;
			    this->data						= static_cast<Base<T>*>(copiedFuseTVParameter);
				break;
			case WHAT_AM_I::EMPTY:
				//std::cout<<"=EMPTY"<<std::endl;

				break;
			default:
				std::cout << "=ERROR!!!!" << std::endl;
		}
		return *this;
	}

	FuseTContainer(Base<T>* obj) : data(obj) {}	// !!! TAKES OWNERSHIP OF POINTER
												// should use shared_ptr

	// Returns type identity
	WHAT_AM_I what() const
	{
		if (data) 
			return data->what();
		else
			return WHAT_AM_I::EMPTY;
	}

	Base<T>* get() const
	{
		return data;
	}

	void set(Base<T>* p)
	{
		data = p;
	}

	~FuseTContainer() 
	{	///////////////////////////////////////////////////////	
		// should be memory leak.
		// please check this out
		//if (data != 0)
		//	delet data;
	}//delete data;}	// what about this?

	template <typename Archive>
	static void do_serialize(const Archive& ar, FuseTContainer& w, bool deserialize)
	{
		int t = static_cast<int>(w.what());
		ar & t;
		if (deserialize) w.allocate(static_cast<WHAT_AM_I>(t));

		if (w.what() == WHAT_AM_I::FuseT_CoeffT)			ar & *static_cast<FuseT_CoeffT<T>*>(w.data);
		else if (w.what() == WHAT_AM_I::FuseT_VCoeffT)		ar & *static_cast<FuseT_VCoeffT<T>*>(w.data);
		else if (w.what() == WHAT_AM_I::FuseT_VParameter)	ar & *static_cast<FuseT_VParameter<T>*>(w.data);
		else if (w.what() == WHAT_AM_I::FuseT_Type)			ar & *static_cast<FuseT_Type<T>*>(w.data);
		else if (w.what() == WHAT_AM_I::EMPTY)				ar & *static_cast<Base<T>*>(w.data);
		else	
			std::cout<<__func__<<" wft"<<std::endl;
	}
};

namespace madness
{
	namespace archive
	{
		template <class Archive, typename T>
		struct ArchiveStoreImpl<Archive, FuseTContainer<T> >
		{
			static void store(const Archive& ar, const FuseTContainer<T>& w)
			{
				w.do_serialize(ar, const_cast<FuseTContainer<T>&>(w), false);
			}
		};

		template <class Archive, typename T>
		struct ArchiveLoadImpl<Archive, FuseTContainer<T> >
		{
			static void load(const Archive& ar, FuseTContainer<T>& w)
			{
				w.do_serialize(ar, w, true);
			}
		};
	}
}

#endif //
