/*
 * indexed_list.hxx
 *
 *  Created on: Nov 29, 2016
 *      Author: hooshi
 *
 *  Automated reindexing of objects makes life waaaay easier.
 */

#ifndef INCLUDE_INDEXED_LIST_HXX_
#define INCLUDE_INDEXED_LIST_HXX_

#include <list>
#include <set>
#include <unordered_map>
#include <cassert>

/*
 * All indexed items like vertices, triangles, or half edges
 * must inherit from this type.
 */
template<typename T> class IndexedList;

template<typename DerivedType>
class IndexedItem
{
private:
	friend class IndexedList<DerivedType>;
	int _idx;
	int& pidx() {return _idx;}
	const int& pidx() const {return _idx;}
protected:
	IndexedItem():
		_idx(-1),
		is_marked(false){}
public:

	bool is_marked;

	std::vector<void *> misc;

	// get the index
	int idx() const {return _idx;}

	// sort according to data
	template <typename TYPE, int IDX>
	struct cmp_functor
	{
		bool operator()(const DerivedType* e1, const DerivedType* e2) const
		{
			assert(IDX < int(e1->misc.size()) );
			assert(IDX < int(e2->misc.size()) );
			return (*static_cast<TYPE*>(e1->misc[IDX])) < (*static_cast<TYPE*>(e2->misc[IDX]));
		}
	};
};

/*
 * This container automates renaming indexed items.
 */
template<typename DerivedType>
class IndexedList
{
	static_assert(std::is_base_of<IndexedItem<DerivedType>, DerivedType>::value, "CAN ONLY STORE DESCENDENTS OF INDEXED ITEM.");
	typedef std::list<DerivedType*> LIST;

	LIST _pool;
	std::vector< typename LIST::iterator > _address;

public:

	int n_mems() const;

	// returns the index of the member
	int add_member(DerivedType* member);

	// delete member
	void remove_member(const int);
	void remove_member(DerivedType*);

	// change index
	// no need to support for now.
	// bool change_member_index(const int, const int);
	// bool change_member_index(const Contained*, const int);

	// get members
	bool is_member(DerivedType*) const;
	int get_member_idx(DerivedType*) const;
	DerivedType* get_member_ptr(const int) const;

	// debuggin stuff
	void print_info() const;

	// clear memory
	void clear_memory();

};


/*******************************************************
 *                  Inline functions
 *******************************************************/

template <typename T>
inline int IndexedList<T>::n_mems () const
{
  return _pool.size();
}

template <typename T>
inline int IndexedList<T>::add_member (T* ptr)
{
	assert(!this->is_member(ptr));
	const int index = this->n_mems();
	static_cast< IndexedItem<T>* >(ptr)->pidx() = index;
	typename LIST::iterator it = _pool.insert(_pool.end(), ptr);
	if( index >= int(_address.size())) _address.resize( std::max( index+1, int(_address.size()*2) ));
	_address[index] = it;
	return index;
}

template <typename T>
inline bool IndexedList<T>::is_member(T* ptr) const
{
	if ( static_cast<IndexedItem<T>* >(ptr)->pidx() < 0 ) return false;
	if ( static_cast<IndexedItem<T>* >(ptr)->pidx() >= this->n_mems() ) return false;
	if( *_address[ static_cast<IndexedItem<T>* >(ptr)->pidx() ] == ptr) return true;
	else return false;
}

template <typename T>
inline int IndexedList<T>::get_member_idx(T* ptr) const
{
	assert(this->is_member(ptr));
	return static_cast<IndexedItem<T>* >(ptr)->pidx();
}

template <typename T>
inline T* IndexedList<T>::get_member_ptr(const int idx) const
{
	assert( idx < int(_pool.size()) );
	assert( _address[ idx ] != _pool.end());
	return *_address[ idx ];
}


template <typename T>
inline void IndexedList<T>::remove_member(T* ptr)
{
	assert(this->is_member(ptr));

	const int idx_last = this->n_mems()-1.;
	typename LIST::iterator it_last = _address[idx_last];
	const int idx = static_cast<IndexedItem<T>* >(ptr)->pidx();
	typename LIST::iterator it = _address[idx];

	delete *it;
	_pool.erase(it);
	if(idx != idx_last)
	{
		_address[idx] = it_last;
		static_cast<IndexedItem<T>* >(*it_last)->pidx() = idx;
	}
}

template <typename T>
inline void IndexedList<T>::remove_member(const int idx)
{
	this->remove_member( this->get_member_ptr(idx) );
}

template <typename T>
inline void IndexedList<T>::clear_memory()
{
	for(auto it = _pool.begin() ; it!= _pool.end() ; ++it) delete *it;
	_address.clear();
	_pool.clear();
}

template <typename T>
inline void IndexedList<T>::print_info() const
{
	for(int i = 0 ; i < n_mems() ; i++)
	{
		std::cout << i << " " << *_address[i] << std::endl;
	}
}

#endif /* INCLUDE_INDEXED_LIST_HXX_ */
