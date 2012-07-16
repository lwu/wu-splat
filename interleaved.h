#ifndef INTERLEAVED_H
#define INTERLEAVED_H

//
// Quasi-STL-compliant interleaved array (random access container)
//
// Leslie Wu
//

// Accesses element of type T from aggregate structure
// STRUCT, where OFFSET is offsetof(STRUCT, <T's name>)
template <typename T, typename STRUCT>
T& TfromSTRUCT(STRUCT& elem, size_t OFFSET) 
{	
	void* data = (void*)((size_t)&elem + OFFSET);
	T* pt = (T*)data;
	return (T&)(*pt);
}

template <typename T, typename S, size_t OFFSET>
class InterleavedArray 
{
public:
	explicit InterleavedArray(S* an_array, S* the_end) 
	{
		array_begin = an_array;
		array_end = the_end;
	}

	T& operator[](int index) {		
		return TfromSTRUCT<T,S>( array_begin[index], OFFSET );
	}

	// STL container traits
	typedef T value_type;
	typedef T& reference;
	typedef const T& const_reference;
	typedef T* pointer;
	typedef size_t difference_type;
	typedef size_t size_type_type;

	class iterator {
	public:
		explicit iterator(S* s) : ptr(s) {}

		bool operator==(const iterator& rhs) {
			return (rhs.ptr == ptr);
		}

		bool operator!=(const iterator& rhs) {
			return !(*this == rhs);
		}


		iterator& operator++() {
			++ptr;
			return (*this);
		}

		iterator& operator--() {
			--ptr;
			return (*this);
		}

		reference operator*() const {
			return TfromSTRUCT<T,S>(*ptr, OFFSET);
		}
	private:
		S* ptr;
	};

	// STL container functions
	iterator begin() {
		return iterator(array_begin);
	}

	iterator end() {
		return iterator(array_end);
	}

	size_t size() {
		return distance(array_begin, array_end);
	}
private:
	S* array_begin;
	S* array_end;
};

#endif // INTERLEAVED_H
