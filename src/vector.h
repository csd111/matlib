//
//  vector.h
//
//  Copyright (c) 2015 Clement DOIRE. All rights reserved.
//
// This program is free software; you can redistribute it and/or modify it under
// the terms of the GNU General Public License as published by the Free Software
// Foundation; either version 2 of the License, or (at your option) any later 
// version.
//
// This program is distributed in the hope that it will be useful, but WITHOUT 
// ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
// FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details
//
// You can obtain a copy of the GNU General Public License from
// http://www.gnu.org/copyleft/gpl.html or by writing to Free Software 
// Foundation, Inc.,675 Mass Ave, Cambridge, MA 02139, USA.

#ifndef MCL_VECTOR_H_
#define MCL_VECTOR_H_

#include <iostream>
#include <memory>
#include <random>
#include <algorithm>
#include <fstream>
#include <ctime>

namespace matlib {

template <class T> class Vector {
public:
    typedef T* iterator;
    typedef const T* const_iterator;
    typedef size_t size_type;
    typedef T value_type;
    
    // Constructors
    Vector() {Create();};
    explicit Vector(size_type n, const T& val = T()) {Create(n, val);};
    // Copy constructor
    Vector(const Vector& v){Create(v.begin(),v.end());};
    // Destructor
    ~Vector() {Uncreate();};
    // Size member function
    size_type size() const {return avail - data;};
    // Overloaded operators -- use of i-1 to have the same index referencing 
    // as Matlab (i.e. starts at 1)
    T& operator[](size_type i) { return data[i-1]; };
    const T& operator[](size_type i) const {return data[i-1];};
    // Assignement operator
    Vector& operator=(const Vector&);
    // Begin iterator
    iterator begin() {return data;};
    const_iterator begin() const {return data;};
    // End iterator
    iterator end() {return avail;};
    const_iterator end() const {return avail;};
    // Dynamic allocation
    void push_back(const T& val) {
        if (avail == limit) { // get space if needed
            Grow(2);
        }
        Unchecked_append(val) ;
    };
    // Clear a vector
    inline void clear() {Uncreate();};
    // Remove element(s) from the vector
    void erase(iterator, iterator);
    void erase(iterator it) {erase(it,it);};
    void erase(size_type st_idx, size_type end_idx) {
        erase(data+st_idx-1,data+end_idx-1);
    };
    void erase(size_type indx) {erase(data+indx-1,data+indx-1);};
    // Remove the last element
    void pop_back() {
        alloc.destroy(--avail);
        if (avail - data < (limit - data)/4 + 1) {// reduce space if needed
            Shrink();
        }
    };
    // Insert elements at chosen place
    void insert(iterator,iterator,size_type);
    void insert(size_type indx,iterator it,size_type size) {
        insert(data + indx,it,size);
    };
    void insert(size_type indx,value_type val) {insert(data+indx,&val,1);};
    // Resize the vector
    void resize(size_type newSize, const value_type& value);
    void resize(size_type newSize){resize(newSize,(value_type)0);};
    // Print the elements
    void print(std::ostream& out = std::cout);
    // Create random vector --- !!!! Complex types not permitted !!!!
    void rand(size_type size);  // uniformly distributed
    void randn(size_type size); // normally distributed
    // Returns the sum of all the elements of the vector
    T sum();
    // Returns the mean value of all the elements of the vector
    T mean();
private:
    iterator data;  // first element in the vector
    iterator avail; // one past the last element in the vector
    iterator limit; // one past the allocated memory
    
    // facilities for memory allocation
    std::allocator<T> alloc; // object to handle memory allocation
    // Allocate and initialize the underlying array
    void Create() ;
    void Create(size_type, const T&);
    void Create(const_iterator, const_iterator);
    // Initialization without copy
    void Create(size_type, iterator);
    // Destroy the elements in the array and free the memory
    void Uncreate();
    // Support functions for push_back
    void Grow(size_type);
    void Unchecked_append(const T&);
    // Support function for Erase
    void Shrink();
};


/* -------------------------------------------------------------------------- */
/* -------------------------------------------------------------------------- */
/* ----------------- Definitions of the member functions -------------------- */
/* -------------------------------------------------------------------------- */
/* -------------------------------------------------------------------------- */

/* ------------------------------- Create ----------------------------------- */
template <class T>
void Vector<T>::Create() {
    data = avail = limit = 0;
}

template <class T>
void Vector<T>::Create(size_type n, const T& val) {
    data = alloc.allocate(n);
    limit = avail = data + n;
    std::uninitialized_fill(data, limit, val);
}

template <class T>
void Vector<T>::Create(const_iterator i, const_iterator j) {
    data = alloc.allocate(j - i);
    limit = avail = std::uninitialized_copy(i, j, data);
}

/* ------------------------------- Uncreate --------------------------------- */
template <class T>
void Vector<T>::Uncreate() {
    if (data) {
        // Destroy (in reverse order) the elements that were constructed
        iterator it = avail;
        while (it != data)
            alloc.destroy(--it);
        // Release all the space that was allocated
        alloc.deallocate(data, limit - data);
    }
    // Reset pointers to indicate that the Vec is empty again
    data = limit = avail = 0;
}

/* --------------------------------- Grow ----------------------------------- */
template <class T>
void Vector<T>::Grow(size_type factor) {
    // When growing, allocate 'factor' as much space as currently in use
    size_type new_size = 
        std::max<size_type>(factor * (limit - data), ptrdiff_t(1));
    // Allocate new space and copy existing elements to the new space
    iterator new_data = alloc.allocate(new_size);
    iterator new_avail = std::uninitialized_copy(data, avail, new_data);
    // Release the old space
    Uncreate();
    // Reset pointers to point to the newly allocated space
    data = new_data;
    avail = new_avail;
    limit = data + new_size;
}

/* --------------------------- Unchecked_append ----------------------------- */
// This asssumes 'avail' points at allocated, but uninitialized space
template <class T>
void Vector<T>::Unchecked_append(const T& val) {
    alloc.construct(avail++, val);
}

/* ----------------------------- Assignement -------------------------------- */
template <class T>
Vector<T>& Vector<T>::operator=(const Vector& rhs) {
    // Check for self-assignment
    if (&rhs != this) {
        // Free the array in the left-hand side
        Uncreate();
        // Copy elements from the right-hand to the left-hand side
        Create(rhs.begin(), rhs.end());
    }
    return *this;
}

/* -------------------------------- erase ----------------------------------- */
// If the user inputs wrong indexes, an exception is thrown from the 
// <memory> thread
template <class T>
void Vector<T>::erase(iterator st_bd, iterator end_bd) {
    // Get the new sizes
    size_type rm_length = std::max(end_bd - st_bd, ptrdiff_t(1));
    size_type new_size = avail - data - rm_length;
    // Copy the elements between end_bd and avail at the correct place
    iterator new_avail = std::uninitialized_copy(end_bd+1, avail, st_bd);
    // Destroy backwards until we reach the new avail
    iterator it = avail;
    while (it != new_avail) {
        alloc.destroy(--it);
    }
    avail = new_avail;
    // Check if the sized reduced a lot, if so then shrink the allocated memory
    if (new_size < (limit - data)/4 + 1) {
        Shrink();
    }
}

/* -------------------------------- shrink ---------------------------------- */
template <class T>
void Vector<T>::Shrink() {
    // When shrinking, allocate twice less space as currently in use
    size_type new_size = std::max((limit - data)/2, ptrdiff_t(1));
    // Allocate new space and copy existing elements to the new space
    iterator new_data = alloc.allocate(new_size);
    iterator new_avail = std::uninitialized_copy(data, avail, new_data);
    // Release the old space
    Uncreate();
    // Reset pointers to point to the newly allocated space
    data = new_data;
    avail = new_avail;
    limit = data + new_size;
}

/* --------------------------------- insert --------------------------------- */
template <class T>
void Vector<T>::insert(iterator place, iterator insert_data, size_type size) {
    //'place' points to where the vector should be augmented (new elements are 
    // added just before 'place')
    //'insert_data' points to the data to be inserted
    //'size' is the length of the data to be inserted
    
    // Get the new size
    size_type total_length = size + avail - data;
    // Allocate new space and copy existing elements to the new space
    iterator new_data   = alloc.allocate(total_length);
    iterator avail_tmp  = std::uninitialized_copy(data, --place, new_data);
    iterator avail_tmp2 = std::uninitialized_copy(insert_data, 
                                                  insert_data + size, 
                                                  avail_tmp);
    iterator new_avail  = std::uninitialized_copy(place, avail, avail_tmp2);
    // Return the old space
    Uncreate();
    // Reset pointers to point to the newly allocated space
    data = new_data;
    avail = new_avail;
    limit = data + total_length;
}

/* -------------------------------- resize ---------------------------------- */
template <class T>
void Vector<T>::resize(size_type newSize, const value_type& value) {
    size_type currentSize = avail - data;
    if (newSize < currentSize) {
        // Allocate new space and copy existing elements to the new space
        iterator new_data = alloc.allocate(newSize);
        iterator new_avail = std::uninitialized_copy(data, 
                                                     data + newSize, 
                                                     new_data);
        // Return the old space
        Uncreate();
        // Reset pointers to point to the newly allocated space
        data = new_data;
        avail = new_avail;
        limit = data + newSize;
    }
    else if (newSize > currentSize) {
        // Allocate new space and copy existing elements to the new space
        iterator new_data = alloc.allocate(newSize);
        iterator tmp_avail = std::uninitialized_copy(data, avail, new_data);
        iterator new_avail = std::uninitialized_fill_n(tmp_avail, 
                                                       newSize - currentSize, 
                                                       value);
        // Return the old space
        Uncreate();
        // Reset pointers to point to the newly allocated space
        data = new_data;
        avail = new_avail;
        limit = data + newSize;
    }
    else {
        // do nothing, the vector stays the same
    }
}

/* -------------------------------- rand ------------------------------------ */
// Uniformly distributed - !!!! Complex types not permitted !!!!
template <class T>
void Vector<T>::rand(size_type size){
    // Allocate new space
    iterator new_data = alloc.allocate(size);
    // Initialize the random number generator
    std::default_random_engine generator(time(NULL));
    std::uniform_real_distribution<T> distribution(T(0.0),T(1.0));
    // Sample
    for (int i=0;i<size;++i) {
        new_data[i] = distribution(generator);
    }
    // Release the old space
    Uncreate();
    // Reset pointers to point to the newly allocated space
    data = new_data;
    limit = avail = new_data + size;
}

/* -------------------------------- randn ----------------------------------- */
// Normally distributed - !!!! Complex types not permitted !!!!
template <class T>
void Vector<T>::randn(size_type size){
    // Allocate new space
    iterator new_data = alloc.allocate(size);
    // Initialize the random number generator
    std::default_random_engine generator(time(NULL));
    std::normal_distribution<T> distribution(T(0.0),T(1.0));
    // Sample
    for (int i=0;i<size;++i) {
        new_data[i] = distribution(generator);
    }
    // Release the old space
    Uncreate();
    // Reset pointers to point to the newly allocated space
    data = new_data;
    limit = avail = new_data + size;
}

/* -------------------------------- print ----------------------------------- */
template <class T>
void Vector<T>::print(std::ostream &out) {
    size_type length = avail - data;
    for (int i = 0;i<length;++i) {
        out << *(data + i) << "\n";
    }
}

/* -------------------------------- sum ------------------------------------- */
template <class T>
T Vector<T>::sum() {
    T Sum = 0;
    for (size_type i = 0; i < avail-data; ++i) {
        Sum += data[i];
    }
    return Sum;
}

/* ------------------------------- mean ------------------------------------- */
template <class T>
T Vector<T>::mean() {
    T Sum = 0;
    for (size_type i = 0; i < avail-data; ++i) {
        Sum += data[i];
    }
    return Sum / (avail-data);
}

}

#endif /* MCL_VECTOR_H_ */