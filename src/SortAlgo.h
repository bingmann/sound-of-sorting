/******************************************************************************
 * src/SortAlgo.h
 *
 * Implementations of many sorting algorithms.
 *
 * Note that these implementations may not be as good/fast as possible. Some
 * are modified so that the visualization is more instructive.
 *
 * Futhermore, some algorithms are annotated using the mark() and watch()
 * functions from SortArray. These functions add colors to the illustratation
 * and thereby makes the algorithm's visualization easier to explain.
 *
 ******************************************************************************
 * Copyright (C) 2013-2014 Timo Bingmann <tb@panthema.net>
 *
 * This program is free software: you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by the Free
 * Software Foundation, either version 3 of the License, or (at your option)
 * any later version.
 *
 * This program is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
 * FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for
 * more details.
 *
 * You should have received a copy of the GNU General Public License along with
 * this program.  If not, see <http://www.gnu.org/licenses/>.
 *****************************************************************************/

#ifndef SORTALGO_H
#define SORTALGO_H

#include <wx/string.h>
#include "SortArray.h"

// *** List of Sorting Algorithms

struct AlgoEntry
{
    wxString name;
    void (*func)(class SortArray&);
    // maximum item count for test runs
    unsigned int max_testsize;
    // count inversions if n <= limit
    unsigned int inversion_count_limit;
    wxString text;
};

extern const struct AlgoEntry g_algolist[];
extern const size_t g_algolist_size;
extern const struct AlgoEntry* g_algolist_end;

// *** Sorting Algorithms

void SelectionSort(class SortArray& a);
void InsertionSort(class SortArray& a);
void BinaryInsertionSort(class SortArray& a);

void MergeSort(class SortArray& a);
void MergeSortIterative(class SortArray& a);

wxArrayString QuickSortPivotText();

enum QuickSortPivotType { PIVOT_FIRST, PIVOT_LAST, PIVOT_MID, PIVOT_RANDOM, PIVOT_MEDIAN3 };
extern QuickSortPivotType g_quicksort_pivot;

void QuickSortLR(class SortArray& a);
void QuickSortLL(class SortArray& a);
void QuickSortTernaryLR(class SortArray& a);
void QuickSortTernaryLL(class SortArray& a);
void QuickSortDualPivot(class SortArray& a);

void BubbleSort(class SortArray& a);
void CocktailShakerSort(class SortArray& a);
void CombSort(class SortArray& a);
void GnomeSort(class SortArray& a);
void OddEvenSort(class SortArray& a);

void ShellSort(SortArray& a);
void HeapSort(class SortArray& a);
void SmoothSort(class SortArray& a);

void BitonicSort(SortArray& a);
void BitonicSortNetwork(SortArray& a);
void BatcherSortNetwork(SortArray& a);

void RadixSortLSD(class SortArray& a);
void RadixSortMSD(class SortArray& a);

void StlSort(class SortArray& a);
void StlStableSort(class SortArray& a);
void StlHeapSort(class SortArray& a);

void TimSort(class SortArray& a);
void WikiSort(class SortArray& a);

void BogoSort(class SortArray& a);
void BozoSort(class SortArray& a);
void StoogeSort(class SortArray& a);
void SlowSort(class SortArray& a);

void CycleSort(class SortArray& a);

// ****************************************************************************
// *** Iterator Adapter

// iterator based on http://zotu.blogspot.de/2010/01/creating-random-access-iterator.html

class MyIterator : public std::iterator<std::random_access_iterator_tag, ArrayItem>
{
protected:
    SortArray*  m_array;
    size_t      m_pos;

public:
    typedef std::iterator<std::random_access_iterator_tag, ArrayItem> base_type;

    typedef std::random_access_iterator_tag iterator_category;

    typedef base_type::value_type value_type;
    typedef base_type::difference_type difference_type;
    typedef base_type::reference reference;
    typedef base_type::pointer pointer;

    MyIterator() : m_array(NULL), m_pos(0) {}

    MyIterator(SortArray* A, size_t p) : m_array(A), m_pos(p) {}

    MyIterator(const MyIterator& r) : m_array(r.m_array), m_pos(r.m_pos) {}

    MyIterator& operator=(const MyIterator& r)
    { m_array = r.m_array, m_pos = r.m_pos; return *this; }

    MyIterator& operator++()
    { ++m_pos; return *this; }

    MyIterator& operator--()
    { --m_pos; return *this; }

    MyIterator operator++(int)
    { return MyIterator(m_array, m_pos++); }

    MyIterator operator--(int)
    { return MyIterator(m_array, m_pos--); }

    MyIterator operator+(const difference_type& n) const
    { return MyIterator(m_array, m_pos + n); }

    MyIterator& operator+=(const difference_type& n)
    { m_pos += n; return *this; }

    MyIterator operator-(const difference_type& n) const
    { return MyIterator(m_array, m_pos - n); }

    MyIterator& operator-=(const difference_type& n)
    { m_pos -= n; return *this; }

    reference operator*() const
    { return m_array->get_mutable(m_pos); }

    pointer operator->() const
    { return &(m_array->get_mutable(m_pos)); }

    reference operator[](const difference_type& n) const
    { return m_array->get_mutable(n); }

    bool operator==(const MyIterator& r)
    { return (m_array == r.m_array) && (m_pos == r.m_pos); }

    bool operator!=(const MyIterator& r)
    { return (m_array != r.m_array) || (m_pos != r.m_pos); }

    bool operator<(const MyIterator& r)
    { return (m_array == r.m_array ? (m_pos < r.m_pos) : (m_array < r.m_array)); }

    bool operator>(const MyIterator& r)
    { return (m_array == r.m_array ? (m_pos > r.m_pos) : (m_array > r.m_array)); }

    bool operator<=(const MyIterator& r)
    { return (m_array == r.m_array ? (m_pos <= r.m_pos) : (m_array <= r.m_array)); }

    bool operator>=(const MyIterator& r)
    { return (m_array == r.m_array ? (m_pos >= r.m_pos) : (m_array >= r.m_array)); }

    difference_type operator+(const MyIterator& r2) const
    { ASSERT(m_array == r2.m_array); return (m_pos + r2.m_pos); }

    difference_type operator-(const MyIterator& r2) const
    { ASSERT(m_array == r2.m_array); return (m_pos - r2.m_pos); }
};


#endif // SORTALGO_H
