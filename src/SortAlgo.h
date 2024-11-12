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
void DoubleSelectionSort(class SortArray& a);
void SandpaperSort(class SortArray& a);
void DoubleSandpaperSort(class SortArray& a);
void InsertionSort(class SortArray& a);
void BinaryInsertionSort(class SortArray& a);
void BinaryInsertSort(class SortArray& a, size_t start, size_t end);

void MergeSort(class SortArray& a);
void MergeSortIterative(class SortArray& a);
void PairwiseSort(class SortArray& a);
void PairwiseIterativeSort(class SortArray& a);
void WeaveMergeSort(class SortArray& a);
void StrandSort(class SortArray& a);
void NewShuffleMergeSort(class SortArray& a);
void AndreyMergeSort(class SortArray& a);
void ProportionMergeSort(class SortArray& a);
void BufferPartitionMergeSort(class SortArray& a);

wxArrayString QuickSortPivotText();

enum QuickSortPivotType { PIVOT_FIRST, PIVOT_LAST, PIVOT_MID, PIVOT_RANDOM, PIVOT_MEDIAN3, PIVOT_MEDIAN3RANDOM, PIVOT_MEDIAN5, PIVOT_MEDIAN5RANDOM, PIVOT_MEDIAN7, PIVOT_MEDIAN7RANDOM, PIVOT_NINTHER, PIVOT_RANDOMNINTHER, PIVOT_MEDIAN9, PIVOT_MEDIAN9RANDOM, PIVOT_MEDIAN11, PIVOT_MEDIAN11RANDOM, PIVOT_MEDIAN15, PIVOT_MEDIAN21, PIVOT_THREENINTHER };
extern QuickSortPivotType g_quicksort_pivot;

void QuickSortLR(class SortArray& a);
void QuickSortLL(class SortArray& a);
void QuickSortTernaryLR(class SortArray& a);
void QuickSortTernaryLL(class SortArray& a);
void QuickSortDualPivot(class SortArray& a);
void QuickLibrarySort(SortArray& a);

void BubbleSort(class SortArray& a);
void OptimizedBubbleSort(class SortArray& a);
void CocktailShakerSort(class SortArray& a);
void DualCocktailShakerSort(class SortArray& a);
void CombSort(class SortArray& a);
void GnomeSort(class SortArray& a);
void OptimizedGnomeSort(class SortArray& a);
void OddEvenSort(class SortArray& a);
void TargetedBubbleSort(class SortArray& a);
void CircleSort(class SortArray& a);
void CircleSort2(class SortArray& a);
void IntroCircleSort(class SortArray& a);
void IntroIteCircleSort(class SortArray& a);

void ShellSort(SortArray& a);
void HeapSort(class SortArray& a);
void SmoothSort(class SortArray& a);

void BitonicSort(SortArray& a);
void BitonicSortNetwork(SortArray& a);
void BatcherSortNetwork(SortArray& a);

void InPlaceRadixSortLSD(class SortArray& a);
void RadixSortLSD(class SortArray& a);
void RadixSortMSD(class SortArray& a);
void RotateRadixSortLSD(class SortArray& a);
void RotateRadixSortMSD(class SortArray& a);
void AmericanFlagSort(class SortArray& a);

void StlSort(class SortArray& a);
void StlStableSort(class SortArray& a);
void StlHeapSort(class SortArray& a);

void TimSort(class SortArray& a);
void WikiSort(class SortArray& a);
void GrailSort(class SortArray& a);
void AuxGrailSort(class SortArray& a);
void PDQSort(class SortArray& a);
void PDQSortBranchless(class SortArray& a);

void BadSort(class SortArray& a);
void BogoSort(class SortArray& a);
void BozoSort(class SortArray& a);
void StoogeSort(class SortArray& a);
void SlowSort(class SortArray& a);
void PancakeSort(class SortArray& a);
void OptimizedPancakeSort(class SortArray& a);
void AdjacencyPancakeSort(class SortArray& a);
void BeadSort(class SortArray& a);
void GravitySort(class SortArray& a);

void CycleSort(class SortArray& a);

// ****************************************************************************
// *** Iterator Adapter

// iterator based on http://zotu.blogspot.de/2010/01/creating-random-access-iterator.html

class MyIterator {
protected:
    SortArray* m_array;
    size_t m_pos;

public:
    using iterator_category = std::random_access_iterator_tag;
    using value_type = ArrayItem;
    using difference_type = std::ptrdiff_t;
    using reference = ArrayItem&;
    using pointer = ArrayItem*;

    MyIterator() : m_array(nullptr), m_pos(0) {}

    MyIterator(SortArray* A, size_t p) : m_array(A), m_pos(p) {}

    MyIterator(const MyIterator& r) : m_array(r.m_array), m_pos(r.m_pos) {}

    MyIterator& operator=(const MyIterator& r) 
    { m_array = r.m_array; m_pos = r.m_pos; return *this; }

    MyIterator& operator++() 
    { ++m_pos; return *this; }

    MyIterator& operator--() 
    { --m_pos; return *this; }

    MyIterator operator++(int) 
    { return MyIterator(m_array, m_pos++); }

    MyIterator operator--(int) 
    { return MyIterator(m_array, m_pos--); }

    MyIterator operator+(difference_type n) const 
    { return MyIterator(m_array, m_pos + n); }

    MyIterator& operator+=(difference_type n) 
    { m_pos += n; return *this; }

    MyIterator operator-(difference_type n) const 
    { return MyIterator(m_array, m_pos - n); }

    MyIterator& operator-=(difference_type n)
    { m_pos -= n; return *this; }

    reference operator*() const 
    { return m_array->get_mutable(m_pos); }

    pointer operator->() const 
    { return &(m_array->get_mutable(m_pos)); }

    reference operator[](difference_type n) const 
    { return m_array->get_mutable(m_pos + n); }

    bool operator==(const MyIterator& r) const 
    { return (m_array == r.m_array) && (m_pos == r.m_pos); }

    bool operator!=(const MyIterator& r) const 
    { return (m_array != r.m_array) || (m_pos != r.m_pos); }

    bool operator<(const MyIterator& r) const 
    { return (m_array == r.m_array ? (m_pos < r.m_pos) : (m_array < r.m_array)); }

    bool operator>(const MyIterator& r) const 
    { return (m_array == r.m_array ? (m_pos > r.m_pos) : (m_array > r.m_array)); }

    bool operator<=(const MyIterator& r) const 
    { return (m_array == r.m_array ? (m_pos <= r.m_pos) : (m_array <= r.m_array)); }

    bool operator>=(const MyIterator& r) const 
    { return (m_array == r.m_array ? (m_pos >= r.m_pos) : (m_array >= r.m_array)); }

    difference_type operator+(const MyIterator& r2) const
    { ASSERT(m_array == r2.m_array); return (m_pos + r2.m_pos); }

    difference_type operator-(const MyIterator& r2) const
    { ASSERT(m_array == r2.m_array); return (m_pos - r2.m_pos); }
};

#endif // SORTALGO_H
