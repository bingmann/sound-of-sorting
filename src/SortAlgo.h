/******************************************************************************
 * src/SortAlgo.h
 *
 * Implementations is many sorting algorithms.
 *
 * Note that these implementations may not be as good/fast as possible. Some
 * are modified so that the visualization is more instructive.
 *
 * Futhermore, some algorithms are annotated using the mark() and watch()
 * functions from WSortView. These functions add colors to the illustratation
 * and thereby makes the algorithm's visualization easier to explain.
 *
 ******************************************************************************
 * Copyright (C) 2013 Timo Bingmann <tb@panthema.net>
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

// *** List of Sorting Algorithms

struct AlgoEntry
{
    const wxChar* name;
    void        (*func)(class WSortView&);
    const wxChar* text;
};

extern const struct AlgoEntry g_algolist[];
extern const size_t g_algolist_size;

// *** Sorting Algorithms

void SelectionSort(class WSortView& a);
void InsertionSort(class WSortView& a);

void MergeSort(class WSortView& a);

extern const wxChar* g_quicksort_pivot_text[];
enum QuickSortPivotType { PIVOT_FIRST, PIVOT_LAST, PIVOT_MID, PIVOT_RANDOM, PIVOT_MEDIAN3 };
extern QuickSortPivotType g_quicksort_pivot;

void QuickSortLR(class WSortView& a);
void QuickSortLL(class WSortView& a);
void QuickSortTernaryLR(class WSortView& a);
void QuickSortTernaryLL(class WSortView& a);

void BubbleSort(class WSortView& a);
void CocktailShakerSort(class WSortView& a);
void CombSort(class WSortView& a);
void GnomeSort(class WSortView& a);
void OddEvenSort(class WSortView& a);

void ShellSort(WSortView& a);
void HeapSort(class WSortView& a);
void SmoothSort(class WSortView& a);

void BitonicSort(WSortView& a);

void RadixSortLSD(class WSortView& a);
void RadixSortMSD(class WSortView& a);

void StlSort(class WSortView& a);
void StlStableSort(class WSortView& a);
void StlHeapSort(class WSortView& a);

void TimSort(class WSortView& a);

void BogoSort(class WSortView& a);
void BozoSort(class WSortView& a);
void StoogeSort(class WSortView& a);
void SlowSort(class WSortView& a);

#endif // SORTALGO_H
