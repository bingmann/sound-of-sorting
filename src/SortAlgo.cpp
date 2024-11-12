/******************************************************************************
 * src/SortAlgo.cpp
 *
 * Implementations is many sorting algorithms.
 *
 * Note that these implementations may not be as good/fast as possible. Some
 * are modified so that the visualization is more instructive.
 *
 * Futhermore, some algorithms are annotated using the mark() and watch()
 * functions from SortArray. These functions add colors to the illustratation
 * and thereby makes the algorithm's visualization easier to explain.
 *
 ******************************************************************************
 * The algorithms in this file are copyrighted by the original authors. All
 * code is freely available.
 *
 * The source code added by myself (Timo Bingmann) and all modifications are
 * copyright (C) 2013-2014 Timo Bingmann <tb@panthema.net>
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

#include "SortAlgo.h"

#include <algorithm>
#include <numeric>
#include <limits>
#include <inttypes.h>
#include <random>
#include <vector>
#include <array>
#include <cmath>

typedef ArrayItem value_type;

// inversion count limit for iterator instrumented algorithms
const unsigned int inversion_count_instrumented = 512;

const struct AlgoEntry g_algolist[] =
{
    { _("Selection Sort"), &SelectionSort, UINT_MAX, UINT_MAX,
      wxEmptyString },
    { _("Double Selection Sort"), &DoubleSelectionSort, UINT_MAX, UINT_MAX,
      wxEmptyString },
    { _("Sandpaper Sort"), &SandpaperSort, UINT_MAX, UINT_MAX,
      _("Also known as Exchange Sort.") },
    { _("Double Sandpaper Sort"), &DoubleSandpaperSort, UINT_MAX, UINT_MAX,
      _("A variant of Exchange Sort that sorts the array bidirectionally.") },
    { _("Insertion Sort"), &InsertionSort, UINT_MAX, UINT_MAX,
      wxEmptyString },
    { _("Binary Insertion Sort"), &BinaryInsertionSort, UINT_MAX, UINT_MAX,
      wxEmptyString },
    { _("Merge Sort"), &MergeSort, UINT_MAX, 512,
      _("Merge sort which merges two sorted sequences into a shadow array, and then copies it back to the shown array.") },
    { _("Merge Sort (iterative)"), &MergeSortIterative, UINT_MAX, 512,
      _("Merge sort variant which iteratively merges "
        "subarrays of sizes of powers of two.") },
    { _("Pairwise Merge Sort (Recursive)"), &PairwiseSort, UINT_MAX, 512,
      wxEmptyString },
    { _("Pairwise Merge Sort (Iterative)"), &PairwiseIterativeSort, UINT_MAX, 512,
      wxEmptyString },
    { _("Weave Merge Sort"), &WeaveMergeSort, UINT_MAX, 512,
      _("An in-place merge sort variant that interleaves 2 halves of the input array, and then uses Insertion Sort to sort the array.")},
    { _("New Shuffle Merge Sort"), &NewShuffleMergeSort, UINT_MAX, 512,
      _("An improvement upon Weave Merge Sort, with faster weave time, and inserting now makes comparisons with a worst-case similar to Merge Sort.") },
    { _("Andrey's In-Place Merge Sort"), &AndreyMergeSort, UINT_MAX, 512,
      wxEmptyString },
    { _("Proportion Extend Merge Sort"), &ProportionMergeSort, UINT_MAX, 512,
      wxEmptyString },
    { _("Buffer Partition Merge Sort"), &BufferPartitionMergeSort, UINT_MAX, 512,
      wxEmptyString },
    { _("Strand Sort"), &StrandSort, UINT_MAX, 512,
      wxEmptyString },
    { _("Quick Sort (LR ptrs)"), &QuickSortLR, UINT_MAX, UINT_MAX,
      _("Quick sort variant with left and right pointers.") },
    { _("Quick Sort (LL ptrs)"), &QuickSortLL, UINT_MAX, UINT_MAX,
      _("Quick sort variant from 3rd edition of CLRS: two pointers on left.") },
    { _("Quick Sort (ternary, LR ptrs)"), &QuickSortTernaryLR, UINT_MAX, UINT_MAX,
      _("Ternary-split quick sort variant, adapted from multikey quicksort by "
        "Bentley & Sedgewick: partitions \"=<?>=\" using two pairs of pointers "
        "at left and right, then copied to middle.") },
    { _("Quick Sort (ternary, LL ptrs)"), &QuickSortTernaryLL, UINT_MAX, UINT_MAX,
      _("Ternary-split quick sort variant: partitions \"<>?=\" using two "
        "pointers at left and one at right. Afterwards copies the \"=\" to middle.") },
    { _("Quick Sort (dual pivot)"), &QuickSortDualPivot, UINT_MAX, UINT_MAX,
      _("Dual pivot quick sort variant: partitions \"<1<2?>\" using three pointers, "
        "two at left and one at right.") },
    { _("PDQ Sort"), &PDQSort, UINT_MAX, inversion_count_instrumented,
      _("Pattern-defeating Quick Sort (pdqsort) is a novel sorting algorithm that combines the fast average case of randomized quicksort"
          " with the fast worst case of heapsort, while achieving linear time on inputs with certain patterns.") },
    { _("Branchless PDQ Sort"), &PDQSortBranchless, UINT_MAX, inversion_count_instrumented,
      _("Provides potential speedup over default Pattern-Defeating Quick Sort for arithmetic data.") },
    { _("Flan Sort (Quick Library Sort)"), &QuickLibrarySort, UINT_MAX, inversion_count_instrumented,
      _("An in-place Library Sort variant that uses Quick Sort partitioning to make space.") },
    { _("Bubble Sort"), &BubbleSort, UINT_MAX, UINT_MAX,
      wxEmptyString },
    { _("Optimized Bubble Sort"), &OptimizedBubbleSort, UINT_MAX, UINT_MAX,
      _("This variant terminates early if the array is already sorted.") },
    { _("Targeted Bubble Sort"), &TargetedBubbleSort, UINT_MAX, 1024,
      _("This variant of Bubble Sort is capable of adjusting the sorting boundaries based on the input array.") },
    { _("Cocktail Shaker Sort"), &CocktailShakerSort, UINT_MAX, UINT_MAX,
      wxEmptyString },
    { _("Dual Cocktail Shaker Sort"), &DualCocktailShakerSort, UINT_MAX, UINT_MAX,
      _("This variant sorts from both directions of the array simultaneously.") },
    { _("Circle Sort"), &CircleSort, UINT_MAX, UINT_MAX,
      _("Circle Sort is a recursive sorting algorithm that works by comparing and swapping elements in a circular manner.") },
    { _("Iterative Circle Sort"), &CircleSort2, UINT_MAX, UINT_MAX,
      _("A variant of Circle Sort that has less recursion overhead.") },
    { _("Introspective Circle Sort"), &IntroCircleSort, UINT_MAX, UINT_MAX,
      _("A variant of Circle Sort that switches to Insertion Sort when a certain threshold has been reached.") },
    { _("Introspective Iterative Circle Sort"), &IntroIteCircleSort, UINT_MAX, UINT_MAX,
      _("A variant of Iterative Circle Sort that switches to Insertion Sort when a certain threshold has been reached.") },
    { _("Gnome Sort"), &GnomeSort, UINT_MAX, UINT_MAX,
      wxEmptyString },
    { _("Optimized Gnome Sort"), &OptimizedGnomeSort, UINT_MAX, UINT_MAX,
      _("This variant avoids scanning through sorted portions of the array after a number has been placed in its correct spot") },
    { _("Comb Sort"), &CombSort, UINT_MAX, UINT_MAX,
      wxEmptyString },
    { _("Shell Sort"), &ShellSort, UINT_MAX, 1024,
      wxEmptyString },
    { _("Heap Sort"), &HeapSort, UINT_MAX, inversion_count_instrumented,
      wxEmptyString },
    { _("Smooth Sort"), &SmoothSort, UINT_MAX, 1024,
      wxEmptyString },
    { _("Odd-Even Sort"), &OddEvenSort, UINT_MAX, 1024,
      wxEmptyString },
    // older sequential implementation, which really makes little sense to do
    //{ _("Bitonic Sort"), &BitonicSort, UINT_MAX, UINT_MAX, wxEmptyString },
    { _("Batcher's Bitonic Sort"), &BitonicSortNetwork, UINT_MAX, UINT_MAX,
      wxEmptyString },
    { _("Batcher's Odd-Even Merge Sort"), &BatcherSortNetwork, UINT_MAX, UINT_MAX,
      wxEmptyString },
    { _("Cycle Sort"), &CycleSort, 512, UINT_MAX,
      wxEmptyString },
    { _("Radix Sort (LSD)"), &RadixSortLSD, UINT_MAX, 512,
      _("Least significant digit radix sort, which copies item into a shadow "
        "array during counting.") },
    { _("In-Place Radix Sort (LSD)"), &InPlaceRadixSortLSD, UINT_MAX, UINT_MAX,
      _("Least significant digit radix sort, performed in O(1) space.") },
    { _("Radix Sort (MSD)"), &RadixSortMSD, UINT_MAX, 512,
      _("Most significant digit radix sort, which permutes items in-place by walking cycles.") },
    { _("Rotate Radix Sort (MSD)"), &RotateRadixSortMSD, UINT_MAX, 512,
      wxEmptyString },
    { _("Rotate Radix Sort (LSD)"), &RotateRadixSortLSD, UINT_MAX, 512,
      wxEmptyString },
    { _("American Flag Sort"), &AmericanFlagSort, UINT_MAX, inversion_count_instrumented,
      _("American Flag Sort is an efficient, in-place variant of radix sort that distributes items into hundreds of buckets.") },
    { _("std::sort (gcc)"), &StlSort, UINT_MAX, inversion_count_instrumented,
      wxEmptyString },
    { _("std::stable_sort (gcc)"), &StlStableSort, UINT_MAX, inversion_count_instrumented,
      wxEmptyString },
    { _("std::sort_heap (gcc)"), &StlHeapSort, UINT_MAX, inversion_count_instrumented,
      wxEmptyString },
    { _("Tim Sort"), &TimSort, UINT_MAX, inversion_count_instrumented,
      wxEmptyString },
    { _("Block Merge Sort (WikiSort)"), &WikiSort, UINT_MAX, inversion_count_instrumented,
      _("An O(1) place O(n log n) time stable merge sort.") },
    { _("Grail Sort (O(1) buffer)"), &GrailSort, UINT_MAX, inversion_count_instrumented,
      _("Grail Sort is a stable, in-place sorting algorithm that efficiently organizes an array by using a block-based merging technique.") },
    { _("Grail Sort (external buffer)"), &AuxGrailSort, UINT_MAX, inversion_count_instrumented,
      _("A variant of Grail Sort that uses an external buffer for a potential speedup.") },
    { _("Bead Sort"), &BeadSort, UINT_MAX, UINT_MAX,
      _("This is a non-comparison based sorting algorithm that uses the concept of stacked beads to sort the elements.") },
    { _("Gravity Sort"), &GravitySort, UINT_MAX, UINT_MAX,
      _("A non-comparison based sorting algorithm that uses the concept of gravitational fall to sort the elements.") },
    { _("Pancake Sort"), &PancakeSort, UINT_MAX, UINT_MAX,
      _("Sorts the array by performing a series of 'flips' to push the maximum element to the correct spot.") },
    { _("Optimized Pancake Sort"), &OptimizedPancakeSort, UINT_MAX, UINT_MAX,
      _("An optimized variant of Pancake Sort that performs 1/2 as many flips.") },
    { _("Adjacency Pancake Sort"), &AdjacencyPancakeSort, UINT_MAX, UINT_MAX,
      _("An improvement upon Pancake Sort, which performs only 5/3 N + O(1) flips.") },
    { _("Bogo Sort"), &BogoSort, 10, UINT_MAX,
      wxEmptyString },
    { _("Bozo Sort"), &BozoSort, 10, UINT_MAX,
      wxEmptyString },
    { _("Stooge Sort"), &StoogeSort, 256, inversion_count_instrumented,
      wxEmptyString },
    { _("Slow Sort"), &SlowSort, 128, inversion_count_instrumented,
      wxEmptyString },
    { _("Bad Sort"), &BadSort, 128, inversion_count_instrumented,
      _("A humorous sorting algorithm with a time complexity of O(n^3).") }
};

const size_t g_algolist_size = sizeof(g_algolist) / sizeof(g_algolist[0]);

const struct AlgoEntry* g_algolist_end = g_algolist + g_algolist_size;

std::random_device rd1;
std::mt19937 gen1(rd1());

// ****************************************************************************
// *** Selection Sort

void DoubleSelectionSort(SortArray& A)
{
    size_t left = 0;
    size_t right = A.size() - 1, n = right;
    std::atomic<ssize_t> max_i{ 0 }, low_i{ 0 };
    A.watch(&max_i, 4);
    A.watch(&low_i, 5);
    while (left < right)
    {
        max_i.store(right);
        low_i.store(left);
        for (size_t i = left; i <= right; ++i)
        {
            if (A[i] < A[low_i.load()])
            { 
                A.mark_swap(i, low_i.load()); 
                low_i.store(i);
            }
            else if (A[i] > A[max_i.load()])
            { 
                A.mark_swap(i, max_i.load());
                max_i.store(i);
            }
        }
        A.swap(left, low_i.load());
        ssize_t l = left;  // This removes comparison warning
        if (max_i.load() == l) { max_i.store(low_i.load()); }
        A.swap(right, max_i.load());
        if (left > 0) { A.unmark(left - 1); }
        if (right < n) { A.unmark(right + 1); }
        A.mark(left);
        A.mark(right);
        ++left; --right;
    }
    A.unwatch_all();
}

void SelectionSort(SortArray& A)
{
    std::atomic<ssize_t> kMin{ 0 };
    A.watch(&kMin, 3);

    for (size_t i = 0; i < A.size()-1; ++i)
    {
        kMin.store(i);

        for (size_t j = i+1; j < A.size(); ++j)
        {
            if (A[j] < A[kMin.load()]) {
                A.mark_swap(j, kMin.load());
                kMin.store(j);
            }
        }

        A.swap(i, kMin.load());

        // mark the last good element
        if (i > 0) A.unmark(i-1);
        A.mark(i);
    }
    A.unwatch_all();
}

// ****************************************************************************
// *** Sandpaper and Double Sandpaper Sort
// Double Sandpaper Sort by Taihennami

void DoubleSandpaperSort(SortArray& A)
{
    for (size_t left = 0, right = A.size() - 1; left < right; ++left, --right)
    {
        if (A[left] > A[right])
        { 
            A.swap(left, right); 
        }
        for (size_t i = left + 1; i < right; ++i)
        {
            if (A[left] > A[i])
            {
               A.swap(i, left);
            }
            else if (A[i] > A[right])
            {
               A.swap(i, right);
            }
        }
    }
}

void SandpaperSort(SortArray& A)
{
    size_t n = A.size();
    for (size_t i = 0; i < n; ++i)
    {
        for (size_t j = i + 1; j < n; ++j)
        {
            if (A[i] > A[j]) 
            {
                A.swap(i, j);
            }
        }
    }
}

// ****************************************************************************
// *** Insertion Sort

void InsertSort(SortArray& A, size_t start, size_t end)
{
    int begin = static_cast<int>(start);
    for (size_t i = start; i < end; ++i)
    {
        A.mark(i);
        int j = i - 1;
        value_type key = A[i];
        while (j >= begin && A[j] > key)
        {
            A.set(j + 1, A[j]);
            --j;
        }
        int index = j + 1;
        A.set(index, key);

        A.unmark(i);
    }
}

// with extra item on stack
void InsertionSort(SortArray& A)
{
    InsertSort(A, 0, A.size());
}

// swaps every time (keeps all values visible)
void InsertionSort2(SortArray& A)
{
    for (size_t i = 1; i < A.size(); ++i)
    {
        value_type key = A[i];
        A.mark(i);

        ssize_t j = i - 1;
        while (j >= 0 && A[j] > key)
        {
            A.swap(j, j + 1);
            j--;
        }

        A.unmark(i);
    }
}

// swaps every time (keeps all values visible)
void BinaryInsertionSort2(SortArray& A)
{
    for (size_t i = 1; i < A.size(); ++i)
    {
        value_type key = A[i];
        A.mark(i);

        int lo = 0, hi = i;
        while (lo < hi) {
            int mid = (lo + hi) / 2;
            if (key < A[mid])
                hi = mid;
            else
                lo = mid + 1;
        }

        // item has to go into position lo

        ssize_t j = i - 1;
        while (j >= lo)
        {
            A.swap(j, j + 1);
            j--;
        }

        A.unmark(i);
    }
}

void BinaryInsertSort(SortArray& A, size_t start, size_t end)
{
    for (size_t i = start; i < end; ++i)
    {
        value_type key = A[i];
        A.mark(i);

        size_t lo = start, hi = i;
        while (lo < hi) {
            size_t mid = (lo + hi) / 2;
            if (key < A[mid])
                hi = mid;
            else
                lo = mid + 1;
        }

        size_t j = i;
        while (j > lo)
        {
            A.set(j, A[j - 1]);
            j--;
        }
        A.set(lo, key);
        A.unmark(i);
    }
}

void BinaryInsertionSort(SortArray& A)
{
    BinaryInsertSort(A, 0, A.size());
}

// ****************************************************************************
// *** Merge Sort (out-of-place with sentinels)

// by myself (Timo Bingmann)

void Merge(SortArray& A, size_t lo, size_t mid, size_t hi)
{
    // mark merge boundaries
    A.mark(lo);
    A.mark(mid,3);
    A.mark(hi-1);

    // allocate output
    std::vector<value_type> out(hi-lo);

    // merge
    size_t i = lo, j = mid, o = 0; // first and second halves
    while (i < mid && j < hi)
    {
        // copy out for fewer time steps
        value_type ai = A[i], aj = A[j];

        out[o++] = (ai < aj ? (++i, ai) : (++j, aj));
    }

    // copy rest
    while (i < mid) out[o++] = A[i++];
    while (j < hi) out[o++] = A[j++];

    ASSERT(o == hi-lo);

    A.unmark(mid);

    // copy back
    for (i = 0; i < hi-lo; ++i)
        A.set(lo + i, out[i]);

    A.unmark(lo);
    A.unmark(hi-1);
}

void MergeSort(SortArray& A, size_t lo, size_t hi)
{
    if (lo + 1 < hi)
    {
        size_t mid = (lo + hi) / 2;

        MergeSort(A, lo, mid);
        MergeSort(A, mid, hi);

        Merge(A, lo, mid, hi);
    }
}

void MergeSort(SortArray& A)
{
    return MergeSort(A, 0, A.size());
}

void MergeSortIterative(SortArray& A)
{
    for (size_t s = 1; s < A.size(); s *= 2)
    {
        for (size_t i = 0; i + s < A.size(); i += 2 * s)
        {
            Merge(A, i, i + s,
                  std::min(i + 2 * s, A.size()));
        }
    }
}

/*
    MIT License

    Copyright (c) 2020 aphitorite/2021 EmeraldBlock

    Permission is hereby granted, free of charge, to any person obtaining a copy
    of this software and associated documentation files (the "Software"), to deal
    in the Software without restriction, including without limitation the rights
    to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
    copies of the Software, and to permit persons to whom the Software is
    furnished to do so, subject to the following conditions:

    The above copyright notice and this permission notice shall be included in all
    copies or substantial portions of the Software.

    THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
    IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
    FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
    AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
    LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
    OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
    SOFTWARE.
*/

void weaveMerge2(SortArray& A, std::vector<value_type>& tmp, size_t len, size_t residue, size_t modulus)
{
    if (residue + modulus >= len) { return; }

    size_t low = residue, high = residue + modulus, dmodulus = modulus << 1;

    weaveMerge2(A, tmp, len, low, dmodulus);
    weaveMerge2(A, tmp, len, high, dmodulus);

    size_t next = residue;
    for (; low < len && high < len; next += modulus)
    { 
        int cmp = 0;
        if (A[low] > A[high]) { cmp = 1; }
        else if (A[low] < A[high]) { cmp = -1; }

        if (cmp == 1 || (cmp == 0 && low > high))
        {
            tmp[next] = A[high];
            high += dmodulus;
        }
        else
        {
            tmp[next] = A[low];
            low += dmodulus;
        }
    }

    if (low >= len)
    {
        while (high < len)
        {
            tmp[next] = A[high];
            next += modulus;
            high += dmodulus;
            
        }
    }
    else
    {
        while (low < len)
        {
            tmp[next] = A[low];
            next += modulus;
            low += dmodulus;
        }
    }

    for (size_t i = residue; i < len; i += modulus)
    {
        A.set(i, tmp[i]);
        A[i].get();
    }
}

void WeaveMergeSort2(SortArray& A)
{
    size_t len = A.size();
    std::vector<value_type> tmp(len);
    weaveMerge2(A, tmp, len, 0, 1);
}

void insertTo(SortArray& A, size_t a, size_t b)
{
    value_type temp = A[a];
    while (a > b) { A.set(a, A[a - 1]); --a; }
    A.set(b, temp);
}

void shiftValue(SortArray& A, size_t a, size_t b, size_t len)
{
    for (size_t i = 0; i < len; ++i)
    {
        A[a + i].get();
        A.swap(a + i, b + i);
    }
}

void rotate(SortArray& A, int a, int m, int b)
{
    int l = m - a, r = b - m;
    while (l > 0 && r > 0)
    {
        if (r < l)
        {
            shiftValue(A, static_cast<size_t>(m - r), static_cast<size_t>(m), static_cast<size_t>(r));
            b -= r;
            m -= r;
            l -= r;
        }

        else
        {
            shiftValue(A, static_cast<size_t>(a), static_cast<size_t>(m), static_cast<size_t>(l));
            a += l;
            m += l;
            r -= l;
        }
    }
}

void bitReversal(SortArray& A, size_t a, size_t b)
{
    size_t len = b - a, m = 0;
    size_t d1 = len >> 1, d2 = d1 + (d1 >> 1);
    for (size_t i = 1; i < len - 1; ++i)
    {
        size_t j = d1;
        for (size_t k = i, n = d2; (k & 1) == 0; j -= n, k >>= 1, n >>= 1) {}
        m += j;
        if (m > i)
        {
            A.swap(a + i, a + m);
        }
    }
}

void weaveInsert(SortArray& A, size_t a, size_t b, bool right)
{
    size_t i = a, j = i + 1;
    while (j < b)
    {
        while (i < j && ((right == true && A[i] <= A[j]) || (right == false && A[i] < A[j]))) { ++i; }
        if (i == j)
        {
            right = !right;
            ++j;
        }
        else
        {
            insertTo(A, j, i++);
            j += 2;
        }
    }
}

void weaveMerge(SortArray& A, size_t a, size_t m, size_t b)
{
    if (b - a < 2) { return; }
    size_t a1 = a, b1 = b;
    bool right = true;
    if ((b - a) % 2 == 1)
    {
        if (m - a < b - m)
        {
            --a1;
            right = false;
        }
        else { ++b1; }
    }

    for (size_t e = b1, f; e - a1 > 2; e = f)
    {
        m = (a1 + e) / 2;
        size_t p = 1 << static_cast<size_t>(log(m - a1) / log(2));
        rotate(A, static_cast<int>(m - p), static_cast<int>(m), static_cast<int>(e - p));

        m = e - p;
        f = m - p;

        bitReversal(A, f, m);
        bitReversal(A, m, e);
        bitReversal(A, f, e);
    }
    weaveInsert(A, a, b, right);
}

void WeaveMergeSort(SortArray& A)
{
    size_t n = A.size(), d = 1 << static_cast<size_t>(log(n - 1) / log(2) + 1);
    while (d > 1)
    {
        size_t i = 0, dec = 0;
        while (i < n)
        {
            size_t j = i;
            dec += n;
            while (dec >= d)
            {
                dec -= d;
                ++j;
            }
            size_t k = j;
            dec += n;
            while (dec >= d)
            {
                dec -= d;
                ++k;
            }
            weaveMerge(A, i, j, k);
            i = k;
        }
        d /= 2;
    }
}
// ****************************************************************************
// *** Quick Sort Pivot Selection

QuickSortPivotType g_quicksort_pivot = PIVOT_FIRST;

static ssize_t SingleMedianOfThree(SortArray& A, ssize_t lo, ssize_t mid, ssize_t hi)
{
    // cases if two are equal
    if (A[lo] == A[mid]) return lo;
    if (A[lo] == A[hi - 1] || A[mid] == A[hi - 1]) return hi - 1;

    // cases if three are different
    return A[lo] < A[mid]
        ? (A[mid] < A[hi - 1] ? mid : (A[lo] < A[hi - 1] ? hi - 1 : lo))
        : (A[mid] > A[hi - 1] ? mid : (A[lo] < A[hi - 1] ? lo : hi - 1));
}

template <size_t N>
static void PivotInsertionSort(SortArray& A, std::array<ssize_t, N>& arr)
{
    for (size_t i = 1; i < N; ++i)
    {
        ssize_t key = arr[i];
        size_t j = i;
        while (j > 0 && A[arr[j - 1]] > A[key])
        {
            arr[j] = arr[j - 1];
            --j;
        }
        arr[j] = key;
    }
}

static std::array<size_t, 3> gaps = { 8, 3, 1 };
template<size_t N>
static void PivotShellSort(SortArray& A, std::array<ssize_t, N>& arr)
{
    for (size_t k = 0; k < 3; ++k)
    {
        for (size_t h = gaps[k], i = h; i < N; ++i)
        {
            ssize_t key = arr[i];
            size_t j = i;
            while (j >= h && A[arr[j - h]] > A[key])
            {
                arr[j] = arr[j - h];
                j -= h;
            }
            arr[j] = key;
        }
    }
}

static ssize_t MedianOfFive(SortArray& A, ssize_t lo, ssize_t hi)
{
    if (hi - lo < 25) { return SingleMedianOfThree(A, lo, (lo + hi) / 2, hi); }
    ssize_t segment = (hi - lo) / 5;
    std::array<ssize_t, 5> nums = { lo, lo + segment,  lo + 2 * segment, lo + 3 * segment, lo + 4 * segment };
    PivotInsertionSort(A, nums);
    return nums[2];
}

static ssize_t MedianOfSeven(SortArray& A, ssize_t lo, ssize_t hi)
{
    if (hi - lo < 49) { return SingleMedianOfThree(A, lo, (lo + hi) / 2, hi); }
    ssize_t segment = (hi - lo) / 7;
    std::array<ssize_t, 7> samples = { lo, lo + segment, lo + segment * 2, lo + segment * 3, lo + segment * 4, lo + segment * 5, lo + segment * 6 };
    PivotInsertionSort(A, samples);
    return samples[3];
}

static ssize_t NintherPivot(SortArray& A, ssize_t lo, ssize_t hi)
{
    if (hi - lo < 81) { return SingleMedianOfThree(A, lo, (lo + hi) / 2, hi); }
    ssize_t segment_size = (hi - lo) / 9;
    ssize_t g1 = SingleMedianOfThree(A, lo, lo + segment_size, (lo + segment_size * 2) + 1);
    ssize_t g2 = SingleMedianOfThree(A, lo + 3 * segment_size, lo + 4 * segment_size, (lo + 5 * segment_size) + 1);
    ssize_t g3 = SingleMedianOfThree(A, lo + 6 * segment_size, lo + 7 * segment_size, (lo + 8 * segment_size) + 1);
    return SingleMedianOfThree(A, g1, g2, g3 + 1);
}

// some quicksort variants use hi inclusive and some exclusive, we require it
// to be _exclusive_. hi == array.end()!
ssize_t QuickSortSelectPivot(SortArray& A, ssize_t lo, ssize_t hi)
{
    switch (g_quicksort_pivot)
    {
        case PIVOT_FIRST:
        {
            return lo;
        }
        case PIVOT_LAST:
        {
            return hi - 1;
        }
        case PIVOT_MID:
        {
            return (lo + hi) / 2;
        }
        case PIVOT_RANDOM:
        {
            return lo + (rand() % (hi - lo));
        }
        case PIVOT_MEDIAN3:
        {
            ssize_t mid = (lo + hi) / 2;
            return SingleMedianOfThree(A, lo, mid, hi);
        }
        case PIVOT_MEDIAN3RANDOM:
        {
            if (lo > hi)
            {
                ssize_t temp = lo;
                lo = hi;
                hi = temp;
            }
            std::uniform_int_distribution<ssize_t> dist(lo, hi - 1);
            return SingleMedianOfThree(A, dist(gen1), dist(gen1), dist(gen1) + 1);
        }
        case PIVOT_MEDIAN5:
        {
            return MedianOfFive(A, lo, hi);
        }
        case PIVOT_MEDIAN5RANDOM:
        {
            if (lo > hi)
            {
                ssize_t temp = lo;
                lo = hi;
                hi = temp;
            }
            std::uniform_int_distribution<ssize_t> dist(lo, hi - 1);
            if (hi - lo < 25)
            {
                return SingleMedianOfThree(A, dist(gen1), dist(gen1), dist(gen1) + 1);
            }
            std::array<ssize_t, 5> samples = { dist(gen1), dist(gen1), dist(gen1), dist(gen1), dist(gen1) };
            PivotInsertionSort(A, samples);
            return samples[2];
        }
        case PIVOT_MEDIAN7:
        {
            return MedianOfSeven(A, lo, hi);
        }
        case PIVOT_MEDIAN7RANDOM:
        {
            if (lo > hi)
            {
                ssize_t temp = lo;
                lo = hi;
                hi = temp;
            }
            std::uniform_int_distribution<ssize_t> dist(lo, hi - 1);
            if (hi - lo < 49)
            {
                return SingleMedianOfThree(A, dist(gen1), dist(gen1), dist(gen1) + 1);
            }
            std::array<ssize_t, 7> samples = { dist(gen1), dist(gen1), dist(gen1), dist(gen1), dist(gen1), dist(gen1), dist(gen1) };
            PivotInsertionSort(A, samples);
            return samples[3];
        }
        case PIVOT_NINTHER:
        {
            return NintherPivot(A, lo, hi);
        }
        case PIVOT_RANDOMNINTHER:
        {
            if (lo > hi)
            {
                ssize_t temp = lo;
                lo = hi;
                hi = temp;
            }
            std::uniform_int_distribution<ssize_t> dist(lo, hi - 1);
            if (hi - lo < 81)
            {
                return SingleMedianOfThree(A, dist(gen1), dist(gen1), dist(gen1) + 1);
            }
            ssize_t lo1 = dist(gen1), lo2 = dist(gen1), lo3 = dist(gen1), lo4 = dist(gen1), lo5 = dist(gen1);
            ssize_t lo6 = dist(gen1), lo7 = dist(gen1), lo8 = dist(gen1), lo9 = dist(gen1);
            ssize_t g1 = SingleMedianOfThree(A, lo1, lo2, lo3 + 1);
            ssize_t g2 = SingleMedianOfThree(A, lo4, lo5, lo6 + 1);
            ssize_t g3 = SingleMedianOfThree(A, lo7, lo8, lo9 + 1);
            return SingleMedianOfThree(A, g1, g2, g3 + 1);
        }
        case PIVOT_MEDIAN9:
        {
            if (hi - lo < 81) { return SingleMedianOfThree(A, lo, (lo + hi) / 2, hi); }
            ssize_t segment = (hi - lo) / 9;
            std::array<ssize_t, 9> samples = { lo, lo + segment, lo + segment * 2, lo + segment * 3, lo + segment * 4, lo + segment * 5, lo + segment * 6, lo + segment * 7, lo + segment * 8 };
            PivotInsertionSort(A, samples);
            return samples[4];
        }
        case PIVOT_MEDIAN9RANDOM:
        {
            if (lo > hi)
            {
                ssize_t temp = lo;
                lo = hi;
                hi = temp;
            }
            std::uniform_int_distribution<ssize_t> dist(lo, hi - 1);
            if (hi - lo < 81)
            {
                return SingleMedianOfThree(A, dist(gen1), dist(gen1), dist(gen1) + 1);
            }
            std::array<ssize_t, 9> samples = { dist(gen1), dist(gen1), dist(gen1), dist(gen1), dist(gen1), dist(gen1), dist(gen1), dist(gen1), dist(gen1) };
            PivotInsertionSort(A, samples);
            return samples[4];
        }
        case PIVOT_MEDIAN11:
        {
            if (hi - lo < 121) { return SingleMedianOfThree(A, lo, (lo + hi) / 2, hi); }
            ssize_t segment = (hi - lo) / 11;
            std::array<ssize_t, 11> samples = {lo, lo + segment, lo + segment * 2, lo + segment * 3, lo + segment * 4, lo + segment * 5, lo + segment * 6, lo + segment * 7, lo + segment * 8, lo + segment * 9, lo + segment * 10 };
            PivotShellSort(A, samples);
            return samples[5];
            break;
        }
        case PIVOT_MEDIAN11RANDOM:
        {
            if (lo > hi)
            {
                ssize_t temp = lo;
                lo = hi;
                hi = temp;
            }
            std::uniform_int_distribution<ssize_t> dist(lo, hi - 1);
            if (hi - lo < 121)
            {
                return SingleMedianOfThree(A, dist(gen1), dist(gen1), dist(gen1) + 1);
            }
            std::array<ssize_t, 11> samples = { dist(gen1), dist(gen1), dist(gen1), dist(gen1), dist(gen1), dist(gen1), dist(gen1), dist(gen1), dist(gen1), dist(gen1), dist(gen1) };
            PivotShellSort(A, samples);
            return samples[5];
            break;
        }
        case PIVOT_MEDIAN15:
        {
            if (hi - lo < 75) { return MedianOfFive(A, lo, hi); }
            ssize_t segment = (hi - lo) / 3;
            ssize_t g1 = MedianOfFive(A, lo, lo + segment);
            ssize_t g2 = MedianOfFive(A, lo + segment, lo + segment * 2);
            ssize_t g3 = MedianOfFive(A, lo + segment * 2, hi);
            return SingleMedianOfThree(A, g1, g2, g3 + 1);
        }
        case PIVOT_MEDIAN21:
        {
            if (hi - lo < 143) { return MedianOfSeven(A, lo, hi); }
            ssize_t segment = (hi - lo) / 3;
            ssize_t g1 = MedianOfSeven(A, lo, lo + segment);
            ssize_t g2 = MedianOfSeven(A, lo + segment, lo + segment * 2);
            ssize_t g3 = MedianOfSeven(A, lo + segment * 2, hi);
            return SingleMedianOfThree(A, g1, g2, g3 + 1);
        }
        case PIVOT_THREENINTHER:
        {
            if (hi - lo < 243) { return NintherPivot(A, lo, hi); }
            ssize_t segment = (hi - lo) / 3;
            ssize_t g1 = NintherPivot(A, lo, lo + segment);
            ssize_t g2 = NintherPivot(A, lo + segment, lo + segment * 2);
            ssize_t g3 = NintherPivot(A, lo + segment * 2, hi);
            return SingleMedianOfThree(A, g1, g2, g3 + 1);
        }
        default:
        {
            return lo;
        }
    }        
}

wxArrayString QuickSortPivotText()
{
    wxArrayString sl;

    sl.Add(_("First Item"));
    sl.Add(_("Last Item"));
    sl.Add(_("Middle Item"));
    sl.Add(_("Random Item"));
    sl.Add(_("Median of Three"));
    sl.Add(_("Random Median of Three"));
    sl.Add(_("Median of Five"));
    sl.Add(_("Random Median of Five"));
    sl.Add(_("Median of Seven"));
    sl.Add(_("Random Median of Seven"));
    sl.Add(_("Ninther"));
    sl.Add(_("Random Ninther"));
    sl.Add(_("Median of Nine"));
    sl.Add(_("Random Median of Nine"));
    sl.Add(_("Median of Eleven"));
    sl.Add(_("Random Median of Eleven"));
    sl.Add(_("Median of Fifteen"));
    sl.Add(_("Median of Twenty-One"));
    sl.Add(_("Median of Three Ninther"));
    return sl;
}

// ****************************************************************************
// *** Quick Sort LR (in-place, pointers at left and right, pivot is middle element)

// by myself (Timo Bingmann), based on Hoare's original code

void QuickSortLR(SortArray& A, ssize_t lo, ssize_t hi)
{
    // pick pivot and watch
    std::atomic<ssize_t> p{ QuickSortSelectPivot(A, lo, hi + 1) };
    value_type pivot = A[p.load()];
    A.watch(&p, 2);

    std::atomic<ssize_t> i{ lo }, j{ hi };
    A.watch(&i, 3);
    A.watch(&j, 3);

    while (i.load() <= j.load())
    {
        while (A[i.load()] < pivot) { i++; }
            
        while (A[j.load()] > pivot) { j--; }

        if (i.load() <= j.load())
        {
            A.swap(i.load(), j.load());

            // follow pivot if it is swapped
            if (p.load() == i.load()) p.store(j.load());
            else if (p.load() == j.load()) p.store(i.load());

            i++, j--;
        }
    }

    A.unwatch_all();

    if (lo < j)
        QuickSortLR(A, lo, j.load());

    if (i < hi)
        QuickSortLR(A, i.load(), hi);
}

void QuickSortLR(SortArray& A)
{
    return QuickSortLR(A, 0, A.size()-1);
}

// ****************************************************************************
// *** Quick Sort LL (in-place, two pointers at left, pivot is first element and moved to right)

// by myself (Timo Bingmann), based on CLRS' 3rd edition

size_t PartitionLL(SortArray& A, size_t lo, size_t hi)
{
    // pick pivot and move to back
    size_t p = QuickSortSelectPivot(A, lo, hi);

    value_type pivot = A[p];
    A.swap(p, hi - 1);
    A.mark(hi - 1);

    std::atomic<ssize_t> i{static_cast<ssize_t>(lo)};
    A.watch(&i, 3);

    for (size_t j = lo; j < hi - 1; ++j)
    {
        if (A[j] <= pivot) {
            A.swap(i.load(), j);
            ++i;
        }
    }

    A.swap(i.load(), hi - 1);
    A.unmark(hi - 1);
    A.unwatch_all();

    return i.load();
}

void QuickSortLL(SortArray& A, size_t lo, size_t hi)
{
    if (lo + 1 < hi)
    {
        size_t mid = PartitionLL(A, lo, hi);

        QuickSortLL(A, lo, mid);
        QuickSortLL(A, mid+1, hi);
    }
}

void QuickSortLL(SortArray& A)
{
    return QuickSortLL(A, 0, A.size());
}

// ****************************************************************************
// *** Quick Sort Ternary (in-place, two pointers at left, pivot is first element and moved to right)

// by myself (Timo Bingmann), loosely based on multikey quicksort by B&S

void QuickSortTernaryLR(SortArray& A, ssize_t lo, ssize_t hi)
{
    if (hi <= lo) return;

    int cmp;

    // pick pivot and swap to back
    ssize_t piv = QuickSortSelectPivot(A, lo, hi + 1);
    A.swap(piv, hi);
    A.mark(hi);
    
    const value_type& pivot = A[hi];

    // schema: |p ===  |i <<< | ??? |j >>> |q === |piv
    std::atomic<ssize_t> i{ lo }, j{hi - 1};
    ssize_t p = lo, q = hi - 1;
    
    A.watch(&i, 3);
    A.watch(&j, 3);

    for (;;)
    {
        // partition on left
        while (i.load() <= j.load() && (cmp = A[i.load()].cmp(pivot)) <= 0)
        {
            if (cmp == 0) {
                A.mark(p, 4);
                A.swap(i.load(), p++);
            }
            ++i;
        }

        // partition on right
        while (i.load() <= j.load() && (cmp = A[j.load()].cmp(pivot)) >= 0)
        {
            if (cmp == 0) {
                A.mark(q, 4);
                A.swap(j.load(), q--);
            }
            --j;
        }

        if (i.load() > j.load()) break;

        // swap item between < > regions
        A.swap(i++, j--);
    }

    // swap pivot to right place
    A.swap(i.load(), hi);
    A.mark_swap(i.load(), hi);

    ssize_t num_less = i.load() - p;
    ssize_t num_greater = q - j.load();

    // swap equal ranges into center, but avoid swapping equal elements
    j.store(i.load() - 1); i.store(i.load() + 1);

    ssize_t pe = lo + std::min(p - lo, num_less);
    for (ssize_t k = lo; k < pe; k++, j--) {
        A.swap(k, j.load());
        A.mark_swap(k, j.load());
    }

    ssize_t qe = hi - 1 - std::min(hi - 1 - q, num_greater - 1); // one already greater at end
    for (ssize_t k = hi - 1; k > qe; k--, i++) {
        A.swap(i.load(), k);
        A.mark_swap(i.load(), k);
    }

    A.unwatch_all();
    A.unmark_all();

    QuickSortTernaryLR(A, lo, lo + num_less - 1);
    QuickSortTernaryLR(A, hi - num_greater + 1, hi);
}

void QuickSortTernaryLR(SortArray& A)
{
    return QuickSortTernaryLR(A, 0, A.size()-1);
}

// ****************************************************************************
// *** Quick Sort LL (in-place, two pointers at left, pivot is first element and moved to right)

// by myself (Timo Bingmann)

std::pair<ssize_t, ssize_t> PartitionTernaryLL(SortArray& A, ssize_t lo, ssize_t hi)
{
    // pick pivot and swap to back
    ssize_t p = QuickSortSelectPivot(A, lo, hi);

    value_type pivot = A[p];
    A.swap(p, hi - 1);
    A.mark(hi - 1);

    std::atomic<ssize_t> i{lo};
    ssize_t k = hi - 1;
    
    A.watch(&i, 3);

    for (ssize_t j = lo; j < k; ++j)
    {
        int cmp = A[j].cmp(pivot); // ternary comparison
        if (cmp == 0) {
            A.swap(--k, j);
            --j; // reclassify A[j]
            A.mark(k, 4);
        }
        else if (cmp < 0) {
            A.swap(i++, j);
        }
    }

    // unwatch i, because the pivot is swapped there
    // in the first step of the following swap loop.
    A.unwatch_all();

    ssize_t j = i.load() + (hi - k);

    for (ssize_t s = 0; s < hi - k; ++s) {
        A.swap(i.load() + s, hi - 1 - s);
        A.mark_swap(i.load() + s, hi - 1 - s);
    }
    A.unmark_all();

    return std::make_pair(i.load(), j);
}

void QuickSortTernaryLL(SortArray& A, size_t lo, size_t hi)
{
    if (lo + 1 < hi)
    {
        std::pair<ssize_t,ssize_t> mid = PartitionTernaryLL(A, lo, hi);

        QuickSortTernaryLL(A, lo, mid.first);
        QuickSortTernaryLL(A, mid.second, hi);
    }
}

void QuickSortTernaryLL(SortArray& A)
{
    return QuickSortTernaryLL(A, 0, A.size());
}

// ****************************************************************************
// *** Dual-Pivot Quick Sort

// by Sebastian Wild

void dualPivotYaroslavskiy(class SortArray& a, int left, int right)
{
    if (right > left)
    {
        if (a[left] > a[right]) {
            a.swap(left, right);
        }

        const value_type p = a[left];
        const value_type q = a[right];

        a.mark(left);
        a.mark(right);

        std::atomic<ssize_t> l{ left + 1 }, g{ right - 1 }, k{ l.load() };
       
        a.watch(&l, 3);
        a.watch(&g, 3);
        a.watch(&k, 3);

        while (k.load() <= g)
        {
            if (a[k.load()] < p) {
                a.swap(k.load(), l.load());
                ++l;
            }
            else if (a[k.load()] >= q) {
                while (a[g.load()] > q && k.load() < g) { --g; }
                a.swap(k.load(), g.load());
                --g;           
                if (a[k.load()] < p) {
                    a.swap(k.load(), l.load());
                    ++l;
                }
            }
            ++k;
        }
        --l;
        ++g;

        a.swap(left, l.load());
        a.swap(right, g.load());

        a.unmark_all();
        a.unwatch_all();

        dualPivotYaroslavskiy(a, left, l.load() - 1);
        dualPivotYaroslavskiy(a, l.load() + 1, g.load() - 1);
        dualPivotYaroslavskiy(a, g.load() + 1, right);
    }
}

void QuickSortDualPivot(class SortArray& a)
{
    return dualPivotYaroslavskiy(a, 0, a.size()-1);
}

// ****************************************************************************
// *** Bubble Sort

void BubbleSort(SortArray& A)
{
    for (size_t i = 0; i < A.size()-1; ++i)
    {
        for (size_t j = 0; j < A.size()-1 - i; ++j)
        {
            if (A[j] > A[j + 1])
            {
                A.swap(j, j+1);
            }
        }
    }
}

void OptimizedBubbleSort(SortArray& A)
{
    for (size_t i = 0; i < A.size() - 1; ++i)
    {
        bool sorted = true;
        for (size_t j = 0; j < A.size() - 1 - i; ++j)
        {
            if (A[j] > A[j + 1])
            {
                A.swap(j, j + 1);
                sorted = false;
            }
        }
        if (sorted == true) { break; }
    }
}

void TargetedBubbleSort(SortArray& A)
{
    bool sorted = false;
    size_t target = A.size() - 1, lastSwapped = 0;
    while (!sorted)
    {
        sorted = true;
        for (size_t i = 0; i < target; ++i)
        {
            if (A[i] > A[i + 1])
            {
                A.swap(i, i + 1);
                lastSwapped = i;
                sorted = false;
            }
        }
        target = lastSwapped;
    }
}

// ****************************************************************************
// *** Cocktail Shaker Sort

// from http://de.wikibooks.org/wiki/Algorithmen_und_Datenstrukturen_in_C/_Shakersort

void CocktailShakerSort(SortArray& A)
{
    size_t lo = 0, hi = A.size()-1, mov = lo;

    while (lo < hi)
    {
        for (size_t i = hi; i > lo; --i)
        {
            if (A[i-1] > A[i])
            {
                A.swap(i-1, i);
                mov = i;
            }
        }

        lo = mov;

        for (size_t i = lo; i < hi; ++i)
        {
            if (A[i] > A[i+1])
            {
                A.swap(i, i+1);
                mov = i;
            }
        }

        hi = mov;
    }
}

void DualCocktailShakerSort(SortArray& A)
{
    size_t lo = 0, hi = A.size() - 1;
    while (lo < hi)
    {
        size_t lo_mov = 0, hi_mov = 0;
        for (size_t i = lo + 1, j = hi - 1; i <= hi; ++i, --j)
        {
            if (A[i - 1] > A[i])
            {
                A.swap(i - 1, i);
                lo_mov = i;
            }
            if (A[j + 1] < A[j])
            {
                A.swap(j + 1, j);
                hi_mov = j;
            }
        }
        lo = hi_mov;
        hi = lo_mov;
    }
}

// ****************************************************************************
// *** Circle Sort

bool CircleSortRec(SortArray& A, size_t low, size_t high, size_t len)
{
    bool swapped = false;
    if (low == high) { return false; }
    size_t lo = low, hi = high;
    while (lo < hi)
    {
        if (hi < len && A[lo] > A[hi])
        {
            A.swap(lo, hi);
            swapped = true;
        }
        ++lo; --hi;
    }
    size_t mid = (high - low) / 2;
    bool firstHalf = CircleSortRec(A, low, low + mid, len);
    bool secondHalf = false;
    if (low + mid + 1 < len)
    {
        secondHalf = CircleSortRec(A, low + mid + 1, high, len);
    }
    return swapped || firstHalf || secondHalf;
}

void CircleSort(SortArray& A)
{
    size_t len = A.size(), n = 1;
    for (; n < len; n *= 2) {}
    while (CircleSortRec(A, 0, n - 1, len)) {}
}

bool CircleSortIte(SortArray& A, size_t length, size_t arr_len)
{
    bool swapped = false;
    for (size_t gap = length / 2; gap > 0; gap /= 2)
    {
        for (size_t start = 0; start + gap < arr_len; start += 2 * gap)
        {
            size_t high = start + 2 * gap - 1, low = start;
            while (low < high)
            {
                if (high < arr_len && A[low] > A[high])
                {
                    A.swap(low, high);
                    swapped = true;
                }
                ++low; --high;
            }
        }
    }
    return swapped;
}

void CircleSort2(SortArray& A)
{
    size_t len = A.size(), n = 1;
    for (; n < len; n *= 2) {}
    while (CircleSortIte(A, n, len)) {}
}

void IntroCircleSort(SortArray& A)
{
    size_t len = A.size(), threshold = 0, n = 1, iterations = 0;
    for (; n < len; n *= 2, ++threshold) {}
    threshold /= 2;
    do
    {
        ++iterations;
        if (iterations >= threshold)
        {
            InsertSort(A, 0, len);
            break;
        }
    } 
    while (CircleSortRec(A, 0, n - 1, len));
}

void IntroIteCircleSort(SortArray& A)
{
    size_t len = A.size(), threshold = 0, n = 1, iterations = 0;
    for (; n < len; n *= 2, ++threshold) {}
    threshold /= 2;
    do
    {
        ++iterations;
        if (iterations >= threshold)
        {
            InsertSort(A, 0, len);
            break;
        }
    } while (CircleSortIte(A, n, len));
}

// ****************************************************************************
// *** Gnome Sort

// from http://en.wikipediA.org/wiki/Gnome_sort

void GnomeSort(SortArray& A)
{
    for (size_t i = 1; i < A.size(); )
    {
        if (A[i] >= A[i-1])
        {
            ++i;
        }
        else
        {
            A.swap(i, i-1);
            if (i > 1) --i;
        }
    }
}

void OptimizedGnomeSort(SortArray& A)
{
    size_t prev = 0;
    for (size_t i = 1; i < A.size(); )
    {
        if (i == 0 || A[i] >= A[i - 1])
        {
            if (prev != 0) { i += prev; prev = 0; }
            ++i;
        }
        else
        {
            A.swap(i, i - 1);
            --i; ++prev;
        }
    }
}

// ****************************************************************************
// *** Comb Sort

// from http://en.wikipediA.org/wiki/Comb_sort

void CombSort(SortArray& A)
{
    const double shrink = 1.3;

    bool swapped = false;
    size_t gap = A.size();

    while ((gap > 1) || swapped)
    {
        if (gap > 1) {
            gap = (size_t)((float)gap / shrink);
        }

        swapped = false;

        for (size_t i = 0; gap + i < A.size(); ++i)
        {
            if (A[i] > A[i + gap])
            {
                A.swap(i, i+gap);
                swapped = true;
            }
        }
    }
}

// ****************************************************************************
// *** Odd-Even Sort

// from http://en.wikipediA.org/wiki/Odd%E2%80%93even_sort

void OddEvenSort(SortArray& A)
{
    bool sorted = false;

    while (!sorted)
    {
        sorted = true;

        for (size_t i = 1; i < A.size()-1; i += 2)
        {
            if(A[i] > A[i+1])
            {
                A.swap(i, i+1);
                sorted = false;
            }
        }

        for (size_t i = 0; i < A.size()-1; i += 2)
        {
            if(A[i] > A[i+1])
            {
                A.swap(i, i+1);
                sorted = false;
            }
        }
    }
}

// ****************************************************************************
// *** Shell Sort

// with gaps by Robert Sedgewick from http://www.cs.princeton.edu/~rs/shell/shell.c

void ShellSort(SortArray& A)
{
    size_t incs[16] = { 1391376, 463792, 198768, 86961, 33936,
                        13776, 4592, 1968, 861, 336,
                        112, 48, 21, 7, 3, 1 };

    for (size_t k = 0; k < 16; k++)
    {
        for (size_t h = incs[k], i = h; i < A.size(); i++)
        {
            value_type v = A[i];
            size_t j = i;

            while (j >= h && A[j-h] > v)
            {
                A.set(j, A[j-h]);
                j -= h;
            }

            A.set(j, v);
        }
    }
}

// ****************************************************************************
// *** Heap Sort

// heavily adapted from http://www.codecodex.com/wiki/Heapsort

bool isPowerOfTwo(size_t x)
{
    return ((x != 0) && !(x & (x - 1)));
}

uint32_t prevPowerOfTwo(uint32_t x)
{
    x |= x >> 1; x |= x >> 2; x |= x >> 4;
    x |= x >> 8; x |= x >> 16;
    return x - (x >> 1);
}

int largestPowerOfTwoLessThan(int n)
{
    int k = 1;
    while (k < n) k = k << 1;
    return k >> 1;
}

void HeapSort(SortArray& A)
{
    size_t n = A.size(), i = n / 2;

    // mark heap levels with different colors
    for (size_t j = i; j < n; ++j)
        A.mark(j, log(prevPowerOfTwo(j+1)) / log(2) + 4);

    while (1)
    {
        if (i > 0) {
            // build heap, sift A[i] down the heap
            i--;
        }
        else {
            // pop largest element from heap: swap front to back, and sift
            // front A[0] down the heap
            n--;
            if (n == 0) return;
            A.swap(0,n);

            A.mark(n);
            if (n+1 < A.size()) A.unmark(n+1);
        }

        size_t parent = i;
        size_t child = i*2 + 1;

        // sift operation - push the value of A[i] down the heap
        while (child < n)
        {
            if (child + 1 < n && A[child + 1] > A[child]) {
                child++;
            }
            if (A[child] > A[parent]) {
                A.swap(parent, child);
                parent = child;
                child = parent*2+1;
            }
            else {
                break;
            }
        }

        // mark heap levels with different colors
        A.mark(i, log(prevPowerOfTwo(i+1)) / log(2) + 4);
    }

}

// ****************************************************************************
// *** Radix Sort (counting sort, most significant digit (MSD) first, in-place redistribute)

// by myself (Timo Bingmann)

void RadixSortMSD2(SortArray& A, size_t lo, size_t hi, size_t depth)
{
    A.mark(lo); A.mark(hi-1);

    // radix and base calculations
    const unsigned int RADIX = 4;

    unsigned int pmax = floor( log(A.array_max()+1) / log(RADIX) );
    ASSERT(depth <= pmax);

    size_t base = pow(RADIX, pmax - depth);

    // count digits
    std::vector<size_t> count(RADIX, 0);

    for (size_t i = lo; i < hi; ++i)
    {
        size_t r = A[i].get() / base % RADIX;
        ASSERT(r < RADIX);
        count[r]++;
    }

    // inclusive prefix sum
    std::vector<size_t> bkt(RADIX, 0);
    std::partial_sum(count.begin(), count.end(), bkt.begin());

    // mark bucket boundaries
    for (size_t i = 0; i < bkt.size(); ++i) {
        if (bkt[i] == 0) continue;
        A.mark(lo + bkt[i]-1, 3);
    }

    // reorder items in-place by walking cycles
    for (size_t i=0, j; i < (hi-lo); )
    {
        while ( (j = --bkt[ (A[lo+i].get() / base % RADIX) ]) > i )
        {
            A.swap(lo + i, lo + j);
        }
        i += count[ (A[lo+i].get() / base % RADIX) ];
    }

    A.unmark_all();

    // no more depth to sort?
    if (depth+1 > pmax) return;

    // recurse on buckets
    size_t sum = lo;
    for (size_t i = 0; i < RADIX; ++i)
    {
        if (count[i] > 1)
            RadixSortMSD2(A, sum, sum+count[i], depth+1);
        sum += count[i];
    }
}

void RadixSortMSD2(SortArray& A)
{
    return RadixSortMSD2(A, 0, A.size(), 0);
}

/*
    MIT License

    Copyright (c) 2019 w0rthy

    Permission is hereby granted, free of charge, to any person obtaining a copy
    of this software and associated documentation files (the "Software"), to deal
    in the Software without restriction, including without limitation the rights
    to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
    copies of the Software, and to permit persons to whom the Software is
    furnished to do so, subject to the following conditions:

    The above copyright notice and this permission notice shall be included in all
    copies or substantial portions of the Software.

    THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
    IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
    FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
    AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
    LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
    OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
    SOFTWARE.
*/

size_t maxLog(SortArray& A, size_t n, size_t base)
{
    int max = A[0];
    for (size_t i = 1, j = n - 1; i <= j; ++i, --j)
    {
        int ele1 = A[i].get(), ele2 = A[j];
        if (ele1 > max) { max = ele1; }
        if (ele2 > max) { max = ele2; }
    }
    size_t digit = static_cast<size_t>(log(max) / log(base));
    return digit;
}

int getDigit2(int a, double power, int radix)
{
    double digit = (a / static_cast<int>(pow(radix, power)) % radix);
    return static_cast<int>(digit);
}

void transcribeMSD(SortArray& A, std::vector<std::vector<value_type>>& registers, size_t start, size_t min)
{
    size_t total = start, temp = 0;
    for (const std::vector<value_type>& arr : registers)
    {
        total += arr.size();
    }

    for (int i = registers.size() - 1; i >= 0; --i)
    {
        for (int j = registers[i].size() - 1; j >= 0; --j)
        {
            size_t loc = total + min - temp - 1;
            A.set(loc, registers[i].at(j));
            A[loc].get();
            ++temp;
        }
    }
}

void radixMSD(SortArray& A, size_t len, size_t min, size_t max, size_t radix, double pow)
{
    if (min >= max || pow < 0) { return; }
    A.mark(min); A.mark(max - 1);

    std::vector<std::vector<value_type>> registers(radix, std::vector<value_type>());

    for (size_t i = min; i < max; ++i)
    {
        int ele = A[i].get();
        int digit = getDigit2(ele, pow, radix);
        registers[digit].push_back(A[i]);
    }

    transcribeMSD(A, registers, 0, min);
    A.unmark_all();

    size_t sum = 0;
    for (size_t i = 0; i < registers.size(); ++i)
    {
        radixMSD(A, len, sum + min, sum + min + registers[i].size(), radix, pow - 1);
        sum += registers[i].size();
        registers[i].clear();
    }
}

void RadixSortMSD(SortArray& A)
{
    size_t n = A.size(), maxPow = maxLog(A, n, 4);
    radixMSD(A, n, 0, n, 4, (double)maxPow);
}

// ****************************************************************************
// *** Radix Sort (counting sort, least significant digit (LSD) first, out-of-place redistribute)

// by myself (Timo Bingmann)

void RadixSortLSD(SortArray& A)
{
    // radix and base calculations
    const unsigned int RADIX = 4;

    unsigned int pmax = ceil( log(A.array_max()+1) / log(RADIX) );

    for (unsigned int p = 0; p < pmax; ++p)
    {
        size_t base = pow(RADIX, p);

        // count digits and copy data
        std::vector<size_t> count(RADIX, 0);
        std::vector<value_type> copy(A.size());

        for (size_t i = 0; i < A.size(); ++i)
        {
            size_t r = (copy[i] = A[i]).get() / base % RADIX;
            ASSERT(r < RADIX);
            count[r]++;
        }

        // exclusive prefix sum
        std::vector<size_t> bkt(RADIX+1, 0);
        std::partial_sum(count.begin(), count.end(), bkt.begin()+1);

        // mark bucket boundaries
        for (size_t i = 0; i < bkt.size()-1; ++i) {
            if (bkt[i] >= A.size()) continue;
            A.mark(bkt[i], 3);
        }

        // redistribute items back into array (stable)
        for (size_t i=0; i < A.size(); ++i)
        {
            size_t r = copy[i].get() / base % RADIX;
            A.set( bkt[r]++, copy[i] );
        }

        A.unmark_all();
    }
}

// ****************************************************************************
// *** In-Place Radix Sort LSD
/*
    Copyright (c) 2019 w0rthy
    Permission is hereby granted, free of charge, to any person obtaining a copy
    of this software and associated documentation files (the "Software"), to deal
    in the Software without restriction, including without limitation the rights
    to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
    copies of the Software, and to permit persons to whom the Software is
    furnished to do so, subject to the following conditions:
    The above copyright notice and this permission notice shall be included in all
    copies or substantial portions of the Software.
    THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
    IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
    FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
    AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
    LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
    OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
    SOFTWARE.
*/

void shiftElement(SortArray& A, size_t start, size_t end)
{
    if (start < end)
    {
        while (start < end)
        {
            A[start].get();
            A.swap(start, start + 1);
            ++start;
        }
    }
    else
    {
        while (start > end)
        {
            A[start].get();
            A.swap(start, start - 1);
            --start;
        }
    }
}

void InPlaceRadixSortLSD(SortArray& A)
{
    size_t pos = 0, n = A.size();
    int bucket = 20;
    std::vector<int> buckets(bucket - 1, 0);
    size_t maxPow = maxLog(A, n, bucket);
    for (size_t p = 0; p <= maxPow; ++p)
    {
        for (size_t i = 0; i < buckets.size(); ++i) { buckets[i] = n - 1; }
        pos = 0;
        for (size_t i = 0; i < n; ++i)
        {
            int ele = A[pos].get();
            int digit = getDigit2(ele, (double)p, bucket);
            if (digit == 0) { ++pos; }
            else
            {
                shiftElement(A, pos, buckets[digit - 1]);
                for (size_t j = digit - 1; j > 0; --j)
                {
                    buckets[j - 1] = buckets[j - 1] - 1;
                }
            }
        }
    }
}

// ****************************************************************************
// *** Rotate Radix MSD/LSD Sort

/*
    Copyright (c) 2020-2021 aphitorite

    Permission is hereby granted, free of charge, to any person obtaining a copy
    of this software and associated documentation files (the "Software"), to deal
    in the Software without restriction, including without limitation the rights
    to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
    copies of the Software, and to permit persons to whom the Software is
    furnished to do so, subject to the following conditions:

    The above copyright notice and this permission notice shall be included in all
    copies or substantial portions of the Software.

    THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
    IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
    FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
    AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
    LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
    OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
    SOFTWARE.
*/

static const int base = 4;

int shift(int n, int q)
{
    while (q > 0)
    {
        n /= base;
        --q;
    }
    return n;
}

int binSearch(SortArray& A, int a, int b, int d, int p)
{
    while (a < b)
    {
        int m = (a + b) / 2;
        int ele = A[m].get();
        int result = static_cast<int>(getDigit2(static_cast<size_t>(ele), static_cast<size_t>(p), static_cast<size_t>(base)));
        if (result >= d) { b = m; }
        else { a = m + 1; }
    }
    return a;
}

void rotateMerge(SortArray& A, int a, int m, int b, int da, int db, int p)
{
    if (b - a < 2 || db - da < 2) { return; }
    int dm = (da + db) / 2;
    int m1 = binSearch(A, a, m, dm, p);
    int m2 = binSearch(A, m, b, dm, p);
    rotate(A, m1, m, m2);
    m = m1 + (m2 - m);
    rotateMerge(A, m, m2, b, dm, db, p);
    rotateMerge(A, a, m1, m, da, dm, p);
}

void rotateMergeSort(SortArray& A, int a, int b, int p)
{
    if (b - a < 2) { return; }
    int m = (a + b) / 2;
    rotateMergeSort(A, a, m, p);
    rotateMergeSort(A, m, b, p);
    rotateMerge(A, a, m, b, 0, base, p);
}

int dist(SortArray& A, int a, int b, int p)
{
    rotateMergeSort(A, a, b, p);
    return binSearch(A, a, b, 1, p);
}

void RotateRadixSortMSD(SortArray& A)
{
    int len = static_cast<int>(A.size());
    int q = static_cast<int>(maxLog(A, static_cast<size_t>(len), static_cast<size_t>(base)));
    int m = 0, i = 0, b = len;
    while (i < len)
    {
        int p = 0;
        if (b - i < 1) { p = i; }
        else { p = dist(A, i, b, q); }
        if (q == 0)
        {
            m += base;
            int t = m / base;
            while (t % base == 0)
            {
                t /= base;
                ++q;
            }
            i = b;
            while (b < len)
            {
                int ele = A[b].get();
                if (shift(ele, q + 1) == shift(m, q + 1)) { ++b; }
                else { break; }
            }
        }
        else
        {
            b = p;
            --q;
        }
    }
}

void RotateRadixSortLSD(SortArray& A)
{
    int len = static_cast<int>(A.size());
    int max = static_cast<int>(maxLog(A, static_cast<size_t>(len), static_cast<size_t>(base)));
    for (int i = 0; i <= max; ++i)
    {
        rotateMergeSort(A, 0, len, i);
    }
}

// ****************************************************************************
// *** Use STL Sorts via Iterator Adapters

void siftDown(SortArray& A, size_t root, size_t dist, size_t start, bool isMax)
{
    int compareVal = 0;
    if (isMax) { compareVal = -1; }
    else { compareVal = 1; }
    while (root <= dist / 2)
    {
        size_t leaf = 2 * root;
        int compVal = 0;

        if (leaf < dist) 
        { 
            if (A[start + leaf - 1] < A[start + leaf]) { compVal = -1; }
            else if (A[start + leaf - 1] > A[start + leaf]) { compVal = 1; }
            if (compVal == compareVal) { ++leaf; }
        }
        
        if (A[start + root - 1] < A[start + leaf -  1]) { compVal = -1; }
        else if (A[start + root - 1] > A[start + leaf - 1]) { compVal = 1; }
        else { compVal = 0; }

        if (compVal == compareVal)
        {
            A.swap(start + root - 1, start + leaf - 1);
            root = leaf;
        }
        else { break; }
    }
}

void heapifyArr(SortArray& A, size_t low, size_t high, bool isMax)
{
    size_t len = high - low;
    for (size_t i = len / 2; i >= 1; --i)
    {
        siftDown(A, i, len, low, isMax);
    }
}

void reverseArr(SortArray& A, size_t start, size_t len)
{
    for (size_t i = start; i < start + ((len - start + 1) / 2); ++i)
    {
        A.swap(i, start + len - 1);
    }
}

void HeapSort2(SortArray& A, size_t start, size_t len, bool isMax)
{
    heapifyArr(A, start, len, isMax);
    for (size_t i = len - start; i > 1; --i)
    {
        A.swap(start, start + i - 1);
        siftDown(A, 1, i - 1, start, isMax);
    }

    if (!isMax)
    {
        reverseArr(A, start, start + len - 1);
    }
}

size_t floorLogBaseTwo(size_t a)
{
    return static_cast<size_t>(floor(log(a) / log(2)));
}

value_type gccmedianof3(SortArray& A, size_t left, size_t mid, size_t right)
{
    if (A[left] < A[mid])
    {
        if (A[mid] < A[right])
        {
            A.swap(left, mid);
        }
        else if (A[left] < A[right])
        {
            A.swap(left, right);
        }
    }
    else if (A[left] < A[right]) { return A[left]; }
    else if (A[mid] < A[right])
    {
        A.swap(left, right);
    }
    else
    {
        A.swap(mid, right);
    }
    return A[left];
}

value_type medianof3(SortArray& A, size_t left, size_t mid, size_t right)
{
    if (A[right] < A[left]) { A.swap(left, right); }
    if (A[mid] < A[left]) { A.swap(mid, left); }
    if (A[right] < A[mid]) { A.swap(right, mid); }
    return A[mid];
}

size_t partitionArr(SortArray& A, size_t lo, size_t hi, value_type x)
{
    size_t i = lo, j = hi;
    while (true)
    {
        while (A[i] < x) { ++i; }
        --j;
        while (x < A[j]) { --j; }

        if (!(i < j)) { return i; }
        A.swap(i, j);
        ++i;
    }
}

void introsortLoop(SortArray& A, size_t lo, size_t hi, size_t depth)
{
    size_t threshold = 16;
    while (hi - lo > threshold)
    {
        if (depth == 0)
        {
            HeapSort2(A, lo, hi, true);
            return;
        }
        --depth;
        size_t p = partitionArr(A, lo, hi, medianof3(A, lo, lo + ((hi - lo) / 2), hi - 1));
        introsortLoop(A, p, hi, depth);
        hi = p;
    }
    return;
}

void StlSort(SortArray& A)
{
    size_t n = A.size();
    introsortLoop(A, 0, n, 2 * floorLogBaseTwo(n));
    InsertionSort(A);
}

void StlSort2(SortArray& A)
{
    std::sort(MyIterator(&A, 0), MyIterator(&A, A.size()));
}

void StlStableSort(SortArray& A)
{
    std::stable_sort(MyIterator(&A, 0), MyIterator(&A, A.size()));
}

void StlHeapSort2(SortArray& A)
{
    std::make_heap(MyIterator(&A, 0), MyIterator(&A, A.size()));
    std::sort_heap(MyIterator(&A, 0), MyIterator(&A, A.size()));
}

void StlHeapSort(SortArray& A)
{
    HeapSort2(A, 0, A.size(), true);
}

// ****************************************************************************
// *** BogoSort and more slow sorts

// by myself (Timo Bingmann)

bool BogoCheckSorted(SortArray& A)
{
    size_t i;
    A.mark(0);
    for (i = 1; i < A.size(); ++i)
    {
        if (A[i - 1] > A[i]) break;
        A.mark(i);
    }

    if (i == A.size()) {
        // this is amazing.
        return true;
    }

    // unmark
    while (i > 0) A.unmark(i--);
    A.unmark(0);

    return false;
}

void BogoSort(SortArray& A)
{
    // keep a permutation of [0,size)
    std::vector<size_t> perm(A.size());

    for (size_t i = 0; i < A.size(); ++i)
        perm[i] = i;

    while (1)
    {
        // check if array is sorted
        if (BogoCheckSorted(A)) break;

        // pick a random permutation of indexes
        std::shuffle(perm.begin(), perm.end(), gen1);

        // permute array in-place
        std::vector<char> pmark(A.size(), 0);

        for (size_t i = 0; i < A.size(); ++i)
        {
            if (pmark[i]) continue;

            // walk a cycle
            size_t j = i;

            //std::cout << "cycle start " << j << " -> " << perm[j] << "\n";

            while ( perm[j] != i )
            {
                ASSERT(!pmark[j]);
                A.swap(j, perm[j]);
                pmark[j] = 1;

                j = perm[j];
                //std::cout << "cycle step " << j << " -> " << perm[j] << "\n";
            }
            //std::cout << "cycle end\n";

            ASSERT(!pmark[j]);
            pmark[j] = 1;
        }

        //std::cout << "permute end\n";

        for (size_t i = 0; i < A.size(); ++i)
            ASSERT(pmark[i]);
    }
}

void BozoSort(SortArray& A)
{
    srand(time(nullptr));

    while (1)
    {
        // check if array is sorted
        if (BogoCheckSorted(A)) break;

        // swap two random items
        A.swap(rand() % A.size(), rand() % A.size());
    }
}

void flip(SortArray& A, size_t high)
{
    size_t low = 0;
    while (low < high)
    {
        A[low].get();
        A.swap(low, high);
        ++low; --high;
    }
}

size_t find_max(SortArray& A, size_t n)  // Optimized find_max method, searching the max element in n/2 time
{
    size_t max = 0;
    for (size_t low = 1, hi = n; low <= hi; ++low, --hi)
    {
        if (A[low] > A[max])
        { 
            max = low; 
        }
        if (A[hi] > A[max])
        {
            max = hi;
        }
    }
    return max;
}

void PancakeSort(SortArray& A)
{
    size_t n = A.size() - 1;
    for (size_t cur_size = n; cur_size >= 1; --cur_size)
    {
        size_t max_idx = find_max(A, cur_size);
        if (max_idx != cur_size)
        {
            if (max_idx > 0) { flip(A, max_idx); }
            flip(A, cur_size);
        }
    }
}


// ****************************************************************************
// *** Optimized Pancake Sort
/*
    The MIT License (MIT)

    Copyright (c) 2021-2023 aphitorite

    Permission is hereby granted, free of charge, to any person obtaining a copy of
    this software and associated documentation files (the "Software"), to deal in
    the Software without restriction, including without limitation the rights to
    use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of
    the Software, and to permit persons to whom the Software is furnished to do so,
    subject to the following conditions:

    The above copyright notice and this permission notice shall be included in all
    copies or substantial portions of the Software.

    THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
    IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS
    FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR
    COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER
    IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
    CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*/

void flip2(SortArray& A, size_t high)
{
    if (high > 0) { --high; }
    size_t low = 0; 
    while (low < high)
    {
        A[low].get();
        A.swap(low, high);
        ++low; --high;
    }
}


bool mergeFlip(SortArray& A, size_t h1, size_t h2)
{
    if (h1 == 1 && h2 == 1)
    {
        if (A[0] > A[1]) { flip2(A, 2); }
        return true;
    }
    size_t n = h1 + h2, m = n / 2;
    if (h2 < h1)
    {
        if (h2 < 1) { return true; }
        size_t i = 0, j = h2;
        while (i < j)
        {
            size_t k = (i + j) / 2, loc = n - 1 - k;
            A[loc].get();
            if (A[n - 1 - k - m] > A[loc]) { i = k + 1; }
            else { j = k; }
        }
        flip2(A, n - m - i);
        flip2(A, n - i);
        if (mergeFlip(A, h2 - i, i + m - h2)) { flip2(A, m); }
        flip2(A, n);
        if (!mergeFlip(A, i, n - m - i)) { flip2(A, n - m); }
    }
    else
    {
        if (h1 < 1) { return false; }
        size_t i = 0, j = h1;
        while (i < j)
        {
            size_t k = (i + j) / 2;
            A[k].get();
            if (A[k] < A[k + m]) { i = k + 1; }
            else { j = k; }
        }
        flip2(A, i);
        flip2(A, i + m);
        if (mergeFlip(A, i + m - h1, h1 - i)) { flip2(A, m); }
        flip2(A, n);
        if (!mergeFlip(A, n - m - i, i)) { flip2(A, n - m); }
    }
    return true;
}

void sortFlip(SortArray& A, size_t n)
{
    if (n < 2) { return; }
    size_t h = n / 2;
    sortFlip(A, h);
    flip2(A, n);
    sortFlip(A, n - h);
    mergeFlip(A, n - h, h);
}

void OptimizedPancakeSort(SortArray& A)
{
    sortFlip(A, A.size());
}

void BeadSort(SortArray& A)
{  
    int max = A[0];
    int len = A.size();
    for (int n = 1; n < len; ++n)
    {
        int m = A[n].get();
        if (m > max)
        { max = m; }
    }

    std::vector<std::vector<int>> beads(len, std::vector<int>(max, 0));

    for (int i = 0; i < len; ++i) 
    {
        int n = A[i].get();
        for (int j = 0; j < n; ++j) 
        { beads[i][j] = 1; }
    }

    for (int j = 0; j < max; ++j) 
    {
        int sum = 0;
        for (int i = 0; i < len; ++i) 
        {
            sum += beads[i][j];
            beads[i][j] = 0;
        }
        for (int i = len - 1; i >= len - sum; --i) 
        {
            size_t k = static_cast<size_t>(i);
            A.set(k, ArrayItem(j + 1));
            A[k].get();
        }
    }
}

/*
    MIT License

    Copyright (c) 2020 aphitorite

    Permission is hereby granted, free of charge, to any person obtaining a copy
    of this software and associated documentation files (the "Software"), to deal
    in the Software without restriction, including without limitation the rights
    to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
    copies of the Software, and to permit persons to whom the Software is
    furnished to do so, subject to the following conditions:

    The above copyright notice and this permission notice shall be included in all
    copies or substantial portions of the Software.

    THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
    IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
    FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
    AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
    LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
    OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
    SOFTWARE.
*/

void GravitySort(SortArray& A)
{
    int n = A.size();
    int min = A[0], max = A[0];
    for (int i = 1; i < n; ++i)
    {
        int ele = A[i].get();
        if (ele < min) { min = ele; }
        if (ele > max) { max = ele; }
    }
    std::vector<int> x(n, 0);
    std::vector<int> y(max - min + 1, 0);
    for (int i = 0; i < n; ++i)
    {
        int ele = A[i].get();
        x[i] = ele - min;
        y[ele - min] = y[ele - min] + 1;
    }

    int y_size = static_cast<int>(y.size() - 1);
    for (int i = y_size; i > 0; --i)
    {
        y[i - 1] = y[i - 1] += y[i];
    }

    for (int j = y_size; j >= 0; --j)
    {
        for (int i = 0; i < n; ++i)
        {
            int val = 0, val2 = 0;
            if (i >= n - y[j]) { val = 1; }
            if (x[i] >= j) { val2 = 1; }
            int inc = val - val2, ele = A[i].get();
            A.set(i, ArrayItem(ele + inc));
        }
    }
}

/*
    MIT License

    Copyright (c) 2024 aphitorite

    Permission is hereby granted, free of charge, to any person obtaining a copy
    of this software and associated documentation files (the "Software"), to deal
    in the Software without restriction, including without limitation the rights
    to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
    copies of the Software, and to permit persons to whom the Software is
    furnished to do so, subject to the following conditions:

    The above copyright notice and this permission notice shall be included in all
    copies or substantial portions of the Software.

    THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
    IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
    FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
    AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
    LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
    OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
    SOFTWARE.

*/

void dualSwap(SortArray& A, std::vector<value_type>& keys, int a, int b)
{
    size_t z = static_cast<size_t>(a);
    A[z].get();
    A.swap(z, static_cast<size_t>(b));
    value_type temp = keys[a];
    keys[a] = keys[b];
    keys[b] = temp;
}

void reversal(SortArray& A, std::vector<value_type>& keys, int a, int b)
{
    while (b - a > 1) { --b; dualSwap(A, keys, a, b); ++a; }
}

bool isAdjacent(std::vector<value_type>& keys, int a, int b, int N)
{
    return (keys[a] + 1) % N == keys[b] || (keys[b] + 1) % N == keys[a];
}

int findAdjacent(std::vector<value_type>& keys, int e, int a, int N)
{
    while (!isAdjacent(keys, a, e, N)) { ++a; }
    return a;
}

void AdjacencyPancakeSort(SortArray& A)
{
    int a = 0, N = static_cast<int>(A.size()), b = N;
    if (N == 2)
    {
        reverseArr(A, a, a + 1);
        return;
    }
    std::vector<value_type> keys(N);

    for (int j = a; j < b; ++j)
    {
        int c = 0;
        for (int i = a; i < b; ++i)
        {
            if (i == j) { continue; }
            if (A[i] < A[j] || (A[i] == A[j] && i < j)) { ++c; }
        }
        keys[j - a] = ArrayItem(c);
    }

    while (true)
    {
        int i = a;
        while (i < b - 1 && isAdjacent(keys, i, i + 1, N)) { ++i; }
        if (i == b - 1) { break; }
        if (i == a)
        {
            int j = findAdjacent(keys, a, a + 2, N);
            if (!isAdjacent(keys, j - 1, j, N))
            { reversal(A, keys, a, j); }
            else
            {
                int k = findAdjacent(keys, a, j + 1, N);
                if (!isAdjacent(keys, k - 1, k, N))
                { reversal(A, keys, a, k); }
                else
                {
                    reversal(A, keys, a, j + 1);
                    reversal(A, keys, a, j);
                    reversal(A, keys, a, k + 1);
                    reversal(A, keys, a, a + k - j);
                }
            }
        }
        else
        {
            int j = findAdjacent(keys, a, i + 1, N);
            if (!isAdjacent(keys, j - 1, j, N))
            { reversal(A, keys, a, j); }
            else
            {
                int k = findAdjacent(keys, i, i + 2, N);
                if (k + 1 < b && isAdjacent(keys, k + 1, k, N))
                {
                    reversal(A, keys, a, i + 1);
                    reversal(A, keys, a, k + 1);
                }
                else
                {
                    reversal(A, keys, a, k + 1);
                    reversal(A, keys, a, a + k - i);
                    if (!isAdjacent(keys, k - 1, k, N))
                    {
                        if (j < k)
                        {
                            reversal(A, keys, a, k + 1);
                            reversal(A, keys, a, i + k - j + 1);
                        }
                        else
                        {
                            reversal(A, keys, a, j + 1);
                            reversal(A, keys, a, a + j - k);
                        }
                    }
                }
            }
        }
    }

    int i = a;
    while (keys[i] != 0 && keys[i] != N - 1) { ++i; }
    if (keys[i] == 0)
    {
        if (i == a) { return; }
        reversal(A, keys, a, b);
        i = b - 2 - (i - a);
    }
    else if (i == a)
    {
        reversal(A, keys, a, b);
        return;
    }
    ++i;
    reversal(A, keys, a, i);
    reversal(A, keys, a, b);
    reversal(A, keys, a, b - (i - a));
}

// ****************************************************************************
// *** Bitonic Sort

// from http://www.iti.fh-flensburg.de/lang/algorithmen/sortieren/bitonic/oddn.htm

namespace BitonicSortNS {

static const bool ASCENDING = true;    // sorting direction

static void compare(SortArray& A, int i, int j, bool dir)
{
    if (dir == (A[i] > A[j]))
        A.swap(i, j);
}

static void bitonicMerge(SortArray& A, int lo, int n, bool dir)
{
    if (n > 1)
    {
        int m = largestPowerOfTwoLessThan(n);

        for (int i = lo; i < lo + n - m; i++)
            compare(A, i, i+m, dir);

        bitonicMerge(A, lo, m, dir);
        bitonicMerge(A, lo + m, n - m, dir);
    }
}

static void bitonicSort(SortArray& A, int lo, int n, bool dir)
{
    if (n > 1)
    {
        int m = n / 2;
        bitonicSort(A, lo, m, !dir);
        bitonicSort(A, lo + m, n - m, dir);
        bitonicMerge(A, lo, n, dir);
    }
}

} // namespace BitonicSortNS

void BitonicSort(SortArray& A)
{
    BitonicSortNS::bitonicSort(A, 0, A.size(), BitonicSortNS::ASCENDING);
}

// ****************************************************************************
// *** Bitonic Sort as "Parallel" Sorting Network

// from http://www.iti.fh-flensburg.de/lang/algorithmen/sortieren/bitonic/oddn.htm

// modified to first record the recursively generated swap sequence, and then
// sort it back into the order a parallel sorting network would perform the
// swaps in

namespace BitonicSortNetworkNS {

struct swappair_type
{
    // swapped positions
    unsigned int i,j;

    // depth of recursions: sort / merge
    unsigned int sort_depth, merge_depth;

    swappair_type(unsigned int _i, unsigned int _j,
                  unsigned int _sort_depth, unsigned int _merge_depth)
        : i(_i), j(_j),
          sort_depth(_sort_depth), merge_depth(_merge_depth)
    { }

    // order relation for sorting swaps
    bool operator < (const swappair_type& b) const
    {
        if (sort_depth != b.sort_depth)
            return sort_depth > b.sort_depth;

        if (merge_depth != b.merge_depth)
            return merge_depth < b.merge_depth;

        return i < b.i;
    }
};

typedef std::vector<swappair_type> sequence_type;
std::vector<swappair_type> sequence;

void replay(SortArray& A)
{
    for (sequence_type::const_iterator si = sequence.begin();
         si != sequence.end(); ++si)
    {
        if (A[si->i] > A[si->j])
            A.swap(si->i, si->j);
    }
}

static const bool ASCENDING = true; // sorting direction

static void compare(SortArray& /* A */, unsigned int i, unsigned int j, bool dir,
                    unsigned int sort_depth, unsigned int merge_depth)
{
    // if (dir == (A[i] > A[j])) A.swap(i, j);

    if (dir)
        sequence.push_back( swappair_type(i,j, sort_depth, merge_depth) );
    else
        sequence.push_back( swappair_type(j,i, sort_depth, merge_depth) );
}

static void bitonicMerge(SortArray& A, unsigned int lo, unsigned int n, bool dir,
                         unsigned int sort_depth, unsigned int merge_depth)
{
    if (n > 1)
    {
        unsigned int m = largestPowerOfTwoLessThan(n);

        for (unsigned int i = lo; i < lo + n - m; i++)
            compare(A, i, i + m, dir, sort_depth, merge_depth);

        bitonicMerge(A, lo, m, dir, sort_depth, merge_depth+1);
        bitonicMerge(A, lo + m, n - m, dir, sort_depth, merge_depth+1);
    }
}

static void bitonicSort(SortArray& A, unsigned int lo, unsigned int n, bool dir,
                        unsigned int sort_depth)
{
    if (n > 1)
    {
        unsigned int m = n / 2;
        bitonicSort(A, lo, m, !dir, sort_depth+1);
        bitonicSort(A, lo + m, n - m, dir, sort_depth+1);
        bitonicMerge(A, lo, n, dir, sort_depth, 0);
    }
}

void sort(SortArray& A)
{
    sequence.clear();
    bitonicSort(A, 0, A.size(), BitonicSortNS::ASCENDING, 0);
    std::sort(sequence.begin(), sequence.end());
    replay(A);
    sequence.clear();
}

} // namespace BitonicSortNS

void BitonicSortNetwork(SortArray& A)
{
    BitonicSortNetworkNS::sort(A);
}

// ****************************************************************************
// *** Batcher's Odd-Even Merge Sort as "Parallel" Sorting Network

// from http://www.iti.fh-flensburg.de/lang/algorithmen/sortieren/networks/oemen.htm

// modified to first record the recursively generated swap sequence, and then
// sort it back into the order a parallel sorting network would perform the
// swaps in

namespace BatcherSortNetworkNS {

struct swappair_type
{
    // swapped positions
    unsigned int i,j;

    // depth of recursions: sort / merge
    unsigned int sort_depth, merge_depth;

    swappair_type(unsigned int _i, unsigned int _j,
                  unsigned int _sort_depth, unsigned int _merge_depth)
        : i(_i), j(_j),
          sort_depth(_sort_depth), merge_depth(_merge_depth)
    { }

    // order relation for sorting swaps
    bool operator < (const swappair_type& b) const
    {
        if (sort_depth != b.sort_depth)
            return sort_depth > b.sort_depth;

        if (merge_depth != b.merge_depth)
            return merge_depth > b.merge_depth;

        return i < b.i;
    }
};

typedef std::vector<swappair_type> sequence_type;
std::vector<swappair_type> sequence;

void replay(SortArray& A)
{
    for (sequence_type::const_iterator si = sequence.begin();
         si != sequence.end(); ++si)
    {
        if (A[si->i] > A[si->j])
            A.swap(si->i, si->j);
    }
}

static void compare(SortArray& A, unsigned int i, unsigned int j,
                    unsigned int sort_depth, unsigned int merge_depth)
{
    // skip all swaps beyond end of array
    ASSERT(i < j);
    if (j >= A.size()) return;

    sequence.push_back( swappair_type(i,j, sort_depth, merge_depth) );

    //if (A[i] > A[j]) A.swap(i, j);
}

// lo is the starting position and n is the length of the piece to be merged, r
// is the distance of the elements to be compared
static void oddEvenMerge(SortArray& A, unsigned int lo, unsigned int n, unsigned int r,
                         unsigned int sort_depth, unsigned int merge_depth)
{
    unsigned int m = r * 2;
    if (m < n)
    {
        // even subsequence
        oddEvenMerge(A, lo, n, m, sort_depth, merge_depth+1);
        // odd subsequence
        oddEvenMerge(A, lo + r, n, m, sort_depth, merge_depth+1);

        for (unsigned int i = lo + r; i + r < lo + n; i += m)
            compare(A, i, i + r, sort_depth, merge_depth);
    }
    else {
        compare(A, lo, lo + r, sort_depth, merge_depth);
    }
}

// sorts a piece of length n of the array starting at position lo
static void oddEvenMergeSort(SortArray& A, unsigned int lo, unsigned int n,
                             unsigned int sort_depth)
{
    if (n > 1)
    {
        unsigned int m = n / 2;
        oddEvenMergeSort(A, lo, m, sort_depth+1);
        oddEvenMergeSort(A, lo + m, m, sort_depth+1);
        oddEvenMerge(A, lo, n, 1, sort_depth, 0);
    }
}

void sort(SortArray& A)
{
    sequence.clear();

    unsigned int n = largestPowerOfTwoLessThan(A.size());
    if (n != A.size()) n *= 2;

    oddEvenMergeSort(A, 0, n, 0);
    std::sort(sequence.begin(), sequence.end());
    replay(A);
    sequence.clear();
}

} // namespace BatcherSortNetworkNS

void BatcherSortNetwork(SortArray& A)
{
    BatcherSortNetworkNS::sort(A);
}

// ****************************************************************************
// *** Smooth Sort

// from http://en.wikipediA.org/wiki/Smoothsort

namespace SmoothSortNS {

static const int LP[] = {
    1, 1, 3, 5, 9, 15, 25, 41, 67, 109,
    177, 287, 465, 753, 1219, 1973, 3193, 5167, 8361, 13529, 21891,
    35421, 57313, 92735, 150049, 242785, 392835, 635621, 1028457,
    1664079, 2692537, 4356617, 7049155, 11405773, 18454929, 29860703,
    48315633, 78176337, 126491971, 204668309, 331160281, 535828591,
    866988873 // the next number is > 31 bits.
};

static void sift(SortArray& A, int pshift, int head)
{
    // we do not use Floyd's improvements to the heapsort sift, because we
    // are not doing what heapsort does - always moving nodes from near
    // the bottom of the tree to the root.

    value_type val = A[head];

    while (pshift > 1)
    {
        int rt = head - 1;
        int lf = head - 1 - LP[pshift - 2];

        if (val.cmp(A[lf]) >= 0 && val.cmp(A[rt]) >= 0)
            break;

        if (A[lf].cmp(A[rt]) >= 0) {
            A.set(head, A[lf]);
            head = lf;
            pshift -= 1;
        }
        else {
            A.set(head, A[rt]);
            head = rt;
            pshift -= 2;
        }
    }

    A.set(head, val);
}

static void trinkle(SortArray& A, int p, int pshift, int head, bool isTrusty)
{
    value_type val = A[head];

    while (p != 1)
    {
        int stepson = head - LP[pshift];

        if (A[stepson].cmp(val) <= 0)
            break; // current node is greater than head. sift.

        // no need to check this if we know the current node is trusty,
        // because we just checked the head (which is val, in the first
        // iteration)
        if (!isTrusty && pshift > 1) {
            int rt = head - 1;
            int lf = head - 1 - LP[pshift - 2];
            if (A[rt].cmp(A[stepson]) >= 0 ||
                A[lf].cmp(A[stepson]) >= 0)
                break;
        }

        A.set(head, A[stepson]);

        head = stepson;
        //int trail = Integer.numberOfTrailingZeros(p & ~1);
        int trail = __builtin_ctz(p & ~1);
        p >>= trail;
        pshift += trail;
        isTrusty = false;
    }

    if (!isTrusty) {
        A.set(head, val);
        sift(A, pshift, head);
    }
}

void sort(SortArray& A, int lo, int hi)
{
    int head = lo; // the offset of the first element of the prefix into m

    // These variables need a little explaining. If our string of heaps
    // is of length 38, then the heaps will be of size 25+9+3+1, which are
    // Leonardo numbers 6, 4, 2, 1.
    // Turning this into a binary number, we get b01010110 = 0x56. We represent
    // this number as a pair of numbers by right-shifting all the zeros and
    // storing the mantissa and exponent as "p" and "pshift".
    // This is handy, because the exponent is the index into L[] giving the
    // size of the rightmost heap, and because we can instantly find out if
    // the rightmost two heaps are consecutive Leonardo numbers by checking
    // (p&3)==3

    int p = 1; // the bitmap of the current standard concatenation >> pshift
    int pshift = 1;

    while (head < hi)
    {
        if ((p & 3) == 3) {
            // Add 1 by merging the first two blocks into a larger one.
            // The next Leonardo number is one bigger.
            sift(A, pshift, head);
            p >>= 2;
            pshift += 2;
        }
        else {
            // adding a new block of length 1
            if (LP[pshift - 1] >= hi - head) {
                // this block is its final size.
                trinkle(A, p, pshift, head, false);
            } else {
                // this block will get merged. Just make it trusty.
                sift(A, pshift, head);
            }

            if (pshift == 1) {
                // LP[1] is being used, so we add use LP[0]
                p <<= 1;
                pshift--;
            } else {
                // shift out to position 1, add LP[1]
                p <<= (pshift - 1);
                pshift = 1;
            }
        }
        p |= 1;
        head++;
    }

    trinkle(A, p, pshift, head, false);

    while (pshift != 1 || p != 1)
    {
        if (pshift <= 1) {
            // block of length 1. No fiddling needed
            //int trail = Integer.numberOfTrailingZeros(p & ~1);
            int trail = __builtin_ctz(p & ~1);
            p >>= trail;
            pshift += trail;
        }
        else {
            p <<= 2;
            p ^= 7;
            pshift -= 2;

            // This block gets broken into three bits. The rightmost bit is a
            // block of length 1. The left hand part is split into two, a block
            // of length LP[pshift+1] and one of LP[pshift].  Both these two
            // are appropriately heapified, but the root nodes are not
            // necessarily in order. We therefore semitrinkle both of them

            trinkle(A, p >> 1, pshift + 1, head - LP[pshift] - 1, true);
            trinkle(A, p, pshift, head - 1, true);
        }

        head--;
    }
}

} // namespace SmoothSortNS

void SmoothSort(SortArray& A)
{
    return SmoothSortNS::sort(A, 0, A.size()-1);
}

// ****************************************************************************
// *** Stooge Sort

void StoogeSort(SortArray& A, int i, int j)
{
    if (A[i] > A[j])
    {
        A.swap(i, j);
    }

    if (j - i + 1 >= 3)
    {
        int t = (j - i + 1) / 3;

        A.mark(i, 3);
        A.mark(j, 3);

        StoogeSort(A, i, j-t);
        StoogeSort(A, i+t, j);
        StoogeSort(A, i, j-t);

        A.unmark(i);
        A.unmark(j);
    }
}

void StoogeSort(SortArray& A)
{
    StoogeSort(A, 0, A.size()-1);
}

void BadSort(SortArray& A)
{
    for (size_t i = 0; i < A.size(); i++)
    {
        size_t shortest = i;
        for (size_t j = i; j < A.size(); j++)
        {
            bool isShortest = true;
            for (size_t k = j + 1; k < A.size(); k++)
            {
                if (A[j] > A[k])
                {
                    isShortest = false;
                    break;
                }
            }
            if (isShortest)
            {
                shortest = j;
                break;
            }
        }
        A.swap(i, shortest);
    }
}

// ****************************************************************************
// *** Slow Sort

void SlowSort(SortArray& A, int i, int j)
{
    if (i >= j) return;

    int m = (i + j) / 2;

    SlowSort(A, i, m);
    SlowSort(A, m+1, j);

    if (A[m] > A[j])
        A.swap(m, j);

    A.mark(j, 2);

    SlowSort(A, i, j-1);

    A.unmark(j);
}

void SlowSort(SortArray& A)
{
    SlowSort(A, 0, A.size()-1);
}

// ****************************************************************************
// *** Cycle Sort

// Adapted from http://en.wikipedia.org/wiki/Cycle_sort

void CycleSort(SortArray& vec_arr, ssize_t n)
{
    std::atomic<ssize_t> cycleStart{ 0 }, rank{ 0 };
    vec_arr.watch(&cycleStart, 16);

    vec_arr.watch(&rank, 5);

    // Loop through the array to find cycles to rotate.
    for (; cycleStart.load() < n - 1; ++cycleStart)
    {
        value_type& item = vec_arr.get_mutable(cycleStart.load());
        do {
            // Find where to put the item.
            rank.store(cycleStart.load());
            for (ssize_t i = cycleStart.load() + 1; i < n; ++i)
            {
                if (vec_arr[i] < item) { rank++; }
            }

            // If the item is already there, this is a 1-cycle.
            if (rank.load() == cycleStart.load()) {
                vec_arr.mark(rank.load(), 2);
                break;
            }

            // Otherwise, put the item after any duplicates.
            while (item == vec_arr[rank.load()]) { rank++; }

            // Put item into right place and colorize
            counted_swap(vec_arr.get_mutable(rank.load()), item);
            vec_arr.mark(rank.load(), 2);

            // Continue for rest of the cycle.
        } while (rank.load() != cycleStart.load());
    }

    vec_arr.unwatch_all();
}

void CycleSort(SortArray& A)
{
    CycleSort(A, A.size());
}


// ****************************************************************************
// *** Pairwise Sorting Network (Recursive and Iterative)
/*
    Copyright (c) 2021 aphitorite
    Permission is hereby granted, free of charge, to any person obtaining a copy
    of this software and associated documentation files (the "Software"), to deal
    in the Software without restriction, including without limitation the rights
    to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
    copies of the Software, and to permit persons to whom the Software is
    furnished to do so, subject to the following conditions:
    The above copyright notice and this permission notice shall be included in all
    copies or substantial portions of the Software.
    THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
    IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
    FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
    AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
    LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
    OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
    SOFTWARE.
*/

void compSwap(SortArray& A, size_t a, size_t b)
{
    size_t n = A.size();
    if (b < n && A[a] > A[b]) { A.swap(a, b); }
}
void PairwiseMerge(SortArray& A, size_t a, size_t b)
{
    size_t m = (a + b) / 2, m1 = (a + m) / 2, g = m - m1;
    for (size_t i = 0; m1 + i < m; ++i)
    {
        for (size_t j = m1, k = g; k > 0; k >>= 1, j -= k - (i & k))
        {
            compSwap(A, j + i, j + i + k);
        }
    }
    if (b - a > 4) { PairwiseMerge(A, m, b); }
}

void PairwiseMergeSort(SortArray& A, size_t a, size_t b)
{
    size_t m = (a + b) / 2;
    for (size_t i = a, j = m; i < m; ++i, ++j)
    {
        compSwap(A, i, j);
    }
    if (b - a > 2)
    {
        PairwiseMergeSort(A, a, m);
        PairwiseMergeSort(A, m, b);
        PairwiseMerge(A, a, b);
    }
}

void PairwiseSort(SortArray& A)
{
    size_t end = A.size();
    size_t n = 1;
    for (; n < end; n <<= 1) {}
    PairwiseMergeSort(A, 0, n);
}

void PairwiseIterativeSort(SortArray& A)
{
    size_t end = A.size(), n = 1;
    for (; n < end; n <<= 1) {}
    for (size_t k = n >> 1; k > 0; k >>= 1)
    {
        for (size_t j = 0; j < end; j += k << 1)
        {
            for (size_t i = 0; i < k; ++i)
            {
                compSwap(A, j + i, j + i + k);
            }
        }
    }
    for (size_t k = 2; k < n; k <<= 1)
    {
        for (size_t j = k >> 1; j > 0; j >>= 1)
        {
            for (size_t i = 0; i < end; i += k << 1)
            {
                for (size_t m = j; m < ((k - j) << 1); m += j << 1)
                {
                    for (size_t o = 0; o < j; ++o)
                    {
                        compSwap(A, i + m + o, i + m + j + o);
                    }
                }
            }
        }
    }
}

// ****************************************************************************
// *** American Flag Sort
// Adapted from https://en.wikipedia.org/wiki/American_flag_sort
/*
    Copyright 2017 Justin Wetherell

    Licensed under the Apache License, Version 2.0 (the "License");
    you may not use this file except in compliance with the License.
    You may obtain a copy of the License at

        http://www.apache.org/licenses/LICENSE-2.0

    Unless required by applicable law or agreed to in writing, software
    distributed under the License is distributed on an "AS IS" BASIS,
    WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
    See the License for the specific language governing permissions and
    limitations under the License.
*/

size_t getDigit(size_t num, size_t divisor, size_t buckets) { return (num / divisor) % buckets; }

int getMaxNumberOfDigits(SortArray& A, size_t len, int buckets)
{
    int max = std::numeric_limits<int>::min();
    int temp = 0;
    for (size_t i = 0; i < len; ++i)
    {
        int ele = A[i].get();
        temp = static_cast<int>(log(ele) / log(buckets)) + 1;
        if (temp > max) { max = temp; }
    }
    return max;
}

void sort(SortArray& A, size_t start, size_t len, size_t divisor)
{
    size_t buckets = 128;
    std::vector<int> count(buckets, 0);
    std::vector<int> offset(buckets, 0);
    size_t digit = 0;

    for (size_t i = start; i < len; ++i)
    {
        int d = A[i].get();
        size_t l = d;
        digit = getDigit(l, divisor, buckets);
        count[digit] = count[digit] + 1;
    }
    int s = start;
    offset[0] = s;

    for (size_t i = 1; i < buckets; ++i)
    {
        offset[i] = count[i - 1] + offset[i - 1];
    }

    for (size_t b = 0; b < buckets; ++b)
    {
        while (count[b] > 0)
        {
            size_t origin = offset[b], from = origin;
            int num = A[from];
            do
            {
                size_t m = num;
                digit = getDigit(m, divisor, buckets);
                size_t to = offset[digit];

                offset[digit] = offset[digit] + 1;
                count[digit] = count[digit] - 1;
                
                int temp = A[to].get();
                A.set(to, ArrayItem(num));

                num = temp;
                from = to;
            } 
            while (from != origin);
        }
    }
    if (divisor > 1)
    {
        for (size_t i = 0; i < buckets; ++i)
        {
            size_t begin = 0;
            if (i > 0) { begin = offset[i - 1]; }
            else { begin = start; }
            size_t last = offset[i];

            if (last - begin > 1)
            {
                sort(A, begin, last, divisor / buckets);
            }
        }
    }
}

void AmericanFlagSort(SortArray& A)
{
    size_t len = A.size();
    int buckets = 128, max = 1;
    int numberOfDigits = getMaxNumberOfDigits(A, len, buckets); // Max number of digits

    for (int i = 0; i < numberOfDigits - 1; ++i) { max *= buckets; }

    size_t m = max;
    sort(A, 0, len, m);
}


// ****************************************************************************
// *** Strand Sort

/*
    MIT License

    Copyright (c) 2020-2021 aphitorite

    Permission is hereby granted, free of charge, to any person obtaining a copy
    of this software and associated documentation files (the "Software"), to deal
    in the Software without restriction, including without limitation the rights
    to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
    copies of the Software, and to permit persons to whom the Software is
    furnished to do so, subject to the following conditions:

    The above copyright notice and this permission notice shall be included in all
    copies or substantial portions of the Software.

    THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
    IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
    FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
    AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
    LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
    OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
    SOFTWARE.
*/

void mergeTo(SortArray& A, std::vector<value_type>& subList, size_t a, size_t m, size_t b)
{
    size_t i = 0, s = m - a;
    while (i < s && m < b)
    {
        if (subList[i] < A[m]) 
        { 
            A.set(a++, subList[i++]);
        }
        else
        {
            A.set(a++, A[m++]);
        }
    }
    while (i < s)
    {
        A.set(a++, subList[i++]);
    }
}

void StrandSort(SortArray& A)
{
    size_t n = A.size(), j = n, k = j;
    std::vector<value_type> subList(n);
    while (j > 0)
    {
        subList[0] = A[0];
        --k;
        for (size_t i = 0, p = 0, m = 1; m < j; ++m)
        {
            if (A[m] >= subList[i])
            {
                subList[++i] = A[m];
                --k;
            }
            else
            {
                A.set(p++, A[m]);
            }
        }
        mergeTo(A, subList, k, j, n);
        j = k;
    }
}


// ****************************************************************************
// *** New Shuffle Merge Sort
/*
    MIT License

    Copyright (c) 2021 EmeraldBlock

    Permission is hereby granted, free of charge, to any person obtaining a copy
    of this software and associated documentation files (the "Software"), to deal
    in the Software without restriction, including without limitation the rights
    to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
    copies of the Software, and to permit persons to whom the Software is
    furnished to do so, subject to the following conditions:

    The above copyright notice and this permission notice shall be included in all
    copies or substantial portions of the Software.

    THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
    IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
    FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
    AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
    LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
    OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
    SOFTWARE.
*/

/*
  * Implements https://www.sciencedirect.com/science/article/pii/S1877050910005478.
  *
  * The shuffle algorithm is at https://arxiv.org/abs/0805.1598.
  * Note that the unshuffle algorithm is not the shuffle algorithm in reverse,
  * but rather, it is a variation of the shuffle algorithm.
  *
  * See also a proof of the time complexity at https://arxiv.org/abs/1508.00292.
  * The implementation is based on the pseudocode found in this.
*/

void rotateEqual(SortArray& A, size_t a, size_t b, size_t size)
{
    for (size_t i = 0; i < size; ++i)
    {
        A.swap(a + i, b + i);
    }
}

void rotateArray(SortArray& A, size_t mid, size_t a, size_t b)
{
    while (a > 0 && b > 0)
    {
        if (a > b)
        {
            rotateEqual(A, mid - b, mid, b);
            mid -= b;
            a -= b;
        }
        else
        {
            rotateEqual(A, mid - a, mid, a);
            mid += a;
            b -= a;
        }
    }
}

void shuffleEasy(SortArray& A, size_t start, size_t size)
{
    for (size_t i = 1; i < size; i *= 3)
    {
        value_type val = A[start + i - 1];
        for (size_t j = i * 2 % size; j != i; j = j * 2 % size)
        {
            value_type nval = A[start + j - 1];
            A.set(start + j - 1, val);
            val = nval;
        }
        A.set(start + i - 1, val);
    }
}

void shuffleArray(SortArray& A, size_t start, size_t end)
{
    while (end - start > 1)
    {
        size_t n = (end - start) / 2, l = 1;
        while (l * 3 - 1 <= 2 * n) { l *= 3; }
        size_t m = (l - 1) / 2;
        rotateArray(A, start + n, n - m, m);
        shuffleEasy(A, start, l);
        start += l - 1;
    }
}

void rotateShuffledEqual(SortArray& A, size_t a, size_t b, size_t size)
{
    for (size_t i = 0; i < size; i += 2)
    {
        A.swap(a + i, b + i);
    }
}

void rotateShuffled(SortArray& A, size_t mid, size_t a, size_t b)
{
    while (a > 0 && b > 0)
    {
        if (a > b)
        {
            rotateShuffledEqual(A, mid - b, mid, b);
            mid -= b;
            a -= b;
        }
        else
        {
            rotateShuffledEqual(A, mid - a, mid, a);
            mid += a;
            b -= a;
        }
    }
}

void rotateShuffledOuter(SortArray& A, size_t mid, size_t a, size_t b) 
{
    if (a > b) 
    {
        rotateShuffledEqual(A, mid - b, mid + 1, b);
        mid -= b;
        a -= b;
        rotateShuffled(A, mid, a, b);
    }
    else 
    {
        rotateShuffledEqual(A, mid - a, mid + 1, a);
        mid += a + 1;
        b -= a;
        rotateShuffled(A, mid, a, b);
    }
}

void unshuffleEasy(SortArray& A, size_t start, size_t size) 
{
    for (size_t i = 1; i < size; i *= 3) 
    {
        size_t prev = i;
        value_type val = A[start + i - 1];
        for (size_t j = i * 2 % size; j != i; j = j * 2 % size) {
            A.set(start + prev - 1, A[start + j - 1]);
            prev = j;
        }
        A.set(start + prev - 1, val);
    }
}

void unshuffle(SortArray& A, size_t start, size_t end) 
{
    while (end - start > 1) 
    {
        size_t n = (end - start) / 2, l = 1;
        while (l * 3 - 1 <= 2 * n) { l *= 3; }
        size_t m = (l - 1) / 2;

        rotateShuffledOuter(A, start + 2 * m, 2 * m, 2 * n - 2 * m);
        unshuffleEasy(A, start, l);
        start += l - 1;
    }
}

void mergeUp(SortArray& A, size_t start, size_t end, bool type)
{
    size_t i = start, j = i + 1;
    while (j < end)
    {
        if (A[i] < A[j] || (!type && A[i] == A[j]))
        {
            ++i;
            if (i == j)
            {
                ++j;
                type = !type;
            }
        }
        else if (end - j == 1)
        {
            rotateArray(A, j, j - i, 1);
            break;
        }
        else
        {
            size_t r = 0;
            if (type)
            {
                while (j + 2 * r < end && A[j + 2 * r] <= A[i]) { ++r; }
            }
            else
            {
                while (j + 2 * r < end && A[j + 2 * r] < A[i]) { ++r; }
            }
            --j;
            unshuffle(A, j, j + 2 * r);
            rotateArray(A, j, j - i, r);
            i += r + 1;
            j += 2 * r + 1;
        }
    }
}

void mergeArray(SortArray& A, size_t start, size_t mid, size_t end)
{
    if (mid - start <= end - mid)
    {
        shuffleArray(A, start, end);
        mergeUp(A, start, end, true);
    }
    else
    {
        shuffleArray(A, start + 1, end);
        mergeUp(A, start, end, false);
    }
}

size_t ceilPowerOfTwo(size_t x)
{
    --x;
    for (size_t i = 16; i > 0; i >>= 1) { x |= x >> i; }
    return ++x;
}

void sortLarge(SortArray& A, size_t len)
{
    for (size_t subarrayCount = ceilPowerOfTwo(len), wholeI = len / subarrayCount, fracI = len % subarrayCount; subarrayCount > 1; )
    {
        for (size_t whole = 0, frac = 0; whole < len; )
        {
            size_t start = whole;
            whole += wholeI;
            frac += fracI;
            if (frac >= subarrayCount)
            {
                ++whole;
                frac -= subarrayCount;
            }
            size_t mid = whole;
            whole += wholeI;
            frac += fracI;
            if (frac >= subarrayCount)
            {
                ++whole;
                frac -= subarrayCount;
            }
            mergeArray(A, start, mid, whole);
        }
        subarrayCount >>= 1;
        wholeI <<= 1;
        if (fracI >= subarrayCount)
        {
            ++wholeI;
            fracI -= subarrayCount;
        }
    }
}

void mergeSortArray(SortArray& A, size_t len)
{
    if (len < 1 << 15)
    {
        for (size_t subarrayCount = ceilPowerOfTwo(len); subarrayCount > 1; subarrayCount >>= 1)
        {
            for (size_t i = 0; i < subarrayCount; i += 2)
            {
                mergeArray(A, len * i / subarrayCount, len * (i + 1) / subarrayCount, len * (i + 2) / subarrayCount);
            }
        }
    }
    else
    {
        sortLarge(A, len);
    }
}

void NewShuffleMergeSort(SortArray& A) 
{
    size_t len = A.size();
    mergeSortArray(A, len);
}

// ****************************************************************************
// *** Andrey Astrelin's In-Place Merge Sort

void sortVector(SortArray& A, size_t a, size_t b)
{
    while (b > 1)
    {
        size_t k = 0;
        for (size_t i = 1; i < b; ++i)
        {
            if (A[a + k] > A[a + i]) { k = i; }
        }
        A.swap(a, a + k);
        ++a; --b;
    }
}

void aswap(SortArray& A, size_t arr1, size_t arr2, size_t l)
{
    while (l-- > 0)
    {
        A.swap(arr1, arr2);
        ++arr1; ++arr2;
    }
}

int backMerge(SortArray& A, size_t arr1, size_t l1, size_t arr2, size_t l2)
{
    size_t arr0 = arr2 + l1;
    for (;;)
    {
        if (A[arr1] > A[arr2])
        {
            A.swap(arr1, arr0);
            --arr1; --arr0;
            if (--l1 == 0) { return 0; }
        }
        else
        {
            A.swap(arr2, arr0);
            --arr2; --arr0;
            if (--l2 == 0) { break; }
        }
    }
    size_t res = l1;
    do
    {
        A.swap(arr1, arr0);
        --arr1; --arr0;
    } 
    while (--l1 != 0);
    return res;
}

void rMerge(SortArray& A, size_t a, size_t l, size_t r)
{
    for (size_t i = 0; i < l; i += r)
    {
        size_t q = i;
        for (size_t j = i + r; j < l; j += r)
        {
            if (A[a + q] > A[a + j]) { q = j; }
        }
        if (q != i) { aswap(A, a + i, a + q, r); }
        if (i != 0)
        {
            aswap(A, a + l, a + i, r);
            backMerge(A, a + (l + r - 1), r, a + (i - 1), r);
        }
    }
}

size_t rbnd(size_t len)
{
    len = len / 2;
    size_t k = 0;
    for (size_t i = 1; i < len; i *= 2) { ++k; }
    len /= k;
    for (k = 1; k <= len; k *= 2) {}
    return k;
}

void msort(SortArray& A, size_t a, size_t len)
{
    if (len < 12) { sortVector(A, a, len); return; }
    size_t r = rbnd(len), lr = (len / r - 1) * r;
    for (size_t p = 2; p <= lr; p += 2)
    {
        if (A[a + (p - 2)] > A[a + (p - 1)]) 
        { A.swap(a + (p - 2), a + (p - 1)); }
        if ((p & 2) != 0) { continue; }
        aswap(A, a + (p - 2), a + p, 2);
        size_t m = len - p, q = 2;
        for (;;)
        {
            size_t q0 = 2 * q;
            if (q0 > m || (p & q0) != 0) { break; }
            backMerge(A, a + (p - q - 1), q, a + (p + q - 1), q);
            q = q0;
        }
        backMerge(A, a + (p + q - 1), q, a + (p - q - 1), q);
        size_t q1 = q;
        q *= 2;
        while ((q & p) == 0)
        {
            q *= 2;
            rMerge(A, a + (p - q), q, q1);
        }
    }

    size_t q1 = 0;
    for (size_t q = r; q < lr; q *= 2)
    {
        if ((lr & q) != 0)
        {
            q1 += q;
            if (q1 != q)
            {
                rMerge(A, a + (lr - q1), q1, r);
            }
        }
    }

    size_t s = len - lr;
    msort(A, a + lr, s);
    aswap(A, a, a + lr, s);
    s += backMerge(A, a + (s - 1), s, a + (lr - 1), lr - s);
    msort(A, a, s);
}

void AndreyMergeSort(SortArray& A)
{
    msort(A, 0, A.size());
}


// ****************************************************************************
// *** Proportion Extend Merge Sort

/*
    MIT License

    Copyright (c) 2023 aphitorite

    Permission is hereby granted, free of charge, to any person obtaining a copy
    of this software and associated documentation files (the "Software"), to deal
    in the Software without restriction, including without limitation the rights
    to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
    copies of the Software, and to permit persons to whom the Software is
    furnished to do so, subject to the following conditions:

    The above copyright notice and this permission notice shall be included in all
    copies or substantial portions of the Software.

    THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
    IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
    FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
    AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
    LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
    OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
    SOFTWARE.
*/

void blockSwap(SortArray& A, size_t a, size_t b, size_t s)
{
    while (s-- > 0)
    {
        A.swap(a, b);
        ++a; ++b;
    }
}

size_t partition(SortArray& A, size_t a, size_t b, size_t p)
{
    size_t i = a - 1, j = b;
    while (true)
    {
        do { ++i; } while (i < j && A[i] < A[p]);
        do { --j; } while (j >= i && A[j] > A[p]);
        if (i < j) { A.swap(i, j); }
        else { return i; }
    }
}

void mergeFW(SortArray& A, size_t a, size_t m, size_t b, size_t p)
{
    size_t pLen = m - a, i = 0, j = m, k = a;
    blockSwap(A, a, p, pLen);
    while (i < pLen && j < b)
    {
        if (A[p + i] <= A[j]) { A.swap(k, p + i); ++k; ++i; }
        else { A.swap(k, j); ++k; ++j; }
    }
    while (i < pLen) { A.swap(k, p + i); ++k; ++i; }
}

void mergeBW(SortArray& A, size_t a, size_t m, size_t b, size_t p)
{
    size_t pLen = b - m;
    int i = static_cast<int>(pLen - 1), 
        j = static_cast<int>(m - 1), 
        k = static_cast<int>(b - 1), 
        z = static_cast<int>(a);
    blockSwap(A, m, p, pLen);
    while (i >= 0 && j >= z)
    {
        if (A[p + i] >= A[j]) { A.swap(k, p + i); --k; --i; }
        else { A.swap(k, j); --k; --j; }
    }
    while (i >= 0) { A.swap(k, p + i); --k; --i; }
}

void smartMerge(SortArray& A, size_t a, size_t m, size_t b, size_t p)
{
    if (m - a < b - m) { mergeFW(A, a, m, b, p); }
    else { mergeBW(A, a, m, b, p); }
}

void mergeTo(SortArray& A, size_t a, size_t m, size_t b, size_t p)
{
    size_t i = a, j = m;
    while (i < m && j < b)
    {
        if (A[i] <= A[j]) { A.swap(p, i); ++p; ++i; }
        else { A.swap(p, j); ++p; ++j; }
    }
    while (i < m) { A.swap(p, i); ++p; ++i; }
    while (j < b) { A.swap(p, j); ++p; ++j; }
}

void pingPongMerge(SortArray& A, size_t a, size_t m1, size_t m, size_t m2, size_t b, size_t p)
{
    size_t p1 = p + m - a, pEnd = p + b - a;
    mergeTo(A, a, m1, m, p);
    mergeTo(A, m, m2, b, p1);
    mergeTo(A, p, p1, pEnd, a);
}

void mergeSort(SortArray& A, size_t a, size_t b, size_t p)
{
    const size_t min_insert = 8;
    size_t n = b - a, j = n;
    for (; (j + 3) / 4 >= min_insert; j = (j + 3) / 4) {}
    for (size_t i = a; i < b; i += j)
    {
        BinaryInsertSort(A, i, std::min(b, i + j));
    }
    for (size_t i; j < n; j *= 4)
    {
        for (i = a; i + 2 * j < b; i += 4 * j)
        {
            pingPongMerge(A, i, i + j, i + 2 * j, std::min(i + 3 * j, b), std::min(i + 4 * j, b), p);

        }
        if (i + j < b) { mergeBW(A, i, i + j, b, p); }
    }
}

void smartMergeSort(SortArray& A, size_t a, size_t b, size_t p, size_t pb) 
{
    if (b - a <= pb - p) { mergeSort(A, a, b, p); return; }
    size_t m = (a + b) >> 1;
    mergeSort(A, a, m, p);
    mergeSort(A, m, b, p);
    mergeFW(A, a, m, b, p);
}

void peSort(SortArray& A, size_t a, size_t m, size_t b)
{
    size_t n = b - a;
    const size_t min_insert = 8;
    if (n < 4 * min_insert) { BinaryInsertSort(A, a, b); return; }
    if (m - a <= n / 3)
    {
        size_t t = (n + 2) / 3;
        smartMergeSort(A, m, b - t, b - t, b);
        smartMerge(A, a, m, b - t, b - t);
        m = b - t;
    }
    size_t m1 = (a + m) >> 1, m2 = partition(A, m, b, m1);
    size_t i = m, j = m2;
    while (i > m1) { --i; --j; A.swap(i, j); }
    m = m2 - (m - m1);
    if (m - m1 < b - m2)
    {
        mergeSort(A, m1, m, m2);
        smartMerge(A, a, m1, m, m2);
        peSort(A, m + 1, m2, b);
    }
    else
    {
        mergeSort(A, m2, b, m1);
        smartMerge(A, m + 1, m2, b, m1);
        peSort(A, a, m1, m);
    }
}

void ProportionMergeSort(SortArray& A)
{
    peSort(A, 0, 0, A.size());
}
 
// ****************************************************************************
// *** Buffer Partition Merge Sort

/*
 *
    MIT License

    Copyright (c) 2021 yuji, implemented by aphitorite

    Permission is hereby granted, free of charge, to any person obtaining a copy
    of this software and associated documentation files (the "Software"), to deal
    in the Software without restriction, including without limitation the rights
    to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
    copies of the Software, and to permit persons to whom the Software is
    furnished to do so, subject to the following conditions:

    The above copyright notice and this permission notice shall be included in all
    copies or substantial portions of the Software.

    THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
    IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
    FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
    AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
    LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
    OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
    SOFTWARE.
*
*/

void shiftBW2(SortArray& A, int a, int m, int b)
{
    while (m > a) 
    {  
        --b; --m;
        A.swap(static_cast<size_t>(b), static_cast<size_t>(m));
    }
}

void inPlaceMerge(SortArray& A, int a, int m, int b)
{
    int i = a, j = m, k = 0;
    while (i < j && j < b)
    {
        if (A[i] > A[j])
        {
            k = j; ++k;
            while (k < b && A[i] > A[k]) { ++k; }
            rotate(A, static_cast<size_t>(i), static_cast<size_t>(j), static_cast<size_t>(k));
            i += k - j;
            j = k;
        }
        else { ++i; }
    }
}

void medianOfThree(SortArray& A, int a, int b)
{
    int m = a + (b - 1 - a) / 2;
    if (A[a] > A[m]) { A.swap(static_cast<size_t>(a), static_cast<size_t>(m)); }
    if (A[m] > A[b - 1]) 
    { 
        A.swap(static_cast<size_t>(m), static_cast<size_t>(b - 1)); 
        if (A[a] > A[m]) { return; }
    }
    A.swap(static_cast<size_t>(a), static_cast<size_t>(m));
}

void medianofMedians(SortArray& A, int a, int b, int s)
{
    int end = b, start = a, i, j;
    bool ad = true;
    while (end - start > 1)
    {
        j = start;
        for (i = start; i + 2 * s <= end; i += s)
        {
            InsertSort(A, static_cast<size_t>(i), static_cast<size_t>(i + s));
            A.swap(static_cast<size_t>(j), static_cast<size_t>(i + s / 2));
            ++j;
        }
        if (i < end)
        {
            InsertSort(A, static_cast<size_t>(i), static_cast<size_t>(end));
            int val = 0;
            if (ad) { val = 1; }
            A.swap(static_cast<size_t>(j), static_cast<size_t>(i + (end - val - i) / 2));
            if ((end - i) % 2 == 0) { ad = !ad; }
            ++j;
        }
        end = j;
    }
}

int partition(SortArray& A, int a, int b)
{
    int i = a, j = b;
    while (true)
    {
        do { ++i; }
        while (i < j && A[i] > A[a]);
        do { --j; }
        while (j >= i && A[j] < A[a]);
        if (i < j) { A.swap(static_cast<size_t>(i), static_cast<size_t>(j)); }
        else { return j; }
    }
}

int quickSelect(SortArray& A, int a, int b, int m)
{
    bool badPartition = false, mom = false;
    int m1 = (m + b + 1) / 2;
    while (true)
    {
        if (badPartition)
        {
            medianofMedians(A, a, b, 5);
            mom = true;
        }
        else { medianOfThree(A, a, b); }
        int p = partition(A, a, b);
        A.swap(static_cast<size_t>(a), static_cast<size_t>(p));
        int l = std::max(1, p - a), r = std::max(1, b - (p + 1));
        badPartition = !mom && (l / r >= 16 || r / l >= 16);
        if (p >= m && p < m1) { return p; }
        else if (p < m) { a = p + 1; }
        else { b = p; }
    }
}

void mergeVector(SortArray& A, int a, int m, int b, int p)
{
    int i = a, j = m;
    while (i < m && j < b)
    {
        if (A[i] <= A[j])
        {
            A.swap(static_cast<size_t>(p), static_cast<size_t>(i));
            ++p; ++i;
        }
        else
        {
            A.swap(static_cast<size_t>(p), static_cast<size_t>(j));
            ++p; ++j;
        }
    }
    while (i < m) 
    { 
        A.swap(static_cast<size_t>(p), static_cast<size_t>(i));
        ++p; ++i;
    }
    while (j < b)
    {
        A.swap(static_cast<size_t>(p), static_cast<size_t>(j));
        ++p; ++j;
    }
}

int mergeFW2(SortArray& A, int p, int a, int m, int b)
{
    int i = a, j = m;
    while (i < m && j < b)
    {
        if (A[i] <= A[j]) 
        { 
            A.swap(static_cast<size_t>(p), static_cast<size_t>(i)); 
            ++p; ++i;
        }
        else
        {
            A.swap(static_cast<size_t>(p), static_cast<size_t>(j));
            ++p; ++j;
        }
    }
    if (i < m) { return i; }
    else { return j; }
}

int getMinLevel(int n)
{
    while (n >= 32) { n = (n + 3) / 4; }
    return n;
}

void mergeSort2(SortArray& A, int a, int b, int p)
{
    int len = b - a;
    if (len < 2) { return; }
    int i, pos, j = getMinLevel(len);
    for (i = a; i + j <= b; i += j)
    {
        BinaryInsertSort(A, static_cast<size_t>(i), static_cast<size_t>(i + j));
    }
    BinaryInsertSort(A, static_cast<size_t>(i), static_cast<size_t>(b));
    while (j < len)
    {
        pos = p;
        for (i = a; i + 2 * j <= b; i += 2 * j, pos += 2 * j)
        { mergeVector(A, i, i + j, i + 2 * j, pos); }
        if (i + j < b) { mergeVector(A, i, i + j, b, pos); }
        else
        {
            while (i < b) 
            {
                A.swap(static_cast<size_t>(i), static_cast<size_t>(pos));
                ++i; ++pos;
            }
        }
        j *= 2;
        pos = a;
        for (i = p; i + 2 * j <= p + len; i += 2 * j, pos += 2 * j)
        { mergeVector(A, i, i + j, i + 2 * j, pos); }
        if (i + j < p + len) { mergeVector(A, i, i + j, p + len, pos); }
        else
        {
            while (i < p + len)
            {
                A.swap(static_cast<size_t>(i), static_cast<size_t>(pos));
                ++i; ++pos;
            }
        }
        j *= 2;
    }
}

void sortArray(SortArray& A, int a, int b)
{
    int minLvl = static_cast<int>(sqrt(b - a));
    int m = (a + b + 1) / 2;
    mergeSort2(A, m, b, a);
    while (m - a > minLvl)
    {
        int m1 = (a + m + 1) / 2;
        m1 = quickSelect(A, a, m, m1);
        mergeSort2(A, m1, m, a);
        int bSize = m1 - a;
        int m2 = std::min(m1 + bSize, b);
        m1 = mergeFW2(A, a, m1, m, m2);
        while (m1 < m)
        {
            shiftBW2(A, m1, m, m2);
            m1 = m2 - (m - m1);
            a = m1 - bSize;
            m = m2;
            if (m == b) { break; }
            m2 = std::min(m2 + bSize, b);
            m1 = mergeFW2(A, a, m1, m, m2);
        }
        m = m1;
        a = m1 - bSize;
    }
    BinaryInsertSort(A, static_cast<size_t>(a), static_cast<size_t>(m));
    inPlaceMerge(A, a, m, b);
}

void BufferPartitionMergeSort(SortArray& A)
{
    sortArray(A, 0, A.size());
}