/******************************************************************************
 * src/SortAlgo.cpp
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
 * The algorithms in this file are copyrighted by the original authors. All
 * code is freely available.
 *
 * The source code added by myself (Timo Bingmann) and all modifications are
 * copyright (C) 2013 Timo Bingmann <tb@panthema.net>
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
#include "WSortView.h"

#include <algorithm>
#include <numeric>
#include <limits>
#include <inttypes.h>

typedef ArrayItem value_type;

const struct AlgoEntry g_algolist[] =
{
    { _("Selection Sort"), &SelectionSort, NULL },
    { _("Insertion Sort"), &InsertionSort, NULL },
    { _("Merge Sort"), &MergeSort,
      _("Merge sort which merges two sorted sequences into a shadow array, and then copies it back to the shown array.") },
    { _("Quick Sort (LR ptrs)"), &QuickSortLR,
      _("Quick sort variant with left and right pointers, the middle element is taken as pivot.") },
    { _("Quick Sort (LL ptrs)"), &QuickSortLL,
      _("Quick sort variant from 3rd edition of CLRS: two pointers on left, always picks first element as pivot.") },
    { _("Quick Sort (ternary, LR ptrs)"), &QuickSortTernaryLR,
      _("Ternary-split quick sort variant, adapted from multikey quicksort by Bentley & Sedgewick: partitions \"=<?>=\" using two pairs of pointers at left and right, then copied to middle.") },
    { _("Quick Sort (ternary, LL ptrs)"), &QuickSortTernaryLL,
      _("Ternary-split quick sort variant: partitions \"<>?=\" using two pointers at left and one at right. Afterwards copies the \"=\" to middle.") },
    { _("Bubble Sort"), &BubbleSort, NULL },
    { _("Cocktail Shaker Sort"), &CocktailShakerSort, NULL },
    { _("Gnome Sort"), &GnomeSort, NULL },
    { _("Comb Sort"), &CombSort, NULL },
    { _("Shell Sort"), &ShellSort, NULL },
    { _("Heap Sort"), &HeapSort, NULL },
    { _("Smooth Sort"), &SmoothSort, NULL },
    { _("Odd-Even Sort"), &OddEvenSort, NULL },
    { _("Bitonic Sort"), &BitonicSort, NULL },
    { _("Radix Sort (LSD)"), &RadixSortLSD,
      _("Least significant digit radix sort, which copies item into a shadow array during counting.") },
    { _("Radix Sort (MSD)"), &RadixSortMSD,
      _("Most significant digit radix sort, which permutes items in-place by walking cycles.") },
    { _("std::sort (gcc)"), &StlSort, NULL },
    { _("std::stable_sort (gcc)"), &StlStableSort, NULL },
    { _("std::sort_heap (gcc)"), &StlHeapSort, NULL },
    { _("Bogo Sort"), &BogoSort, NULL },
    { _("Bozo Sort"), &BozoSort, NULL },
    { NULL, NULL, NULL },
};

const size_t g_algolist_size = sizeof(g_algolist) / sizeof(g_algolist[0]) - 1;

// ****************************************************************************
// *** Selection Sort

void SelectionSort(WSortView& a)
{
    volatile ssize_t jMin = 0;
    a.watch(&jMin,2);

    for (size_t i = 0; i < a.size()-1; ++i)
    {
        jMin = i;

        for (size_t j = i+1; j < a.size(); ++j)
        {
            if (a[j] < a[jMin]) {
                a.mark_swap(j, jMin);
                jMin = j;
            }
        }

        a.swap(i, jMin);

        // mark the last good element
        if (i > 0) a.unmark(i-1);
        a.mark(i);
    }
    a.unwatch_all();
}

// ****************************************************************************
// *** Insertion Sort

// swaps every time (keeps all values visible)
void InsertionSort(WSortView& a)
{
    for (size_t i = 1; i < a.size(); ++i)
    {
        value_type key = a[i];
        a.mark(i);

        ssize_t j = i - 1;
        while (j >= 0 && a[j] > key)
        {
            a.swap(j, j+1);
            j--;
        }

        a.unmark(i);
    }
}

// with extra item on stack
void InsertionSort2(WSortView& a)
{
    for (size_t i = 1; i < a.size(); ++i)
    {
        value_type tmp, key = a[i];
        a.mark(i);

        ssize_t j = i - 1;
        while (j >= 0 && (tmp = a[j]) > key)
        {
            a.set(j + 1, tmp);
            j--;
        }
        a.set(j + 1, key);

        a.unmark(i);
    }
}

// ****************************************************************************
// *** Merge Sort (out-of-place with sentinels)

// by myself (Timo Bingmann)

void Merge(WSortView& a, size_t lo, size_t mid, size_t hi)
{
    // mark merge boundaries
    a.mark(lo);
    a.mark(mid,2);
    a.mark(hi-1);

    // allocate output
    std::vector<value_type> out(hi-lo);

    // merge
    size_t i = lo, j = mid, o = 0; // first and second halves
    while (i < mid && j < hi)
    {
        // copy out for fewer time steps
        value_type ai = a[i], aj = a[j];

        out[o++] = (ai < aj ? (++i, ai) : (++j, aj));
    }

    // copy rest
    while (i < mid) out[o++] = a[i++];
    while (j < hi) out[o++] = a[j++];

    ASSERT(o == hi-lo);

    a.unmark(mid);

    // copy back
    for (i = 0; i < hi-lo; ++i)
        a.set(lo + i, out[i]);

    a.unmark(lo);
    a.unmark(hi-1);
}

void MergeSort(WSortView& a, size_t lo, size_t hi)
{
    if (lo + 1 < hi)
    {
        size_t mid = (lo + hi) / 2;

        MergeSort(a, lo, mid);
        MergeSort(a, mid, hi);

        Merge(a, lo, mid, hi);
    }
}

void MergeSort(WSortView& a)
{
    return MergeSort(a, 0, a.size());
}

// ****************************************************************************
// *** Quick Sort Pivot Selection

QuickSortPivotType g_quicksort_pivot = PIVOT_FIRST;

// some quicksort variants use hi inclusive and some exclusive, we require it
// to be _exclusive_. hi == array.end()!
ssize_t QuickSortSelectPivot(WSortView& a, ssize_t lo, ssize_t hi)
{
    if (g_quicksort_pivot == PIVOT_FIRST)
        return lo;

    if (g_quicksort_pivot == PIVOT_LAST)
        return hi-1;

    if (g_quicksort_pivot == PIVOT_MID)
        return (lo + hi) / 2;

    if (g_quicksort_pivot == PIVOT_RANDOM)
        return lo + (rand() % (hi - lo));

    if (g_quicksort_pivot == PIVOT_MEDIAN3)
    {
        ssize_t mid = (lo + hi) / 2;

        // cases if two are equal
        if (a[lo] == a[mid]) return lo;
        if (a[lo] == a[hi-1] || a[mid] == a[hi-1]) return hi-1;

        // cases if three are different
        return a[lo] < a[mid]
            ? (a[mid] < a[hi-1] ? mid : (a[lo] < a[hi-1] ? hi-1 : lo))
            : (a[mid] > a[hi-1] ? mid : (a[lo] < a[hi-1] ? lo : hi-1));
    }

    return lo;
}

const wxChar* g_quicksort_pivot_text[] = {
    _("First Item"),
    _("Last Item"),
    _("Middle Item"),
    _("Random Item"),
    _("Median of Three"),
    NULL
};

// ****************************************************************************
// *** Quick Sort LR (in-place, pointers at left and right, pivot is middle element)

// by myself (Timo Bingmann), based on Hoare's original code

void QuickSortLR(WSortView& a, ssize_t lo, ssize_t hi)
{
    // pick pivot and watch
    volatile ssize_t p = QuickSortSelectPivot(a, lo, hi+1);

    value_type pivot = a[p];
    a.watch(&p,1);

    volatile ssize_t i = lo, j = hi;
    a.watch(&i,2);
    a.watch(&j,2);

    while (i <= j)
    {
        while (a[i] < pivot)
            i++;

        while (a[j] > pivot)
            j--;

        if (i <= j)
        {
            a.swap(i,j);

            // follow pivot if it is swapped
            if (p == i) p = j;
            else if (p == j) p = i;

            i++, j--;
        }
    }

    a.unwatch_all();

    if (lo < j)
        QuickSortLR(a, lo, j);

    if (i < hi)
        QuickSortLR(a, i, hi);
}

void QuickSortLR(WSortView& a)
{
    return QuickSortLR(a, 0, a.size()-1);
}

// ****************************************************************************
// *** Quick Sort LL (in-place, two pointers at left, pivot is first element and moved to right)

// by myself (Timo Bingmann), based on CLRS' 3rd edition

size_t PartitionLL(WSortView& a, size_t lo, size_t hi)
{
    // pick pivot and move to back
    size_t p = QuickSortSelectPivot(a, lo, hi);

    value_type pivot = a[p];
    a.swap(p, hi-1);
    a.mark(hi-1);

    volatile ssize_t i = lo;
    a.watch(&i,2);

    for (size_t j = lo; j < hi-1; ++j)
    {
        if (a[j] <= pivot) {
            a.swap(i, j);
            ++i;
        }
    }

    a.swap(i, hi-1);
    a.unmark(hi-1);
    a.unwatch_all();

    return i;
}

void QuickSortLL(WSortView& a, size_t lo, size_t hi)
{
    if (lo + 1 < hi)
    {
        size_t mid = PartitionLL(a, lo, hi);

        QuickSortLL(a, lo, mid);
        QuickSortLL(a, mid+1, hi);
    }
}

void QuickSortLL(WSortView& a)
{
    return QuickSortLL(a, 0, a.size());
}

// ****************************************************************************
// *** Quick Sort Ternary (in-place, two pointers at left, pivot is first element and moved to right)

// by myself (Timo Bingmann), loosely based on multikey quicksort by B&S

void QuickSortTernaryLR(WSortView& a, ssize_t lo, ssize_t hi)
{
    if (hi <= lo) return;

    int cmp;

    // pick pivot and swap to back
    ssize_t piv = QuickSortSelectPivot(a, lo, hi+1);
    a.swap(piv, hi);
    a.mark(hi);

    const value_type& pivot = a[hi];

    // schema: |p ===  |i <<< | ??? |j >>> |q === |piv
    volatile ssize_t i = lo, j = hi-1;
    volatile ssize_t p = lo, q = hi-1;

    a.watch(&i,2);
    a.watch(&j,2);

    for (;;)
    {
        // partition on left
        while (i <= j && (cmp = a[i].cmp(pivot)) <= 0)
        {
            if (cmp == 0) {
                a.mark(p,3);
                a.swap(i, p++);
            }
            ++i;
        }

        // partition on right
        while (i <= j && (cmp = a[j].cmp(pivot)) >= 0)
        {
            if (cmp == 0) {
                a.mark(q,3);
                a.swap(j, q--);
            }
            --j;
        }

        if (i > j) break;

        // swap item between < > regions
        a.swap(i++, j--);
    }

    // swap pivot to right place
    a.swap(i,hi);
    a.mark_swap(i,hi);

    ssize_t num_less = i - p;
    ssize_t num_greater = q - j;

    // swap equal ranges into center, but avoid swapping equal elements
    j = i-1; i = i+1;

    ssize_t pe = lo + std::min(p-lo, num_less);
    for (ssize_t k = lo; k < pe; k++, j--) {
        a.swap(k,j);
        a.mark_swap(k,j);
    }

    ssize_t qe = hi-1 - std::min(hi-1-q, num_greater-1); // one already greater at end
    for (ssize_t k = hi-1; k > qe; k--, i++) {
        a.swap(i,k);
        a.mark_swap(i,k);
    }

    a.unwatch_all();
    a.unmark_all();

    QuickSortTernaryLR(a, lo, lo + num_less - 1);
    QuickSortTernaryLR(a, hi - num_greater + 1, hi);
}

void QuickSortTernaryLR(WSortView& a)
{
    return QuickSortTernaryLR(a, 0, a.size()-1);
}

// ****************************************************************************
// *** Quick Sort LL (in-place, two pointers at left, pivot is first element and moved to right)

// by myself (Timo Bingmann)

std::pair<ssize_t,ssize_t> PartitionTernaryLL(WSortView& a, ssize_t lo, ssize_t hi)
{
    // pick pivot and swap to back
    ssize_t p = QuickSortSelectPivot(a, lo, hi);

    value_type pivot = a[p];
    a.swap(p, hi-1);
    a.mark(hi-1);

    volatile ssize_t i = lo, k = hi-1;
    a.watch(&i,2);

    for (ssize_t j = lo; j < k; ++j)
    {
        int cmp = a[j].cmp(pivot); // ternary comparison
        if (cmp == 0) {
            a.swap(--k, j);
            --j; // reclassify a[j]
            a.mark(k,3);
        }
        else if (cmp < 0) {
            a.swap(i++, j);
        }
    }

    // unwatch i, because the pivot is swapped there
    // in the first step of the following swap loop.
    a.unwatch_all();

    ssize_t j = i + (hi-k);

    for (ssize_t s = 0; s < hi-k; ++s) {
        a.swap(i+s, hi-1-s);
        a.mark_swap(i+s, hi-1-s);
    }
    a.unmark_all();

    return std::make_pair(i,j);
}

void QuickSortTernaryLL(WSortView& a, size_t lo, size_t hi)
{
    if (lo + 1 < hi)
    {
        std::pair<ssize_t,ssize_t> mid = PartitionTernaryLL(a, lo, hi);

        QuickSortTernaryLL(a, lo, mid.first);
        QuickSortTernaryLL(a, mid.second, hi);
    }
}

void QuickSortTernaryLL(WSortView& a)
{
    return QuickSortTernaryLL(a, 0, a.size());
}

// ****************************************************************************
// *** Bubble Sort

void BubbleSort(WSortView& a)
{
    for (size_t i = 0; i < a.size()-1; ++i)
    {
        for (size_t j = 0; j < a.size()-1 - i; ++j)
        {
            if (a[j] > a[j + 1])
            {
                a.swap(j, j+1);
            }
        }
    }
}

// ****************************************************************************
// *** Cocktail Shaker Sort

// from http://de.wikibooks.org/wiki/Algorithmen_und_Datenstrukturen_in_C/_Shakersort

void CocktailShakerSort(WSortView& a)
{
    size_t lo = 0, hi = a.size()-1, mov = lo;

    while (lo < hi)
    {
        for (size_t i = hi; i > lo; --i)
        {
            if (a[i-1] > a[i])
            {
                a.swap(i-1, i);
                mov = i;
            }
        }

        lo = mov;

        for (size_t i = lo; i < hi; ++i)
        {
            if (a[i] > a[i+1])
            {
                a.swap(i, i+1);
                mov = i;
            }
        }

        hi = mov;
    }
}

// ****************************************************************************
// *** Gnome Sort

// from http://en.wikipedia.org/wiki/Gnome_sort

void GnomeSort(WSortView& a)
{
    for (size_t i = 1; i < a.size(); )
    {
        if (a[i] >= a[i-1])
        {
            ++i;
        }
        else
        {
            a.swap(i, i-1);
            if (i > 1) --i;
        }
    }
}

// ****************************************************************************
// *** Comb Sort

// from http://en.wikipedia.org/wiki/Comb_sort

void CombSort(WSortView& a)
{
    const double shrink = 1.3;

    bool swapped = false;
    size_t gap = a.size();

    while ((gap > 1) || swapped)
    {
        if (gap > 1) {
            gap = (size_t)((float)gap / shrink);
        }

        swapped = false;

        for (size_t i = 0; gap + i < a.size(); ++i)
        {
            if (a[i] > a[i + gap])
            {
                a.swap(i, i+gap);
                swapped = true;
            }
        }
    }
}

// ****************************************************************************
// *** Odd-Even Sort

// from http://en.wikipedia.org/wiki/Odd%E2%80%93even_sort

void OddEvenSort(WSortView& a)
{
    bool sorted = false;

    while (!sorted)
    {
        sorted = true;

        for (size_t i = 1; i < a.size()-1; i += 2)
        {
            if(a[i] > a[i+1])
            {
                a.swap(i, i+1);
                sorted = false;
            }
        }

        for (size_t i = 0; i < a.size()-1; i += 2)
        {
            if(a[i] > a[i+1])
            {
                a.swap(i, i+1);
                sorted = false;
            }
        }
    }
}

// ****************************************************************************
// *** Shell Sort

// with gaps by Robert Sedgewick from http://www.cs.princeton.edu/~rs/shell/shell.c

void ShellSort(WSortView& a)
{
    size_t incs[16] = { 1391376, 463792, 198768, 86961, 33936,
                        13776, 4592, 1968, 861, 336,
                        112, 48, 21, 7, 3, 1 };

    for (size_t k = 0; k < 16; k++)
    {
        for (size_t h = incs[k], i = h; i < a.size(); i++)
        {
            value_type v = a[i];
            size_t j = i;

            while (j >= h && a[j-h] > v)
            {
                a.set(j, a[j-h]);
                j -= h;
            }

            a[j] = v;
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

void HeapSort(WSortView& a)
{
    size_t n = a.size(), i = n / 2;

    // mark heap levels with different colors
    for (size_t j = i; j < n; ++j)
        a.mark(j, log(prevPowerOfTwo(j+1)) / log(2) + 3);

    while (1)
    {
        if (i > 0) {
            // build heap, sift a[i] down the heap
            i--;
        }
        else {
            // pop largest element from heap: swap front to back, and sift
            // front a[0] down the heap
            n--;
            if (n == 0) return;
            a.swap(0,n);

            a.mark(n);
            if (n+1 < a.size()) a.unmark(n+1);
        }

        size_t parent = i;
        size_t child = i*2 + 1;

        // sift operation - push the value of a[i] down the heap
        while (child < n)
        {
            if (child + 1 < n && a[child + 1] > a[child]) {
                child++;
            }
            if (a[child] > a[parent]) {
                a.swap(parent, child);
                parent = child;
                child = parent*2+1;
            }
            else {
                break;
            }
        }

        // mark heap levels with different colors
        a.mark(i, log(prevPowerOfTwo(i+1)) / log(2) + 3);
    }

}

// ****************************************************************************
// *** Radix Sort (counting sort, most significant digit (MSD) first, in-place redistribute)

// by myself (Timo Bingmann)

void RadixSortMSD(WSortView& a, size_t lo, size_t hi, size_t depth)
{
    a.mark(lo); a.mark(hi-1);

    // radix and base calculations
    const unsigned int RADIX = 4;

    unsigned int pmax = floor( log(a.array_max()) / log(RADIX) );
    ASSERT(depth <= pmax);

    size_t base = pow(RADIX, pmax - depth);

    // count digits
    std::vector<size_t> count(RADIX, 0);

    for (size_t i = lo; i < hi; ++i)
    {
        size_t r = a[i].get() / base % RADIX;
        ASSERT(r < RADIX);
        count[r]++;
    }

    // inclusive prefix sum
    std::vector<size_t> bkt(RADIX, 0);
    std::partial_sum(count.begin(), count.end(), bkt.begin());

    // mark bucket boundaries
    for (size_t i = 0; i < bkt.size(); ++i) {
        if (bkt[i] == 0) continue;
        a.mark(lo + bkt[i]-1, 2);
    }

    // reorder items in-place by walking cycles
    for (size_t i=0, j; i < (hi-lo); )
    {
        while ( (j = --bkt[ (a[lo+i].get() / base % RADIX) ]) > i )
        {
            a.swap(lo + i, lo + j);
        }
        i += count[ (a[lo+i].get() / base % RADIX) ];
    }

    a.unmark_all();

    // no more depth to sort?
    if (depth+1 > pmax) return;

    // recurse on buckets
    size_t sum = lo;
    for (size_t i = 0; i < RADIX; ++i)
    {
        if (count[i] <= 1) continue;
        RadixSortMSD(a, sum, sum+count[i], depth+1);
        sum += count[i];
    }
}

void RadixSortMSD(WSortView& a)
{
    return RadixSortMSD(a, 0, a.size(), 0);
}

// ****************************************************************************
// *** Radix Sort (counting sort, least significant digit (LSD) first, out-of-place redistribute)

// by myself (Timo Bingmann)

void RadixSortLSD(WSortView& a)
{
    // radix and base calculations
    const unsigned int RADIX = 4;

    unsigned int pmax = floor( log(a.array_max()) / log(RADIX) );

    for (unsigned int p = 0; p <= pmax; ++p)
    {
        size_t base = pow(RADIX, p);

        // count digits and copy data
        std::vector<size_t> count(RADIX, 0);
        std::vector<value_type> copy(a.size());

        for (size_t i = 0; i < a.size(); ++i)
        {
            size_t r = (copy[i] = a[i]).get() / base % RADIX;
            ASSERT(r < RADIX);
            count[r]++;
        }

        // exclusive prefix sum
        std::vector<size_t> bkt(RADIX+1, 0);
        std::partial_sum(count.begin(), count.end(), bkt.begin()+1);

        // mark bucket boundaries
        for (size_t i = 0; i < bkt.size()-1; ++i) {
            if (bkt[i] >= a.size()) continue;
            a.mark(bkt[i], 2);
        }

        // redistribute items back into array (stable)
        for (size_t i=0; i < a.size(); ++i)
        {
            size_t r = copy[i].get() / base % RADIX;
            a[ bkt[r]++ ] = copy[i];
        }

        a.unmark_all();
    }
}

// ****************************************************************************
// *** Use STL Sorts via Iterator Adapters

// iterator based on http://zotu.blogspot.de/2010/01/creating-random-access-iterator.html

class MyIterator : public std::iterator<std::random_access_iterator_tag, value_type>
{
protected:
    WSortView*  m_array;
    size_t      m_pos;

public:
    typedef typename std::iterator<std::random_access_iterator_tag, value_type> base_type;

    typedef std::random_access_iterator_tag iterator_category;

    typedef typename base_type::value_type value_type;
    typedef typename base_type::difference_type difference_type;
    typedef typename base_type::reference reference;
    typedef typename base_type::pointer pointer;

    MyIterator() : m_array(NULL), m_pos(0) {}

    MyIterator(WSortView* a, size_t p) : m_array(a), m_pos(p) {}

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
    { return (*m_array)[m_pos]; }

    pointer operator->() const
    { return &(*m_array)[m_pos]; }

    reference operator[](const difference_type& n) const
    { return (*m_array)[n]; }

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

void StlSort(WSortView& a)
{
    std::sort(MyIterator(&a,0), MyIterator(&a,a.size()));
}

void StlStableSort(WSortView& a)
{
    std::stable_sort(MyIterator(&a,0), MyIterator(&a,a.size()));
}

void StlHeapSort(WSortView& a)
{
    std::make_heap(MyIterator(&a,0), MyIterator(&a,a.size()));
    std::sort_heap(MyIterator(&a,0), MyIterator(&a,a.size()));
}

// ****************************************************************************
// *** BogoSort and more slow sorts

// by myself (Timo Bingmann)

bool BogoCheckSorted(WSortView& a)
{
    size_t i;
    value_type prev = a[0];
    a.mark(0);
    for (i = 1; i < a.size(); ++i)
    {
        value_type val = a[i];
        if (prev > val) break;
        prev = val;
        a.mark(i);
    }

    if (i == a.size()) {
        // this is amazing.
        return true;
    }

    // unmark
    while (i > 0) a.unmark(i--);
    a.unmark(0);

    return false;
}

void BogoSort(WSortView& a)
{
    // keep a permutation of [0,size)
    std::vector<size_t> perm(a.size());

    for (size_t i = 0; i < a.size(); ++i)
        perm[i] = i;

    while (1)
    {
        // check if array is sorted
        if (BogoCheckSorted(a)) break;

        // pick a random permutation of indexes
        std::random_shuffle(perm.begin(), perm.end());

        // permute array in-place
        std::vector<char> pmark(a.size(), 0);

        for (size_t i = 0; i < a.size(); ++i)
        {
            if (pmark[i]) continue;

            // walk a cycle
            size_t j = i;

            //std::cout << "cycle start " << j << " -> " << perm[j] << "\n";

            while ( perm[j] != i )
            {
                ASSERT(!pmark[j]);
                a.swap(j, perm[j]);
                pmark[j] = 1;

                j = perm[j];
                //std::cout << "cycle step " << j << " -> " << perm[j] << "\n";
            }
            //std::cout << "cycle end\n";

            ASSERT(!pmark[j]);
            pmark[j] = 1;
        }

        //std::cout << "permute end\n";

        for (size_t i = 0; i < a.size(); ++i)
            ASSERT(pmark[i]);
    }
}

void BozoSort(WSortView& a)
{
    srand(time(NULL));

    while (1)
    {
        // check if array is sorted
        if (BogoCheckSorted(a)) break;

        // swap two random items
        a.swap(rand() % a.size(), rand() % a.size());
    }
}

// ****************************************************************************
// *** Bitonic Sort

// from http://www.iti.fh-flensburg.de/lang/algorithmen/sortieren/bitonic/oddn.htm

namespace BitonicSortNS {

static const bool ASCENDING = true;    // sorting direction

static void compare(WSortView& a, int i, int j, bool dir)
{
    if (dir == (a[i] > a[j]))
        a.swap(i, j);
}

static int greatestPowerOfTwoLessThan(int n)
{
    int k = 1;
    while (k < n) k = k << 1;
    return k >> 1;
}

static void bitonicMerge(WSortView& a, int lo, int n, bool dir)
{
    if (n > 1)
    {
        int m = greatestPowerOfTwoLessThan(n);

        for (int i = lo; i < lo + n - m; i++)
            compare(a, i, i+m, dir);

        bitonicMerge(a, lo, m, dir);
        bitonicMerge(a, lo + m, n - m, dir);
    }
}

static void bitonicSort(WSortView& a, int lo, int n, bool dir)
{
    if (n > 1)
    {
        int m = n / 2;
        bitonicSort(a, lo, m, !dir);
        bitonicSort(a, lo + m, n - m, dir);
        bitonicMerge(a, lo, n, dir);
    }
}

} // namespace BitonicSortNS

void BitonicSort(WSortView& a)
{
    BitonicSortNS::bitonicSort(a, 0, a.size(), BitonicSortNS::ASCENDING);
}

// ****************************************************************************
// *** Smooth Sort

// from http://en.wikipedia.org/wiki/Smoothsort

namespace SmoothSortNS {

static const int LP[] = {
    1, 1, 3, 5, 9, 15, 25, 41, 67, 109,
    177, 287, 465, 753, 1219, 1973, 3193, 5167, 8361, 13529, 21891,
    35421, 57313, 92735, 150049, 242785, 392835, 635621, 1028457,
    1664079, 2692537, 4356617, 7049155, 11405773, 18454929, 29860703,
    48315633, 78176337, 126491971, 204668309, 331160281, 535828591,
    866988873 // the next number is > 31 bits.
};

static void sift(WSortView& a, int pshift, int head)
{
    // we do not use Floyd's improvements to the heapsort sift, because we
    // are not doing what heapsort does - always moving nodes from near
    // the bottom of the tree to the root.

    value_type val = a[head];

    while (pshift > 1)
    {
        int rt = head - 1;
        int lf = head - 1 - LP[pshift - 2];

        if (val.compareTo(a[lf]) >= 0 && val.compareTo(a[rt]) >= 0)
            break;

        if (a[lf].compareTo(a[rt]) >= 0) {
            a[head] = a[lf];
            head = lf;
            pshift -= 1;
        }
        else {
            a[head] = a[rt];
            head = rt;
            pshift -= 2;
        }
    }

    a[head] = val;
}

static void trinkle(WSortView& a, int p, int pshift, int head, bool isTrusty)
{
    value_type val = a[head];

    while (p != 1)
    {
        int stepson = head - LP[pshift];

        if (a[stepson].compareTo(val) <= 0)
            break; // current node is greater than head. sift.

        // no need to check this if we know the current node is trusty,
        // because we just checked the head (which is val, in the first
        // iteration)
        if (!isTrusty && pshift > 1) {
            int rt = head - 1;
            int lf = head - 1 - LP[pshift - 2];
            if (a[rt].compareTo(a[stepson]) >= 0
                || a[lf].compareTo(a[stepson]) >= 0)
                break;
        }

        a[head] = a[stepson];

        head = stepson;
        //int trail = Integer.numberOfTrailingZeros(p & ~1);
        int trail = __builtin_ctz(p & ~1);
        p >>= trail;
        pshift += trail;
        isTrusty = false;
    }

    if (!isTrusty) {
        a[head] = val;
        sift(a, pshift, head);
    }
}

void sort(WSortView& a, int lo, int hi)
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
            sift(a, pshift, head);
            p >>= 2;
            pshift += 2;
        }
        else {
            // adding a new block of length 1
            if (LP[pshift - 1] >= hi - head) {
                // this block is its final size.
                trinkle(a, p, pshift, head, false);
            } else {
                // this block will get merged. Just make it trusty.
                sift(a, pshift, head);
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

    trinkle(a, p, pshift, head, false);

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

            trinkle(a, p >> 1, pshift + 1, head - LP[pshift] - 1, true);
            trinkle(a, p, pshift, head - 1, true);
        }

        head--;
    }
}

} // namespace SmoothSortNS

void SmoothSort(WSortView& a)
{
    return SmoothSortNS::sort(a, 0, a.size()-1);
}

// ****************************************************************************
