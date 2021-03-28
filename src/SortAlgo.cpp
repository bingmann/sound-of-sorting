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

typedef ArrayItem value_type;

// inversion count limit for iterator instrumented algorithms
const unsigned int inversion_count_instrumented = 512;

const struct AlgoEntry g_algolist[] =
{
    { _("Selection Sort"), &SelectionSort, UINT_MAX, UINT_MAX,
      wxEmptyString },
    { _("Insertion Sort"), &InsertionSort, UINT_MAX, UINT_MAX,
      wxEmptyString },
    { _("Binary Insertion Sort"), &BinaryInsertionSort, UINT_MAX, UINT_MAX,
      wxEmptyString },
    { _("Adaptive Binary Insertion Sort"), &AdBinaryInsertionSort, UINT_MAX, UINT_MAX,
      wxEmptyString },
    { _("Merge Sort"), &MergeSort, UINT_MAX, UINT_MAX,
      _("Merge sort which merges two sorted sequences into a shadow array, "
        "and then copies it back to the shown array.") },
    { _("Merge Sort (iterative)"), &MergeSortIterative, UINT_MAX, UINT_MAX,
      _("Merge sort variant which iteratively merges "
        "subarrays of sizes of powers of two.") },
    { _("Merge Sort 2"), &MergeSort2, UINT_MAX, UINT_MAX,
      _("Merge sort which uses binary insertion sort up to 256 inversions.") },
    { _("Adaptive Merge Sort"), &AdaptiveMergeSort, UINT_MAX, UINT_MAX,
      _("Uses adaptive binary insertion sort in it with an adaptive limit for insertion runs") },
    { _("Linked Merge Sort"), &LinkedMergeSort, UINT_MAX, UINT_MAX,
      _("Merge sort which relies on links to determine the locations of the items.") },
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
    { _("Bubble Sort"), &BubbleSort, UINT_MAX, UINT_MAX,
      wxEmptyString },
    { _("Bubble Sort 2"), &BubbleSort2, UINT_MAX, UINT_MAX,
      wxEmptyString },
    { _("Cocktail Shaker Sort"), &CocktailShakerSort, UINT_MAX, UINT_MAX,
      wxEmptyString },
    { _("Cashew Sort"), &CashewSort, UINT_MAX, UINT_MAX,
      wxEmptyString },
    { _("Gnome Sort"), &GnomeSort, UINT_MAX, UINT_MAX,
      wxEmptyString },
    { _("Comb Sort"), &CombSort, UINT_MAX, UINT_MAX,
      wxEmptyString },
    { _("Comb Sort 11"), &CombSort11, UINT_MAX, UINT_MAX,
      wxEmptyString },
    { _("Brick Sort"), &BrickSort, UINT_MAX, UINT_MAX,
      wxEmptyString },
    { _("Shell Sort"), &ShellSort, UINT_MAX, UINT_MAX,
      _("Uses 1391376, 463792, 198768, 86961, 33936, 13776, 4592, 1968, 861, 336, 112, 48, 21, 7, 3, 1") },
    { _("Shell Sort 2"), &ShellSort2, UINT_MAX, UINT_MAX,
      _("Uses division by 2 but adds 1 if even") },
    { _("Shell Sort 3"), &ShellSort3, UINT_MAX, UINT_MAX,
      _("Uses division by 3 but to the nearest gap that is true mod 3") },
    { _("Heap Sort"), &HeapSort, UINT_MAX, UINT_MAX,
      wxEmptyString },
    { _("Smooth Sort"), &SmoothSort, UINT_MAX, UINT_MAX,
      wxEmptyString },
    { _("Odd-Even Sort"), &OddEvenSort, UINT_MAX, UINT_MAX,
      wxEmptyString },
    // older sequential implementation, which really makes little sense to do
    { _("Bitonic Sort"), &BitonicSort, UINT_MAX, UINT_MAX,
      wxEmptyString },
    { _("Batcher's Bitonic Sort"), &BitonicSortNetwork, UINT_MAX, UINT_MAX,
      wxEmptyString },
    { _("Batcher's Odd-Even Merge Sort"), &BatcherSortNetwork, UINT_MAX, UINT_MAX,
      wxEmptyString },
    { _("Iterative Bitonic Sort"), &ItBitonic, UINT_MAX, UINT_MAX,
      wxEmptyString },
    { _("Iterative Odd Even Merge Sort"), &ItOddEvenMerge, UINT_MAX, UINT_MAX,
      wxEmptyString },
    { _("Iterative Pairwise Sorting Network"), &ItPairwise, UINT_MAX, UINT_MAX,
      wxEmptyString },
    { _("Cycle Sort"), &CycleSort, 512, UINT_MAX,
      wxEmptyString },
    { _("Radix Sort (LSD)"), &RadixSortLSD, UINT_MAX, UINT_MAX,
      _("Least significant digit radix sort, which copies item into a shadow "
        "array during counting.") },
    { _("Radix Sort (MSD)"), &RadixSortMSD, UINT_MAX, UINT_MAX,
      _("Most significant digit radix sort, which permutes items in-place by walking cycles.") },
    { _("std::sort (gcc)"), &StlSort, UINT_MAX, inversion_count_instrumented,
      wxEmptyString },
    { _("std::stable_sort (gcc)"), &StlStableSort, UINT_MAX, inversion_count_instrumented,
      wxEmptyString },
    { _("std::stable_sort (gcc) in-place"), &StlInPlaceStableSort, UINT_MAX, inversion_count_instrumented,
      wxEmptyString },
    { _("std::sort_heap (gcc)"), &StlHeapSort, UINT_MAX, inversion_count_instrumented,
      wxEmptyString },
    { _("std::__insertion_sort (gcc)"), &StlInsertion, UINT_MAX, inversion_count_instrumented,
      wxEmptyString },
    { _("std::rotate (insertion) (gcc)"), &StlRotateInsert, UINT_MAX, inversion_count_instrumented,
      wxEmptyString },
    { _("Tim Sort"), &TimSort, UINT_MAX, inversion_count_instrumented,
      wxEmptyString },
    { _("Block Merge Sort (WikiSort)"), &WikiSort, UINT_MAX, inversion_count_instrumented,
      _("An O(1) place O(n log n) time stable merge sort.") },
    { _("Block Merge Sort (Grail Sort)"), &GrailSort, UINT_MAX, UINT_MAX,
      wxEmptyString },
    { _("Grail Lazy Stable Sort"), &LazyStableSort, UINT_MAX, UINT_MAX,
      wxEmptyString },
    { _("Grail Rotate Merge Sort"), &RotateMergeSort, UINT_MAX, UINT_MAX,
      wxEmptyString },
    { _("Sugar Coat 5991 sort"), &SugarCoat5991sort, UINT_MAX, UINT_MAX,
      wxEmptyString },
    { _("Bogo Sort"), &BogoSort, 10, UINT_MAX,
      wxEmptyString },
    { _("Bozo Sort"), &BozoSort, 10, UINT_MAX,
      wxEmptyString },
    { _("Zvaray Sort"), &ZvaraySort, 10, UINT_MAX,
      wxEmptyString },
    { _("Bubblegum Hill Sort"), &BubblegumHillSort, 10, UINT_MAX,
      wxEmptyString },
    { _("Bubblegum Hill Sort II"), &BubblegumHillSortII, 256, UINT_MAX,
      wxEmptyString },
    { _("Rainbow Reef Sort"), &RainbowReefSort, 10, UINT_MAX,
      wxEmptyString },
    { _("Stupid Sort"), &StupidSort, 256, UINT_MAX,
      wxEmptyString },
    { _("Stooge Sort"), &StoogeSort, 256, UINT_MAX,
      wxEmptyString },
    { _("Slow Sort"), &SlowSort, 128, UINT_MAX,
      wxEmptyString },
    { _("Parallel Bitonic Sort"), &ParBitonic, UINT_MAX, UINT_MAX,
      wxEmptyString },
    { _("Parallel Odd Even Merge Sort"), &ParOddEvenMerge, UINT_MAX, UINT_MAX,
      wxEmptyString },
    { _("Parallel Pairwise Sorting Network"), &ParPairwise, UINT_MAX, UINT_MAX,
      wxEmptyString },
    { _("Parallel Odd-Even Sort"), &ParOddEvenSort, UINT_MAX, UINT_MAX,
      wxEmptyString },
    { _("Parallel Brick Sort"), &ParBrickSort, UINT_MAX, UINT_MAX,
      wxEmptyString },
};

const size_t g_algolist_size = sizeof(g_algolist) / sizeof(g_algolist[0]);

const struct AlgoEntry* g_algolist_end = g_algolist + g_algolist_size;

// ****************************************************************************
// *** Selection Sort

void SelectionSort(SortArray& A)
{
    volatile ssize_t jMin = 0;
    A.watch(&jMin, 3);

    for (size_t i = 0; i < A.size()-1; ++i)
    {
        jMin = i;

        for (size_t j = i+1; j < A.size(); ++j)
        {
            if (A[j] < A[jMin]) {
                //A.mark_swap(j, jMin);
                jMin = j;
            }
        }

        A.swap(i, jMin);

        // mark the last good element
        if (i > 0) A.unmark(i-1);
        A.mark(i);
    }
    A.unwatch_all();
}

// ****************************************************************************
// *** Insertion Sort

// swaps every time (keeps all values visible)
void InsertionSort(SortArray& A)
{
    for (size_t i = 1; i < A.size(); ++i)
    {
        value_type key = A[i];
        A.mark(i);

        ssize_t j = i - 1;
        while (j >= 0 && A[j] > key)
        {
            A.swap(j, j+1);
            j--;
        }

        A.unmark(i);
    }
}

// with extra item on stack
void InsertionSort2(SortArray& A)
{
    for (size_t i = 1; i < A.size(); ++i)
    {
        value_type tmp, key = A[i];
        A.mark(i);

        ssize_t j = i - 1;
        while (j >= 0 && (tmp = A[j]) > key)
        {
            A.set(j + 1, tmp);
            j--;
        }
        A.set(j + 1, key);

        A.unmark(i);
    }
}

// swaps every time (keeps all values visible)
void BinaryInsertionSort(SortArray& A)
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
            A.swap(j, j+1);
            j--;
        }

        A.unmark(i);
    }
}

void AdBinaryInsertionSort(SortArray& a)
{
	int length = a.size();
        int count = 0;
        int start = 0;
        int start2 = start;
        int end = length;
        int flag = 0;
        count = 0;
        flag = 0;
        for (int i = start + 1; i < end; i++) {
            int v = (2*count/(i-start2))+1;
            int lo = std::max(i-v, start2), hi = i;
            while((lo>=start2)&&((a[i]<a[lo]))){
		lo-=v; hi-=v;
            }
            lo++;
            if(lo<start2)lo=start2;
            while (lo < hi) {
                int mid = lo + ((hi - lo) / 2); // avoid int overflow!
                if ((a[mid] > a[i])) { // do NOT move equal elements to right of inserted element; this maintains stability!
                    hi = mid;
                }
                else {
                    lo = mid + 1;
                }
            }
// item has to go into position lo
            count += (i - lo);
            int j = i - 1;

            while (j >= lo)
            {
			a.swap(j+1, j);
                j--;
            }
        }
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

        out[o++] = (ai <= aj ? (++i, ai) : (++j, aj));
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

void LinkedMerge(SortArray& A, size_t lo, size_t mid, size_t hi, std::vector<size_t>& buffer)
{
    // mark merge boundaries
    A.mark(lo);
    A.mark(mid,3);
    A.mark(hi-1);

    // allocate output
    //std::vector<value_type> out(hi-lo);
    for(int i=0; i<mid-lo; i++){
	buffer[i]=i; // what merge location each array item is
	buffer[i+(mid-lo)]=i; // where in the array each merge item is
    }

    // merge
    size_t i = lo, j = mid, o = 0; // first and second halves
    while (i < mid && j < hi)
    {
        // copy out for fewer time steps
        value_type ai = A[buffer[i-lo+(mid-lo)]+lo], aj = A[j];
	if(ai<=aj){
		A.swap(o+lo, buffer[i-lo+(mid-lo)]+lo);
		buffer[buffer[i-lo+(mid-lo)]%(mid-lo)]=buffer[(o)%(mid-lo)];
		buffer[buffer[(o)%(mid-lo)]+(mid-lo)]=buffer[i-lo+(mid-lo)];
		++i;
	}
	else{
		A.swap(o+lo, j);
		buffer[(j-lo)%(mid-lo)]=buffer[(o)%(mid-lo)];
		buffer[buffer[(o)%(mid-lo)]+(mid-lo)]=j-lo;
		++j;

	}
	o++;
    }

    // copy rest
    while (i < mid){
		A.swap(o+lo, buffer[i-lo+(mid-lo)]+lo);
		buffer[buffer[i-lo+(mid-lo)]%(mid-lo)]=buffer[(o)%(mid-lo)];
		buffer[buffer[(o)%(mid-lo)]+(mid-lo)]=buffer[i-lo+(mid-lo)];
		++i;
		++o;
    }
    //while (j < hi) out[o++] = A[j++];

    //ASSERT(o == hi-lo);

    A.unmark(mid);

    // copy back
    /*for (i = 0; i < hi-lo; ++i)
        A.set(lo + i, out[i]);*/

    A.unmark(lo);
    A.unmark(hi-1);
}

void LinkedMergeSort(SortArray& A, size_t lo, size_t hi, std::vector<size_t>& buffer)
{
    if (lo + 1 < hi)
    {
        size_t mid = (lo + hi) / 2;

        LinkedMergeSort(A, lo, mid, buffer);
        LinkedMergeSort(A, mid, hi, buffer);

        LinkedMerge(A, lo, mid, hi, buffer);
    }
}

void LinkedMergeSort(SortArray& A)
{
	std::vector<size_t> buffer(A.size());
    return LinkedMergeSort(A, 0, A.size(), buffer);
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

uint8_t numberOfTrailingZeros(uint64_t a){
	uint8_t i;
	for(i=0; i<64; i++){
		if(a%2){
			return i;
		}
		a /= 2;
	}
	return i;
}

void merge2(SortArray& array, value_type* merge, int min1, int max1, int min2, int max2) {
      if ((min2 - min1) <= (max2 - max1)){
        for (int g = 0; g < (min2 - min1); g++){
		merge[g]=array[max1-g];
        }
        int point1 = min2;
        int point2 = (min2 - min1) - 1;
        int point3 = min1;
        while ((point1 <= max2) && (point2 >= 0)){
            if (array[point1] < merge[point2]){
		array.set(point3, array[point1]);
                point3++;
                point1++;
            }
            else{
		array.set(point3, merge[point2]);
                point2--;
                point3++;
            }
        }
        while ((point2 >= 0)){
		array.set(point3, merge[point2]);
                point2--;
                point3++;
        }
      }
      else {
        for (int g = 0; g < (max2 - max1); g++){
		merge[g]=array[min2+g];
        }
        int point1 = max1;
        int point2 = (max2 - max1) - 1;
        int point3 = max2;
        while ((point1 >= min1) && (point2 >= 0)){
            if (array[point1] > merge[point2]){
		array.set(point3, array[point1]);
                point3--;
                point1--;
            }
            else{
		array.set(point3, merge[point2]);
                point2--;
                point3--;
            }
        }
        while ((point2 >= 0)){
		array.set(point3, merge[point2]);
                point2--;
                point3--;
        }
      }
    }

void MergeSort2(SortArray& a)
{
	int length = a.size();
	int limit = 256;
        int sequence = 0;
        int* stack = (int*)malloc((int)floor(log(length+0.0)/log(2.0)+2.0)*sizeof(int));
        stack[0]=0;
        int stackend = 0;
        value_type* merge = (value_type*)malloc(length/2*sizeof(value_type));
        int count = 0;
        int i = 0;
        int start = 0;
        int test = start;
        int start2 = start;
        int end = length;
        int flag = 0;
        while (test < (end - 1)) {
        count = 0;
        i = test + 1;
        start2 = test;
        flag = 0;
        while (i < end && flag == 0) {
            int lo = start2, hi = i;

            while (lo < hi) {
                int mid = lo + ((hi - lo) / 2); // avoid int overflow!
                if (a[i] < a[mid]) { // do NOT move equal elements to right of inserted element; this maintains stability!
                    hi = mid;
                }
                else {
                    lo = mid + 1;
                }
            }

            // item has to go into position lo
            count += (i - lo);
            if (count > limit){
            flag = 1;
            }
            else {
            if (i > lo){
            int j = i - 1;

            while (j >= lo)
            {
			a.swap(j+1, j);
                j--;
            }

            }
            i++;
            }
            test = i;
        }
        sequence++;
        stackend++;
        stack[stackend] = test;
        for (int r = 0; r < numberOfTrailingZeros(sequence); r++){
            merge2(a, merge, stack[stackend - 2], stack[stackend - 1] - 1, stack[stackend - 1], stack[stackend] - 1);
            stackend--;
            stack[stackend] = stack[stackend+1];
        }
        }
        if (stack[stackend] == (end - 1)){
            stackend++;
            stack[stackend]=end;
        }
        while (stackend > 1){
            merge2(a, merge, stack[stackend - 2], stack[stackend - 1] - 1, stack[stackend - 1], stack[stackend] - 1);
            stackend--;
            stack[stackend] = stack[stackend+1];
        }
        free(stack);
}

void AdaptiveMergeSort(SortArray& a)
{
	int length = a.size();
        int sequence = 0;
        int* stack = (int*)malloc((int)floor(log(length+0.0)/log(2.0)+2.0)*sizeof(int));
        stack[0]=0;
        int stackend = 0;
        value_type* merge = (value_type*)malloc(length/2*sizeof(value_type));
        int count = 0;
        int i = 0;
        int start = 0;
        int test = start;
        int start2 = start;
        int end = length;
        int flag = 0;
        while (test < (end - 1)) {
        count = 0;
        i = test + 1;
        start2 = test;
        flag = 0;
        bool dir = 0;
        while (i < end && flag == 0) {
	if((i-start2)==16&&count>90){
		dir=1;
		for(int q=0; q<(i-start2)/2; q++){
			a.swap(start2+q, i-1-q);
		}
		count = 120-count;
	}
            int v = (2*count/(i-start2))+1;
            int lo = std::max(i-v, start2), hi = i;
            while((lo>=start2)&&((a[i]<a[lo])^dir)){
		lo-=v; hi-=v;
            }
            lo++;
            if(lo<start2)lo=start2;
            while (lo < hi) {
                int mid = lo + ((hi - lo) / 2); // avoid int overflow!
                if ((a[mid] > a[i])^dir) { // do NOT move equal elements to right of inserted element; this maintains stability!
                    hi = mid;
                }
                else {
                    lo = mid + 1;
                }
            }
            ArrayItem num = a[i];
	int limit = (int)((i-start2)*(log2((double)(i-start2)))*1.5);
            // item has to go into position lo
            count += (i - lo);
            if (count > limit && (i-start2)>=16){
            flag = 1;
            }
            else {
            if (i > lo){
            int j = i - 1;

            while (j >= lo)
            {
			a.swap(j+1, j);
                j--;
            }

            }
            i++;
            }
            test = i;
        }
        if(dir){
		for(int q=0; q<(test-start2)/2; q++){
			a.swap(start2+q, test-1-q);
		}
        }
        sequence++;
        stackend++;
        stack[stackend] = test;
        for (int r = 0; r < numberOfTrailingZeros(sequence); r++){
            merge2(a, merge, stack[stackend - 2], stack[stackend - 1] - 1, stack[stackend - 1], stack[stackend] - 1);
            stackend--;
            stack[stackend] = stack[stackend+1];
        }
        }
        if (stack[stackend] == (end - 1)){
            stackend++;
            stack[stackend]=end;
        }
        while (stackend > 1){
            merge2(a, merge, stack[stackend - 2], stack[stackend - 1] - 1, stack[stackend - 1], stack[stackend] - 1);
            stackend--;
            stack[stackend] = stack[stackend+1];
        }
        free(stack);
}

// ****************************************************************************
// *** Quick Sort Pivot Selection

QuickSortPivotType g_quicksort_pivot = PIVOT_FIRST;
QuickSortDualPivotType g_quicksort_dualpivot = DUALPIVOT_FIRSTLAST;

// some quicksort variants use hi inclusive and some exclusive, we require it
// to be _exclusive_. hi == array.end()!
ssize_t QuickSortSelectPivot(SortArray& A, ssize_t lo, ssize_t hi)
{
    if (g_quicksort_pivot == PIVOT_FIRST)
        return lo;

    if (g_quicksort_pivot == PIVOT_LAST)
        return hi-1;

    if (g_quicksort_pivot == PIVOT_MID)
        return (lo + hi) / 2;

    if (g_quicksort_pivot == PIVOT_RANDOM)
        return lo + (next() % (hi - lo));

    if (g_quicksort_pivot == PIVOT_MEDIAN3)
    {
        ssize_t mid = (lo + hi) / 2;
    	if(hi-lo<=2)return lo; // skip all comparisons if less than 3 items

        // cases if two are equal - skip to make 3 comparisons worst case
        /*if (A[lo] == A[mid]) return lo;
        if (A[lo] == A[hi-1] || A[mid] == A[hi-1]) return hi-1;*/

        // cases if three are different
        return A[lo] < A[mid]
            ? (A[mid] < A[hi-1] ? mid : (A[lo] < A[hi-1] ? hi-1 : lo))
            : (A[mid] > A[hi-1] ? mid : (A[lo] < A[hi-1] ? lo : hi-1));
    }

    if (g_quicksort_pivot == PIVOT_NINTHER)
    {
        ssize_t mid = (lo + hi) / 2;
    	if(hi-lo<=2)return lo;
    	if(hi-lo<=8) // median of 3 if less than 9 items
        return A[lo] < A[mid]
            ? (A[mid] < A[hi-1] ? mid : (A[lo] < A[hi-1] ? hi-1 : lo))
            : (A[mid] > A[hi-1] ? mid : (A[lo] < A[hi-1] ? lo : hi-1));
        ssize_t a = (lo);
        ssize_t c = (lo + lo + (hi-1) + 1) / 3;
        ssize_t g = (hi) - (c - a);
        ssize_t k = (hi);
        ssize_t b = (a + c) / 2;
        ssize_t h = (g + k) / 2;
        ssize_t d = (c) ;
        ssize_t f = (g) ;
        ssize_t e = (d + f) / 2;

        ssize_t x;
        ssize_t y;
        ssize_t z;

        x = A[a] < A[b]
            ? (A[b] < A[c-1] ? b : (A[a] < A[c-1] ? c-1 : a))
            : (A[b] > A[c-1] ? b : (A[a] < A[c-1] ? a : c-1));

        y = A[d] < A[e]
            ? (A[e] < A[f-1] ? e : (A[d] < A[f-1] ? f-1 : d))
            : (A[e] > A[f-1] ? e : (A[d] < A[f-1] ? d : f-1));

        z = A[g] < A[h]
            ? (A[h] < A[k-1] ? h : (A[g] < A[k-1] ? k-1 : g))
            : (A[h] > A[k-1] ? h : (A[g] < A[k-1] ? g : k-1));

        return A[x] < A[y]
            ? (A[y] < A[z-1] ? y : (A[x] < A[z-1] ? z-1 : x))
            : (A[y] > A[z-1] ? y : (A[x] < A[z-1] ? x : z-1));
    }

    return lo;
}

ssize_t* QuickSortSelectDualPivot(SortArray& A, ssize_t lo, ssize_t hi, ssize_t* pivot)
{
    if (g_quicksort_dualpivot == DUALPIVOT_FIRSTLAST){
        pivot[0]=lo; pivot[1]=hi-1; return pivot;
    }
    if (g_quicksort_dualpivot == DUALPIVOT_THIRDS){
        pivot[0]=(lo+lo+(hi-1)+1)/3; pivot[1]=(lo+(hi-1)+(hi-1)+1)/3; return pivot;
    }
    if (g_quicksort_dualpivot == DUALPIVOT_RANDOM){
        pivot[0]=next()%(hi-lo)+lo; pivot[1]=next()%(hi-lo-1)+lo; if(pivot[1]==pivot[0])pivot[1]=hi-1; return pivot;
    }
    if (g_quicksort_dualpivot == DUALPIVOT_MEDIAN4){
		hi-=1;
    	ssize_t temp;
        ssize_t mid1 = (lo+lo+hi+1)/3;
        ssize_t mid2 = (lo+hi+hi+1)/3;
    	if(hi-lo<=2){pivot[0]=lo; pivot[1]=hi; return pivot;} // skip all comparisons if less than 4 items

    	if(A[lo] > A[mid1]){temp=lo;lo=mid1;mid1=temp;}
    	if(A[mid2] > A[hi]){temp=mid2;mid2=hi;hi=temp;}
    	if(A[lo] > A[mid2]){temp=lo;lo=mid2;mid2=temp;}
    	if(A[mid1] > A[hi]){temp=mid1;mid1=hi;hi=temp;}
    	pivot[0]=mid1; pivot[1]=mid2; return pivot;
    }
    pivot[0]=lo; pivot[1]=hi-1; return pivot;
}

wxArrayString QuickSortPivotText()
{
    wxArrayString sl;

    sl.Add( _("First Item") );
    sl.Add( _("Last Item") );
    sl.Add( _("Middle Item") );
    sl.Add( _("Random Item") );
    sl.Add( _("Median of Three") );
    sl.Add( _("Ninther") );

    return sl;
}

wxArrayString QuickSortDualPivotText()
{
    wxArrayString sl;

    sl.Add( _("First and Last Item") );
    sl.Add( _("Thirds Items") );
    sl.Add( _("Random Items") );
    sl.Add( _("Medians of Four") );

    return sl;
}

// ****************************************************************************
// *** Quick Sort LR (in-place, pointers at left and right, pivot is middle element)

// by myself (Timo Bingmann), based on Hoare's original code

void QuickSortLR(SortArray& A, ssize_t lo, ssize_t hi)
{
    // pick pivot and watch
    volatile ssize_t p = QuickSortSelectPivot(A, lo, hi+1);

    value_type pivot = A[p];
    A.watch(&p, 2);

    volatile ssize_t i = lo, j = hi;
    A.watch(&i, 3);
    A.watch(&j, 3);

    while (i <= j)
    {
        while (A[i] < pivot)
            i++;

        while (A[j] > pivot)
            j--;

        if (i <= j)
        {
            A.swap(i,j);

            // follow pivot if it is swapped
            if (p == i) p = j;
            else if (p == j) p = i;

            i++, j--;
        }
    }

    A.unwatch_all();

    if (lo < j)
        QuickSortLR(A, lo, j);

    if (i < hi)
        QuickSortLR(A, i, hi);
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
    A.swap(p, hi-1);
    A.mark(hi-1);

    volatile ssize_t i = lo;
    A.watch(&i, 3);

    for (size_t j = lo; j < hi-1; ++j)
    {
        if (A[j] <= pivot) {
            A.swap(i, j);
            ++i;
        }
    }

    A.swap(i, hi-1);
    A.unmark(hi-1);
    A.unwatch_all();

    return i;
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
    ssize_t piv = QuickSortSelectPivot(A, lo, hi+1);
    A.swap(piv, hi);
    A.mark(hi);

    const value_type& pivot = A[hi];

    // schema: |p ===  |i <<< | ??? |j >>> |q === |piv
    volatile ssize_t i = lo, j = hi-1;
    volatile ssize_t p = lo, q = hi-1;

    A.watch(&i, 3);
    A.watch(&j, 3);

    for (;;)
    {
        // partition on left
        while (i <= j && (cmp = A[i].cmp(pivot)) <= 0)
        {
            if (cmp == 0) {
                A.mark(p,4);
                A.swap(i, p++);
            }
            ++i;
        }

        // partition on right
        while (i <= j && (cmp = A[j].cmp(pivot)) >= 0)
        {
            if (cmp == 0) {
                A.mark(q,4);
                A.swap(j, q--);
            }
            --j;
        }

        if (i > j) break;

        // swap item between < > regions
        A.swap(i++, j--);
    }

    // swap pivot to right place
    A.swap(i,hi);
    A.mark_swap(i,hi);

    ssize_t num_less = i - p;
    ssize_t num_greater = q - j;

    // swap equal ranges into center, but avoid swapping equal elements
    j = i-1; i = i+1;

    ssize_t pe = lo + std::min(p-lo, num_less);
    for (ssize_t k = lo; k < pe; k++, j--) {
        A.swap(k,j);
        A.mark_swap(k,j);
    }

    ssize_t qe = hi-1 - std::min(hi-1-q, num_greater-1); // one already greater at end
    for (ssize_t k = hi-1; k > qe; k--, i++) {
        A.swap(i,k);
        A.mark_swap(i,k);
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

std::pair<ssize_t,ssize_t> PartitionTernaryLL(SortArray& A, ssize_t lo, ssize_t hi)
{
    // pick pivot and swap to back
    ssize_t p = QuickSortSelectPivot(A, lo, hi);

    value_type pivot = A[p];
    A.swap(p, hi-1);
    A.mark(hi-1);

    volatile ssize_t i = lo, k = hi-1;
    A.watch(&i, 3);

    for (ssize_t j = lo; j < k; ++j)
    {
        int cmp = A[j].cmp(pivot); // ternary comparison
        if (cmp == 0) {
            A.swap(--k, j);
            --j; // reclassify A[j]
            A.mark(k,4);
        }
        else if (cmp < 0) {
            A.swap(i++, j);
        }
    }

    // unwatch i, because the pivot is swapped there
    // in the first step of the following swap loop.
    A.unwatch_all();

    ssize_t j = i + (hi-k);

    for (ssize_t s = 0; s < hi-k; ++s) {
        A.swap(i+s, hi-1-s);
        A.mark_swap(i+s, hi-1-s);
    }
    A.unmark_all();

    return std::make_pair(i,j);
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
    	ssize_t pivot[2];
    	QuickSortSelectDualPivot(a, left, right+1, pivot);
    	if(pivot[1]<pivot[0]){ssize_t a = pivot[1];pivot[1]=pivot[0];pivot[0]=a;}
    	if(pivot[0]!=left){
		a.swap(left, pivot[0]);
    	}
    	if(pivot[1]!=right){
		a.swap(right, pivot[1]);
    	}
        if (a[left] > a[right]) {
            a.swap(left, right);
        }

        const value_type p = a[left];
        const value_type q = a[right];

        a.mark(left);
        a.mark(right);

        volatile ssize_t l = left + 1;
        volatile ssize_t g = right - 1;
        volatile ssize_t k = l;

        a.watch(&l, 3);
        a.watch(&g, 3);
        a.watch(&k, 3);

        while (k <= g)
        {
            if (a[k] < p) {
                a.swap(k, l);
                ++l;
            }
            else if (a[k] >= q) {
                while (a[g] > q && k < g)  --g;
                a.swap(k, g);
                --g;

                if (a[k] < p) {
                    a.swap(k, l);
                    ++l;
                }
            }
            ++k;
        }
        --l;
        ++g;
        a.swap(left, l);
        a.swap(right, g);

        a.unmark_all();
        a.unwatch_all();

        dualPivotYaroslavskiy(a, left, l - 1);
        dualPivotYaroslavskiy(a, l + 1, g - 1);
        dualPivotYaroslavskiy(a, g + 1, right);
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
            if (A[j + 1] < A[j])
            {
                A.swap(j, j+1);
            }
        }
    }
}

void BubbleSort2(SortArray& A)
{
    size_t lo = 0, hi = A.size()-1, mov = lo;

    while (lo < hi)
    {
    	mov = lo;
        for (size_t i = lo; i < hi; ++i)
        {
            if (A[i+1] < A[i])
            {
                A.swap(i, i+1);
                mov = i;
            }
        }

        hi = mov;
    }
}

// ****************************************************************************
// *** Cocktail Shaker Sort

// from http://de.wikibooks.org/wiki/Algorithmen_und_Datenstrukturen_in_C/_Shakersort

void CocktailShakerSort(SortArray& A)
{
    size_t lo = 0, hi = A.size()-1, mov = hi;

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
            if (A[i+1] < A[i])
            {
                A.swap(i, i+1);
                mov = i;
            }
        }

        hi = mov;
    }
}

void CashewSort(SortArray& A)
{
    size_t lo = 0, hi = A.size()-1, mov = lo;

    while (lo < hi)
    {
        for (size_t i = lo; i < hi; ++i)
        {
            if (A[i+1] < A[i])
            {
                A.swap(i, i+1);
                mov = i;
            }
        }

        hi = mov;

        for (size_t i = lo; i < hi; ++i)
        {
            if (A[i+1] < A[i])
            {
                A.swap(i, i+1);
                mov = i;
            }
        }

        hi = mov;

        for (size_t i = hi; i > lo; --i)
        {
            if (A[i-1] > A[i])
            {
                A.swap(i-1, i);
                mov = i;
            }
        }

        lo = mov;

        for (size_t i = hi; i > lo; --i)
        {
            if (A[i-1] > A[i])
            {
                A.swap(i-1, i);
                mov = i;
            }
        }

        lo = mov;
    }
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

void StupidSort(SortArray& A)
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
            i = 1;
        }
    }
}

// ****************************************************************************
// *** Comb Sort

// from http://en.wikipediA.org/wiki/Comb_sort

void CombSort(SortArray& A)
{
    const unsigned long long shrinknumr = 13;
    const unsigned long long shrinkdnom = 10;

    bool swapped = false;
    size_t gap = A.size();

    while ((gap > 1) || swapped)
    {
        if (gap > 1) {
            gap = (size_t)(gap * shrinkdnom / shrinknumr);
        }

        swapped = false;

        for (size_t i = 0; gap + i < A.size(); ++i)
        {
            if (A[i + gap] < A[i])
            {
                A.swap(i, i+gap);
                swapped = true;
            }
        }
    }
}

void CombSort11(SortArray& A)
{
    const unsigned long long shrinknumr = 13;
    const unsigned long long shrinkdnom = 10;

    bool swapped = false;
    size_t gap = A.size();

    while ((gap > 1) || swapped)
    {
        if (gap > 1) {
            gap = (size_t)(gap * shrinkdnom / shrinknumr);
        }
        if(gap==9||gap==10)gap=11;

        swapped = false;

        for (size_t i = 0; gap + i < A.size(); ++i)
        {
            if (A[i + gap] < A[i])
            {
                A.swap(i, i+gap);
                swapped = true;
            }
        }
    }
}

void BrickSort(SortArray& A)
{
    const unsigned long long shrinknumr = 6;
    const unsigned long long shrinkdnom = 5;

    bool swapped = false;
    size_t gap = A.size();

    while ((gap > 1) || swapped)
    {
        if (gap > 1) {
            gap = (size_t)(gap * shrinkdnom / shrinknumr);
        }

        swapped = false;

        for (size_t i = 0; gap + i < A.size(); ++i)
        {
            if (!((i/gap)&1) && A[i + gap] < A[i])
            {
                A.swap(i, i+gap);
                swapped = true;
            }
        }
        for (size_t i = 0; gap + i < A.size(); ++i)
        {
            if ( ((i/gap)&1) && A[i + gap] < A[i])
            {
                A.swap(i, i+gap);
                swapped = true;
            }
        }
    }
}

void ParBrickSort(SortArray& A)
{
    const unsigned long long shrinknumr = 6;
    const unsigned long long shrinkdnom = 5;

    bool swapped = false;
    size_t gap = A.size();

    while ((gap > 1) || swapped)
    {
        if (gap > 1) {
            gap = (size_t)(gap * shrinkdnom / shrinknumr);
        }

        swapped = false;

        #pragma omp parallel for
        for (size_t i = 0; i < A.size()-gap; ++i)
        {
            if (!((i/gap)&1) && A[i + gap] < A[i])
            {
                A.swap(i, i+gap);
                swapped = true;
            }
        }
        #pragma omp parallel for
        for (size_t i = 0; i < A.size()-gap; ++i)
        {
            if ( ((i/gap)&1) && A[i + gap] < A[i])
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
//#pragma omp parallel for
        for (size_t i = 1; i < A.size()-1; i += 2)
        {
            if(A[i] > A[i+1])
            {
                A.swap(i, i+1);
                sorted = false;
            }
        }

//#pragma omp parallel for
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

void ParOddEvenSort(SortArray& A)
{
    bool sorted = false;

    while (!sorted)
    {
        sorted = true;
#pragma omp parallel for
        for (size_t i = 1; i < A.size()-1; i += 2)
        {
            if(A[i] > A[i+1])
            {
                A.swap(i, i+1);
                sorted = false;
            }
        }

#pragma omp parallel for
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
   /* for (size_t j = 0; j < A.size(); ++j)
       {A.mark(j, log((j+1)) / log(2) + 4);
       	A.mark(j, 0);
       }*/
    size_t incs[16] = { 1391376, 463792, 198768, 86961, 33936,
                        13776, 4592, 1968, 861, 336,
                        112, 48, 21, 7, 3, 1 };

    for (size_t k = 0; k < 16; k++)
    {
        for (size_t h = incs[k], i = h; i < A.size(); i++)
        {
            size_t j = i;

            while (j >= h && A[j-h] > A[j])
            {
                A.swap(j, j-h);
                j -= h;
            }
        }
    }
}

void ShellSort2(SortArray& A)
{
	size_t gap = A.size();

    while(gap>1)
    {
    	gap = (gap>>1)|1;
        for (size_t h = gap, i = h; i < A.size(); i++)
        {
            size_t j = i;

            while (j >= h && A[j-h] > A[j])
            {
                A.swap(j, j-h);
                j -= h;
            }
        }
    }
}

void ShellSort3(SortArray& A)
{
	size_t gap = A.size();

    while(gap>1)
    {
    	gap = (gap/9*3)+(gap%9/5+1);
        for (size_t h = gap, i = h; i < A.size(); i++)
        {
            size_t j = i;

            while (j >= h && A[j-h] > A[j])
            {
                A.swap(j, j-h);
                j -= h;
            }
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
            A.mark(i, log(prevPowerOfTwo(i+1)) / log(2) + 4);
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
    }

}

// ****************************************************************************
// *** Radix Sort (counting sort, most significant digit (MSD) first, in-place redistribute)

// by myself (Timo Bingmann)

void RadixSortMSD(SortArray& A, size_t lo, size_t hi, size_t depth)
{
    A.mark(lo); A.mark(hi-1);

    // radix and base calculations
    const unsigned int RADIX = 4;

    unsigned int pmax = floor( log(A.array_max()+1) / log(RADIX) );
    ASSERT(depth <= pmax);

    size_t base = /*pow(RADIX, pmax - depth)*/1<<((pmax-depth)*2);

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
            RadixSortMSD(A, sum, sum+count[i], depth+1);
        sum += count[i];
    }
}

void RadixSortMSD(SortArray& A)
{
    return RadixSortMSD(A, 0, A.size(), 0);
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
        size_t base = /*pow(RADIX, p);*/1<<(p*2);

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
// *** Use STL Sorts via Iterator Adapters

void StlSort(SortArray& A)
{
    std::sort(MyIterator(&A,0), MyIterator(&A,A.size()));
}

void StlStableSort(SortArray& A)
{
    std::stable_sort(MyIterator(&A,0), MyIterator(&A,A.size()));
}

void StlInPlaceStableSort(SortArray& A)
{
      // concept requirements
      __glibcxx_function_requires(_Mutable_RandomAccessIteratorConcept<
	    _RandomAccessIterator>)
      __glibcxx_function_requires(_LessThanComparableConcept<
	    typename iterator_traits<_RandomAccessIterator>::value_type>)
      __glibcxx_requires_valid_range(MyIterator(&A,0), MyIterator(&A,A.size()));
      __glibcxx_requires_irreflexive(MyIterator(&A,0), MyIterator(&A,A.size()));

    std::__inplace_stable_sort(MyIterator(&A,0), MyIterator(&A,A.size()), __gnu_cxx::__ops::__iter_less_iter());
}

void StlHeapSort(SortArray& A)
{
    std::make_heap(MyIterator(&A,0), MyIterator(&A,A.size()));
    std::sort_heap(MyIterator(&A,0), MyIterator(&A,A.size()));
}

void StlInsertion(SortArray& A){
      __glibcxx_function_requires(_Mutable_RandomAccessIteratorConcept<
	    _RandomAccessIterator>)
      __glibcxx_function_requires(_LessThanComparableConcept<
	    typename iterator_traits<_RandomAccessIterator>::value_type>)
      __glibcxx_requires_valid_range(MyIterator(&A,0), MyIterator(&A,A.size()));
      __glibcxx_requires_irreflexive(MyIterator(&A,0), MyIterator(&A,A.size()));
    std::__insertion_sort(MyIterator(&A,0), MyIterator(&A,A.size()), __gnu_cxx::__ops::__iter_less_iter());
}

void StlRotateInsert(SortArray& A){      // concept requirements
      __glibcxx_function_requires(_Mutable_RandomAccessIteratorConcept<
	    _RandomAccessIterator>)
      __glibcxx_function_requires(_LessThanComparableConcept<
	    typename iterator_traits<_RandomAccessIterator>::value_type>)
      __glibcxx_requires_valid_range(MyIterator(&A,0), MyIterator(&A,A.size()));
      __glibcxx_requires_irreflexive(MyIterator(&A,0), MyIterator(&A,A.size()));

    for (MyIterator it = MyIterator(&A,0); it != MyIterator(&A,A.size()); ++it) {
        std::rotate(std::upper_bound(MyIterator(&A,0), it, *it, __gnu_cxx::__ops::__iter_less_iter()), it, it+1);
    }
}



void SugarCoat5991sort(SortArray& A, size_t lo, size_t hi, std::vector<value_type>& buffer, bool mode, bool dir)
{
	if(mode){
	if(hi-lo<=0)return;
	if(hi-lo==1){A.set(lo, buffer[lo]);return;}
		if(hi-lo==2){
		if(dir ? buffer[lo]>=buffer[lo+1] : buffer[lo]>buffer[lo+1]){A.set(lo, buffer[lo+1]);A.set(lo+1, buffer[lo]);return;}
		else {A.set(lo, buffer[lo]);A.set(lo+1, buffer[lo+1]);}
				return;}
		if(hi-lo==3){
				if(dir){
			if(buffer[lo]>=buffer[lo+1]){
				if(buffer[lo+1]>=buffer[lo+2]){A.set(lo, buffer[lo+2]);A.set(lo+1, buffer[lo+1]);A.set(lo+2, buffer[lo]);}
				else{
					if(buffer[lo]>=buffer[lo+2]){A.set(lo, buffer[lo+1]);A.set(lo+1, buffer[lo+2]);A.set(lo+2, buffer[lo]);}
					else {A.set(lo, buffer[lo+1]);A.set(lo+1, buffer[lo]);A.set(lo+2, buffer[lo+2]);};}
			}
			else{
				if(buffer[lo+1]>=buffer[lo+2])
					if(buffer[lo]>=buffer[lo+2]){A.set(lo, buffer[lo+2]);A.set(lo+1, buffer[lo]);A.set(lo+2, buffer[lo+1]);}
					else {A.set(lo, buffer[lo]);A.set(lo+1, buffer[lo+2]);A.set(lo+2, buffer[lo+1]);}
				else {A.set(lo, buffer[lo]);A.set(lo+1, buffer[lo+1]);A.set(lo+2, buffer[lo+2]);} return;
			}
				}
				else{
			if(buffer[lo]>buffer[lo+1]){
				if(buffer[lo+1]>buffer[lo+2]){A.set(lo, buffer[lo+2]);A.set(lo+1, buffer[lo+1]);A.set(lo+2, buffer[lo]);}
				else{
					if(buffer[lo]>buffer[lo+2]){A.set(lo, buffer[lo+1]);A.set(lo+1, buffer[lo+2]);A.set(lo+2, buffer[lo]);}
					else {A.set(lo, buffer[lo+1]);A.set(lo+1, buffer[lo]);A.set(lo+2, buffer[lo+2]);};}
			}
			else{
				if(buffer[lo+1]>buffer[lo+2])
					if(buffer[lo]>buffer[lo+2]){A.set(lo, buffer[lo+2]);A.set(lo+1, buffer[lo]);A.set(lo+2, buffer[lo+1]);}
					else {A.set(lo, buffer[lo]);A.set(lo+1, buffer[lo+2]);A.set(lo+2, buffer[lo+1]);}
				else {A.set(lo, buffer[lo]);A.set(lo+1, buffer[lo+1]);A.set(lo+2, buffer[lo+2]);} return;
			}
				}
				return;
		}
		if(hi-lo==4){
				std::array<size_t, 4> q;
				if(dir) q = {lo+3, lo+2, lo+1, lo+0} ; else q = {lo, lo+1, lo+2, lo+3};
				size_t temp;
				if(buffer[q[0]]>buffer[q[1]]){temp=q[0];q[0]=q[1];q[1]=temp;}
				if(buffer[q[1]]>buffer[q[2]]){temp=q[1];q[1]=q[2];q[2]=temp;}
				if(q[2]!=lo+(2-dir)&&buffer[q[0]]>buffer[q[1]]){temp=q[0];q[0]=q[1];q[1]=temp;}
				if(buffer[q[1]]>buffer[q[3]]){
					temp=q[3];q[3]=q[2];q[2]=q[1];q[1]=temp;
					if(buffer[q[0]]>buffer[q[1]]){temp=q[0];q[0]=q[1];q[1]=temp;}
				}
				else{
					if(buffer[q[2]]>buffer[q[3]]){temp=q[2];q[2]=q[3];q[3]=temp;}
				}
				for(int i=0; i<4; i++){
					A.set(lo+i, buffer[q[i]]);
				}
				return;
		}
		value_type pivotvalue;
		if(hi-lo<16){
			size_t pivot = next()%(hi-lo)+lo;
			pivotvalue = buffer[pivot];
		}
		else{
			size_t pivot1 = next()%(hi-lo)+lo;
			size_t pivot2 = next()%(hi-lo)+lo;
			size_t pivot3 = next()%(hi-lo)+lo;
			size_t pivot = buffer[pivot1]>buffer[pivot2] ? (buffer[pivot2]>buffer[pivot3] ? pivot2 : (buffer[pivot1]>buffer[pivot3] ? pivot3 : pivot1)) : (buffer[pivot2]>buffer[pivot3] ? (buffer[pivot1]>buffer[pivot3] ? pivot1 : pivot3) : pivot2);
			pivotvalue = buffer[pivot];
		}
		int p = lo;
		int q = hi-1;
		int i;
		if(dir)i=hi-1;else i=lo;
		size_t point = next()%(hi-lo-1)+lo+1;
		while(i<hi&&i>=lo){
			int a = pivotvalue.cmp(buffer[i]);
			if(a+(dir^(i<point))>0){
				A.set(p, buffer[i]);
				p++;
			}
			else{
				A.set(q, buffer[i]);
				q--;
			}
			dir ? i-- : i++;
		}
		SugarCoat5991sort(A, lo, p, buffer, 0, 0);
		SugarCoat5991sort(A, p, hi, buffer, 0, 1);
	}
	else{
	if(hi-lo<=1)return;
		if(hi-lo==2){ if(dir ? A[lo]>=A[lo+1] : A[lo]>A[lo+1])A.swap(lo, lo+1); return;}
		if(hi-lo==3){ if(dir) {if(A[lo]>=A[lo+1])if(A[lo+1]>=A[lo+2])A.swap(lo,lo+2);else
			if(A[lo]>=A[lo+2]){A.swap(lo, lo+1);A.swap(lo+1, lo+2);}else A.swap(lo, lo+1);
			else if(A[lo+1]>=A[lo+2]) if(A[lo]>=A[lo+2]){A.swap(lo, lo+2);A.swap(lo+1, lo+2);}
			else A.swap(lo+1, lo+2); else ; return;} else
				{if(A[lo]>A[lo+1])if(A[lo+1]>A[lo+2])A.swap(lo,lo+2);else
			if(A[lo]>A[lo+2]){A.swap(lo, lo+1);A.swap(lo+1, lo+2);}else A.swap(lo, lo+1);
			else if(A[lo+1]>A[lo+2]) if(A[lo]>A[lo+2]){A.swap(lo, lo+2);A.swap(lo+1, lo+2);}
			else A.swap(lo+1, lo+2); else ; }
				return;}
		if(hi-lo==4){
				std::array<size_t, 4> q;
				if(dir) q = {lo+3, lo+2, lo+1, lo+0} ; else q = {lo, lo+1, lo+2, lo+3};
				size_t temp;
				if(A[q[0]]>A[q[1]]){temp=q[0];q[0]=q[1];q[1]=temp;}
				if(A[q[1]]>A[q[2]]){temp=q[1];q[1]=q[2];q[2]=temp;}
				if(q[2]!=lo+(2-dir)&&A[q[0]]>A[q[1]]){temp=q[0];q[0]=q[1];q[1]=temp;}
				if(A[q[1]]>A[q[3]]){
					temp=q[3];q[3]=q[2];q[2]=q[1];q[1]=temp;
					if(A[q[0]]>A[q[1]]){temp=q[0];q[0]=q[1];q[1]=temp;}
				}
				else{
					if(A[q[2]]>A[q[3]]){temp=q[2];q[2]=q[3];q[3]=temp;}
				}
				for(int i=0; i<4; i++){
					if(q[i]!=lo+i){
						A.swap(lo+i, q[i]);
						for(int j=i+1; j<4; j++){
							if(q[j]==lo+i)q[j]=q[i];
						}
					}
					/*while(q[i]!=lo+i){
						A.swap(lo+i, q[i]);
						temp=q[q[i]-lo];q[q[i]-lo]=q[i];q[i]=temp;
					}*/
				}
				return;
		}
		value_type pivotvalue;
		if(hi-lo<16){
			size_t pivot = next()%(hi-lo)+lo;
			pivotvalue = A[pivot];
		}
		else{
			size_t pivot1 = next()%(hi-lo)+lo;
			size_t pivot2 = next()%(hi-lo)+lo;
			size_t pivot3 = next()%(hi-lo)+lo;
			size_t pivot = A[pivot1]>A[pivot2] ? (A[pivot2]>A[pivot3] ? pivot2 : (A[pivot1]>A[pivot3] ? pivot3 : pivot1)) : (A[pivot2]>A[pivot3] ? (A[pivot1]>A[pivot3] ? pivot1 : pivot3) : pivot2);
			pivotvalue = A[pivot];
		}
		int p = lo;
		int q = hi-1;
		int i;
		if(dir)i=hi-1;else i=lo;
		size_t point = next()%(hi-lo-1)+lo+1;
		while(i<hi&&i>=lo){
			int a = pivotvalue.cmp(A[i]);
			if(a+(dir^(i<point))>0){
				buffer[p] = A[i];
				p++;
			}
			else{
				buffer[q] = A[i];
				q--;
			}
			dir ? i-- : i++;
		}
		SugarCoat5991sort(A, lo, p, buffer, 1, 0);
		SugarCoat5991sort(A, p, hi, buffer, 1, 1);
	}
}

void SugarCoat5991sort(SortArray& A)
{
	std::vector<value_type> buffer(A.size());
    return SugarCoat5991sort(A, 0, A.size(), buffer, 0, 0);
}

// ****************************************************************************
// *** BogoSort and more slow sorts

// by myself (Timo Bingmann)

bool BogoCheckSorted(SortArray& A)
{
    size_t i;
    value_type prev = A[0];
    A.mark(0);
    for (i = 1; i < A.size(); ++i)
    {
        value_type val = A[i];
        if (prev > val) break;
        prev = val;
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
        std::random_shuffle(perm.begin(), perm.end());

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
    while (1)
    {
        // check if array is sorted
        if (BogoCheckSorted(A)) break;

        // swap two random items
        A.swap(next() % A.size(), next() % A.size());
    }
}

void ZvaraySort(SortArray& A){
    while (1)
    {
        // check if array is sorted
        if (BogoCheckSorted(A)) break;

        // swap six random items
        size_t r[6];
        for(int i=0; i<6; i++)r[i]=next()%A.size();
        for(int i=6; i>1; i--){
		A.swap(r[i-1], r[next()%i]);
        }
    }
}

void BubblegumHillSort(SortArray& A){
    // keep a permutation of [0,size)
    std::vector<size_t> perm(A.size());

    for (size_t i = 0; i < A.size(); ++i)
        perm[i] = i;
    while (1)
    {
        // check if array is sorted
        if (BogoCheckSorted(A)) break;

        // swap eleven random items
        size_t r[11];
        for(int i=0; i<11; i++)r[i]=next()%A.size();
        for(int i=11; i>1; i--){
		A.swap(r[i-1], r[next()%(i-1)]);

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

        }
    }
}

void BubblegumHillSortII(SortArray& A){
        bool swapped = 1;
    while (1)
    {
        // check if array is sorted
        if (swapped && BogoCheckSorted(A)) break;
	swapped = 0;
        // swap eleven random items
        size_t r[11];
        for(int i=0; i<11; i++)r[i]=next()%A.size();
        for(int i=11; i>1; i--){
		int q = next()%(i-1);
		if(r[q]>r[i-1] ? A[r[q]]<A[r[i-1]] : A[r[q]]>A[r[i-1]]){
			A.swap(r[i-1], r[q]);
			swapped = 1;
		}
        }
    }
}

void RainbowReefSort(SortArray& A){
    while (1)
    {
        // check if array is sorted
        if (BogoCheckSorted(A)) break;

        // swap nine random items
        size_t r[9];
        for(int i=0; i<9; i++){r[i]=next()%A.size();A.mark(r[i], 3);}
        for(int i=9; i>1; i--){
		A.swap(r[i-1], r[next()%i]);
        }
        A.unmark_all();
    }
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
        // odd subsequence
    //#pragma omp parallel for
    for(int q=0; q<2; q++){
        oddEvenMerge(A, lo + q*r, n, m, sort_depth, merge_depth+1);
    }

    //#pragma omp parallel for
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
    //#pragma omp parallel for
    for(int q=0; q<2; q++){
        oddEvenMergeSort(A, lo + q*m, m, sort_depth+1);
    }
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

void ItBitonic(SortArray& array){
        int i, j, k;
        for(k = 2; k < array.size()*2; k = 2 * k) {
		int m = ((array.size()+(k-1))/k)%2;
            for(j = k >> 1; j > 0; j = j >> 1) {
                for(i = 0; i < array.size(); i++) {
                    int ij = i ^ j;
                    if((ij) > i && ij < array.size()) {
                        if((( !(i & k))^(!m)) && array[i] > array[ij])
                            array.swap(i, ij);
                        if(((!!(i & k))^(!m)) && array[i] < array[ij])
                            array.swap(i, ij);
                    }
                }
            }
        }
}

void ParBitonic(SortArray& array){
        int i, j, k;
        for(k = 2; k < array.size()*2; k = 2 * k) {
		int m = ((array.size()+(k-1))/k)%2;
            for(j = k >> 1; j > 0; j = j >> 1) {
#pragma omp parallel for
                for(i = 0; i < array.size(); i++) {
                    int ij = i ^ j;
                    if((ij) > i && ij < array.size()) {
                        if((( !(i & k))^(!m)) && array[i] > array[ij])
                            array.swap(i, ij);
                        if(((!!(i & k))^(!m)) && array[i] < array[ij])
                            array.swap(i, ij);
                    }
                }
            }
        }
}

void ItOddEvenMerge(SortArray& array){
        for (int p = 1; p < array.size(); p += p)
          for (int k = p; k > 0; k /= 2)
            for (int j = k % p; j + k < array.size(); j += k + k)
              for (int i = 0; i < k; i++)
                if ((i + j)/(p + p) == (i + j + k)/(p + p)) {
			if(i+j+k < array.size()){
                   if(array[i+j] > array[i+j+k])
			array.swap(i+j, i+j+k);
			}
                }
}

void ParOddEvenMerge(SortArray& array){
        for (int p = 1; p < array.size(); p += p)
          for (int k = p; k > 0; k /= 2)
#pragma omp parallel for
            for (int j = k % p; j < array.size(); j ++ )
              if (!((j-k%p)/k%2))
                if ((j)/(p + p) == (j + k)/(p + p)) {
			if(j+k < array.size()){
                   if(array[j] > array[j+k])
			array.swap(j, j+k);
			}
                }
}

void ItPairwise(SortArray& a)
{
	int start = 0;
	int end = a.size();
        int a2 = 1;
        int b = 0;
        int c = 0;
        int d = 0;
        int e = 0;
        while (a2 < end){
            b = start+a2;
            c = 0;
            while (b < end){
                if(a[b-a2]>a[b]) a.swap(b-a2, b);
                c = (c + 1) % a2;
                b++;
                if (c == 0){
                    b += a2;
                }
            }
            a2 *= 2;
        }
        a2 /= 4;
        e = 1;
        while (a2 > 0){
            d = e;
            while (d > 0){
                b = start+((d + 1) * a2);
                c = 0;
                while (b < end){
                if(a[b - (d * a2)]>a[b]) a.swap(b - (d * a2), b);
                    c = (c + 1) % a2;
                    b++;
                    if (c == 0){
                        b += a2;
                    }
                }
                d /= 2;
            }
            a2 /= 2;
            e = (e * 2) + 1;
        }
}

void ParPairwise(SortArray& a)
{
	int start = 0;
	int end = a.size();
        int a2 = 1;
        int b = 0;
        int c = 0;
        int d = 0;
        int e = 0;
        while (a2 < end){
            c = 0;
#pragma omp parallel for
            for (b=0; b < end-start-a2; b++){
		if(!(b/a2%2))
                if(a[b+start]>a[b+start+a2]) a.swap(b+start, b+start+a2);
            }
            a2 *= 2;
        }
        a2 /= 4;
        e = 1;
        while (a2 > 0){
            d = e;
            while (d > 0){
                c = 0;
#pragma omp parallel for
                for (b=start+((d + 1) * a2); b < end; b++){
			if(!((b-start+((d + 1) * a2))/a2%2))
			if(a[b - (d * a2)]>a[b]) a.swap(b - (d * a2), b);

                }
                d /= 2;
            }
            a2 /= 2;
            e = (e * 2) + 1;
        }
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

void CycleSort(SortArray& array, ssize_t n)
{
    volatile ssize_t cycleStart = 0;
    array.watch(&cycleStart, 16);

    volatile ssize_t rank = 0;
    array.watch(&rank, 3);

    // Loop through the array to find cycles to rotate.
    for (cycleStart = 0; cycleStart < n - 1; ++cycleStart)
    {
        value_type& item = array.get_mutable(cycleStart);

        do {
            // Find where to put the item.
            rank = cycleStart;
            for (ssize_t i = cycleStart + 1; i < n; ++i)
            {
                if (array[i] < item)
                    rank++;
            }

            // If the item is already there, this is a 1-cycle.
            if (rank == cycleStart) {
                array.mark(rank, 2);
                break;
            }

            // Otherwise, put the item after any duplicates.
            while (item == array[rank])
                rank++;

            // Put item into right place and colorize
            std::swap(array.get_mutable(rank), item);
            array.mark(rank, 2);

            // Continue for rest of the cycle.
        }
        while (rank != cycleStart);
    }

    array.unwatch_all();
}

void CycleSort(SortArray& A)
{
    CycleSort(A, A.size());
}

void grailMergeLeft(SortArray& arr, int pos, int leftLen, int rightLen, int dist);

    int leftOverLen;
    int leftOverFrag;

    struct GrailState{
        int leftOverLen;
        int leftOverFrag;
	int getLeftOverLen() {
		return leftOverLen;
	}
	int getLeftOverFrag() {
		return leftOverFrag;
	}
    };

    GrailState grailSmartMergeWithBuffer(SortArray& arr, int pos, int leftOverLen, int leftOverFrag, int blockLen);
    GrailState grailSmartMergeWithoutBuffer(SortArray& arr, int pos, int leftOverLen, int leftOverFrag, int regBlockLen);

    int grailStaticBufferLen = 32; //Buffer length changed due to less numbers in this program being sorted than what Mr. Astrelin used for testing.

    int getStaticBuffer() {
        return grailStaticBufferLen;
    }

    void grailSwap(SortArray& arr, int a, int b) {
        arr.swap(a, b);
    }

    void grailMultiSwap(SortArray& arr, int a, int b, int swapsLeft) {
        while(swapsLeft != 0) {
            grailSwap(arr, a++, b++);
            swapsLeft--;
        }
    }

    void grailRotate(SortArray& array, int pos, int lenA, int lenB) {
        while(lenA != 0 && lenB != 0) {
            if(lenA <= lenB) {
                grailMultiSwap(array, pos, pos + lenA, lenA);
                pos += lenA;
                lenB -= lenA;
            }
            else {
                grailMultiSwap(array, pos + (lenA - lenB), pos + lenA, lenB);
                lenA -= lenB;
            }
        }
    }

    void grailInsertSort(SortArray& A, int pos, int len) {
	    for (size_t i = pos+1; i < pos+len; ++i)
	    {
		ssize_t j = i - 1;
		while (j >= pos && A[j] > A[j+1])
		{
		    A.swap(j, j+1);
		    j--;
		}
	    }
    }

    //boolean argument determines direction
    int grailBinSearch(SortArray& arr, int pos, int len, int keyPos, boolean isLeft) {
        int left = -1, right = len;
        while(left < right - 1) {
            int mid = left + ((right - left) >> 1);
            if(isLeft) {
                if(arr[pos + mid] >= arr[keyPos]) {
                    right = mid;
                } else {
                    left = mid;
                }
            } else {
                if(arr[pos + mid] > arr[keyPos]) {
                    right = mid;
                } else left = mid;
            }
        }
        return right;
    }

    // cost: 2 * len + numKeys^2 / 2
    int grailFindKeys(SortArray& arr, int pos, int len, int numKeys) {
        int dist = 1, foundKeys = 1, firstKey = 0;  // first key is always here

        while(dist < len && foundKeys < numKeys) {
            //Binary Search left
            int loc = grailBinSearch(arr, pos + firstKey, foundKeys, pos + dist, true);
            if(loc == foundKeys || arr[pos + dist] != arr[pos + (firstKey + loc)]) {
                grailRotate(arr, pos + firstKey, foundKeys, dist - (firstKey + foundKeys));
                firstKey = dist - foundKeys;
                grailRotate(arr, pos + (firstKey + loc), foundKeys - loc, 1);
                foundKeys++;
            }
            dist++;
        }
        grailRotate(arr, pos, firstKey, foundKeys);

        return foundKeys;
    }

    // cost: min(len1, len2)^2 + max(len1, len2)
    void grailMergeWithoutBuffer(SortArray& arr, int pos, int len1, int len2) {
        if(len1 < len2) {
            while(len1 != 0) {
                //Binary Search left
                int loc = grailBinSearch(arr, pos + len1, len2, pos, true);
                if(loc != 0) {
                    grailRotate(arr, pos, len1, loc);
                    pos += loc;
                    len2 -= loc;
                }
                if(len2 == 0) break;
                do {
                    pos++;
                    len1--;
                } while(len1 != 0 && arr[pos] <= arr[pos + len1]);
            }
        } else {
            while(len2 != 0) {
                //Binary Search right
                int loc = grailBinSearch(arr, pos, len1, pos + (len1 + len2 - 1), false);
                if(loc != len1) {
                    grailRotate(arr, pos + loc, len1 - loc, len2);
                    len1 = loc;
                }
                if(len1 == 0) break;
                do {
                    len2--;
                } while(len2 != 0 && arr[pos + len1 - 1] <= arr[pos + len1 + len2 - 1]);
            }
        }
    }

    // arr - starting array. arr[0 - regBlockLen..-1] - buffer (if havebuf).
    // regBlockLen - length of regular blocks. First blockCount blocks are stable sorted by 1st elements and key-coded
    // keysPos - arrays of keys, in same order as blocks. keysPos < midkey means stream A
    // aBlockCount are regular blocks from stream A.
    // lastLen is length of last (irregular) block from stream B, that should go before nblock2 blocks.
    // lastLen = 0 requires aBlockCount = 0 (no irregular blocks). lastLen > 0, aBlockCount = 0 is possible.
    void grailMergeBuffersLeft(SortArray& arr, int keysPos, int midkey, int pos,
            int blockCount, int blockLen, boolean havebuf, int aBlockCount,
            int lastLen) {

        if(blockCount == 0) {
            int aBlocksLen = aBlockCount * blockLen;
            if(havebuf) grailMergeLeft(arr, pos, aBlocksLen, lastLen, 0 - blockLen);
            else grailMergeWithoutBuffer(arr, pos, aBlocksLen, lastLen);
            return;
        }

        int leftOverLen = blockLen;
        int leftOverFrag = arr[keysPos] >= arr[midkey];
        int processIndex = blockLen;
        int restToProcess;

        for(int keyIndex = 1; keyIndex < blockCount; keyIndex++, processIndex += blockLen) {
            restToProcess = processIndex - leftOverLen;
            int nextFrag = arr[keysPos + keyIndex] >= arr[midkey];

            if(nextFrag == leftOverFrag) {
                if(havebuf) grailMultiSwap(arr, pos + restToProcess - blockLen, pos + restToProcess, leftOverLen);
                restToProcess = processIndex;
                leftOverLen = blockLen;
            } else {
                if(havebuf) {
                    GrailState results = grailSmartMergeWithBuffer(arr, pos + restToProcess, leftOverLen, leftOverFrag, blockLen);
                    leftOverLen = results.getLeftOverLen();
                    leftOverFrag = results.getLeftOverFrag();
                } else {
                    GrailState results = grailSmartMergeWithoutBuffer(arr, pos + restToProcess, leftOverLen, leftOverFrag, blockLen);
                    leftOverLen = results.getLeftOverLen();
                    leftOverFrag = results.getLeftOverFrag();
                }
            }
        }
        restToProcess = processIndex - leftOverLen;

        if(lastLen != 0) {
            if(leftOverFrag != 0) {
                if(havebuf) {
                    grailMultiSwap(arr, pos + restToProcess - blockLen, pos + restToProcess, leftOverLen);
                }
                restToProcess = processIndex;
                leftOverLen = blockLen * aBlockCount;
                leftOverFrag = 0;
            } else {
                leftOverLen += blockLen * aBlockCount;
            }
            if(havebuf) {
                grailMergeLeft(arr, pos + restToProcess, leftOverLen, lastLen, -blockLen);
            }
            else {
                grailMergeWithoutBuffer(arr, pos + restToProcess, leftOverLen, lastLen);
            }
        } else {
            if(havebuf) {
                grailMultiSwap(arr, pos + restToProcess, pos + (restToProcess - blockLen), leftOverLen);
            }
        }
    }

    // arr[dist..-1] - buffer, arr[0, leftLen - 1] ++ arr[leftLen, leftLen + rightLen - 1]
    // -> arr[dist, dist + leftLen + rightLen - 1]
    void grailMergeLeft(SortArray& arr, int pos, int leftLen, int rightLen, int dist) {
        int left = 0;
        int right = leftLen;

        rightLen += leftLen;

        while(right < rightLen) {
            if(left == leftLen || arr[pos + left] > arr[pos + right]) {
                grailSwap(arr, pos + (dist++), pos + (right++));
            }
            else grailSwap(arr, pos + (dist++), pos + (left++));
        }

        if(dist != left) grailMultiSwap(arr, pos + dist, pos + left, leftLen - left);
    }
    void grailMergeRight(SortArray& arr, int pos, int leftLen, int rightLen, int dist) {
        int mergedPos = leftLen + rightLen + dist - 1;
        int right = leftLen + rightLen - 1;
        int left = leftLen - 1;

        while(left >= 0) {
            if(right < leftLen || arr[pos + left] > arr[pos + right]) {
                grailSwap(arr, pos + (mergedPos--), pos + (left--));
            }
            else grailSwap(arr, pos + (mergedPos--), pos + (right--));
        }

        if(right != mergedPos) {
            while(right >= leftLen) grailSwap(arr, pos + (mergedPos--), pos + (right--));
        }
    }

    //returns the leftover length, then the leftover fragment
    GrailState grailSmartMergeWithoutBuffer(SortArray& arr, int pos, int leftOverLen, int leftOverFrag, int regBlockLen) {
        if(regBlockLen == 0){	GrailState a = {leftOverLen, leftOverFrag}; return a; }

        int len1 = leftOverLen;
        int len2 = regBlockLen;
        int typeFrag = 1 - leftOverFrag; //1 if inverted

        if(len1 != 0 && arr[pos + (len1 - 1)].cmp(arr[pos + len1]) - typeFrag >= 0) {

            while(len1 != 0) {
                int foundLen;
                if (typeFrag != 0) {
                    //Binary Search left
                    foundLen = grailBinSearch(arr, pos + len1, len2, pos, true);
                } else {
                    //Binary Search right
                    foundLen = grailBinSearch(arr, pos + len1, len2, pos, false);
                }
                if(foundLen != 0) {
                    grailRotate(arr, pos, len1, foundLen);
                    pos += foundLen;
                    len2 -= foundLen;
                }
                if(len2 == 0) {
			GrailState a = {len1, leftOverFrag};
			return a;
                }
                do {
                    pos++;
                    len1--;
                } while(len1 != 0 && arr[pos].cmp(arr[pos + len1]) - typeFrag < 0);
            }
        }
        GrailState a = {len2, typeFrag};
        return a;
    }

    //returns the leftover length, then the leftover fragment
    GrailState grailSmartMergeWithBuffer(SortArray& arr, int pos, int leftOverLen, int leftOverFrag, int blockLen) {
        int dist = 0 - blockLen, left = 0, right = leftOverLen, leftEnd = right, rightEnd = right + blockLen;
        int typeFrag = 1 - leftOverFrag;  // 1 if inverted

        while(left < leftEnd && right < rightEnd) {
            if(arr[pos + left].cmp(arr[pos + right]) - typeFrag < 0) {
                grailSwap(arr, pos + (dist++), pos + (left++));
            }
            else grailSwap(arr, pos + (dist++), pos + (right++));
        }

        int length, fragment = leftOverFrag;
        if(left < leftEnd) {
            length = leftEnd - left;
            while(left < leftEnd) grailSwap(arr, pos + (--leftEnd), pos + (--rightEnd));
        } else {
            length = rightEnd - right;
            fragment = typeFrag;
        }
        GrailState a = {length, fragment};
        return a;
    }


    /***** Sort With Extra Buffer *****/

    //returns the leftover length, then the leftover fragment
    GrailState grailSmartMergeWithXBuf(SortArray& arr, int pos, int leftOverLen, int leftOverFrag, int blockLen) {
        int dist = 0 - blockLen, left = 0, right = leftOverLen, leftEnd = right, rightEnd = right + blockLen;
        int typeFrag = 1 - leftOverFrag;  // 1 if inverted

        while(left < leftEnd && right < rightEnd) {
            if(arr[pos + left].cmp(arr[pos + right]) - typeFrag < 0) {
                arr.set(pos + dist++, arr[pos + left++]);
            }
            else arr.set(pos + dist++, arr[pos + right++]);
        }

        int length, fragment = leftOverFrag;
        if(left < leftEnd) {
            length = leftEnd - left;
            while(left < leftEnd) arr.set(pos + --rightEnd, arr[pos + --leftEnd]);
        } else {
            length = rightEnd - right;
            fragment = typeFrag;
        }
        GrailState a = {length, fragment};
        return a;
    }

    // arr[dist..-1] - free, arr[0, leftEnd - 1] ++ arr[leftEnd, leftEnd + rightEnd - 1]
    // -> arr[dist, dist + leftEnd + rightEnd - 1]
    void grailMergeLeftWithXBuf(SortArray& arr, int pos, int leftEnd, int rightEnd, int dist) {
        int left = 0;
        int right = leftEnd;
        rightEnd += leftEnd;

        while(right < rightEnd) {
            if(left == leftEnd || arr[pos + left] > arr[pos + right]) {
                arr.set(pos + dist++, arr[pos + right++]);
            }
            else arr.set(pos + dist++, arr[pos + left++]);
        }

        if(dist != left) {
            while(left < leftEnd) arr.set(pos + dist++, arr[pos + left++]);
        }
    }

    // arr - starting array. arr[0 - regBlockLen..-1] - buffer (if havebuf).
    // regBlockLen - length of regular blocks. First blockCount blocks are stable sorted by 1st elements and key-coded
    // keysPos - where keys are in array, in same order as blocks. keysPos < midkey means stream A
    // aBlockCount are regular blocks from stream A.
    // lastLen is length of last (irregular) block from stream B, that should go before aCountBlock blocks.
    // lastLen = 0 requires aBlockCount = 0 (no irregular blocks). lastLen > 0, aBlockCount = 0 is possible.
    void grailMergeBuffersLeftWithXBuf(SortArray& arr, int keysPos, int midkey, int pos,
            int blockCount, int regBlockLen, int aBlockCount, int lastLen) {

        if(blockCount == 0) {
            int aBlocksLen = aBlockCount * regBlockLen;
            grailMergeLeftWithXBuf(arr, pos, aBlocksLen, lastLen, 0 - regBlockLen);
            return;
        }

        int leftOverLen = regBlockLen;
        int leftOverFrag = arr[keysPos] >= arr[midkey];
        int processIndex = regBlockLen;

        int restToProcess;
        for(int keyIndex = 1; keyIndex < blockCount; keyIndex++, processIndex += regBlockLen) {
            restToProcess = processIndex - leftOverLen;
            int nextFrag = arr[keysPos + keyIndex] >= arr[midkey];

            if(nextFrag == leftOverFrag) {
		for(int i=0; i<leftOverLen; i++){
			arr.set(pos+restToProcess-regBlockLen+i, arr[pos+restToProcess+i]);
		}

                restToProcess = processIndex;
                leftOverLen = regBlockLen;
            } else {
                GrailState results = grailSmartMergeWithXBuf(arr, pos + restToProcess, leftOverLen, leftOverFrag, regBlockLen);
                leftOverLen = results.getLeftOverLen();
                leftOverFrag = results.getLeftOverFrag();
            }
        }
        restToProcess = processIndex - leftOverLen;

        if(lastLen != 0) {
            if(leftOverFrag != 0) {
		for(int i=0; i<leftOverLen; i++){
			arr.set(pos+restToProcess-regBlockLen+i, arr[pos+restToProcess+i]);
		}

                restToProcess = processIndex;
                leftOverLen = regBlockLen * aBlockCount;
                leftOverFrag = 0;
            } else {
                leftOverLen += regBlockLen * aBlockCount;
            }
            grailMergeLeftWithXBuf(arr, pos + restToProcess, leftOverLen, lastLen, 0 - regBlockLen);
        } else {
		for(int i=0; i<leftOverLen; i++){
			arr.set(pos+restToProcess-regBlockLen+i, arr[pos+restToProcess+i]);
		}
        }
    }

    /***** End Sort With Extra Buffer *****/

    // build blocks of length buildLen
    // input: [-buildLen, -1] elements are buffer
    // output: first buildLen elements are buffer, blocks 2 * buildLen and last subblock sorted
    void grailBuildBlocks(SortArray& arr, int pos, int len, int buildLen,
            SortArray& extbuf, int bufferPos, int extBufLen) {

        int buildBuf = buildLen < extBufLen ? buildLen : extBufLen;
        while((buildBuf & (buildBuf - 1)) != 0) buildBuf &= buildBuf - 1;  // max power or 2 - just in case

        int extraDist, part;
        if(buildBuf != 0) {
		for(int i=0; i<buildBuf; i++){
			extbuf.set(bufferPos+i, arr[pos-buildBuf+i]);
		}

            for(int dist = 1; dist < len; dist += 2) {
                extraDist = 0;
                if(arr[pos + (dist - 1)] > arr[pos + dist]) extraDist = 1;
                arr.set(pos + dist - 3, arr[pos + dist - 1 + extraDist]);
                arr.set(pos + dist - 2, arr[pos + dist - extraDist]);
            }
            if(len % 2 != 0) arr.set(pos + len - 3, arr[pos + len - 1]);
            pos -= 2;

            for(part = 2; part < buildBuf; part *= 2) {
                int left = 0;
                int right = len - 2 * part;
                while(left <= right) {
                    grailMergeLeftWithXBuf(arr, pos + left, part, part, 0 - part);
                    left += 2 * part;
                }
                int rest = len - left;

                if(rest > part) {
                    grailMergeLeftWithXBuf(arr, pos + left, part, rest - part, 0 - part);
                } else {
                    for(; left < len; left++) arr.set(pos + left - part, arr[pos + left]);
                }
                pos -= part;
            }
            for(int i=0; i<buildBuf; i++){
		arr.set(pos+len+i, extbuf[bufferPos+i]);
            }
        }
        else {
            for(int dist = 1; dist < len; dist += 2) {
                extraDist = 0;
                if(arr[pos + (dist - 1)] > arr[pos + dist]) extraDist = 1;
                grailSwap(arr, pos + (dist - 3), pos + (dist - 1 + extraDist));
                grailSwap(arr, pos + (dist - 2), pos + (dist - extraDist));
            }
            if(len % 2 != 0) grailSwap(arr, pos + (len - 1), pos + (len - 3));
            pos -= 2;
            part = 2;
        }

        for(; part < buildLen; part *= 2) {
            int left = 0;
            int right = len - 2 * part;
            while(left <= right) {
                grailMergeLeft(arr, pos + left, part, part, 0 - part);
                left += 2 * part;
            }
            int rest = len - left;
            if(rest > part) {
                grailMergeLeft(arr, pos + left, part, rest - part, 0 - part);
            } else {
                grailRotate(arr, pos + left - part, part, rest);
            }
            pos -= part;
        }
        int restToBuild = len % (2 * buildLen);
        int leftOverPos = len - restToBuild;

        if(restToBuild <= buildLen) grailRotate(arr, pos + leftOverPos, restToBuild, buildLen);
        else grailMergeRight(arr, pos + leftOverPos, buildLen, restToBuild - buildLen, buildLen);

        while(leftOverPos > 0) {
            leftOverPos -= 2 * buildLen;
            grailMergeRight(arr, pos + leftOverPos, buildLen, buildLen, buildLen);
        }
    }

    // keys are on the left of arr. Blocks of length buildLen combined. We'll combine them in pairs
    // buildLen and nkeys are powers of 2. (2 * buildLen / regBlockLen) keys are guaranteed
    void grailCombineBlocks(SortArray& arr, int keyPos, int pos, int len, int buildLen,
            int regBlockLen, boolean havebuf, SortArray& buffer, int bufferPos) {

        int combineLen = len / (2 * buildLen);
        int leftOver = len % (2 * buildLen);
        if(leftOver <= buildLen) {
            len -= leftOver;
            leftOver = 0;
        }

        if(buffer.size()) {
		for(int i=0; i<regBlockLen; i++){
			buffer.set(bufferPos+i, arr[pos-regBlockLen+i]);
		}
	}

        for(int i = 0; i <= combineLen; i++) {
            if(i == combineLen && leftOver == 0) break;

            int blockPos = pos + i * 2 * buildLen;
            int blockCount = (i == combineLen ? leftOver : 2 * buildLen) / regBlockLen;

            grailInsertSort(arr, keyPos, blockCount + (i == combineLen ? 1 : 0));

            int midkey = buildLen / regBlockLen;

            for(int index = 1; index < blockCount; index++) {
                int leftIndex = index - 1;

                for(int rightIndex = index; rightIndex < blockCount; rightIndex++) {
                    int rightComp = arr[blockPos + leftIndex * regBlockLen].cmp(arr[blockPos + rightIndex * regBlockLen]);
                    if(rightComp > 0 || (rightComp == 0 && arr[keyPos + leftIndex] > arr[keyPos + rightIndex])) {
                        leftIndex = rightIndex;
                    }
                }
                if(leftIndex != index - 1) {
                    grailMultiSwap(arr, blockPos + (index - 1) * regBlockLen, blockPos + leftIndex * regBlockLen, regBlockLen);
                    grailSwap(arr, keyPos + (index - 1), keyPos + leftIndex);
                    if(midkey == index - 1 || midkey == leftIndex) {
                        midkey ^= (index - 1) ^ leftIndex;
                    }
                }
            }

            int aBlockCount = 0;
            int lastLen = 0;
            if(i == combineLen) lastLen = leftOver % regBlockLen;

            if(lastLen != 0) {
                while(aBlockCount < blockCount && arr[blockPos + blockCount * regBlockLen]
                         < arr[blockPos + (blockCount - aBlockCount - 1) * regBlockLen]) {
                    aBlockCount++;
                }
            }

            if(buffer.size()) {
                grailMergeBuffersLeftWithXBuf(arr, keyPos, keyPos + midkey, blockPos,
                        blockCount - aBlockCount, regBlockLen, aBlockCount, lastLen);
            }
            else grailMergeBuffersLeft(arr, keyPos, keyPos + midkey, blockPos,
                    blockCount - aBlockCount, regBlockLen, havebuf, aBlockCount, lastLen);
        }
        if(buffer.size()) {
            for(int i = len; --i >= 0;) arr.set(pos + i, arr[pos + i - regBlockLen]);
            for(int i=0; i<regBlockLen; i++){
		arr.set(pos-regBlockLen+i, buffer[bufferPos+i]);
            }
        }
        else if(havebuf) {
            while(--len >= 0) {
                grailSwap(arr, pos + len, pos + len - regBlockLen);
            }
        }
    }

    void grailLazyStableSort(SortArray& arr, int pos, int len) {
        for(int dist = 1; dist < len; dist += 2) {
            if(arr[pos + dist - 1] > arr[pos + dist]) {
                grailSwap(arr, pos + (dist - 1), pos + dist);
            }
        }

        for(int part = 2; part < len; part *= 2) {
            int left = 0;
            int right = len - 2 * part;

            while(left <= right) {
                grailMergeWithoutBuffer(arr, pos + left, part, part);
                left += 2 * part;
            }

            int rest = len - left;
            if(rest > part) {
                grailMergeWithoutBuffer(arr, pos + left, part, rest - part);
            }
        }
    }

    void grailCommonSort(SortArray& arr, int pos, int len, SortArray& buffer, int bufferPos, int bufferLen) {

        if(len <= 16) {
            grailInsertSort(arr, pos, len);
            return;
        }

        int blockLen = 1;
        while(blockLen * blockLen < len) blockLen *= 2;

        int numKeys = (len - 1) / blockLen + 1;
        int keysFound = grailFindKeys(arr, pos, len, numKeys + blockLen);

        boolean bufferEnabled = true;

        if(keysFound < numKeys + blockLen) {
            if(keysFound < 4) {
                grailLazyStableSort(arr, pos, len);
                return;
            }
            numKeys = blockLen;
            while(numKeys > keysFound) numKeys /= 2;
            bufferEnabled = false;
            blockLen = 0;
        }

        int dist = blockLen + numKeys;
        int buildLen = bufferEnabled ? blockLen : numKeys;

        if(bufferEnabled) {
            grailBuildBlocks(arr, pos + dist, len - dist, buildLen, buffer, bufferPos, bufferLen);
        }
        else {
            grailBuildBlocks(arr, pos + dist, len - dist, buildLen, buffer, bufferPos, 0);
        }

        // 2 * buildLen are built
        while(len - dist > (buildLen *= 2)) {
            int regBlockLen = blockLen;
            boolean buildBufEnabled = bufferEnabled;

            if(!bufferEnabled) {
                if(numKeys > 4 && numKeys / 8 * numKeys >= buildLen) {
                    regBlockLen = numKeys / 2;
                    buildBufEnabled = true;
                } else {
                    int calcKeys = 1;
                    int i = buildLen * keysFound / 2;
                    while(calcKeys < numKeys && i != 0) {
                        calcKeys *= 2;
                        i /= 8;
                    }
                    regBlockLen = (2 * buildLen) / calcKeys;
                }
            }
            grailCombineBlocks(arr, pos, pos + dist, len - dist, buildLen, regBlockLen, buildBufEnabled,
                    buildBufEnabled && regBlockLen <= bufferLen ? buffer : buffer, bufferPos);

        }

        grailInsertSort(arr, pos, dist);
        grailMergeWithoutBuffer(arr, pos, dist, len - dist);
    }

    void grailInPlaceMerge(SortArray& arr, int pos, int len1, int len2) {
        if(len1 < 3 || len2 < 3) {
            grailMergeWithoutBuffer(arr, pos, len1, len2);
            return;
        }

        int midpoint;
        if(len1 < len2) midpoint = len1 + len2 / 2;
        else midpoint = len1 / 2;

        //Left binary search
        int len1Left, len1Right;
        len1Left = len1Right = grailBinSearch(arr, pos, len1, pos + midpoint, true);

        //Right binary search
        if(len1Right < len1 && arr[pos + len1Right] == arr[pos + midpoint]) {
            len1Right = grailBinSearch(arr, pos + len1Left, len1 - len1Left, pos + midpoint, false) + len1Left;
        }

        int len2Left, len2Right;
        len2Left = len2Right = grailBinSearch(arr, pos + len1, len2, pos + midpoint, true);

        if(len2Right < len2 && arr[pos + len1 + len2Right] == arr[pos + midpoint]) {
            len2Right = grailBinSearch(arr, pos + len1 + len2Left, len2 - len2Left, pos + midpoint, false) + len2Left;
        }

        if(len1Left == len1Right) grailRotate(arr, pos + len1Right, len1 - len1Right, len2Right);
        else {
            grailRotate(arr, pos + len1Left, len1 - len1Left, len2Left);

            if(len2Right != len2Left) {
                grailRotate(arr, pos + (len1Right + len2Left), len1 - len1Right, len2Right - len2Left);
            }
        }

        grailInPlaceMerge(arr, pos + (len1Right + len2Right), len1 - len1Right, len2 - len2Right);
        grailInPlaceMerge(arr, pos, len1Left, len2Left);
    }
    void grailInPlaceMergeSort(SortArray& arr, int start, int len) {
        for(int dist = start + 1; dist < len; dist += 2) {
            if(arr[dist - 1] > arr[dist]) grailSwap(arr, dist - 1, dist);
        }
        for(int part = 2; part < len; part *= 2) {
            int left = start, right = len - 2 * part;

            while(left <= right) {
                grailInPlaceMerge(arr, left, part, part);
                left += 2 * part;
            }

            int rest = len - left;
            if(rest > part) grailInPlaceMerge(arr, left, part, rest - part);
        }
    }

void GrailSort(class SortArray& a){
	SortArray buffer = {};
	grailStaticBufferLen = buffer.size();
    grailCommonSort(a, 0, a.size(), buffer, 0, grailStaticBufferLen);
}

void LazyStableSort(class SortArray& a){
    grailLazyStableSort(a, 0, a.size());
}

void RotateMergeSort(class SortArray& a){
    grailInPlaceMergeSort(a, 0, a.size());
}
