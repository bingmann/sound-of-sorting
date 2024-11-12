/*
 *
	MIT License

	Copyright (c) 2021-2024 aphitorite, edited by Flanlaina (a.k.a. Ayako-chan)

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

/**
	ultimate sorting algorithm:

	unstable and probabilistic sorting algorithm performing an average of
	n log n + ~4.4 n comparisons and O(n) moves (~12.1 n) in O(1) space

	makes the fewest comparisons among the sorts of its kind
	achieves the ultimate goal of a n log n + O(n) comps in place O(n) moves sort

	implements a modified library sort where the original algorithm is proven to be O(n log n)
	from http://www.cs.sunysb.edu/~bender/newpub/BenderFaMo06-librarysort.pdf

	@author aphitorite
*/

#include "../SortAlgo.h"
#include <random>

const int MIN_INSERT = 32;
const int G = 7;
const int R = 3;

extern std::random_device rd1;
extern std::mt19937 gen1;

void shiftBW(SortArray& A, int a, int m, int b)
{
	while (m > a) { --b; --m; A[b].get(); A.swap(b, m); }
}

int randGapSearch(SortArray& A, int a, int b, int val)
{
	int s = G + 1, randCnt = 0;
	while (a < b)
	{
		int m = a + (((b - a) / s) / 2) * s;
		int ele = A[m].get();
		if (val < ele) { b = m; }
		else if (val > ele) { a = m + s; }
		else if (randCnt++ < 1)
		{
			std::uniform_int_distribution<> distr(0, ((b - a) / s) - 1);
			m = a + distr(gen1) * s;
			ele = A[m].get();
			if (val < ele) { b = m; }
			else if (val > ele) { a = m + s; }
			else { a = m; break; }
		}
		else { a = m; break; }
	}
	return a;
}

int rightBinSearch(SortArray& A, int a, int b, int val, bool bw)
{
	int cmp = -1;
	if (bw == true) { cmp = 1; }
	while (a < b)
	{
		int m = a + (b - a) / 2;
		int ele = A[m].get(), result = 0;
		if (val < ele) { result = -1; }
		else if (val > ele) { result = 1; }

		if (result == cmp) { b = m; }
		else { a = m + 1; }
	}
	return a;
}

void insertTo(SortArray& A, int tmp, int a, int b)
{
	while (a > b) 
	{ 
		size_t z = static_cast<size_t>(a);
		A[z].get();
		A.set(z, A[z - 1]); 
		--a; 
	}
	size_t z = static_cast<size_t>(b);
	A.set(z, ArrayItem(tmp));
	A[z].get();
}

void binaryInsert(SortArray& A, int a, int b)
{
	for (int i = a + 1; i < b; ++i)
	{
		int ele = A[i].get();
		insertTo(A, ele, i, rightBinSearch(A, a, i, ele, false));
	}
}

void retrieve(SortArray& A, int b, int p, int pEnd, int bsv, bool bw)
{
	int j = b - 1, m;
	for (int k = pEnd - (G + 1); k > p + G; )
	{
		m = rightBinSearch(A, k - G, k, bsv, bw) - 1;
		k -= G + 1;
		while (m >= k) 
		{ 
			size_t x = static_cast<size_t>(j), y = static_cast<size_t>(m);
			A[x].get();
			A.swap(x, y); 
			--j; --m; 
		}
	}
	m = rightBinSearch(A, p, p + G, bsv, bw) - 1;
	while (m >= p) 
	{
		size_t x = static_cast<size_t>(j), y = static_cast<size_t>(m);
		A[x].get();
		A.swap(x, y); 
		--j; --m; 
	}
}

void rescatter(SortArray & A, int i, int p, int pEnd, int bsv, bool bw)
{
	int j = i - 2 * (G + 1), m;
	for (int k = pEnd - (G + 1); k > p + G; )
	{
		m = rightBinSearch(A, k - G, k, bsv, bw) - 1;
		k -= G + 1;
		while (m >= k)
		{
			size_t x = static_cast<size_t>(j), y = static_cast<size_t>(m);
			A[x].get();
			A.swap(x, y);
			j -= G + 1; --m;
		}
	}
	m = rightBinSearch(A, p, p + G, bsv, bw) - 1;
	while (m >= p)
	{
		size_t x = static_cast<size_t>(j), y = static_cast<size_t>(m);
		A[x].get();
		A.swap(x, y);
		j -= G + 1; --m;
	}
}

void rebalance(SortArray& A, int a, int b, int p, int pEnd, int pb, int bsv, bool bw)
{
	retrieve(A, b, p, pEnd, bsv, bw);
	int gCnt = (pb - p) / (G + 1);
	int gapElems = (b - a) - gCnt + 1;
	int baseCnt = gapElems / gCnt;
	int extra = gapElems - gCnt * baseCnt;
	for (int k = p + G;; k += G + 1)
	{
		int val = 0;
		if (extra > 0) { val = 1; }
		int iter = baseCnt + val;
		--extra;
		for (int j = 0; j < iter; ++j)
		{
			size_t x = static_cast<size_t>(a), y = static_cast<size_t>(k - G + j);
			A[x].get();
			A.swap(x, y);
			++a;
		}
		if (k < pb - (G + 1))
		{
			size_t x = static_cast<size_t>(a), y = static_cast<size_t>(k);
			A[x].get();
			A.swap(x, y);
			++a;
		}
		else { break; }
	}
}

void librarySort(SortArray& A, int a, int b, int p, int pb, int bsv, bool bw)
{
	int len = b - a;
	if (len <= MIN_INSERT)
	{
		binaryInsert(A, a, b);
		return;
	}
	int s = len;
	while (s > MIN_INSERT) { s = (s - 1) / R + 1; }
	int i = a + s, j = a + R * s, pEnd = p + (s + 1) * (G + 1) + G;
	binaryInsert(A, a, i);
	for (int k = 0; k < s; ++k)
	{
		size_t x = static_cast<size_t>(a + k);
		size_t y = static_cast<size_t>(p + k * (G + 1) + G);
		A[x].get();
		A.swap(x, y);
	}

	while (i < b)
	{
		if (i == j)
		{
			s = i - a;
			int pEndNew = p + (s + 1) * (G + 1) + G;
			if (pEndNew > pb)
			{
				rebalance(A, a, i, p, pEnd, pb, bsv, bw);
				pEnd = pb;
				j = a;
			}
			else
			{
				rescatter(A, pEndNew, p, pEnd, bsv, bw);
				pEnd = pEndNew;
				j = a + (j - a) * R;
			}
		}
		int ele = A[i].get();
		int bLoc = randGapSearch(A, p + G, pEnd - (G + 1), ele);
		int loc = rightBinSearch(A, bLoc - G, bLoc, bsv, bw);
		if (loc == bLoc)
		{
			int rotP = -1;
			do { bLoc += G + 1; } 
			while (bLoc < pEnd && (rotP = rightBinSearch(A, bLoc - G, bLoc, bsv, bw)) == bLoc);
			if (bLoc == pb) { rebalance(A, a, i, p, pEnd, pb, bsv, bw); }
			else if (bLoc == pEnd)
			{
				int rotS = G / 2 + 1;
				shiftBW(A, loc - rotS, bLoc - (G + 1), bLoc - (G + 1) + rotS);
				pEnd += G + 1;
			}
			else
			{
				int rotS = bLoc - std::max(rotP, bLoc - (G + 1) / 2);
				shiftBW(A, loc - rotS, bLoc - rotS, bLoc);
			}
		}
		else
		{
			size_t k = static_cast<size_t>(i), l = static_cast<size_t>(loc);
			int ele = A[k].get();
			A.set(k, A[l]);
			++i;
			insertTo(A, ele, loc, rightBinSearch(A, bLoc - G, loc, ele, false));
		}
	}
	retrieve(A, b, p, pEnd, bsv, bw);
}

int medianOfThree(SortArray& A, int a, int m, int b)
{
	if (A[m] > A[a])
	{
		if (A[m] < A[b]) { return m; }
		if (A[a] > A[b]) { return a; }
		else { return b; }
	}
	else
	{
		if (A[m] > A[b]) { return m; }
		if (A[a] < A[b]) { return a; }
		else { return b; }
	}
}

int ninther(SortArray& A, int a, int b)
{
	int s = (b - a) / 9;
	int a1 = medianOfThree(A, a, a + s, a + 2 * s);
	int m1 = medianOfThree(A, a + 3 * s, a + 4 * s, a + 5 * s);
	int b1 = medianOfThree(A, a + 6 * s, a + 7 * s, a + 8 * s);
	return medianOfThree(A, a1, m1, b1);
}

int medianOfThreeNinthers(SortArray& A, int a, int b)
{
	int s = (b - a) / 3;
	int a1 = ninther(A, a, a + s);
	int m1 = ninther(A, a + s, a + 2 * s);
	int b1 = ninther(A, a + 2 * s, b);
	return medianOfThree(A, a1, m1, b1);
}

void quickLibrary(SortArray& A, int a, int b, int p, int pb, int minSize, int bsv, bool bw)
{
	if (b - a <= minSize)
	{
		librarySort(A, a, b, p, pb, bsv, bw);
		return;
	}
	int piv = A[medianOfThreeNinthers(A, a, b)];
	int i = a - 1, j = b;
	do
	{
		int val = A[0];
		do
		{
			++i;
			val = A[i].get();
		} 
		while (i < j && val < piv);
		val = A[0];
		do
		{
			--j;
			val = A[j].get();
		} 
		while (j >= i && val > piv);
		if (i < j) 
		{ 
			size_t x = static_cast<size_t>(i), y = static_cast<size_t>(j);
			A.swap(x, y); 
		}
		else { break; }
	} 
	while (true);
	quickLibrary(A, a, i, p, pb, minSize, bsv, bw);
	quickLibrary(A, i, b, p, pb, minSize, bsv, bw);
}

void QuickLibrarySort(SortArray& A)
{
	int length = static_cast<int>(A.size()), a = 0, b = length;
	while (b - a > MIN_INSERT)
	{
		int piv = A[medianOfThreeNinthers(A, a, b)];
		int i1 = a, i = a - 1, j = b, j1 = b;
		for (;;)
		{
			++i;
			while (i < j)
			{
				int ele = A[i].get();
				if (ele == piv) 
				{ 
					size_t x = static_cast<size_t>(i1), y = static_cast<size_t>(i);
					A.swap(x, y); ++i1; 
				}
				else if (ele > piv) { break; }
				++i;
			}
			--j;
			while (j > i)
			{
				int ele = A[j].get();
				if (ele == piv)
				{ 
					--j1;
					size_t x = static_cast<size_t>(j1), y = static_cast<size_t>(j);
					A.swap(x, y); 
				}
				else if (ele < piv) { break; }
				--j;
			}
			if (i < j)
			{
				size_t x = static_cast<size_t>(i), y = static_cast<size_t>(j);
				A.swap(x, y);
			}
			else
			{
				if (i1 == b) { return; }
				else if (j < i) { ++j; }
				while (i1 > a)
				{
					--i; --i1;
					size_t x = static_cast<size_t>(i), y = static_cast<size_t>(i1);
					A[x].get();
					A.swap(x, y);
				}
				while (j1 < b)
				{
					size_t x = static_cast<size_t>(j), y = static_cast<size_t>(j1);
					A[x].get();
					A.swap(x, y);
					++j; ++j1;
				}
				break;
			}
		}
		int left = i - a, right = b - j;
		if (left <= right)
		{
			right -= (right + 1) % (G + 1) + (G + 1);
			left = std::max(right / (G + 1) * R, MIN_INSERT);
			quickLibrary(A, a, i, j, j + right, left, piv, false);
			a = j;
		}
		else
		{
			left -= (left + 1) % (G + 1) + (G + 1);
			right = std::max(left / (G + 1) * R, MIN_INSERT);
			quickLibrary(A, j, b, a, a + left, right, piv, true);
			b = i;
		}
	}
	binaryInsert(A, a, b);
}
