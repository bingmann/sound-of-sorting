/*
 * MIT License
 *
 * Copyright (c) 2013 Andrey Astrelin
 * Copyright (c) 2020-2021 The Holy Grail Sort Project
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */

 /*
  * The Holy Grail Sort Project
  * Project Manager:      Summer Dragonfly
  * Project Contributors: 666666t
  *                       Anonymous0726
  *                       aphitorite
  *                       Control
  *                       dani_dlg
  *                       DeveloperSort
  *                       EilrahcF
  *                       Enver
  *                       Gaming32
  *                       lovebuny
  *                       Morwenn
  *                       MP
  *                       phoenixbound
  *                       Spex_guy
  *                       thatsOven
  *                       _fluffyy
  *
  * Special thanks to "The Studio" Discord community!
  */

#include "../SortAlgo.h"
#include <algorithm>
#include <functional>
#include <iterator>
#include <type_traits>
#include <utility>

namespace grailsort_detail
{
    // Credit to phoenixbound for this clever idea
    enum struct Subarray
    {
        LEFT,
        RIGHT
    };

    template<typename Compare>
    struct ThreeWayCompare
    {
        Compare compare;

        explicit ThreeWayCompare(Compare&& comp) :
            compare(std::forward<Compare>(comp))
        {}

        template<typename T, typename U>
        int operator()(T&& lhs, U&& rhs)
        {
            if (compare(lhs, rhs)) {
                return -1;
            }
            if (compare(rhs, lhs)) {
                return 1;
            }
            return 0;
        }
    };

    struct GrailSort
    {
        int currBlockLen;
        Subarray currBlockOrigin;

        template<typename RandomAccessIterator>
        static void BlockSwap(RandomAccessIterator array, int a, int b, int blockLen) {
            counted_swap_ranges(array + a, array + (a + blockLen), array + b);
        }

        // Swaps the order of two adjacent blocks whose lengths may or may not be equal.
        // Variant of the Gries-Mills algorithm, which is basically recursive block swaps
        template<typename RandomAccessIterator>
        static void Rotate(RandomAccessIterator array, int start, int leftLen, int rightLen) {
            while (leftLen > 0 && rightLen > 0) {
                if (leftLen <= rightLen) {
                    BlockSwap(array, start, start + leftLen, leftLen);
                    start += leftLen;
                    rightLen -= leftLen;
                }
                else {
                    BlockSwap(array, start + leftLen - rightLen, start + leftLen, rightLen);
                    leftLen -= rightLen;
                }
            }
        }

        // Variant of Insertion Sort that utilizes swaps instead of overwrites.
        // Also known as "Optimized Gnomesort".
        template<typename RandomAccessIterator, typename Compare>
        static void InsertSort(RandomAccessIterator array, int start, int length, Compare comp) {
            for (int item = 1; item < length; item++) {
                int left = start + item - 1;
                int right = start + item;

                while (left >= start && comp(array[left], array[right]) > 0) {
                    counted_iter_swap(array + left, array + right);
                    left--;
                    right--;
                }
            }
        }

        template<typename RandomAccessIterator, typename Compare, typename T>
        static int BinarySearchLeft(RandomAccessIterator array, int start, int length, const T& target, Compare comp) {
            int left = 0;
            int right = length;

            while (left < right) {
                // equivalent to (left + right) / 2 with added overflow protection
                int middle = left + ((right - left) / 2);

                if (comp(array[start + middle], target) < 0) {
                    left = middle + 1;
                }
                else {
                    right = middle;
                }
            }
            return left;
        }

        // Credit to Anonymous0726 for debugging
        template<typename RandomAccessIterator, typename Compare, typename T>
        static int BinarySearchRight(RandomAccessIterator array, int start, int length, const T& target, Compare comp) {
            int left = 0;
            int right = length;

            while (left < right) {
                // equivalent to (left + right) / 2 with added overflow protection
                int middle = left + ((right - left) / 2);
                if (comp(array[start + middle], target) > 0) {
                    right = middle;
                }
                else {
                    left = middle + 1;
                }
            }
            // OFF-BY-ONE BUG FIXED: used to be `return right - 1;`
            return right;
        }

        // cost: 2 * length + idealKeys^2 / 2
        template<typename RandomAccessIterator, typename Compare>
        static int CollectKeys(RandomAccessIterator array, int start, int length, int idealKeys, Compare comp) {
            int keysFound = 1; // by itself, the first item in the array is our first unique key
            int firstKey = 0; // the first item in the array is at the first position in the array
            int currKey = 1; // the index used for finding potentially unique items ("keys") in the array

            while (currKey < length && keysFound < idealKeys) {

                // Find the location in the key-buffer where our current key can be inserted in sorted order.
                // If the key at insertPos is equal to currKey, then currKey isn't unique and we move on.
                int insertPos = BinarySearchLeft(array, start + firstKey, keysFound, array[start + currKey], comp);

                // The second part of this conditional does the equal check we were just talking about; however,
                // if currKey is larger than everything in the key-buffer (meaning insertPos == keysFound),
                // then that also tells us it wasn't *equal* to anything in the key-buffer. Magic! :)
                if (insertPos == keysFound || comp(array[start + currKey], array[start + firstKey + insertPos]) != 0) {

                    // First, rotate the key-buffer over to currKey's immediate left...
                    // (this helps save a TON of swaps/writes!!!)
                    Rotate(array, start + firstKey, keysFound, currKey - (firstKey + keysFound));

                    // Update the new position of firstKey...
                    firstKey = currKey - keysFound;

                    // Then, "insertion sort" currKey to its spot in the key-buffer!
                    Rotate(array, start + firstKey + insertPos, keysFound - insertPos, 1);

                    // One step closer to idealKeys.
                    keysFound++;
                }
                // Move on and test the next key...
                currKey++;
            }

            // Bring however many keys we found back to the beginning of our array,
            // and return the number of keys collected.
            Rotate(array, start, firstKey, keysFound);
            return keysFound;
        }

        template<typename RandomAccessIterator, typename Compare>
        static void PairwiseSwaps(RandomAccessIterator array, int start, int length, Compare comp) {
            int index;
            for (index = 1; index < length; index += 2) {
                int  left = start + index - 1;
                int right = start + index;

                if (comp(array[left], array[right]) > 0) {
                    counted_swap(array[left - 2], array[right]);
                    counted_swap(array[right - 2], array[left]);
                }
                else {
                    counted_swap(array[left - 2], array[left]);
                    counted_swap(array[right - 2], array[right]);
                }
            }

            int left = start + index - 1;
            if (left < start + length) {
                counted_swap(array[left - 2], array[left]);
            }
        }

        template<typename RandomAccessIterator, typename Compare>
        static void PairwiseWrites(RandomAccessIterator array, int start, int length, Compare comp) {
            int index;
            for (index = 1; index < length; index += 2) {
                int  left = start + index - 1;
                int right = start + index;

                if (comp(array[left], array[right]) > 0) {
                    array[left - 2] = std::move(array[right]);
                    array[right - 2] = std::move(array[left]);
                }
                else {
                    array[left - 2] = std::move(array[left]);
                    array[right - 2] = std::move(array[right]);
                }
            }

            int left = start + index - 1;
            if (left < start + length) {
                array[left - 2] = std::move(array[left]);
            }
        }

        // array[buffer .. start - 1] <=> "scrolling buffer"
        //
        // "scrolling buffer" + array[start, middle - 1] + array[middle, end - 1]
        // --> array[buffer, buffer + end - 1] + "scrolling buffer"
        template<typename RandomAccessIterator, typename Compare>
        static void MergeForwards(RandomAccessIterator array, int start, int leftLen, int rightLen, int bufferOffset, Compare comp) {
            auto buffer = array + (start - bufferOffset);
            auto left = array + start;
            auto middle = array + (start + leftLen);
            auto right = middle;
            auto end = middle + rightLen;

            while (right != end) {
                if (left == middle || comp(*left, *right) > 0) {
                    counted_iter_swap(buffer, right);
                    ++right;
                }
                else {
                    counted_iter_swap(buffer, left);
                    ++left;
                }
                ++buffer;
            }

            while (left != middle) {
                counted_iter_swap(buffer, left);
                ++buffer;
                ++left;
            }
        }

        // credit to 666666t for thorough bug-checking/fixing
        template<typename RandomAccessIterator, typename Compare>
        static void MergeBackwards(RandomAccessIterator array, int start, int leftLen, int rightLen, int bufferOffset, Compare comp) {
            // used to be '= start'
            int    end = start - 1;
            // used to be '= start + leftLen - 1'
            int   left = end + leftLen;
            int middle = left;
            // OFF-BY-ONE BUG FIXED: used to be `int  right = middle + rightLen - 1;`
            int  right = middle + rightLen;
            // OFF-BY-ONE BUG FIXED: used to be `int buffer = right  + bufferOffset - 1;`
            int buffer = right + bufferOffset;

            // used to be 'left >= end'
            while (left > end) {
                if (right == middle || comp(array[left], array[right]) > 0) {

                    counted_swap(array[buffer], array[left]);
                    left--;
                }
                else {
                    counted_swap(array[buffer], array[right]);
                    right--;
                }
                buffer--;
            }

            if (right != buffer) {
                while (right > middle) {
                    counted_swap(array[buffer], array[right]);
                    buffer--;
                    right--;
                }
            }
        }

        // array[buffer .. start - 1] <=> "free space"
        //
        // "free space" + array[start, middle - 1] + array[middle, end - 1]
        // --> array[buffer, buffer + end - 1] + "free space"
        //
        // FUNCTION RENAMED: More consistent with "out-of-place" being at the end
        template<typename RandomAccessIterator, typename Compare>
        static void MergeOutOfPlace(RandomAccessIterator array, int start, int leftLen, int rightLen, int bufferOffset, Compare comp) {
            int buffer = start - bufferOffset;
            int   left = start;
            int middle = start + leftLen;
            int  right = middle;
            int    end = middle + rightLen;

            while (right < end) {
                if (left == middle || comp(array[left], array[right]) > 0) {
                    array[buffer] = std::move(array[right]);
                    right++;
                }
                else {
                    array[buffer] = std::move(array[left]);
                    left++;
                }
                buffer++;
            }

            if (buffer != left) {
                while (left < middle) {
                    array[buffer] = std::move(array[left]);
                    buffer++;
                    left++;
                }
            }
        }

        template<typename RandomAccessIterator, typename Compare>
        static void BuildInPlace(RandomAccessIterator array, int start, int length, int currentLen, int bufferLen, Compare comp) {
            for (int mergeLen = currentLen; mergeLen < bufferLen; mergeLen *= 2) {
                int fullMerge = 2 * mergeLen;

                int mergeIndex;
                int mergeEnd = start + length - fullMerge;
                int bufferOffset = mergeLen;

                for (mergeIndex = start; mergeIndex <= mergeEnd; mergeIndex += fullMerge) {
                    MergeForwards(array, mergeIndex, mergeLen, mergeLen, bufferOffset, comp);
                }

                int leftOver = length - (mergeIndex - start);

                if (leftOver > mergeLen) {
                    MergeForwards(array, mergeIndex, mergeLen, leftOver - mergeLen, bufferOffset, comp);
                }
                else {
                    Rotate(array, mergeIndex - mergeLen, mergeLen, leftOver);
                }

                start -= mergeLen;
            }

            int fullMerge = 2 * bufferLen;
            int lastBlock = length % fullMerge;
            int lastOffset = start + length - lastBlock;

            if (lastBlock <= bufferLen) {
                Rotate(array, lastOffset, lastBlock, bufferLen);
            }
            else {
                MergeBackwards(array, lastOffset, bufferLen, lastBlock - bufferLen, bufferLen, comp);
            }

            for (int mergeIndex = lastOffset - fullMerge; mergeIndex >= start; mergeIndex -= fullMerge) {
                MergeBackwards(array, mergeIndex, bufferLen, bufferLen, bufferLen, comp);
            }
        }

        template<typename RandomAccessIterator, typename BufferIterator, typename Compare>
        void BuildOutOfPlace(RandomAccessIterator array, int start, int length, int bufferLen, int extLen,
            BufferIterator extBuffer, Compare comp) {
            std::move(array + (start - extLen), array + start, extBuffer);

            PairwiseWrites(array, start, length, comp);
            start -= 2;

            int mergeLen;
            for (mergeLen = 2; mergeLen < extLen; mergeLen *= 2) {
                int fullMerge = 2 * mergeLen;

                int mergeIndex;
                int mergeEnd = start + length - fullMerge;
                int bufferOffset = mergeLen;

                for (mergeIndex = start; mergeIndex <= mergeEnd; mergeIndex += fullMerge) {
                    MergeOutOfPlace(array, mergeIndex, mergeLen, mergeLen, bufferOffset, comp);
                }

                int leftOver = length - (mergeIndex - start);

                if (leftOver > mergeLen) {
                    MergeOutOfPlace(array, mergeIndex, mergeLen, leftOver - mergeLen, bufferOffset, comp);
                }
                else {
                    // MINOR CHANGE: Used to be a loop; much clearer now
                    std::move(array + mergeIndex, array + (mergeIndex + leftOver), array + (mergeIndex - mergeLen));
                }

                start -= mergeLen;
            }

            std::move(extBuffer, extBuffer + extLen, array + (start + length));
            BuildInPlace(array, start, length, mergeLen, bufferLen, comp);
        }

        // build blocks of length 'bufferLen'
        // input: [start - mergeLen, start - 1] elements are buffer
        // output: first 'bufferLen' elements are buffer, blocks (2 * bufferLen) and last subblock sorted
        template<typename RandomAccessIterator, typename BufferIterator, typename Compare>
        void BuildBlocks(RandomAccessIterator array, int start, int length, int bufferLen,
            BufferIterator extBuffer, int extBufferLen, Compare comp) {
            if (extBufferLen != 0) {
                int extLen;

                if (bufferLen < extBufferLen) {
                    extLen = bufferLen;
                }
                else {
                    // max power of 2 -- just in case
                    extLen = 1;
                    while ((extLen * 2) <= extBufferLen) {
                        extLen *= 2;
                    }
                }

                BuildOutOfPlace(array, start, length, bufferLen, extLen, extBuffer, comp);
            }
            else {
                PairwiseSwaps(array, start, length, comp);
                BuildInPlace(array, start - 2, length, 2, bufferLen, comp);
            }
        }

        // Returns the final position of 'medianKey'
        template<typename RandomAccessIterator, typename Compare>
        static int BlockSelectSort(RandomAccessIterator array, int firstKey, int start, int medianKey, int blockCount, int blockLen, Compare comp) {
            for (int firstBlock = 0; firstBlock < blockCount; firstBlock++) {
                int selectBlock = firstBlock;

                for (int currBlock = firstBlock + 1; currBlock < blockCount; currBlock++) {
                    int compare = comp(array[start + (currBlock * blockLen)],
                        array[start + (selectBlock * blockLen)]);

                    if (compare < 0 || (compare == 0 && comp(array[firstKey + currBlock],
                        array[firstKey + selectBlock]) < 0)) {
                        selectBlock = currBlock;
                    }
                }

                if (selectBlock != firstBlock) {
                    // Swap the left and right selected blocks...
                    BlockSwap(array, start + (firstBlock * blockLen), start + (selectBlock * blockLen), blockLen);

                    // Swap the keys...
                    counted_swap(array[firstKey + firstBlock], array[firstKey + selectBlock]);

                    // ...and follow the 'medianKey' if it was swapped

                    // ORIGINAL LOC: if(midkey==u-1 || midkey==p) midkey^=(u-1)^p;
                    // MASSIVE, MASSIVE credit to lovebuny for figuring this one out!
                    if (medianKey == firstBlock) {
                        medianKey = selectBlock;
                    }
                    else if (medianKey == selectBlock) {
                        medianKey = firstBlock;
                    }
                }
            }

            return medianKey;
        }

        // Swaps Grailsort's "scrolling buffer" from the right side of the array all the way back to 'start'.
        // Costs O(n) swaps (amortized).
        //
        // OFF-BY-ONE BUG FIXED: used to be `int index = start + resetLen`; credit to 666666t for debugging
        // RESTRUCTED, BETTER NAMES: 'resetLen' is now 'length' and 'bufferLen' is now 'bufferOffset'
        template<typename RandomAccessIterator>
        static void InPlaceBufferReset(RandomAccessIterator array, int start, int length, int bufferOffset) {
            int  index = start + length - 1;
            int buffer = index - bufferOffset;
            ArrayItem pos = ArrayItem(index);
            while (index >= start) {
                pos.get();
                counted_swap(array[index], array[buffer]);
                index--; buffer--;
                pos = ArrayItem(index);
            }
        }

        // Shifts entire array over 'bufferOffset' spaces to move the out-of-place merging buffer back to
        // the beginning of the array.
        // Costs O(n) writes (amortized).
        //
        // OFF-BY-ONE BUG FIXED: used to be `int index = start + resetLen`; credit to 666666t for debugging
        // RESTRUCTED, BETTER NAMES: 'resetLen' is now 'length' and 'bufferLen' is now 'bufferOffset'
        template<typename RandomAccessIterator>
        static void OutOfPlaceBufferReset(RandomAccessIterator array, int start, int length, int bufferOffset) {
            int  index = start + length - 1;
            int buffer = index - bufferOffset;
            ArrayItem pos = ArrayItem(index);
            while (index >= start) {
                pos.get();
                array[index] = std::move(array[buffer]);
                index--; buffer--;
                pos = ArrayItem(index);
            }
        }

        // Rewinds Grailsort's "scrolling buffer" to the left any items belonging to the left subarray block
        // left over by a "smart merge". This is used to continue an ongoing merge that has run out of buffer space.
        // Costs O(sqrt n) swaps (amortized) in the *absolute* worst-case.
        //
        // NAMING IMPROVED: the left over items are in the middle of the merge while the buffer is at the end
        // BETTER ORDER-OF-OPERATIONS, NAMING IMPROVED: the left over items (now called 'leftBlock') are in the
        //                                              middle of the merge while the buffer is at the end
        template<typename RandomAccessIterator>
        static void InPlaceBufferRewind(RandomAccessIterator array, int start, int leftBlock, int buffer) {
            while (leftBlock >= start) {
                counted_swap(array[buffer], array[leftBlock]);
                leftBlock--;
                buffer--;
            }
        }

        // Rewinds Grailsort's out-of-place buffer to the left of any items belonging to the left subarray block
        // left over by a "smart merge". This is used to continue an ongoing merge that has run out of buffer space.
        // Costs O(sqrt n) writes (amortized) in the *absolute* worst-case.
        //
        // INCORRECT ORDER OF PARAMETERS BUG FIXED: `leftOvers` should be the middle, and `buffer` should be the end
        // BETTER ORDER, INCORRECT ORDER OF PARAMETERS BUG FIXED: `leftOvers` (now called 'leftBlock') should be
        //                                                        the middle, and `buffer` should be the end
        template<typename RandomAccessIterator>
        static void OutOfPlaceBufferRewind(RandomAccessIterator array, int start, int leftBlock, int buffer) {
            while (leftBlock >= start) {
                array[buffer] = std::move(array[leftBlock]);
                leftBlock--;
                buffer--;
            }
        }

        template<typename RandomAccessIterator, typename Compare>
        static Subarray GetSubarray(RandomAccessIterator array, int currKey, int medianKey, Compare comp) {
            if (comp(array[currKey], array[medianKey]) < 0) {
                return Subarray::LEFT;
            }
            else {
                return Subarray::RIGHT;
            }
        }

        // FUNCTION RE-RENAMED: last/final left blocks are used to calculate the length of the final merge
        template<typename RandomAccessIterator, typename Compare>
        static int CountLastMergeBlocks(RandomAccessIterator array, int offset, int blockCount, int blockLen, Compare comp) {
            int blocksToMerge = 0;

            int lastRightFrag = offset + (blockCount * blockLen);
            int   prevLeftBlock = lastRightFrag - blockLen;

            while (blocksToMerge < blockCount && comp(array[lastRightFrag], array[prevLeftBlock]) < 0) {
                blocksToMerge++;
                prevLeftBlock -= blockLen;
            }

            return blocksToMerge;
        }

        template<typename RandomAccessIterator, typename Compare>
        void SmartMerge(RandomAccessIterator array, int start, int leftLen, Subarray leftOrigin, int rightLen, int bufferOffset, Compare comp) {
            auto buffer = array + (start - bufferOffset);
            auto left = array + start;
            auto middle = left + leftLen;
            auto right = middle;
            auto end = middle + rightLen;

            if (leftOrigin == Subarray::LEFT) {
                while (left != middle && right != end) {
                    if (comp(*left, *right) <= 0) {
                        counted_iter_swap(buffer, left);
                        ++left;
                    }
                    else {
                        counted_iter_swap(buffer, right);
                        ++right;
                    }
                    ++buffer;
                }
            }
            else {
                while (left != middle && right != end) {
                    if (comp(*left, *right) < 0) {
                        counted_iter_swap(buffer, left);
                        ++left;
                    }
                    else {
                        counted_iter_swap(buffer, right);
                        ++right;
                    }
                    ++buffer;
                }
            }

            if (left < middle) {
                this->currBlockLen = middle - left;
                // UPDATED ARGUMENTS: 'middle' and 'end' now 'middle - 1' and 'end - 1'
                InPlaceBufferRewind(array, left - array, (middle - array) - 1, (end - array) - 1);
            }
            else {
                this->currBlockLen = end - right;
                if (leftOrigin == Subarray::LEFT) {
                    this->currBlockOrigin = Subarray::RIGHT;
                }
                else {
                    this->currBlockOrigin = Subarray::LEFT;
                }
            }
        }

        // MINOR CHANGE: better naming -- 'insertPos' is now 'mergeLen' -- and "middle" calculation simplified
        template<typename RandomAccessIterator, typename Compare>
        void SmartLazyMerge(RandomAccessIterator array, int start, int leftLen, Subarray leftOrigin, int rightLen, Compare comp) {
            int middle = start + leftLen;

            if (leftOrigin == Subarray::LEFT) {
                if (comp(array[middle - 1], array[middle]) > 0) {
                    while (leftLen != 0) {
                        int mergeLen = BinarySearchLeft(array, middle, rightLen, array[start], comp);

                        if (mergeLen != 0) {
                            Rotate(array, start, leftLen, mergeLen);
                            start += mergeLen;
                            rightLen -= mergeLen;
                        }

                        middle += mergeLen;

                        if (rightLen == 0) {
                            this->currBlockLen = leftLen;
                            return;
                        }
                        else {
                            do {
                                start++;
                                leftLen--;
                            } while (leftLen != 0 && comp(array[start], array[middle]) <= 0);
                        }
                    }
                }
            }
            else {
                if (comp(array[middle - 1], array[middle]) >= 0) {
                    while (leftLen != 0) {
                        int mergeLen = BinarySearchRight(array, start + leftLen, rightLen, array[start], comp);

                        if (mergeLen != 0) {
                            Rotate(array, start, leftLen, mergeLen);
                            start += mergeLen;
                            rightLen -= mergeLen;
                        }

                        middle += mergeLen;

                        if (rightLen == 0) {
                            this->currBlockLen = leftLen;
                            return;
                        }
                        else {
                            do {
                                start++;
                                leftLen--;
                            } while (leftLen != 0 && comp(array[start], array[middle]) < 0);
                        }
                    }
                }
            }

            this->currBlockLen = rightLen;
            if (leftOrigin == Subarray::LEFT) {
                this->currBlockOrigin = Subarray::RIGHT;
            }
            else {
                this->currBlockOrigin = Subarray::LEFT;
            }
        }

        // FUNCTION RENAMED: more consistent with other "out-of-place" merges
        template<typename RandomAccessIterator, typename Compare>
        void SmartMergeOutOfPlace(RandomAccessIterator array, int start, int leftLen, Subarray leftOrigin, int rightLen, int bufferOffset, Compare comp) {
            int buffer = start - bufferOffset;
            int   left = start;
            int middle = start + leftLen;
            int  right = middle;
            int    end = middle + rightLen;

            if (leftOrigin == Subarray::LEFT) {
                while (left < middle && right < end) {
                    if (comp(array[left], array[right]) <= 0) {
                        array[buffer] = std::move(array[left]);
                        left++;
                    }
                    else {
                        array[buffer] = std::move(array[right]);
                        right++;
                    }
                    buffer++;
                }
            }
            else {
                while (left < middle && right < end) {
                    if (comp(array[left], array[right]) < 0) {
                        array[buffer] = std::move(array[left]);
                        left++;
                    }
                    else {
                        array[buffer] = std::move(array[right]);
                        right++;
                    }
                    buffer++;
                }
            }

            if (left < middle) {
                this->currBlockLen = middle - left;
                // UPDATED ARGUMENTS: 'middle' and 'end' now 'middle - 1' and 'end - 1'
                OutOfPlaceBufferRewind(array, left, middle - 1, end - 1);
            }
            else {
                this->currBlockLen = end - right;
                if (leftOrigin == Subarray::LEFT) {
                    this->currBlockOrigin = Subarray::RIGHT;
                }
                else {
                    this->currBlockOrigin = Subarray::LEFT;
                }
            }
        }

        // Credit to Anonymous0726 for better variable names such as "nextBlock"
        // Also minor change: removed unnecessary "currBlock = nextBlock" lines
        template<typename RandomAccessIterator, typename Compare>
        void MergeBlocks(RandomAccessIterator array, int firstKey, int medianKey, int start, int blockCount, int blockLen,
            int lastMergeBlocks, int lastLen, Compare comp) {
            int buffer;

            int currBlock;
            int nextBlock = start + blockLen;

            this->currBlockLen = blockLen;
            this->currBlockOrigin = GetSubarray(array, firstKey, medianKey, comp);

            Subarray nextBlockOrigin;
            for (int keyIndex = 1; keyIndex < blockCount; keyIndex++, nextBlock += blockLen) {
                currBlock = nextBlock - this->currBlockLen;
                nextBlockOrigin = GetSubarray(array, firstKey + keyIndex, medianKey, comp);

                if (nextBlockOrigin == this->currBlockOrigin) {
                    buffer = currBlock - blockLen;

                    BlockSwap(array, buffer, currBlock, this->currBlockLen);
                    this->currBlockLen = blockLen;
                }
                else {
                    SmartMerge(array, currBlock, this->currBlockLen, this->currBlockOrigin, blockLen, blockLen, comp);
                }
            }

            currBlock = nextBlock - this->currBlockLen;
            buffer = currBlock - blockLen;

            if (lastLen != 0) {
                if (this->currBlockOrigin == Subarray::RIGHT) {
                    BlockSwap(array, buffer, currBlock, this->currBlockLen);

                    currBlock = nextBlock;
                    this->currBlockLen = blockLen * lastMergeBlocks;
                    this->currBlockOrigin = Subarray::LEFT;
                }
                else {
                    this->currBlockLen += blockLen * lastMergeBlocks;
                }

                MergeForwards(array, currBlock, this->currBlockLen, lastLen, blockLen, comp);
            }
            else {
                BlockSwap(array, buffer, currBlock, this->currBlockLen);
            }
        }

        template<typename RandomAccessIterator, typename Compare>
        void LazyMergeBlocks(RandomAccessIterator array, int firstKey, int medianKey, int start, int blockCount, int blockLen,
            int lastMergeBlocks, int lastLen, Compare comp) {
            int currBlock;
            int nextBlock = start + blockLen;

            this->currBlockLen = blockLen;
            this->currBlockOrigin = GetSubarray(array, firstKey, medianKey, comp);

            Subarray nextBlockOrigin;
            for (int keyIndex = 1; keyIndex < blockCount; keyIndex++, nextBlock += blockLen) {
                currBlock = nextBlock - this->currBlockLen;

                nextBlockOrigin = GetSubarray(array, firstKey + keyIndex, medianKey, comp);

                if (nextBlockOrigin == this->currBlockOrigin) {
                    this->currBlockLen = blockLen;
                }
                else {
                    // These checks were included in the original code... but why???
                    if (blockLen != 0 && this->currBlockLen != 0) {
                        SmartLazyMerge(array, currBlock, this->currBlockLen, this->currBlockOrigin, blockLen, comp);
                    }
                }
            }

            currBlock = nextBlock - this->currBlockLen;

            if (lastLen != 0) {
                if (this->currBlockOrigin == Subarray::RIGHT) {
                    currBlock = nextBlock;
                    this->currBlockLen = blockLen * lastMergeBlocks;
                    this->currBlockOrigin = Subarray::LEFT;
                }
                else {
                    this->currBlockLen += blockLen * lastMergeBlocks;
                }

                LazyMerge(array, currBlock, this->currBlockLen, lastLen, comp);
            }
        }

        template<typename RandomAccessIterator, typename Compare>
        void MergeBlocksOutOfPlace(RandomAccessIterator array, int firstKey, int medianKey, int start, int blockCount, int blockLen,
            int lastMergeBlocks, int lastLen, Compare comp) {
            int buffer;
            int currBlock;
            int nextBlock = start + blockLen;

            this->currBlockLen = blockLen;
            this->currBlockOrigin = GetSubarray(array, firstKey, medianKey, comp);

            Subarray nextBlockOrigin;
            for (int keyIndex = 1; keyIndex < blockCount; keyIndex++, nextBlock += blockLen) {
                currBlock = nextBlock - this->currBlockLen;
                nextBlockOrigin = GetSubarray(array, firstKey + keyIndex, medianKey, comp);

                if (nextBlockOrigin == this->currBlockOrigin) {
                    buffer = currBlock - blockLen;
                    std::move(array + currBlock, array + (currBlock + this->currBlockLen), array + buffer);
                    this->currBlockLen = blockLen;
                }
                else {
                    SmartMergeOutOfPlace(array, currBlock, this->currBlockLen, this->currBlockOrigin, blockLen, blockLen, comp);
                }
            }

            currBlock = nextBlock - this->currBlockLen;
            buffer = currBlock - blockLen;

            if (lastLen != 0) {
                if (this->currBlockOrigin == Subarray::RIGHT) {
                    std::move(array + currBlock, array + (currBlock + this->currBlockLen), array + buffer);
                    currBlock = nextBlock;
                    this->currBlockLen = blockLen * lastMergeBlocks;
                    this->currBlockOrigin = Subarray::LEFT;
                }
                else {
                    this->currBlockLen += blockLen * lastMergeBlocks;
                }

                MergeOutOfPlace(array, currBlock, this->currBlockLen, lastLen, blockLen, comp);
            }
            else {
                std::move(array + currBlock, array + (currBlock + this->currBlockLen), array + buffer);
            }
        }

        //TODO: Double-check "Merge Blocks" arguments
        template<typename RandomAccessIterator, typename Compare>
        void CombineInPlace(RandomAccessIterator array, int firstKey, int start, int length, int subarrayLen, int blockLen,
            int mergeCount, int lastSubarrays, bool buffer, Compare comp) {
            int fullMerge = 2 * subarrayLen;
            // SLIGHT OPTIMIZATION: 'blockCount' only needs to be calculated once for regular merges
            int blockCount = fullMerge / blockLen;

            for (int mergeIndex = 0; mergeIndex < mergeCount; mergeIndex++) {
                int offset = start + (mergeIndex * fullMerge);

                InsertSort(array, firstKey, blockCount, comp);

                // INCORRECT PARAMETER BUG FIXED: `block select sort` should be using `offset`, not `start`
                int medianKey = subarrayLen / blockLen;
                medianKey = BlockSelectSort(array, firstKey, offset, medianKey, blockCount, blockLen, comp);

                if (buffer) {
                    MergeBlocks(array, firstKey, firstKey + medianKey, offset, blockCount, blockLen, 0, 0, comp);
                }
                else {
                    LazyMergeBlocks(array, firstKey, firstKey + medianKey, offset, blockCount, blockLen, 0, 0, comp);
                }
            }

            // INCORRECT CONDITIONAL/PARAMETER BUG FIXED: Credit to 666666t for debugging.
            if (lastSubarrays != 0) {
                int offset = start + (mergeCount * fullMerge);
                blockCount = lastSubarrays / blockLen;

                InsertSort(array, firstKey, blockCount + 1, comp);

                // INCORRECT PARAMETER BUG FIXED: `block select sort` should be using `offset`, not `start`
                int medianKey = subarrayLen / blockLen;
                medianKey = BlockSelectSort(array, firstKey, offset, medianKey, blockCount, blockLen, comp);

                // MISSING BOUNDS CHECK BUG FIXED: `lastFragment` *can* be 0 if the last two subarrays are evenly
                //                                 divided into blocks. This prevents Grailsort from going out-of-bounds.
                int lastFragment = lastSubarrays - (blockCount * blockLen);
                int lastMergeBlocks;
                if (lastFragment != 0) {
                    lastMergeBlocks = CountLastMergeBlocks(array, offset, blockCount, blockLen, comp);
                }
                else {
                    lastMergeBlocks = 0;
                }

                int smartMerges = blockCount - lastMergeBlocks;

                //TODO: Double-check if this micro-optimization works correctly like the original
                if (smartMerges == 0) {
                    int leftLen = lastMergeBlocks * blockLen;

                    // INCORRECT PARAMETER BUG FIXED: these merges should be using `offset`, not `start`
                    if (buffer) {
                        MergeForwards(array, offset, leftLen, lastFragment, blockLen, comp);
                    }
                    else {
                        LazyMerge(array, offset, leftLen, lastFragment, comp);
                    }
                }
                else {
                    if (buffer) {
                        MergeBlocks(array, firstKey, firstKey + medianKey, offset,
                            smartMerges, blockLen, lastMergeBlocks, lastFragment, comp);
                    }
                    else {
                        LazyMergeBlocks(array, firstKey, firstKey + medianKey, offset,
                            smartMerges, blockLen, lastMergeBlocks, lastFragment, comp);
                    }
                }
            }

            if (buffer) {
                InPlaceBufferReset(array, start, length, blockLen);
            }
        }

        template<typename RandomAccessIterator, typename BufferIterator, typename Compare>
        void CombineOutOfPlace(RandomAccessIterator array, int firstKey, int start, int length, int subarrayLen, int blockLen,
            int mergeCount, int lastSubarrays, BufferIterator extBuffer, Compare comp) {
            std::move(array + (start - blockLen), array + start, extBuffer);

            int fullMerge = 2 * subarrayLen;
            // SLIGHT OPTIMIZATION: 'blockCount' only needs to be calculated once for regular merges
            int blockCount = fullMerge / blockLen;

            for (int mergeIndex = 0; mergeIndex < mergeCount; mergeIndex++) {
                int offset = start + (mergeIndex * fullMerge);

                InsertSort(array, firstKey, blockCount, comp);

                // INCORRECT PARAMETER BUG FIXED: `block select sort` should be using `offset`, not `start`
                int medianKey = subarrayLen / blockLen;
                medianKey = BlockSelectSort(array, firstKey, offset, medianKey, blockCount, blockLen, comp);

                MergeBlocksOutOfPlace(array, firstKey, firstKey + medianKey, offset,
                    blockCount, blockLen, 0, 0, comp);
            }

            // INCORRECT CONDITIONAL/PARAMETER BUG FIXED: Credit to 666666t for debugging.
            if (lastSubarrays != 0) {
                int offset = start + (mergeCount * fullMerge);
                blockCount = lastSubarrays / blockLen;

                InsertSort(array, firstKey, blockCount + 1, comp);

                // INCORRECT PARAMETER BUG FIXED: `block select sort` should be using `offset`, not `start`
                int medianKey = subarrayLen / blockLen;
                medianKey = BlockSelectSort(array, firstKey, offset, medianKey, blockCount, blockLen, comp);

                // MISSING BOUNDS CHECK BUG FIXED: `lastFragment` *can* be 0 if the last two subarrays are evenly
                //                                 divided into blocks. This prevents Grailsort from going out-of-bounds.
                int lastFragment = lastSubarrays - (blockCount * blockLen);
                int lastMergeBlocks;
                if (lastFragment != 0) {
                    lastMergeBlocks = CountLastMergeBlocks(array, offset, blockCount, blockLen, comp);
                }
                else {
                    lastMergeBlocks = 0;
                }

                int smartMerges = blockCount - lastMergeBlocks;

                if (smartMerges == 0) {
                    // MINOR CHANGE: renamed for consistency (used to be 'leftLength')
                    int leftLen = lastMergeBlocks * blockLen;

                    // INCORRECT PARAMETER BUG FIXED: this merge should be using `offset`, not `start`
                    MergeOutOfPlace(array, offset, leftLen, lastFragment, blockLen, comp);
                }
                else {
                    MergeBlocksOutOfPlace(array, firstKey, firstKey + medianKey, offset,
                        smartMerges, blockLen, lastMergeBlocks, lastFragment, comp);
                }
            }

            OutOfPlaceBufferReset(array, start, length, blockLen);
            std::move(extBuffer, extBuffer + blockLen, array + (start - blockLen));
        }

        // Keys are on the left side of array. Blocks of length 'subarrayLen' combined. We'll combine them in pairs
        // 'subarrayLen' is a power of 2. (2 * subarrayLen / blockLen) keys are guaranteed
        //
        // IMPORTANT RENAMES: 'lastSubarray' is now 'lastSubarrays' because it includes the length of the last left
        //                    subarray AND last right subarray (if there is a right subarray at all).
        //
        //                    *Please also check everything surrounding 'if(lastSubarrays != 0)' inside
        //                    'combine in-/out-of-place' methods for other renames!!*
        template<typename RandomAccessIterator, typename BufferIterator, typename Compare>
        void CombineBlocks(RandomAccessIterator array, int firstKey, int start, int length, int subarrayLen, int blockLen,
            bool buffer, BufferIterator extBuffer, int extBufferLen, Compare comp) {
            int    fullMerge = 2 * subarrayLen;
            int   mergeCount = length / fullMerge;
            int lastSubarrays = length - (fullMerge * mergeCount);

            if (lastSubarrays <= subarrayLen) {
                length -= lastSubarrays;
                lastSubarrays = 0;
            }

            // INCOMPLETE CONDITIONAL BUG FIXED: In order to combine blocks out-of-place, we must check if a full-sized
            //                                   block fits into our external buffer.
            if (buffer && blockLen <= extBufferLen) {
                CombineOutOfPlace(array, firstKey, start, length, subarrayLen, blockLen, mergeCount, lastSubarrays,
                    extBuffer, comp);
            }
            else {
                CombineInPlace(array, firstKey, start, length, subarrayLen, blockLen,
                    mergeCount, lastSubarrays, buffer, comp);
            }
        }

        // "Classic" in-place merge sort using binary searches and rotations
        //
        // cost: min(leftLen, rightLen)^2 + max(leftLen, rightLen)
        // MINOR CHANGES: better naming -- 'insertPos' is now 'mergeLen' -- and "middle"/"end" calculations simplified
        template<typename RandomAccessIterator, typename Compare>
        static void LazyMerge(RandomAccessIterator array, int start, int leftLen, int rightLen, Compare comp) {
            if (leftLen < rightLen) {
                int middle = start + leftLen;

                while (leftLen != 0) {
                    int mergeLen = BinarySearchLeft(array, middle, rightLen, array[start], comp);

                    if (mergeLen != 0) {
                        Rotate(array, start, leftLen, mergeLen);
                        start += mergeLen;
                        rightLen -= mergeLen;
                    }

                    middle += mergeLen;

                    if (rightLen == 0) {
                        break;
                    }
                    else {
                        do {
                            start++;
                            leftLen--;
                        } while (leftLen != 0 && comp(array[start], array[middle]) <= 0);
                    }
                }
            }
            // INDEXING BUG FIXED: Credit to Anonymous0726 for debugging.
            else {
                int end = start + leftLen + rightLen - 1;


                while (rightLen != 0) {
                    int mergeLen = BinarySearchRight(array, start, leftLen, array[end], comp);

                    if (mergeLen != leftLen) {
                        Rotate(array, start + mergeLen, leftLen - mergeLen, rightLen);
                        end -= leftLen - mergeLen;
                        leftLen = mergeLen;
                    }

                    if (leftLen == 0) {
                        break;
                    }
                    else {
                        int middle = start + leftLen;
                        do {
                            rightLen--;
                            end--;
                        } while (rightLen != 0 && comp(array[middle - 1], array[end]) <= 0);
                    }
                }
            }
        }

        template<typename RandomAccessIterator, typename Compare>
        static void LazyStableSort(RandomAccessIterator array, int start, int length, Compare comp) {
            for (int index = 1; index < length; index += 2) {
                int  left = start + index - 1;
                int right = start + index;

                if (comp(array[left], array[right]) > 0) {
                    counted_swap(array[left], array[right]);
                }
            }
            for (int mergeLen = 2; mergeLen < length; mergeLen *= 2) {
                int fullMerge = 2 * mergeLen;

                int mergeIndex;
                int mergeEnd = length - fullMerge;

                for (mergeIndex = 0; mergeIndex <= mergeEnd; mergeIndex += fullMerge) {
                    LazyMerge(array, start + mergeIndex, mergeLen, mergeLen, comp);
                }

                int leftOver = length - mergeIndex;
                if (leftOver > mergeLen) {
                    LazyMerge(array, start + mergeIndex, mergeLen, leftOver - mergeLen, comp);
                }
            }
        }

        // Calculates the minimum between numKeys and cbrt(2 * subarrayLen * keysFound).
        // Math will be further explained later, but just like in CommonSort, this
        // loop is rendered completely useless by the scrolling buffer optimization;
        // minKeys will always equal numKeys.
        //
        // Code still here for preservation purposes.
        /*static int CalcMinKeys(int numKeys, long halfSubarrKeys) {
            int minKeys = 1;
            while(minKeys < numKeys && halfSubarrKeys != 0) {
                minKeys *= 2;
                halfSubarrKeys /= 8;
            }
            return minKeys;
        }*/

        template<typename RandomAccessIterator, typename BufferIterator, typename Compare>
        void CommonSort(RandomAccessIterator array, int start, int length, BufferIterator extBuf, int extBufLen, Compare comp) {
            if (length < 16) {
                InsertSort(array, start, length, comp);
                return;
            }

            BufferIterator extBuffer{};
            int extBufferLen = 0;

            int blockLen = 1;

            // find the smallest power of two greater than or equal to
            // the square root of the input's length
            while ((blockLen * blockLen) < length) {
                blockLen *= 2;
            }

            // '((a - 1) / b) + 1' is actually a clever and very efficient
            // formula for the ceiling of (a / b)
            //
            // credit to Anonymous0726 for figuring this out!
            int keyLen = ((length - 1) / blockLen) + 1;

            // Grailsort is hoping to find `2 * sqrt(n)` unique items
            // throughout the array
            int idealKeys = keyLen + blockLen;

            //TODO: Clean up `start +` offsets
            int keysFound = CollectKeys(array, start, length, idealKeys, comp);

            bool idealBuffer;
            if (keysFound < idealKeys) {
                if (keysFound == 1) return;
                if (keysFound < 4) {
                    // GRAILSORT STRATEGY 3 -- No block swaps or scrolling buffer; resort to Lazy Stable Sort
                    LazyStableSort(array, start, length, comp);
                    return;
                }
                else {
                    // GRAILSORT STRATEGY 2 -- Block swaps with small scrolling buffer and/or lazy merges
                    keyLen = blockLen;
                    blockLen = 0;
                    idealBuffer = false;

                    while (keyLen > keysFound) {
                        keyLen /= 2;
                    }
                }
            }
            else {
                // GRAILSORT STRATEGY 1 -- Block swaps with scrolling buffer
                idealBuffer = true;
            }

            int bufferEnd = blockLen + keyLen;
            int subarrayLen;
            if (idealBuffer) {
                subarrayLen = blockLen;
            }
            else {
                subarrayLen = keyLen;
            }

            if (idealBuffer && extBufLen != 0) {
                // GRAILSORT + EXTRA SPACE
                extBuffer = extBuf;
                extBufferLen = extBufLen;
            }

            BuildBlocks(array, start + bufferEnd, length - bufferEnd, subarrayLen,
                extBuffer, extBufferLen, comp);

            while ((length - bufferEnd) > (2 * subarrayLen)) {
                subarrayLen *= 2;

                int currBlockLen = blockLen;
                bool scrollingBuffer = idealBuffer;

                // Huge credit to Anonymous0726, phoenixbound, and DeveloperSort for their tireless efforts
                // towards deconstructing this math.
                if (!idealBuffer) {
                    int keyBuffer = keyLen / 2;
                    // TODO: Rewrite explanation for this math
                    if (keyBuffer >= (2 * subarrayLen) / keyBuffer) {
                        currBlockLen = keyBuffer;
                        scrollingBuffer = true;
                    }
                    else {
                        // This is a very recent discovery, and the math will be spelled out later, but this
                        // "minKeys" calculation is *completely unnecessary*. "minKeys" would be less than
                        // "keyLen" iff keyBuffer >= (2 * subarrayLen) / keyBuffer... but this situation is
                        // already covered by our scrolling buffer optimization right above!! Consequently,
                        // "minKeys" will *always* be equal to "keyLen" when Grailsort resorts to smart lazy
                        // merges. Removing this loop is by itself a decent optimization, as well!
                        //
                        // Code still here for preservation purposes.
                        /*long halfSubarrKeys = ((long) subarrayLen * keysFound) / 2;
                        int minKeys = CalcMinKeys(keyLen, halfSubarrKeys);*/

                        currBlockLen = (2 * subarrayLen) / keyLen;
                    }
                }

                // WRONG VARIABLE BUG FIXED: 4th argument should be `length - bufferEnd`, was `length - bufferLen` before.
                // Credit to 666666t and Anonymous0726 for debugging.
                CombineBlocks(array, start, start + bufferEnd, length - bufferEnd,
                    subarrayLen, currBlockLen, scrollingBuffer,
                    extBuffer, extBufferLen, comp);
            }

            InsertSort(array, start, bufferEnd, comp);
            LazyMerge(array, start, bufferEnd, length - bufferEnd, comp);
        }
    };
}

template<
    typename RandomAccessIterator,
    typename Compare = std::less<typename std::iterator_traits<RandomAccessIterator>::value_type>
>
void grailsort(RandomAccessIterator first, RandomAccessIterator last, Compare comp = {})
{
    using value_type = typename std::iterator_traits<RandomAccessIterator>::value_type;

    grailsort_detail::GrailSort gsort;
    gsort.CommonSort(first, 0, last - first,
        (value_type*)nullptr, 0,
        grailsort_detail::ThreeWayCompare<Compare>(std::move(comp)));
}

template<
    typename RandomAccessIterator1,
    typename RandomAccessIterator2,
    typename Compare = std::less<typename std::iterator_traits<RandomAccessIterator1>::value_type>
>
void grailsort(RandomAccessIterator1 first, RandomAccessIterator1 last,
    RandomAccessIterator2 buff_first, RandomAccessIterator2 buff_last,
    Compare comp = {})
{
    grailsort_detail::GrailSort gsort;
    gsort.CommonSort(first, 0, last - first,
        buff_first, buff_last - buff_first,
        grailsort_detail::ThreeWayCompare<Compare>(std::move(comp)));
}


void GrailSort(SortArray& A)
{
    grailsort(MyIterator(&A, 0), MyIterator(&A, A.size()));
}

void AuxGrailSort(SortArray& A)
{
    std::vector<ArrayItem> space_buf(512, ArrayItem(0));
    grailsort(MyIterator(&A, 0), MyIterator(&A, A.size()), space_buf.begin(), space_buf.end());
}
// ****************************************************************************
