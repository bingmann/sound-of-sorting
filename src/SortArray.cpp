/*******************************************************************************
 * src/SortArray.cpp
 *
 * SortArray represents a simple array, which is displayed by WSortView to the
 * user is real-time.
 *
 *******************************************************************************
 * Copyright (C) 2014 Timo Bingmann <tb@panthema.net>
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
 ******************************************************************************/

#include "SortArray.h"
#include "SortAlgo.h"

#include <algorithm>

extern void SoundAccess(size_t i);

// *****************************************************************************
// *** Comparisons of ArrayItems

size_t g_compare_count = 0;

size_t g_access_count = 0;

void ArrayItem::OnAccess(const ArrayItem& a)
{
    SoundAccess(a.get_direct());
}

void ArrayItem::OnComparison(const ArrayItem& a, const ArrayItem& b)
{
    ++g_compare_count;

    SoundAccess(a.get_direct());
    SoundAccess(b.get_direct());
}

// *****************************************************************************
// *** SortArray

SortArray::SortArray()
    : m_calc_inversions(false),
      m_delay(NULL)
{
}

void SortArray::OnAlgoLaunch(const AlgoEntry& ae)
{
    if (size() <= ae.inversion_count_limit)
    {
        m_calc_inversions = true;
        RecalcInversions();
    }
    else
    {
        m_calc_inversions = false;
        m_inversions = -1;
    }
}

void SortArray::ResetArray(size_t size)
{
    m_array.resize(size, ArrayItem(0));
    m_mark.resize(size);
}

void SortArray::FinishFill()
{
    ASSERT(m_array.size() > 0);

    // calculate max value in array
    m_array_max = m_array[0].get_direct();
    for (size_t i = 1; i < size(); ++i)
    {
        if (m_array_max < m_array[i].get_direct())
            m_array_max = m_array[i].get_direct();
    }

    // reset access markers
    unmark_all();
    unwatch_all();

    // reset counters and info
    m_is_sorted = false;
    g_access_count = 0;
    g_compare_count = 0;
    m_calc_inversions = true;

    RecalcInversions();
}

void SortArray::FillInputlist(wxArrayString& list)
{
    list.Add(_("Random Shuffle"));
    list.Add(_("Ascending"));
    list.Add(_("Descending"));
    list.Add(_("Shuffled Cubic"));
    list.Add(_("Shuffled Quintic"));
    list.Add(_("Shuffled n-2 Equal"));
}

void SortArray::FillData(unsigned int schema, size_t arraysize)
{
    if (arraysize == 0) arraysize = 1;

    ResetArray(arraysize);

    if (schema == 0) // Shuffle of [1,n]
    {
        for (size_t i = 0; i < m_array.size(); ++i)
            m_array[i] = ArrayItem(i+1);

        std::random_shuffle(m_array.begin(), m_array.end());
    }
    else if (schema == 1) // Ascending [1,n]
    {
        for (size_t i = 0; i < m_array.size(); ++i)
            m_array[i] = ArrayItem(i+1);
    }
    else if (schema == 2) // Descending [1,n]
    {
        for (size_t i = 0; i < m_array.size(); ++i)
            m_array[i] = ArrayItem(m_array.size() - i);
    }
    else if (schema == 3) // Cubic skew of [1,n]
    {
        for (size_t i = 0; i < m_array.size(); ++i)
        {
            // normalize to [-1,+1]
            double x = (2.0 * (double)i / m_array.size()) - 1.0;
            // calculate x^3
            double v = x * x * x;
            // normalize to array size
            double w = (v + 1.0) / 2.0 * arraysize + 1;
            // decrease resolution for more equal values
            w /= 3.0;
            m_array[i] = ArrayItem(w + 1);
        }

        std::random_shuffle(m_array.begin(), m_array.end());
    }
    else if (schema == 4) // Quintic skew of [1,n]
    {
        for (size_t i = 0; i < m_array.size(); ++i)
        {
            // normalize to [-1,+1]
            double x = (2.0 * (double)i / m_array.size()) - 1.0;
            // calculate x^5
            double v = x * x * x * x * x;
            // normalize to array size
            double w = (v + 1.0) / 2.0 * arraysize + 1;
            // decrease resolution for more equal values
            w /= 3.0;
            m_array[i] = ArrayItem(w + 1);
        }

        std::random_shuffle(m_array.begin(), m_array.end());
    }
    else if (schema == 5) // shuffled n-2 equal values in [1,n]
    {
        m_array[0] = ArrayItem(1);
        for (size_t i = 1; i < m_array.size()-1; ++i)
        {
            m_array[i] = ArrayItem( arraysize / 2 + 1 );
        }
        m_array[m_array.size()-1] = ArrayItem(arraysize);

        std::random_shuffle(m_array.begin(), m_array.end());
    }
    else // fallback
    {
        return FillData(0, arraysize);
    }

    FinishFill();
}

void SortArray::OnAccess()
{
    ++g_access_count;

    if (m_delay)
        m_delay->OnAccess();
}

bool SortArray::CheckSorted()
{
    unmark_all();
    // needed because iterator instrumentated algorithms may have changed the array
    RecalcInversions();

    ArrayItem prev = get_nocount(0);
    mark(0);

    bool is_sorted = true;

    for (size_t i = 1; i < size(); ++i)
    {
        ArrayItem key = get_nocount(i);
        g_compare_count--; // dont count the following comparison
        if (!(prev <= key)) {
            wxLogError(_T("Result of sorting algorithm is incorrect!"));
            is_sorted = false;
            break;
        }
        mark(i);
        prev = key;
    }

    unmark_all();

    return (m_is_sorted = is_sorted);
}

void SortArray::SetCalcInversions(bool on)
{
    // toggle boolean
    m_calc_inversions = on;

    if (!m_calc_inversions)
        m_inversions = -1;
}

void SortArray::ToggleCalcInversions()
{
    // toggle boolean
    SetCalcInversions(!m_calc_inversions);
}

void SortArray::RecalcInversions()
{
    if (!m_calc_inversions) {
        m_inversions = -1;
        return;
    }

    unsigned int inversions = 0;

    for (size_t i = 0; i < size(); ++i)
    {
        const ArrayItem& a = direct(i);

        for (size_t j = i+1; j < size(); ++j)
        {
            const ArrayItem& b = direct(j);

            if ( a.greater_direct(b) )
            {
                inversions++;
            }
        }
    }

    m_inversions = inversions;
}

void SortArray::UpdateInversions(size_t i, size_t j)
{
    if (!m_calc_inversions) {
        m_inversions = -1;
        return;
    }
    if (m_inversions < 0) return RecalcInversions();

    if (i == j) return;

    unsigned int lo = i, hi = j;
    if (lo > hi) std::swap(lo, hi);

    const ArrayItem& ilo = m_array[lo];
    const ArrayItem& ihi = m_array[hi];
    int invdelta = 0;

    for (size_t k = lo + 1; k < hi; ++k)
    {
        if (m_array[k].less_direct(ilo))
            invdelta--;
        if (m_array[k].greater_direct(ilo))
            invdelta++;

        if (m_array[k].less_direct(ihi))
            invdelta++;
        if (m_array[k].greater_direct(ihi))
            invdelta--;
    }

    if (ilo.less_direct(ihi))
        invdelta++;
    if (ilo.greater_direct(ihi))
        invdelta--;

    m_inversions += invdelta;
}

size_t SortArray::GetRuns() const
{
    unsigned int runs = 1;

    for (size_t i = 1; i < size(); ++i)
    {
        const ArrayItem& a = direct(i-1);
        const ArrayItem& b = direct(i);

        if ( a.greater_direct(b) )
        {
            runs++;
        }
    }

    return runs;
}

short SortArray::InAccessList(ssize_t idx)
{
    if (idx < 0) return -1;

    signed color = -1;
    signed priority = -1;

    for (std::vector<Access>::iterator it = m_access_list.begin();
         it != m_access_list.end(); )
    {
        if (it->index != (size_t)idx) {
            ++it;
            continue;
        }

        if (it->priority >= priority)
        {
            priority = it->priority;
            color = it->color;
        }

        if (it->sustain == 0) {
            if (it->index == m_access1.index ||
                it->index == m_access2.index)
            {
                ++it;
            }
            else
            {
                it = m_access_list.erase(it);
            }
        }
        else {
            it->sustain--;
            ++it;
        }
    }

    return color;
}

unsigned short SortArray::InWatchList(ssize_t idx) const
{
    for (size_t i = 0; i < m_watch.size(); ++i)
    {
        if (m_watch[i].first == NULL) continue;

        // compare watched value
        if (*m_watch[i].first != idx) continue;

        return m_watch[i].second;
    }
    return 0;
}

int SortArray::GetIndexColor(size_t idx)
{
    int clr;

    // select color
    if (idx == m_access1.index)
    {
        clr = m_access1.color;
    }
    else if (idx == m_access2.index)
    {
        clr = m_access2.color;
    }
    else if ( (clr = InWatchList(idx)) != 0 )
    {
        // clr already set
    }
    else if (m_mark[idx] != 0)
    {
        clr = m_mark[idx];
    }
    else if ( (clr = InAccessList(idx)) >= 0 )
    {
    }
    else
    {
        clr = 0;
    }

    return clr;
}

// *****************************************************************************
