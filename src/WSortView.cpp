/******************************************************************************
 * src/WSortView.cpp
 *
 * WSortView contains both the instrumentated array operators and draws the
 * displayed visualization.
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

#include "WSortView.h"
#include "SortAlgo.h"
#include "WMain.h"

#include <wx/dcbuffer.h>

// ****************************************************************************
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

// ****************************************************************************
// *** WSortView

double g_delay = 0;

bool g_algo_running = false;

wxString g_algo_name;

WSortView::WSortView(wxWindow *parent, int id, class WMain_wxg* wmain)
    : wxPanel(parent, id),
      wmain(reinterpret_cast<WMain*>(wmain))
{
    SetBackgroundStyle(wxBG_STYLE_CUSTOM);

    m_stepwise = false;
}

WSortView::~WSortView()
{
}

void WSortView::OnAlgoLaunch(const AlgoEntry& ae)
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

void WSortView::ResetArray(size_t size)
{
    m_array.resize(size, ArrayItem(0));
    m_mark.resize(size);
}

void WSortView::FinishFill()
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

    RecalcInversions();
}

void WSortView::FillInputlist(wxControlWithItems* list)
{
    list->Append(_("Random Shuffle"));
    list->Append(_("Ascending"));
    list->Append(_("Descending"));
    list->Append(_("Shuffled Cubic"));
    list->Append(_("Shuffled Quintic"));
    list->Append(_("Shuffled n-2 Equal"));
}

void WSortView::FillData(unsigned int schema, size_t arraysize)
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

void WSortView::DoDelay(double delay)
{
    // must be called by the algorithm thread
    ASSERT(wmain->m_thread);
    ASSERT(wxThread::GetCurrentId() == wmain->m_thread->GetId());

    if (wmain->m_thread_terminate)
        wmain->m_thread->Exit();

    // idle until main thread signals a condition
    while (m_stepwise)
    {
        wxSemaError se = m_step_semaphore.WaitTimeout(200);
        if (se == wxSEMA_NO_ERROR)
            break;
        // else timeout, recheck m_stepwise and loop
        wmain->m_thread->TestDestroy();
        wmain->m_thread->Yield();
    }

    wmain->m_thread->TestDestroy();

#if __WXGTK__
    wxMicroSleep(delay * 1000.0);
#else
    // wxMSW does not have a high resolution timer, maybe others do?
    wxMilliSleep(delay);
#endif
}

void WSortView::DoStepwise()
{
    wxSemaError se = m_step_semaphore.Post();
    if (se != wxSEMA_NO_ERROR)
        wxLogError(_T("Error posting to semaphore: %d"), se);
}

void WSortView::OnAccess()
{
    ++g_access_count;

    DoDelay(g_delay);
}

void WSortView::CheckSorted()
{
    unmark_all();
    // needed because iterator instrumentated algorithms may have changed the array
    RecalcInversions();

    ArrayItem prev = get_nocount(0);
    mark(0);

    for (size_t i = 1; i < size(); ++i)
    {
        ArrayItem key = get_nocount(i);
        g_compare_count--; // dont count the following comparison
        if (!(prev <= key)) {
            wxLogError(_T("Result of sorting algorithm is incorrect!"));
            break;
        }
        mark(i);
        prev = key;
    }

    unmark_all();

    m_is_sorted = true;
}

void WSortView::ToggleCalcInversions()
{
    // toggle boolean
    m_calc_inversions = !m_calc_inversions;

    if (!m_calc_inversions)
        m_inversions = -1;
}

void WSortView::RecalcInversions()
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

void WSortView::UpdateInversions(size_t i, size_t j)
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

void WSortView::RepaintNow()
{
    if (!IsShownOnScreen()) return;

    wxClientDC dc(this);
    wxBufferedDC bdc(&dc, GetSize(), wxBUFFER_CLIENT_AREA);
    paint(bdc, GetSize());
}

void WSortView::OnPaint(wxPaintEvent&)
{
    wxAutoBufferedPaintDC dc(this);
    // on wxMSW, wxAutoBufferedPaintDC holds a bitmap. The bitmaps size =
    // wxDC's size may not match the window size, thus we pass the window size
    // explicitly.
    paint(dc, GetSize());
}

void WSortView::OnSize(wxSizeEvent&)
{
    // full redraw on resize
    Refresh(false);
}

short WSortView::InAccessList(ssize_t idx)
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

unsigned short WSortView::InWatchList(ssize_t idx) const
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

void WSortView::paint(wxDC& dc, const wxSize& dcsize)
{
    dc.SetBackground(*wxBLACK_BRUSH);
    dc.Clear();

    if (size() == 0) return;

    size_t fwidth = dcsize.GetWidth();
    size_t fheight = dcsize.GetHeight();

    size_t width = fwidth - 20;
    size_t height = fheight - 20;

    dc.SetDeviceOrigin(10,10);

    // *** draw array element bars

    // draw | | | |  bars: each bar is width w, separation is w/2
    // thus n bars need n * w + (n-1) * w/2 width

    // 1st variant: space is 0.5 of bar size
    //double wbar = width / (size() + (size()-1) / 2.0);
    //double bstep = 1.5 * wbar;

    // 2nd variant: one pixel between bars
    double wbar = (width - (size()-1)) / (double)size();
    if (width <= (size()-1)) wbar = 0.0;

    double bstep = wbar + 1.0;

    // special case for bstep = 2 pixel -> draw 2 pixel bars instead of 1px
    // bar/1px gaps.
    if ( fabs(wbar - 1.0) < 0.1 && fabs(bstep - 2.0) < 0.1 ) wbar = 2, bstep = 2;

    static const wxPen pens[] = {
        *wxWHITE_PEN,
        *wxRED_PEN,
        *wxGREEN_PEN,
        *wxCYAN_PEN,
        wxPen(wxColour(255,255,0)),   //  4 yellow
        wxPen(wxColour(255,0,255)),   //  5 magenta
        wxPen(wxColour(255,192,128)), //  6 orange
        wxPen(wxColour(255,128,192)), //  7 pink
        wxPen(wxColour(128,192,255)), //  8 darker cyan
        wxPen(wxColour(192,255,128)), //  9 darker green
        wxPen(wxColour(192,128,255)), // 10 purple
        wxPen(wxColour(128,255,192)), // 11 light green
        wxPen(wxColour(128,128,255)), // 12 blue
        wxPen(wxColour(192,128,192)), // 13 dark purple
        wxPen(wxColour(128,192,192)), // 14 dark cyan
        wxPen(wxColour(192,192,128)), // 15 dark yellow
        wxPen(wxColour(0,128,255)),   // 16 blue/cyan mix
    };

    static const wxBrush brushes[] = {
        *wxWHITE_BRUSH,
        *wxRED_BRUSH,
        *wxGREEN_BRUSH,
        *wxCYAN_BRUSH,
        wxBrush(wxColour(255,255,0)),   //  4 yellow
        wxBrush(wxColour(255,0,255)),   //  5 magenta
        wxBrush(wxColour(255,192,128)), //  6 orange
        wxBrush(wxColour(255,128,192)), //  7 pink
        wxBrush(wxColour(128,192,255)), //  8 darker cyan
        wxBrush(wxColour(192,255,128)), //  9 darker green
        wxBrush(wxColour(192,128,255)), // 10 purple
        wxBrush(wxColour(128,255,192)), // 11 light green
        wxBrush(wxColour(128,128,255)), // 12 blue
        wxBrush(wxColour(192,128,192)), // 13 dark purple
        wxBrush(wxColour(128,192,192)), // 14 dark cyan
        wxBrush(wxColour(192,192,128)), // 15 dark yellow
        wxBrush(wxColour(0,128,255)),   // 16 blue/cyan mix
    };

    wxMutexLocker lock(m_mutex);
    ASSERT(lock.IsOk());

    for (size_t i = 0; i < size(); ++i)
    {
        int clr;

        // select color
        if (i == m_access1.index)
        {
            clr = m_access1.color;
        }
        else if (i == m_access2.index)
        {
            clr = m_access2.color;
        }
        else if ( (clr = InWatchList(i)) != 0 )
        {
            // clr already set
        }
        else if (m_mark[i] != 0)
        {
            clr = m_mark[i];
        }
        else if ( (clr = InAccessList(i)) >= 0 )
        {
        }
        else
        {
            clr = 0;
        }

        ASSERT(clr < (int)(sizeof(brushes) / sizeof(brushes[0])));
        dc.SetPen( pens[clr] );
        dc.SetBrush( brushes[clr] );

        dc.DrawRectangle(i*bstep, height,
                         wxMax(1, // draw at least 1 pixel
                               (wxCoord((i+1)*bstep) - wxCoord(i*bstep)) // integral gap to next bar
                               - (bstep - wbar)    // space between bars
                             ),
                         -(double)height * m_array[i].get_direct() / m_array_max);
    }
}

BEGIN_EVENT_TABLE(WSortView, wxWindow)

    EVT_PAINT		(WSortView::OnPaint)
    EVT_SIZE            (WSortView::OnSize)

END_EVENT_TABLE()

// ****************************************************************************
// *** Threading

SortAlgoThread::SortAlgoThread(WMain* wmain, class WSortView& array, size_t algo)
    : wxThread(wxTHREAD_JOINABLE),
      m_wmain(wmain),
      m_array(array),
      m_algo(algo)
{
}

void* SortAlgoThread::Entry()
{
    ASSERT(m_algo < g_algolist_size);
    const AlgoEntry& ae = g_algolist[m_algo];

    m_array.OnAlgoLaunch(ae);

    ae.func(m_array);

    m_array.CheckSorted();

    wxCommandEvent evt(wxEVT_COMMAND_BUTTON_CLICKED, WMain::ID_RUN_FINISHED);
    m_wmain->GetEventHandler()->AddPendingEvent(evt);

    return NULL;
}

void SortAlgoThread::Exit()
{
    wxThread::Exit();
}

// ****************************************************************************
