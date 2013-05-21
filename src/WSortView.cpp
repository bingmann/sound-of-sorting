/******************************************************************************
 * src/WSortView.cpp
 *
 * WSortView contains both the instrumentated array operators and draws the
 * displayed visualization.
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
      wmain(reinterpret_cast<WMain*>(wmain)),
      m_step_condition(m_step_mutex)
{
    SetBackgroundStyle(wxBG_STYLE_CUSTOM);

    m_stepwise = false;
    m_step_mutex.Lock();
}

WSortView::~WSortView()
{
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
}

void WSortView::FillInputlist(wxControlWithItems* list)
{
    list->Append(_("Random Shuffle"));
    list->Append(_("Ascending"));
    list->Append(_("Descending"));
    list->Append(_("Shuffled Cubic"));
    list->Append(_("Shuffled Quintic"));
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
        wxCondError ce = m_step_condition.WaitTimeout(100);
        if (ce == wxCOND_NO_ERROR)
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

void WSortView::OnAccess()
{
    ++g_access_count;

    DoDelay(g_delay);
}

void WSortView::CheckSorted()
{
    unmark_all();

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

template <typename Vector>
static inline bool InList(const Vector& vec, const typename Vector::value_type& v)
{
    for (typename Vector::const_iterator it = vec.begin();
         it != vec.end(); ++it)
    {
        if (*it == v) return true;
    }
    return false;
}

unsigned char WSortView::InWatchList(ssize_t idx) const
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
        *wxGREEN_PEN,
        *wxCYAN_PEN,
        wxPen(wxColour(255,255,0)), // yellow
        wxPen(wxColour(255,0,255)), // magenta
        wxPen(wxColour(255,192,128)), // orange
        wxPen(wxColour(255,128,192)), // pink
        wxPen(wxColour(128,192,255)), // darker cyan
        wxPen(wxColour(192,255,128)), // darker green
        wxPen(wxColour(192,128,255)), // purple
        wxPen(wxColour(128,255,192)), // light green
        wxPen(wxColour(128,128,255)), // blue
        wxPen(wxColour(192,128,192)), // dark purple
        wxPen(wxColour(128,192,192)), // dark cyan
        wxPen(wxColour(192,192,128)), // dark yellow
    };

    static const wxBrush brushes[] = {
        *wxWHITE_BRUSH,
        *wxGREEN_BRUSH,
        *wxCYAN_BRUSH,
        wxBrush(wxColour(255,255,0)), // yellow
        wxBrush(wxColour(255,0,255)), // magenta
        wxBrush(wxColour(255,192,128)), // orange
        wxBrush(wxColour(255,128,192)), // pink
        wxBrush(wxColour(128,192,255)), // darker cyan
        wxBrush(wxColour(192,255,128)), // darker green
        wxBrush(wxColour(192,128,255)), // purple
        wxBrush(wxColour(128,255,192)), // light green
        wxBrush(wxColour(128,128,255)), // blue
        wxBrush(wxColour(192,128,192)), // dark purple
        wxBrush(wxColour(128,192,192)), // dark cyan
        wxBrush(wxColour(192,192,128)), // dark yellow
    };

    wxMutexLocker lock(m_mutex);
    ASSERT(lock.IsOk());

    unsigned char clr;

    for (size_t i = 0; i < size(); ++i)
    {
        if (i == m_access1 || i == m_access2)
        {
            dc.SetPen( *wxRED_PEN );
            dc.SetBrush( *wxRED_BRUSH );
        }
        else if ( (clr = InWatchList(i)) != 0 )
        {
            ASSERT(clr < sizeof(brushes) / sizeof(brushes[0]));
            dc.SetPen( pens[clr] );
            dc.SetBrush( brushes[clr] );
        }
        else if (m_mark[i] != 0)
        {
            ASSERT(m_mark[i] < sizeof(brushes) / sizeof(brushes[0]));
            dc.SetPen( pens[m_mark[i]] );
            dc.SetBrush( brushes[m_mark[i]] );
        }
        else if (InList(m_access_list, i))
        {
            dc.SetPen( *wxRED_PEN );
            dc.SetBrush( *wxRED_BRUSH );
        }
        else
        {
            dc.SetPen( *wxWHITE_PEN );
            dc.SetBrush( *wxWHITE_BRUSH );
        }

        dc.DrawRectangle(i*bstep, height,
                         wxMax(1, // draw at least 1 pixel
                               (wxCoord((i+1)*bstep) - wxCoord(i*bstep)) // integral gap to next bar
                               - (bstep - wbar)    // space between bars
                             ),
                         -(double)height * m_array[i].get_direct() / m_array_max);
    }

    m_access_list.clear();
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

    ae.func(m_array);

    m_array.CheckSorted();

    wxCommandEvent evt(wxEVT_COMMAND_BUTTON_CLICKED, WMain::ID_RUN_FINISHED);
    m_wmain->AddPendingEvent(evt);

    return NULL;
}

void SortAlgoThread::Exit()
{
    wxThread::Exit();
}

// ****************************************************************************
