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
    m_array.SetSortDelay(this);
}

WSortView::~WSortView()
{
}

#if MSW_PERFORMANCECOUNTER

void mswMicroSleep(int microseconds)
{
    static LARGE_INTEGER s_liFreq = { 0, 0 };

    if (s_liFreq.QuadPart == 0)
        QueryPerformanceFrequency(&s_liFreq);

    LARGE_INTEGER liStart, liGoal;
    QueryPerformanceCounter(&liStart);

    liGoal.QuadPart = liStart.QuadPart;
    liGoal.QuadPart += (s_liFreq.QuadPart * microseconds) / 1000000;
    LARGE_INTEGER liAct;
    do
    {
        QueryPerformanceCounter(&liAct);
        if (liStart.QuadPart > liAct.QuadPart) break;
    } while(liAct.QuadPart < liGoal.QuadPart);
}

#endif // MSW_PERFORMANCECOUNTER

static double returnFrac(double n) { return n - static_cast<int>(n); }

void WSortView::DoDelay(double delay)
{
    // must be called by the algorithm thread
    ASSERT(wmain->m_thread);
    ASSERT(wxThread::GetCurrentId() == wmain->m_thread->GetId());

    if (wmain->m_thread_terminate)
    {
        wmain->m_thread->Exit();
        return;
    }

    // idle until main thread signals a condition
    while (m_stepwise)
    {
        wxSemaError se = m_step_semaphore.WaitTimeout(200);
        if (se == wxSEMA_NO_ERROR)
            break;
        // else timeout, recheck m_stepwise and loop
        wxMilliSleep(1);
    }
    double secs = 0.0, microDelay = returnFrac(delay) * 1000.0;
    delay = std::trunc(delay);
    #if __WXGTK__
        // wxMicroSleep(delay * 1000.0);
        while (secs < delay)
        {
            wxMilliSleep(1);
            secs += 1.0;
            if (wmain->m_thread_terminate)
            {
                wmain->m_thread->Exit();
                return;
            }
        }
        secs = 0.0;
        while (secs < microDelay)
        {
            wxMicroSleep(1);
            secs += 1.0;
            if (wmain->m_thread_terminate)
            {
                wmain->m_thread->Exit();
                return;
            }
        }
    #elif MSW_PERFORMANCECOUNTER
        // mswMicroSleep(delay * 1000.0);
        while (secs < delay)
        {
            mswMicroSleep(1000);
            secs += 1.0;
            if (wmain->m_thread_terminate)
            {
                wmain->m_thread->Exit();
                return;
            }
        }
        secs = 0.0;
        while (secs < microDelay)
        {
            mswMicroSleep(1);
            secs += 1.0;
            if (wmain->m_thread_terminate)
            {
                wmain->m_thread->Exit();
                return;
            }
        }
    #else
        // wxMSW does not have a high resolution timer, maybe others do?
        // wxMilliSleep(delay);
        if (delay < 1.0) { delay = 1.0; }
        while (secs < delay)
        {
            wxMilliSleep(1);
            secs += 1.0;
            if (wmain->m_thread_terminate)
            {
                wmain->m_thread->Exit();
                return;
            }
        }
    #endif
}

void WSortView::OnAccess()
{
    DoDelay(g_delay);
}

void WSortView::DoStepwise()
{
    wxSemaError se = m_step_semaphore.Post();
    if (se != wxSEMA_NO_ERROR) {
        wxLogError(_T("Error posting to semaphore: %d"), se);
    }
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

void WSortView::paint(wxDC& dc, const wxSize& dcsize)
{
    dc.SetBackground(*wxBLACK_BRUSH);
    dc.Clear();

    if (m_array.size() == 0) return;

    size_t size = m_array.size();

    size_t fwidth = dcsize.GetWidth();
    size_t fheight = dcsize.GetHeight();

    size_t width = fwidth - 20;
    size_t height = fheight - 20;

    dc.SetDeviceOrigin(10, 10);

    // *** draw array element bars

    // draw | | | |  bars: each bar is width w, separation is w/2
    // thus n bars need n * w + (n-1) * w/2 width

    // 1st variant: space is 0.5 of bar size
    //double wbar = width / (size + (size-1) / 2.0);
    //double bstep = 1.5 * wbar;

    // 2nd variant: one pixel between bars
    size_t step = 1;
    double wbar = (width - (size - 1)) / (double)size;
    double bstep = wbar + 1.0;
    if (size > width) {
        step = size / width;
        wbar = wxMax(0.03, width / (double)size);  // Scale down further if too many elements
        bstep = wbar;
    }

    // special case for bstep = 2 pixel -> draw 2 pixel bars instead of 1px
    // bar/1px gaps.
    if (fabs(wbar - 1.0) < 0.1 && fabs(bstep - 2.0) < 0.1) wbar = 2, bstep = 2;

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
        wxPen(wxColour(255,64,64)),   // 17 bright red
        wxPen(wxColour(64,255,64)),   // 18 bright green
        wxPen(wxColour(64,64,255)),   // 19 bright blue
        wxPen(wxColour(255,128,64)),  // 20 coral
        wxPen(wxColour(128,255,64)),  // 21 lime green
        wxPen(wxColour(64,255,255)),  // 22 aqua
        wxPen(wxColour(255,215,0)),   // 23 gold
        wxPen(wxColour(128,0,128)),   // 24 purple
        wxPen(wxColour(0,255,127)),   // 25 spring green
        wxPen(wxColour(255,69,0)),    // 26 orange red
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
        wxBrush(wxColour(255,64,64)),   // 17 bright red
        wxBrush(wxColour(64,255,64)),   // 18 bright green
        wxBrush(wxColour(64,64,255)),   // 19 bright blue
        wxBrush(wxColour(255,128,64)),  // 20 coral
        wxBrush(wxColour(128,255,64)),  // 21 lime green
        wxBrush(wxColour(64,255,255)),  // 22 aqua
        wxBrush(wxColour(255,215,0)),   // 23 bright yellow
        wxBrush(wxColour(128,0,128)),   // 24 purple
        wxBrush(wxColour(0,255,127)),   // 25 spring green
        wxBrush(wxColour(255,69,0)),    // 26 orange red
    };

    wxMutexLocker lock(m_array.m_mutex);
    ASSERT(lock.IsOk());
    const int numBrushes = sizeof(brushes) / sizeof(brushes[0]);
    if (step > 1)
    {
        size_t i_step = 0;
        const size_t last = size - 1;
        for (size_t i = 0; i < size; ++i)
        {
            int clr = m_array.GetIndexColor(i);
            if (i == i_step || i == last)
            {
                clr = clr % numBrushes;
                dc.SetPen(pens[clr]);
                dc.SetBrush(brushes[clr]);
                dc.DrawRectangle(i * bstep, height,
                    wxMax(1, // draw at least 1 pixel
                        (wxCoord((i + 1) * bstep) - wxCoord(i * bstep)) // integral gap to next bar
                        - (bstep - wbar)    // space between bars
                    ),
                    -(double)height * m_array.direct(i).get_direct() / m_array.array_max());
                i_step += step;
            }
        }
    }
    else
    {
        for (size_t i = 0; i < size; ++i)
        {
            int clr = m_array.GetIndexColor(i);
            ASSERT(clr < numBrushes);
            dc.SetPen(pens[clr]);
            dc.SetBrush(brushes[clr]);
            dc.DrawRectangle(i * bstep, height,
                wxMax(1, // draw at least 1 pixel
                    (wxCoord((i + 1) * bstep) - wxCoord(i * bstep)) // integral gap to next bar
                    - (bstep - wbar)    // space between bars
                ),
                -(double)height * m_array.direct(i).get_direct() / m_array.array_max());
        }
    }
};

BEGIN_EVENT_TABLE(WSortView, wxWindow)

    EVT_PAINT(WSortView::OnPaint)
    EVT_SIZE(WSortView::OnSize)

END_EVENT_TABLE();

// ****************************************************************************
// *** Threading

SortAlgoThread::SortAlgoThread(WMain* wmain, WSortView& sortview, size_t algo)
    : wxThread(wxTHREAD_JOINABLE),
      m_wmain(wmain),
      m_sortview(sortview),
      m_algo(algo)
{
}

void* SortAlgoThread::Entry()
{
    ASSERT(m_algo < g_algolist_size);
    const AlgoEntry& ae = g_algolist[m_algo];

    m_sortview.m_array.OnAlgoLaunch(ae);

    ae.func(m_sortview.m_array);

    m_sortview.m_array.CheckSorted();

    wxCommandEvent evt(wxEVT_COMMAND_BUTTON_CLICKED, WMain::ID_RUN_FINISHED);
    m_wmain->GetEventHandler()->AddPendingEvent(evt);

    return nullptr;
}

void SortAlgoThread::Exit()
{
    wxThread::Exit();
}

// ****************************************************************************
