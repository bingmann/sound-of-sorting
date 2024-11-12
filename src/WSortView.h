/******************************************************************************
 * src/WSortView.h
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

#ifndef WSORTVIEW_H
#define WSORTVIEW_H

#include <wx/wx.h>
#include "SortArray.h"

// ----------------------------------------------------------------------------

/// global sound processing on/off
extern bool g_sound_on;

/// multiplying g_delay with this value yields the duration a sound is sustained
extern double g_sound_sustain;

/// the SDL sound callback
void SoundCallback(void *udata, unsigned char* stream, int len);

/// reset internal state of sound generator
void SoundReset();

/// append access to start output in next sound callback
void SoundAccess(size_t i);

// ----------------------------------------------------------------------------

/// global delay for one step
extern double g_delay;

/// global flag when algorithm is running
extern bool g_algo_running;

/// global string of running algorithm
extern wxString g_algo_name;

class WSortView : public wxPanel, public SortDelay
{
public:
    WSortView(wxWindow *parent, int id, class WMain_wxg* wmain);
    ~WSortView();

    /// reference to WMain
    class WMain* wmain;

    /// the instrumentated array to sort
    SortArray m_array;

    void RepaintNow();

    virtual void OnPaint(wxPaintEvent& pe);
    virtual void OnSize(wxSizeEvent& se);

    // paint the visualization, (dcsize is passed along due to wxMSW issues)
    void paint(wxDC& dc, const wxSize& dcsize);

    /// central access function called by each array access of the algorithms
    virtual void OnAccess();

    /// delay algorithm time by this amount
    void DoDelay(double delay);

protected:

    /// flag for step-wise processing (idle while true)
    bool m_stepwise;

    /// semaphore for signaling steps from main thread
    wxSemaphore m_step_semaphore;

public:

    /// set stepwise processing for Step button
    void SetStepwise(bool v) { m_stepwise = v; }

    /// stepwise processing for Step button
    void DoStepwise();

public:
    DECLARE_EVENT_TABLE()
};

// ----------------------------------------------------------------------------

class SortAlgoThread : public wxThread {
protected:
    class WMain*        m_wmain;
    class WSortView&    m_sortview;

    size_t              m_algo;

public:

    SortAlgoThread(class WMain* wmain, class WSortView& sortview, size_t algo);

    virtual ExitCode Entry();

    void    Exit();
};

// ----------------------------------------------------------------------------

#endif // WSORTVIEW_H
