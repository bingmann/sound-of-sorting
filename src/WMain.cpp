/******************************************************************************
 * src/WMain.cpp
 *
 * Implementation of the main window's functions.
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

#include "WMain.h"
#include "SortAlgo.h"
#include <SDL.h>

#include "wxg/WAbout_wxg.h"

WMain::WMain(wxWindow* parent)
    : WMain_wxg(parent, wxID_ANY, wxEmptyString, wxDefaultPosition, wxDefaultSize, wxDEFAULT_FRAME_STYLE)
{
    m_thread = nullptr;
    g_sound_on = false;

    recordButton->Hide();
    panelQuickSortPivot->Hide();
    infoTextctrl->Hide();

    // program icon
    {    
        #include "sos.xpm"
	SetIcon( wxIcon(sos) );
    }

    // program version
    SetTitle(_("The Sound of Sorting " PACKAGE_VERSION " - http://panthema.net/2013/sound-of-sorting"));

    // resize right split window
    splitter_0->SetSashPosition(GetSize().x - 280);
    splitter_0->SetMinimumPaneSize(1);

    // insert list of algorithms into wxListBox
    for (const AlgoEntry* ae = g_algolist; ae != g_algolist_end; ++ae)
        algoList->Append(ae->name);

    algoList->SetSelection(0);

    // insert list of data templates into wxChoice
    wxArrayString inputList;
    sortview->m_array.FillInputlist(inputList);
    inputTypeChoice->Append(inputList);

    sortview->m_array.FillData(0, 100);
    inputTypeChoice->SetSelection(0);
    arraySizeSlider->SetValue(100);
    SetArraySize(100);

    // insert quicksort pivot rules into wxChoice
    pivotRuleChoice->Append( QuickSortPivotText() );
    pivotRuleChoice->SetSelection(0);

    // set default speed
    speedSlider->SetValue(1000);
    SetDelay(1000);

    // set default sound sustain
    soundSustainSlider->SetValue(700);
    SetSoundSustain(700);

    // create audio output
    SDL_AudioSpec sdlaudiospec;

    // Set the audio format
    sdlaudiospec.freq = 44100;
    sdlaudiospec.format = AUDIO_S16SYS;
    sdlaudiospec.channels = 1;    	/* 1 = mono, 2 = stereo */
    sdlaudiospec.samples = 4096;  	/* Good low-latency value for callback */
    sdlaudiospec.callback = SoundCallback;
    sdlaudiospec.userdata = sortview;

    // Open the audio device, forcing the desired format
    if ( SDL_OpenAudio(&sdlaudiospec, nullptr) < 0 ) {
        wxLogError(_("Couldn't open audio: ") + wxString(SDL_GetError(), wxConvISO8859_1));
        soundButton->Disable();
    }
    soundButton->SetValue(false);

    Show();

    // start refreshing timer
    m_reftimer = new RefreshTimer(this);
    m_reftimer->Start(1000 / g_framerate);
}

WMain::~WMain()
{
    AbortAlgorithm();

    SDL_CloseAudio();
}

BEGIN_EVENT_TABLE(WMain, WMain_wxg)

    EVT_TOGGLEBUTTON(ID_RUN_BUTTON, WMain::OnRunButton)
    EVT_BUTTON(ID_RESET_BUTTON, WMain::OnResetButton)
    EVT_BUTTON(ID_STEP_BUTTON, WMain::OnStepButton)
    EVT_TOGGLEBUTTON(ID_SOUND_BUTTON, WMain::OnSoundButton)
    EVT_BUTTON(ID_RANDOM_BUTTON, WMain::OnRandomButton)
    EVT_BUTTON(wxID_ABOUT, WMain::OnAboutButton)
    EVT_SPLITTER_DCLICK(wxID_ANY, WMain::OnDClick)

    EVT_COMMAND_SCROLL(ID_SPEED_SLIDER, WMain::OnSpeedSliderChange)
    EVT_COMMAND_SCROLL(ID_SOUND_SUSTAIN_SLIDER, WMain::OnSoundSustainSliderChange)
    EVT_COMMAND_SCROLL(ID_ARRAY_SIZE_SLIDER, WMain::OnArraySizeSliderChange)
    EVT_LISTBOX(ID_ALGO_LIST, WMain::OnAlgoList)
    EVT_LISTBOX_DCLICK(ID_ALGO_LIST, WMain::OnAlgoListDClick)

    EVT_COMMAND(ID_RUN_FINISHED, wxEVT_COMMAND_BUTTON_CLICKED, WMain::OnRunFinished)

    EVT_COMMAND(ID_INVERSION_LABEL, wxEVT_COMMAND_BUTTON_CLICKED, WMain::OnInversionLabelClick)

END_EVENT_TABLE();

bool WMain::RunAlgorithm()
{
    ASSERT(!m_thread);

    if (algoList->GetSelection() == wxNOT_FOUND)
    {
        wxLogError(_("Please select a sorting algorithm"));
        runButton->SetValue(false);
        return false;
    }
    else
    {
        if (sortview->m_array.IsSorted())
            sortview->m_array.FillData( inputTypeChoice->GetSelection(), m_array_size );

        sortview->SetStepwise(false);

        g_algo_name = algoList->GetStringSelection();
        g_quicksort_pivot = (QuickSortPivotType)pivotRuleChoice->GetSelection();

        m_thread = new SortAlgoThread(this, *sortview, algoList->GetSelection());

        m_thread_terminate = false;
        m_thread->Create();

        g_algo_running = true;
        m_thread->Run();

        runButton->SetValue(true);
        return true;
    }
}

void WMain::OnDClick(wxSplitterEvent& event) { event.Veto(); }

void WMain::AbortAlgorithm()
{
    if (!m_thread) return;
    
    m_thread_terminate = true;
    if (soundButton->GetValue() == true)
    {
        SoundReset();
        g_sound_on = false;
        SDL_PauseAudio(1);
    }

    if (m_thread->IsPaused()) { m_thread->Resume(); }
    sortview->SetStepwise(false);

    while (m_thread->IsAlive()) { wxMilliSleep(1); }
    g_algo_running = false;
    if (soundButton->GetValue() == true)
    {
        g_sound_on = true;
        SDL_PauseAudio(0);
    }

    delete m_thread;
    m_thread = nullptr;
}

void WMain::OnRunButton(wxCommandEvent &event)
{
    // join finished thread
    if (m_thread && !m_thread->IsAlive())
    {
        g_algo_running = false;
        delete m_thread;
        m_thread = nullptr;
    }

    if (event.IsChecked())
    {
        if (!m_thread)
        {
            RunAlgorithm();
        }
        else
        {
            g_algo_running = true;
            sortview->SetStepwise(false);
            if (m_thread->IsPaused()) m_thread->Resume();
        }
    }
    else
    {
        if (m_thread)
        {
            m_thread->Pause();
            g_algo_running = false;
        }
    }
}

void WMain::OnRunFinished(wxCommandEvent&)
{
    // join finished thread
    if (m_thread)
    {
        m_thread->Wait();
        g_algo_running = false;

        delete m_thread;
        m_thread = nullptr;
    }

    runButton->SetValue(false);
}

void WMain::OnResetButton(wxCommandEvent&)
{
    // terminate running algorithm.
    AbortAlgorithm();

    runButton->SetValue(false);

    sortview->m_array.FillData( inputTypeChoice->GetSelection(), m_array_size );
}

void WMain::OnStepButton(wxCommandEvent&)
{
    // lock thread to one step
    if (!m_thread)
    {
        if (!RunAlgorithm()) return;
        sortview->SetStepwise(true);
    }
    else
    {
        if (m_thread->IsPaused()) m_thread->Resume();
        sortview->SetStepwise(true);    // in case not already set
        runButton->SetValue(false);
        g_algo_running = true;
    }

    // signal to perform one step
    sortview->DoStepwise();
}

void WMain::OnSoundButton(wxCommandEvent&)
{
    if (soundButton->GetValue())
    {
        SoundReset();
        SDL_PauseAudio(0);
        g_sound_on = true;
    }
    else
    {
        g_sound_on = false;
        SDL_PauseAudio(1);
    }
}

void WMain::OnRandomButton(wxCommandEvent&)
{
    AbortAlgorithm();

    algoList->SetSelection( rand() % algoList->GetCount() );
    sortview->m_array.FillData( inputTypeChoice->GetSelection(), m_array_size );

    RunAlgorithm();

    algoList->SetSelection(wxNOT_FOUND);
}

class WAbout : public WAbout_wxg
{
public:
    WAbout(wxWindow* parent)
        : WAbout_wxg(parent, wxID_ANY, wxEmptyString)
    {
        labelTitle->SetLabel(_("The Sound of Sorting " PACKAGE_VERSION));
        labelBuildDate->SetLabel(_("Build Date: " __DATE__));

        GetSizer()->Fit(this);
        Layout();
        Centre();
    }
};

void WMain::OnAboutButton(wxCommandEvent&)
{
    WAbout dlg(this);
    dlg.ShowModal();
}

void WMain::OnSpeedSliderChange(wxScrollEvent &event)
{
    SetDelay(event.GetPosition());
}

void WMain::SetDelay(size_t pos)
{
    // scale slider with this base to 10 seconds = 2 * 1000 ms
    const double base = 4;

    // different slider scale for Linux/GTK: (faster)
#if __WXGTK__ || MSW_PERFORMANCECOUNTER
    // 0.001 ms formula, tested on Windows
    // Warning! The slider will never go past 0.005 ms with this formula
    // g_delay = pow(base, pos / 15000.0 * log(2 * 1000.0 * 10.0) / log(base)) / 800.0;
    g_delay = pow(base, pos / 2000.0 * log(2 * 1000.0 * 10.0) / log(base)) / 10.0;
#else
    // other systems probably have sucking real-time performance anyway
    g_delay = pow(base, pos / 2000.0 * log(2 * 1000.0) / log(base));
    
#endif
    if (pos == 0) g_delay = 0.1;

    if (g_delay > 10)
        labelDelayValue->SetLabel(wxString::Format(_("%.0f ms"), g_delay));
    else if (g_delay > 1)
        labelDelayValue->SetLabel(wxString::Format(_("%.2f ms"), g_delay));
    else
        labelDelayValue->SetLabel(wxString::Format(_("%.3f ms"), g_delay));
}

void WMain::OnSoundSustainSliderChange(wxScrollEvent &event)
{
    SetSoundSustain(event.GetPosition());
}

void WMain::OnInversionLabelClick(wxCommandEvent &)
{
    sortview->m_array.ToggleCalcInversions();
}

void WMain::SetSoundSustain(size_t pos)
{
    // scale slider with this base to 50
    const double base = 4;

    g_sound_sustain = pow(base, pos / 2000.0 * log(50) / log(base));

    labelSoundSustainValue->SetLabel(wxString::Format(_("%.1f"), g_sound_sustain));

    labelSoundSustainValue->GetContainingSizer()->Layout();
}

void WMain::OnArraySizeSliderChange(wxScrollEvent &event)
{
    SetArraySize(event.GetPosition());
}

void WMain::SetArraySize(size_t pos)
{
    // scale slider with this base to 2048
    //const double base = 2;
    //m_array_size = pow(base, pos / 10000.0 * log(2048) / log(base));
    m_array_size = pos;

    labelArraySizeValue->SetLabel(wxString::Format(_("%4lu"), (long unsigned)m_array_size));

    labelArraySizeValue->GetContainingSizer()->Layout();
}

void WMain::OnAlgoList(wxCommandEvent&)
{
    int sel = algoList->GetSelection();
    wxString text;

    bool isQuickSort = (algoList->GetStringSelection().Contains(_("Quick Sort")));
    panelQuickSortPivot->Show(isQuickSort);

    if (sel >= 0 && sel < (int)g_algolist_size && !g_algolist[sel].text.IsEmpty())
    {
        text = g_algolist[sel].text;

        infoTextctrl->SetValue(text);

        infoTextctrl->Show();
        infoTextctrl->GetContainingSizer()->Layout();
    }
    else
    {
        infoTextctrl->Hide();
        infoTextctrl->GetContainingSizer()->Layout();
    }
}

void WMain::OnAlgoListDClick(wxCommandEvent&)
{
    // terminate running algorithm.
    if (m_thread)
    {
        AbortAlgorithm();

        sortview->m_array.FillData( inputTypeChoice->GetSelection(), m_array_size );
    }

    // start new one
    RunAlgorithm();
}

// ----------------------------------------------------------------------------

WMain::RefreshTimer::RefreshTimer(WMain* wmain)
    : wm(*wmain)
{
}

void WMain::RefreshTimer::Notify()
{
    // update some labels with current values
    long int accesses = g_access_count;
    wm.labelAccessCount->SetLabel(wxString::Format(_("%ld"), accesses));

    long int compares = g_compare_count;
    wm.labelComparisonsValue->SetLabel(wxString::Format(_("%ld"), compares));

    long int swaps = m_swaps;
    wm.labelSwapsCount->SetLabel(wxString::Format(_("%ld"), swaps));

    long int inversions = wm.sortview->m_array.GetInversions();
    if (inversions >= 0)
        wm.labelInversionCount->SetLabel(wxString::Format(_("%ld"), inversions));
    else
        wm.labelInversionCount->SetLabel(_("skipped"));

    long int runs = wm.sortview->m_array.GetRuns();
    wm.labelRunsCount->SetLabel(wxString::Format(_("%ld"), runs));

    // repaint sortview
    wm.sortview->Refresh(false);
}

// ----------------------------------------------------------------------------

class MyApp : public wxApp
{
public:

    virtual bool OnInit()
    {
        if (SDL_Init(SDL_INIT_AUDIO) < 0)
        {
            wxLogError(wxT("Couldn't initialize SDL: ") + wxString(SDL_GetError(), wxConvISO8859_1));
            wxLog::FlushActive();
            return false;
        }

        srand((int)wxGetLocalTime());

        WMain* wmain = new WMain(nullptr);
        SetTopWindow(wmain);
        wmain->Show();

        return true;
    }

    virtual int OnExit()
    {
        SDL_Quit();

        // return the standard exit code
        return wxApp::OnExit();
    }
};

IMPLEMENT_APP(MyApp)

// ----------------------------------------------------------------------------
