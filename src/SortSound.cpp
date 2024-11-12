/******************************************************************************
 * src/SortSound.cpp
 *
 * This file contains the audio callback which generates all sound.
 *
 * The sound is created by mixing many many oscillators, whose freqencies are
 * defined by the values compared. Each comparison adds two oscillators.
 *
 * In the final version the oscillators generate a triangular waveform and are
 * enveloped using a ADSR function. The total time an oscillator yields sound
 * is defined by the "Sound Sustain" user slider.
 *
 * The callback function also automatically scales the volume when many
 * oscillators overlap. In effect, this is dynamic range compression.
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
#include <SDL.h>
#include <limits>

// *** All time counters in the sound system are in sample units.

static const size_t s_samplerate = 44100;

/// global sound on/off switch
bool g_sound_on = false;

/// multiplying g_delay with this value yields the duration a sound is sustained
double g_sound_sustain = 2.0;

/// limit the number of oscillators to avoid overloading the callback
static const size_t s_max_oscillators = 512;

/// Oscillator generating sine or triangle waves
class Oscillator
{
protected:
    /// frequency of generated wave
    double              m_freq;

    /// start and end of wave in sample time
    size_t              m_tstart, m_tend;

    /// duration of oscillation note
    size_t              m_duration;

public:
    /// construct new oscillator
    Oscillator(double freq, size_t tstart, size_t duration = 44100 / 8)
        : m_freq(freq), m_tstart(tstart),
          m_tend( m_tstart + duration ),
          m_duration(duration)
    {
    }

    // *** Wave Forms

    /// simple sine wave
    static double wave_sin(double x)
    {
        return sin(x * 2*M_PI);
    }

    /// sin^3 wave
    static double wave_sin3(double x)
    {
        double s = sin(x * 2*M_PI);
        return s * s * s;
    }

    /// triangle wave
    static double wave_triangle(double x)
    {
        x = fmod(x, 1.0);

        if (x <= 0.25) return 4.0 * x;
        if (x <= 0.75) return 2.0 - 4.0 * x;
        return 4.0 * x - 4.0;
    }

    /// picking a waveform
    static double wave(double x)
    {
        //return wave_sin(x);
        //return wave_sin3(x);
        return wave_triangle(x);
    }

    // *** Envelope

    /// envelope applied to wave (uses ADSR)
    double envelope(size_t i) const
    {
        double x = (double)i / m_duration;
        if (x > 1.0) x = 1.0;

        // simple envelope functions:

        //return 1.0 - x;
        //return cos(M_PI_2 * x);

        // *** ADSR envelope

        static const double attack = 0.025; // percentage of duration
        static const double decay = 0.1;    // percentage of duration
        static const double sustain = 0.9;  // percentage of amplitude
        static const double release = 0.3;  // percentage of duration

        if (x < attack)
            return 1.0 / attack * x;

        if (x < attack + decay)
            return 1.0 - (x - attack) / decay * (1.0 - sustain);

        if (x < 1.0 - release)
            return sustain;

        return sustain / release * (1.0 - x);
    }

    // *** Generate Wave and Mix

    /// mix in the output of this oscillator on the wave output
    void mix(double* data, int size, size_t p) const
    {
        for (int i = 0; i < size; ++i)
        {
            if (p+i < m_tstart) continue;
            if (p+i >= m_tend) break;

            size_t trel = (p + i - m_tstart);

            data[i] += envelope(trel) * wave((double)trel / s_samplerate * m_freq);
        }
    }

    /// return start time
    size_t tstart() const
    { return m_tstart; }

    /// true if the oscillator is silent at time p
    bool is_done(size_t p) const
    {
        return (p >= m_tend);
    }
};

/// array of oscillators
static std::vector<Oscillator> s_osclist;

/// global timestamp of callback in sound sample units
static size_t s_pos = 0;

/// add an oscillator to the list (reuse finished ones)
static void add_oscillator(double freq, size_t p, size_t pstart, size_t duration)
{
    // find free oscillator
    size_t oldest = 0, toldest = std::numeric_limits<size_t>::max();
    for (size_t i = 0; i < s_osclist.size(); ++i)
    {
        if (s_osclist[i].is_done(p))
        {
            s_osclist[i] = Oscillator(freq, pstart, duration);
            return;
        }

        if (s_osclist[i].tstart() < toldest) {
            oldest = i;
            toldest = s_osclist[i].tstart();
        }
    }

    if (s_osclist.size() < s_max_oscillators)
    {
        // add new one
        s_osclist.push_back( Oscillator(freq, pstart, duration) );
    }
    else
    {
        // replace oldest oscillator
        s_osclist[oldest] = Oscillator(freq, pstart, duration);
    }
}

/// list of array accesses since last callback
static std::vector<unsigned int> s_access_list;

// mutex to s_access_list
static wxMutex s_mutex_access_list;

/// "public" function to add a new array access
void SoundAccess(size_t i)
{
    if (!g_sound_on) return;

    wxMutexLocker lock(s_mutex_access_list);
    ASSERT(lock.IsOk());

    s_access_list.push_back(i);
}

/// function mapping array index (normalized to [0,1]) to frequency
static double arrayindex_to_frequency(double aindex)
{
    return 120 + 1200 * (aindex*aindex);
}

/// reset internal sound data (called from main thread)
void SoundReset()
{
    wxMutexLocker lock(s_mutex_access_list);
    ASSERT(lock.IsOk());

    s_pos = 0;
    s_osclist.clear();
}

/// sound generator callback run by SDL
void SoundCallback(void* udata, unsigned char* stream, int len)
{
    if (!g_sound_on) {
        memset(stream, 0, len);
        return;
    }

    // current sample time (32-bit size_t wraps after 27 hours, 64-bit size_t
    // wraps after 13 million years).
    size_t& p = s_pos;

    // reference to sortview
    WSortView& sv = *reinterpret_cast<WSortView*>(udata);

    // we use 16-bit mono output at 44.1 kHz
    int16_t* data = (int16_t*)stream;
    size_t size = len / sizeof(int16_t);

    // fetch new access list and create oscillators
    {
        wxMutexLocker lock(s_mutex_access_list);
        ASSERT(lock.IsOk());

        // spread out access list over time of one callback
        double pscale = (double)size / s_access_list.size();

        for (size_t i = 0; i < s_access_list.size(); ++i)
        {
            double relindex = s_access_list[i] / (double)sv.m_array.array_max();
            double freq = arrayindex_to_frequency(relindex);

            add_oscillator( freq, p, p + i * pscale,
                            g_delay / 1000.0 * g_sound_sustain * s_samplerate );
        }

        s_access_list.clear();
    }

    // calculate waveform
    std::vector<double> wave(size, 0.0);
    size_t wavecount = 0;

    for (std::vector<Oscillator>::const_iterator it = s_osclist.begin();
         it != s_osclist.end(); ++it)
    {
        if (!it->is_done(p))
        {
            it->mix(wave.data(), size, p);
            ++wavecount;
        }
    }

    // scale output with number of waves mixed

    if (wavecount == 0)
    {
        // set zero, the function below messes up when vol = 0.0
        memset(stream, 0, len);
    }
    else
    {
        // count maximum wave amplitude
        double vol = *std::max_element(wave.begin(), wave.end());

        static double oldvol = 1.0;

        if (vol > oldvol) {
            // immediately ramp down volume
        }
        else {
            // but slowly ramp up volume
            vol = 0.9 * oldvol;
        }

        // convert waveform to samples, with ramping of volume

        for (size_t i = 0; i < size; ++i)
        {
            int32_t v = 24000.0 * wave[i] / (oldvol + (vol - oldvol) * (i / (double)size));

            if (v > 32200) {
                //std::cout << "clip " << p << "\n";
                v = 32200;
            }
            if (v < -32200) {
                //std::cout << "clip " << p << "\n";
                v = -32200;
            }

            data[i] = v;
        }

        oldvol = vol;
    }

    // advance sample timestamp
    p += size;
}
