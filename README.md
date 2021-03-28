# The Sound of Sorting - Visualization and "Audibilization" of Sorting Algorithms

Sorting algorithms are an **essential chapter** in undergraduate computer
science education. Due to their easy to explain nature and fairly
straight-forward analysis, this set of algorithms offers a convenient
introduction to the methods and techniques of theoretical computer science and
algorithm analysis.

This source archive presents my own **demo program for sortings algorithms**,
called **"The Sound of Sorting"**, which both **visualizes** the algorithms
internals and their operations, and **generates sound effects** from the values
being compared.

The demo is implemented using the **cross-platform** toolkits wxWidgets and
SDL, can be executed on Windows, Linux and Mac, and **runs in real time**.

There are many resources on the Internet about sorting algorithms, among them
are <a href="http://en.wikipedia.org/wiki/Sorting_algorithm">Wikipedia</a>, <a
href="http://www.sorting-algorithms.com">animated sorting algorithms</a> by
David R. Martin and various Java applets by many college or university staff.

## Website and License

The current source package and some binaries can be downloaded from
https://github.com/yg8ijvjvjv/sound-of-sorting/

The program and code is published under the GNU General Public License v3
(GPL), which can also be found in the file COPYING. A few of the sorting
algorithms' implementations were written by other authors and may have
different licenses.

## Usage

The Sound of Sorting demo program is very intuitive to use. It contains many
sorting algorithms, which are selectable via the list box on the right. For the
quick sort variants the pivot rule can be selected separately.

When pressing "Run", the algorithm is started on selected input. "Step" will
stop a running algorithm, or start a new one, and halt it after performing one
operation. With "Reset" a running algorithm is stopped, and a new random input
is created. When "Sound" is activated, the program will generate sound effects
on the default audio output. The "Speed" slider changes the algorithms
execution speed by adding a delay for each array access.

*Both Windows and Linux versions have a higher time resolution of 0.001 ms!*

The algorithm's visualization contains mostly white bars representing the value
of the array position corresponding to the x-axis. When the algorithm gets or
sets an array item, the white bar runs red for one algorithmic step. A swap
operation is represented by two bars turning red and their values being
exchanged. Some algorithms specially colorize bars which represent indexes,
pointers or other information about the internal mechanisms of the algorithm.

Both value comparisons and array accesses are counted and shown in real time.
The comparison counter includes ternary comparisons as just one operation. Due
to algorithms often using extra memory or local variables, the array access
counter highly depends on the actual algorithm implementation.

The generated sound effects depend on the values being compared. Only
comparisons yield sound output (except for in distributive sorts)! The length of
each comparison's sound effect can be modified using the "Sound Sustain"
slider. The frequency of the sound is calculated from the compared values. The
sound wave itself is square wave and modulated with a smooth curve. This
yields the "8-bit game tune" feeling. An item value is scaled (with double
precision) to the frequency range 120 Hz - 1300 Hz, which is large but not too
high to be annoying.

## Source Code Overview and Implementation Notes

The demo program uses the wxWidgets toolkit for cross-platform dialogs and
painting. For sound output, the audio component of the cross-platform SDL
library is used.

All sources resides in `src/`. The main window's GUI functions are grouped in
`WMain.h/cpp`. The sorting visualization, including the instrumented array
class and painting methods are in `WSortView.h/cpp`.

`SortAlgo.h/cpp` contains all sorting algorithms. These were mainly modified to
operate on a `WSortView` object, which exposes most of the usual array
operators such as `operator[]`, and many special functions to create nicer
visualizations. Most notable among these, are a special `swap()` function and
`mark()` to colorize bars. There is also `watch()` to do *live tracking* of
indexes stored as local variables (use this with care)!

Comparison counting and sound effects are signaled by the operators of
`ArrayItem`, which is the item class of the instrumented array `WSortView`. As
such, *all comparisons* of the sorting algorithms will be intercepted by this
class.

On each comparison, the values are used to generate sound. All sound generating
methods are in `SortSound.cpp`. The main class here is an `Oscillator`, which
generates an enveloped square waveform of a specific frequency. Oscillators
are mixed together for the output sound. The output volume is scaled
automatically depending on the number of oscillators active.

If you want to contribute a sorting algorithm, please notify me.

## Exits

Written 2013-05-22 by Timo Bingmann
