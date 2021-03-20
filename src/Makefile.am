# $Id: Makefile.am 2 2007-01-05 14:04:27Z tb $

bin_PROGRAMS = sound-of-sorting

COMMON = SortArray.cpp SortArray.h \
	SortSound.cpp \
	SortAlgo.cpp SortAlgo.h \
	algorithms/timsort.cpp \
	algorithms/wikisort.cpp

sound_of_sorting_SOURCES = \
	WMain.cpp WMain.h \
	wxg/WMain_wxg.cpp wxg/WMain_wxg.h \
	wxg/WAbout_wxg.cpp wxg/WAbout_wxg.h \
	WSortView.cpp WSortView.h \
	wxClickText.cpp wxClickText.h \
	$(COMMON)

noinst_PROGRAMS = sorting-test

TESTS = sorting-test

sorting_test_SOURCES = \
	SortTest.cpp \
	$(COMMON)

AM_CXXFLAGS = -W -Wall @WX_CXXFLAGS@ @SDL_CFLAGS@
LDADD = @WX_LIBS@ @SDL_LIBS@

if GOT_RESCOMP
resources.o: resources.rc
	$(WX_RESCOMP) $(srcdir)/resources.rc resources.o

LDADD += resources.o
endif

EXTRA_DIST = \
	resources.rc sos.ico sos.xpm \
	wxg/sos.wxg
