# -*- mode: autoconf -*-

AC_PREREQ(2.59)

AC_INIT(sound-of-sorting, 0.6.5)
AC_CONFIG_SRCDIR(src/WMain.cpp)
AC_CONFIG_AUX_DIR(acscripts)

AC_CANONICAL_TARGET

AM_INIT_AUTOMAKE(foreign subdir-objects)
AM_MAINTAINER_MODE

# check for Win32 system and set compilation flag.

AC_MSG_CHECKING(building for Win32)
case "$target_os" in
    *cygwin* | *mingw32*)
        ON_WIN32=true
        AC_MSG_RESULT(yes)
        ;;
    *)
        ON_WIN32=false
        AC_MSG_RESULT(no)
        ;;
esac
AM_CONDITIONAL(ON_WIN32, $ON_WIN32)

# set debug info flag if no optimization flags are set.

if test "$CFLAGS" == ""; then
    CFLAGS="-g -Wno-unused-local-typedefs"
fi
if test "$CXXFLAGS" == ""; then
    CXXFLAGS="-O2 -Wno-unused-local-typedefs"
fi

# enable full optimization by configure switch

AC_ARG_ENABLE(optimize,
	AS_HELP_STRING([--enable-optimize],
		       [Build with full optimization @<:@default=no@:>@]),
        [ case "${enableval}" in
          yes)
                CFLAGS="$CFLAGS -O3";
                CXXFLAGS="$CXXFLAGS -O3";
          ;;
          no) ;;
          *) AC_MSG_ERROR(bad value ${enableval} for --enable-optimize) ;;
          esac ],
        [ optimize=false ])

# check for programs.

AC_PROG_CC
AC_PROG_CXX
AC_PROG_RANLIB

AC_LANG([C++])

# check for SDL

SDL_VERSION=2
AM_PATH_SDL2($SDL_VERSION,
            :,
	    AC_MSG_ERROR([*** SDL version $SDL_VERSION not found!])
)

# check for wxWidgets 2.8.9 or later

AM_OPTIONS_WXCONFIG

AM_PATH_WXCONFIG(2.8.9, [wxWin=1], [wxWin=0], [adv,core,base])

if test "$wxWin" != 1; then
    AC_MSG_ERROR([
	wxWidgets does not seem to be installed on your system.

        If you think wxWidgets >= 2.8.9 is installed,
	please check that wx-config is in path, the directory
	where wxWidgets libraries are installed (returned by
	'wx-config --libs' command) is in LD_LIBRARY_PATH or
	equivalent variable.
    ])
fi

# Test for function "QueryPerformanceCounter" on Win32

if $ON_WIN32; then

   AC_CHECK_DECL(QueryPerformanceCounter,
        [HavePerformanceCounter=1],
        [HavePerformanceCounter=0],
        [#include <windows.h>])

   if test "$HavePerformanceCounter" != 1; then
      AC_MSG_WARN([
        Could not find QueryPerformanceCounter() function. This
        reduces the precision of delays at array accesses.
      ])
   fi

   AC_DEFINE_UNQUOTED(MSW_PERFORMANCECOUNTER, $HavePerformanceCounter)
fi

# Win32 tools to add icons to the programs

if $ON_WIN32; then
   AM_CONDITIONAL(GOT_RESCOMP, test "$WX_RESCOMP" != "")
else
   AM_CONDITIONAL(GOT_RESCOMP, false)
fi

# transform Makefiles

AC_OUTPUT([Makefile src/Makefile])
