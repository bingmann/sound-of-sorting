###############################################################################
# CMakeLists.txt
#
# Build script for The Sound of Sorting
#
###############################################################################
# Copyright (C) 2024 Timo Bingmann <tb@panthema.net>
#
# This program is free software: you can redistribute it and/or modify it under
# the terms of the GNU General Public License as published by the Free Software
# Foundation, either version 3 of the License, or (at your option) any later
# version.
#
# This program is distributed in the hope that it will be useful, but WITHOUT
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
# FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
# details.
#
# You should have received a copy of the GNU General Public License along with
# this program.  If not, see <http://www.gnu.org/licenses/>.
###############################################################################

include_directories(${SDL2_INCLUDE_DIRS})

add_library(sorting-algorithms
  SortAlgo.cpp
  SortAlgo.h
  SortArray.cpp
  SortArray.h
  SortSound.cpp
  algorithms/timsort.cpp
  algorithms/wikisort.cpp
)

add_executable(sound-of-sorting
  WMain.cpp
  WMain.h
  WSortView.cpp
  WSortView.h
  wxClickText.cpp
  wxClickText.h
  wxg/WAbout_wxg.cpp
  wxg/WAbout_wxg.h
  wxg/WMain_wxg.cpp
  wxg/WMain_wxg.h
  resources.rc
)

target_link_libraries(sound-of-sorting
  sorting-algorithms ${wxWidgets_LIBRARIES} ${SDL2_LIBRARIES})

add_executable(sorting-test
  SortTest.cpp
)

target_link_libraries(sorting-test
  sorting-algorithms ${wxWidgets_LIBRARIES})

add_test(sorting-test sorting-test)
