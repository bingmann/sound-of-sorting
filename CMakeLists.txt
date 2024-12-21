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

cmake_minimum_required(VERSION 3.10)

project(sound-of-sorting)

add_definitions(-DPACKAGE_VERSION="0.7.0")

# prohibit in-source builds
if("${PROJECT_SOURCE_DIR}" STREQUAL "${PROJECT_BINARY_DIR}")
  message(SEND_ERROR "In-source builds are not allowed.")
endif()

# default to Debug building for single-config generators
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  message("Defaulting CMAKE_BUILD_TYPE to Debug")
  set(CMAKE_BUILD_TYPE "Debug")
endif()

# enable warnings
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -W -Wall")

# enable use of "make test"
enable_testing()

### Find wxWidgets Libraries ###

find_package(wxWidgets REQUIRED COMPONENTS adv core base)
if(wxWidgets_USE_FILE)
  include(${wxWidgets_USE_FILE})
endif()

### Find SDL2 Libraries ###

find_package(SDL2 REQUIRED)

### Done with Libraries ###

# descend into source
add_subdirectory(src)

################################################################################
