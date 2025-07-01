cmake_minimum_required(VERSION 3.10)
project(bamcleaner)

# Use C++17
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

find_package(ZLIB REQUIRED)

#------------------------------------------------------------------------------
# Options for SeqLib
option(SEQ_LIB_INCLUDE_DIR "Path to SeqLib include directory" "")
option(SEQ_LIB_LIBRARY     "Path to SeqLib library file"      "")

# Options for HTSlib
option(HTSLIB_INCLUDE_DIR "Path to htslib include directory" "")
option(HTSLIB_LIBRARY     "Path to htslib library file"      "")

#------------------------------------------------------------------------------
# Validate user-provided paths
if(NOT SEQ_LIB_INCLUDE_DIR)
  message(FATAL_ERROR "You must set SEQ_LIB_INCLUDE_DIR to the SeqLib include path.")
endif()
if(NOT SEQ_LIB_LIBRARY)
  message(FATAL_ERROR "You must set SEQ_LIB_LIBRARY to the SeqLib library file.")
endif()
if(NOT HTSLIB_INCLUDE_DIR)
  message(FATAL_ERROR "You must set HTSLIB_INCLUDE_DIR to the htslib include path.")
endif()
if(NOT HTSLIB_LIBRARY)
  message(FATAL_ERROR "You must set HTSLIB_LIBRARY to the htslib library file.")
endif()

#------------------------------------------------------------------------------
# Include headers
include_directories(
  ${SEQ_LIB_INCLUDE_DIR}
  ${HTSLIB_INCLUDE_DIR}
)

#------------------------------------------------------------------------------
# Include robin-map (header-only)
set(ROBIN_MAP_ROOT "$ENV{HOME}/git/robin-map")
add_subdirectory(${ROBIN_MAP_ROOT} ${CMAKE_BINARY_DIR}/robin-map-build)

#------------------------------------------------------------------------------
# Define executable
add_executable(bamcleaner
    bamcleaner.cpp
)

# Include paths for bamcleaner
target_include_directories(bamcleaner PRIVATE
    ${SEQ_LIB_INCLUDE_DIR}
    ${HTSLIB_INCLUDE_DIR}
    ${ROBIN_MAP_ROOT}/include
)

# Link libraries: SeqLib, HTSlib
target_link_libraries(bamcleaner PRIVATE
    ${SEQ_LIB_LIBRARY}
    ${HTSLIB_LIBRARY}
    ZLIB::ZLIB            
	
)

#------------------------------------------------------------------------------
# Installation
install(TARGETS bamcleaner
    RUNTIME DESTINATION bin
)
