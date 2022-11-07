cmake_minimum_required (VERSION 2.6)
set(CMAKE_CXX_COMPILER g++)
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

project (DecayGenerator)


include_directories("${PROJECT_SOURCE_DIR}/include")
add_library(DecayGenerator SHARED "${PROJECT_SOURCE_DIR}/src/DecayGenerator.cxx")

find_package(GSL REQUIRED)

message("Will build python interfaces")
# 
find_package(Python REQUIRED COMPONENTS NumPy)
include_directories(${Python_INCLUDE_DIRS})
if(NOT DEFINED Boost_PYTHON_TYPE)
  set(Boost_PYTHON_TYPE python)
  set(Boost_NUMPY_TYPE numpy)
  # Hack for some old Boost CMake modules
  set(_Boost_PYTHON${Python_VERSION_MAJOR}_HEADERS "boost/python.hpp")
  set(_Boost_PYTHON${Python_VERSION_MAJOR}${Python_VERSION_MINOR}_HEADERS "boost/python.hpp")
  find_package(Boost QUIET COMPONENTS python${Python_VERSION_MAJOR}${Python_VERSION_MINOR})
  if (${Boost_PYTHON${Python_VERSION_MAJOR}${Python_VERSION_MINOR}_FOUND})
    set(Boost_PYTHON_TYPE python${Python_VERSION_MAJOR}${Python_VERSION_MINOR})
    set(Boost_NUMPY_TYPE numpy${Python_VERSION_MAJOR}${Python_VERSION_MINOR})
  else()
  find_package(Boost QUIET COMPONENTS python${Python_VERSION_MAJOR})
  if (${Boost_PYTHON${Python_VERSION_MAJOR}_FOUND})
    set(Boost_PYTHON_TYPE python${Python_VERSION_MAJOR})
    set(Boost_NUMPY_TYPE numpy${Python_VERSION_MAJOR})
    endif()
  endif()
endif()

message("Will search for boost python component : " ${Boost_PYTHON_TYPE})
find_package(Boost REQUIRED COMPONENTS ${Boost_PYTHON_TYPE} system ${Boost_NUMPY_TYPE})
include_directories(${Boost_INCLUDE_DIR} ${PYTHON_INCLUDE_DIRS})

set(LIBRARY_OUTPUT_PATH "lib")
set(INCLUDE_OUTPUT_PATH "include")

target_link_libraries (DecayGenerator ${Boost_LIBRARIES})
target_link_libraries(DecayGenerator ${PYTHON_LIBRARIES} )
target_link_libraries (DecayGenerator GSL::gsl)
