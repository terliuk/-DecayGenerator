cmake_minimum_required (VERSION 2.6)
set(CMAKE_CXX_COMPILER g++)
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

project (DecayGenerator)


include_directories("${PROJECT_SOURCE_DIR}/include")
add_library(DecayGenerator SHARED "${PROJECT_SOURCE_DIR}/src/DecayGenerator.cxx")


find_package(GSL REQUIRED)
find_package(Boost REQUIRED COMPONENTS python system numpy)
find_package(PythonLibs REQUIRED )
find_package(Python REQUIRED COMPONENTS NumPy)


set(LIBRARY_OUTPUT_PATH "lib")
set(INCLUDE_OUTPUT_PATH "include")
message("Boost includes : " ${Boost_INCLUDE_DIR} )
message("Boost libraries : "  ${Boost_LIBRARIES} )
message("Boost library dir : " ${Boost_LIBRARY_DIR} )

message("Python libraries : " ${PYTHON_LIBRARIES})
message("Python includes : " ${PYTHON_INCLUDE_DIRS} )

include_directories(${Boost_INCLUDE_DIR} ${PYTHON_INCLUDE_DIRS})
target_link_libraries (DecayGenerator ${Boost_LIBRARIES})
target_link_libraries(DecayGenerator ${PYTHON_LIBRARIES} )
target_link_libraries (DecayGenerator GSL::gsl)

#link_directories(${Boost_LIBRARY_DIR})

