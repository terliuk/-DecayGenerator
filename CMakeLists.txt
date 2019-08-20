cmake_minimum_required (VERSION 2.6)

project (main)


include_directories("${PROJECT_SOURCE_DIR}/include")
add_library(DecayGenerator "${PROJECT_SOURCE_DIR}/src/DecayGenerator.cxx")

add_executable(main main.cxx)
find_package(GSL REQUIRED)
target_link_libraries (main DecayGenerator)
target_link_libraries (main GSL::gsl)


