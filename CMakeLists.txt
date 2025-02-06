cmake_minimum_required(VERSION 3.31)
project(IAPWS95)
add_library(IAPWS95 STATIC)
target_sources(IAPWS95 PUBLIC ${PROJECT_SOURCE_DIR}/src/LibIAPWS95.cpp
                              ${PROJECT_SOURCE_DIR}/src/InterCSpline.cpp)
target_include_directories(IAPWS95 PUBLIC "${PROJECT_SOURCE_DIR}/include")