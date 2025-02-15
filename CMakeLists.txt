cmake_minimum_required(VERSION 3.31)

project(IAPWS76)
add_library(IAPWS76 STATIC)
target_sources(IAPWS76 PUBLIC ${PROJECT_SOURCE_DIR}/src/BaseIAPWS.cpp
                              ${PROJECT_SOURCE_DIR}/src/LibIAPWS76.cpp)
target_include_directories(IAPWS76 PUBLIC "${PROJECT_SOURCE_DIR}/include")

project(IAPWS95)
add_library(IAPWS95 STATIC)
target_sources(IAPWS95 PUBLIC ${PROJECT_SOURCE_DIR}/src/BaseIAPWS.cpp
                              ${PROJECT_SOURCE_DIR}/src/LibIAPWS95.cpp
                              ${PROJECT_SOURCE_DIR}/src/InterCSpline.cpp)
target_include_directories(IAPWS95 PUBLIC "${PROJECT_SOURCE_DIR}/include")

project(IAPWS06)
add_library(IAPWS06 STATIC)
target_sources(IAPWS06 PUBLIC ${PROJECT_SOURCE_DIR}/src/BaseIAPWS.cpp
                              ${PROJECT_SOURCE_DIR}/src/LibIAPWS06.cpp)
target_include_directories(IAPWS06 PUBLIC "${PROJECT_SOURCE_DIR}/include")

project(make_tab_IAPWS95)
add_executable(make_tab_IAPWS95 ${PROJECT_SOURCE_DIR}/src/make_tab_IAPWS95.cpp)
set_target_properties(make_tab_IAPWS95 PROPERTIES OUTPUT_NAME "make_tab_IAPWS95.exec")
target_link_libraries(make_tab_IAPWS95 PUBLIC IAPWS95)
target_include_directories(make_tab_IAPWS95 PUBLIC "${PROJECT_SOURCE_DIR}/include")

project(test_IAPWS76)
add_executable(test_IAPWS76 ${PROJECT_SOURCE_DIR}/src/test_IAPWS76.cpp)
set_target_properties(test_IAPWS76 PROPERTIES OUTPUT_NAME "test_IAPWS76.exec")
target_link_libraries(test_IAPWS76 PUBLIC IAPWS76)
target_include_directories(test_IAPWS76 PUBLIC "${PROJECT_SOURCE_DIR}/include")

project(test_IAPWS95)
add_executable(test_IAPWS95 ${PROJECT_SOURCE_DIR}/src/test_IAPWS95.cpp)
set_target_properties(test_IAPWS95 PROPERTIES OUTPUT_NAME "test_IAPWS95.exec")
target_link_libraries(test_IAPWS95 PUBLIC IAPWS95)
target_include_directories(test_IAPWS95 PUBLIC "${PROJECT_SOURCE_DIR}/include")

project(test_IAPWS06)
add_executable(test_IAPWS06 ${PROJECT_SOURCE_DIR}/src/test_IAPWS06.cpp)
set_target_properties(test_IAPWS06 PROPERTIES OUTPUT_NAME "test_IAPWS06.exec")
target_link_libraries(test_IAPWS06 PUBLIC IAPWS06)
target_include_directories(test_IAPWS06 PUBLIC "${PROJECT_SOURCE_DIR}/include")
