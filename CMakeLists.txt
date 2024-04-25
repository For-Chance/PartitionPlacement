cmake_minimum_required (VERSION 3.17)
project(libDrawingsProcessing)
# Find package if using VCPKG
find_package(CGAL CONFIG OPTIONAL_COMPONENTS Qt5)
find_package(jsoncpp CONFIG REQUIRED)

# Show packages' version and dir
message("jsoncpp:\t${jsoncpp_VERSION},\t\t${jsoncpp_DIR}")
message("CGAL: \t\t${CGAL_VERSION},\t\t${CGAL_DIR}")
message("Boost: \t\t${Boost_VERSION},\t\t${Boost_DIR}")
message("Qt5: \t\t${Qt5_VERSION},\t${Qt5_DIR}")
string(REGEX MATCH "^[0-9]+" CGAL_Major "${CGAL_VERSION}")  # extract major version of CGAL

# Control building args
set(CMAKE_INSTALL_PREFIX "${PROJECT_BINARY_DIR}/install" CACHE STRING "" FORCE)
set(EXECUTABLE_OUTPUT_PATH "${PROJECT_BINARY_DIR}/bin")
set(CMAKE_CXX_STANDARD 11)
add_definitions(-DCGAL_USE_BASIC_VIEWER -DQT_NO_KEYWORDS)
add_compile_options("$<$<C_COMPILER_ID:MSVC>:/source-charset:utf-8>")

#### define library and test example ####


# 定义一个函数来添加一个库和其依赖
function(add_library_with_test PROJECT_NAME DEPENDENCY_LIBS)
    set(TEST_SUFFIX "Test")
    set(PROJECT_TEST_NAME ${PROJECT_NAME}${TEST_SUFFIX})

    file(GLOB header_impl "include/${PROJECT_NAME}/*.hpp" "include/${PROJECT_NAME}/*.h")
    file(GLOB lib_src "src/${PROJECT_NAME}.cpp")
    file(GLOB test_headers "test/${PROJECT_NAME}/*.hpp" "test/${PROJECT_NAME}/*.h")
    file(GLOB test_src "test/${PROJECT_NAME}/*.cpp")

    add_library(${PROJECT_NAME} STATIC ${lib_src} ${header_impl} )
    target_include_directories(${PROJECT_NAME} PUBLIC include "include/${PROJECT_NAME}")
    target_link_libraries(${PROJECT_NAME} PUBLIC CGAL::CGAL jsoncpp_lib CGAL::CGAL_Qt5)
    foreach(LIB IN LISTS DEPENDENCY_LIBS)
        target_link_libraries(${PROJECT_NAME} PUBLIC ${LIB})
    endforeach()

    if(${Qt5_FOUND})
        add_executable (${PROJECT_TEST_NAME} "test/${PROJECT_NAME}/main.cpp"  ${test_headers} ${test_src})
        target_include_directories(${PROJECT_TEST_NAME} PRIVATE test "test/${PROJECT_NAME}")
        target_link_libraries(${PROJECT_TEST_NAME} PUBLIC CGAL::CGAL_Qt5 PUBLIC ${PROJECT_NAME})
    endif()
endfunction()

add_library_with_test("SimplifyBoundary" "")
add_library_with_test("CenterLine" "SimplifyBoundary")
add_library_with_test("PartitionPlacement" "SimplifyBoundary;CenterLine")
