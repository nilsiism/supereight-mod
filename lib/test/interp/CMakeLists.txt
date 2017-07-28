set(PROJECT_TEST_NAME interp_unittest)
add_executable(${PROJECT_TEST_NAME} interp_unittest.cpp)
target_include_directories(${PROJECT_TEST_NAME} PUBLIC ${GTEST_INCLUDE_DIRS})
target_link_libraries(${PROJECT_TEST_NAME} ${GTEST_BOTH_LIBRARIES} pthread)

GTEST_ADD_TESTS(${PROJECT_TEST_NAME} "" ${PROJECT_TEST_NAME}.cpp)
