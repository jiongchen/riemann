link_directories($ENV{HOME}/usr/lib)

set(test_source test_diff.cc)

foreach(testfile ${test_source})
  string(REPLACE ".cc" "" testname ${testfile})
  add_executable(${testname} ${testfile})
  target_link_libraries(${testname} gtest gtest_main riemann)
endforeach()
