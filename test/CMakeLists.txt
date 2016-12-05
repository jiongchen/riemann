link_directories($ENV{HOME}/usr/lib)

set(test_source test_diff.cc test_ipopt.cc test_ipopt_wrapper.cc)

foreach(testfile ${test_source})
  string(REPLACE ".cc" "" testname ${testfile})
  add_executable(${testname} ${testfile})
  target_link_libraries(${testname} gtest gtest_main riemann ${IPOPT_LIBRARIES})
endforeach()
