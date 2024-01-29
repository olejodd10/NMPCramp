add_subdirectory(src)


add_test(NAME TestNMPCramp_Compiles COMMAND make)


add_test(NAME TestLmpcConstant_BadPath_ShouldFail COMMAND TestLmpcConstant /bad/path /bad/path2 10 100)
set_tests_properties(TestLmpcConstant_BadPath_ShouldFail PROPERTIES WILL_FAIL TRUE)

add_test(NAME TestLmpcConstant_MissingArguments_ShouldFail COMMAND TestLmpcConstant /bad/path /bad/path2 10)
set_tests_properties(TestLmpcConstant_MissingArguments_ShouldFail PROPERTIES WILL_FAIL TRUE)

add_test(NAME TestLmpcConstant_QPrampExample1 COMMAND TestLmpcConstant ${CMAKE_CURRENT_SOURCE_DIR}/data/QPramp/example1 ${CMAKE_CURRENT_SOURCE_DIR}/data/QPramp/example1 10 100)
add_test(NAME TestLmpcConstant_QPrampExample1_Verify_Uout COMMAND python ${CMAKE_CURRENT_SOURCE_DIR}/scripts/compare_csvs.py ${CMAKE_CURRENT_SOURCE_DIR}/data/QPramp/example1/uref.csv ${CMAKE_CURRENT_SOURCE_DIR}/data/QPramp/example1/uout.csv 4e-5)
add_test(NAME TestLmpcConstant_QPrampExample1_Verify_Xout COMMAND python ${CMAKE_CURRENT_SOURCE_DIR}/scripts/compare_csvs.py ${CMAKE_CURRENT_SOURCE_DIR}/data/QPramp/example1/xref.csv ${CMAKE_CURRENT_SOURCE_DIR}/data/QPramp/example1/xout.csv 5e-5)

add_test(NAME TestLmpcConstant_QPrampExample2 COMMAND TestLmpcConstant ${CMAKE_CURRENT_SOURCE_DIR}/data/QPramp/example2 ${CMAKE_CURRENT_SOURCE_DIR}/data/QPramp/example2 10 100)
add_test(NAME TestLmpcConstant_QPrampExample2_Verify_Uout COMMAND python ${CMAKE_CURRENT_SOURCE_DIR}/scripts/compare_csvs.py ${CMAKE_CURRENT_SOURCE_DIR}/data/QPramp/example2/uref.csv ${CMAKE_CURRENT_SOURCE_DIR}/data/QPramp/example2/uout.csv 3e-2)
add_test(NAME TestLmpcConstant_QPrampExample2_Verify_Xout COMMAND python ${CMAKE_CURRENT_SOURCE_DIR}/scripts/compare_csvs.py ${CMAKE_CURRENT_SOURCE_DIR}/data/QPramp/example2/xref.csv ${CMAKE_CURRENT_SOURCE_DIR}/data/QPramp/example2/xout.csv 5e-3)

add_test(NAME TestLmpcConstant_QPrampExample3 COMMAND TestLmpcConstant ${CMAKE_CURRENT_SOURCE_DIR}/data/QPramp/example3 ${CMAKE_CURRENT_SOURCE_DIR}/data/QPramp/example3 10 100)
add_test(NAME TestLmpcConstant_QPrampExample3_Verify_Uout COMMAND python ${CMAKE_CURRENT_SOURCE_DIR}/scripts/compare_csvs.py ${CMAKE_CURRENT_SOURCE_DIR}/data/QPramp/example3/uref.csv ${CMAKE_CURRENT_SOURCE_DIR}/data/QPramp/example3/uout.csv 6e-6)
add_test(NAME TestLmpcConstant_QPrampExample3_Verify_Xout COMMAND python ${CMAKE_CURRENT_SOURCE_DIR}/scripts/compare_csvs.py ${CMAKE_CURRENT_SOURCE_DIR}/data/QPramp/example3/xref.csv ${CMAKE_CURRENT_SOURCE_DIR}/data/QPramp/example3/xout.csv 5e-5)
