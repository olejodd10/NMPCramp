add_subdirectory(src)


add_test(NAME TestLtiMpc_BadPath_ShouldFail COMMAND TestLtiMpc /bad/path /bad/path2 10 100)
set_tests_properties(TestLtiMpc_BadPath_ShouldFail PROPERTIES WILL_FAIL TRUE)

add_test(NAME TestLtiMpc_MissingArguments_ShouldFail COMMAND TestLtiMpc /bad/path /bad/path2 10)
set_tests_properties(TestLtiMpc_MissingArguments_ShouldFail PROPERTIES WILL_FAIL TRUE)


add_test(NAME TestLtiMpc_QPrampExample1 COMMAND TestLtiMpc ${CMAKE_CURRENT_SOURCE_DIR}/data/QPramp/example1 ${CMAKE_CURRENT_SOURCE_DIR}/data/QPramp/example1 10 100)
add_test(NAME TestLtiMpc_QPrampExample1_Verify_Uout COMMAND python ${CMAKE_CURRENT_SOURCE_DIR}/scripts/compare_csvs.py ${CMAKE_CURRENT_SOURCE_DIR}/data/QPramp/example1/uref.csv ${CMAKE_CURRENT_SOURCE_DIR}/data/QPramp/example1/uout.csv 4e-5)
add_test(NAME TestLtiMpc_QPrampExample1_Verify_Xout COMMAND python ${CMAKE_CURRENT_SOURCE_DIR}/scripts/compare_csvs.py ${CMAKE_CURRENT_SOURCE_DIR}/data/QPramp/example1/xref.csv ${CMAKE_CURRENT_SOURCE_DIR}/data/QPramp/example1/xout.csv 5e-5)

add_test(NAME TestLtiMpc_QPrampExample2 COMMAND TestLtiMpc ${CMAKE_CURRENT_SOURCE_DIR}/data/QPramp/example2 ${CMAKE_CURRENT_SOURCE_DIR}/data/QPramp/example2 10 100)
add_test(NAME TestLtiMpc_QPrampExample2_Verify_Uout COMMAND python ${CMAKE_CURRENT_SOURCE_DIR}/scripts/compare_csvs.py ${CMAKE_CURRENT_SOURCE_DIR}/data/QPramp/example2/uref.csv ${CMAKE_CURRENT_SOURCE_DIR}/data/QPramp/example2/uout.csv 3e-2)
add_test(NAME TestLtiMpc_QPrampExample2_Verify_Xout COMMAND python ${CMAKE_CURRENT_SOURCE_DIR}/scripts/compare_csvs.py ${CMAKE_CURRENT_SOURCE_DIR}/data/QPramp/example2/xref.csv ${CMAKE_CURRENT_SOURCE_DIR}/data/QPramp/example2/xout.csv 5e-3)

add_test(NAME TestLtiMpc_QPrampExample3 COMMAND TestLtiMpc ${CMAKE_CURRENT_SOURCE_DIR}/data/QPramp/example3 ${CMAKE_CURRENT_SOURCE_DIR}/data/QPramp/example3 10 100)
add_test(NAME TestLtiMpc_QPrampExample3_Verify_Uout COMMAND python ${CMAKE_CURRENT_SOURCE_DIR}/scripts/compare_csvs.py ${CMAKE_CURRENT_SOURCE_DIR}/data/QPramp/example3/uref.csv ${CMAKE_CURRENT_SOURCE_DIR}/data/QPramp/example3/uout.csv 6e-6)
add_test(NAME TestLtiMpc_QPrampExample3_Verify_Xout COMMAND python ${CMAKE_CURRENT_SOURCE_DIR}/scripts/compare_csvs.py ${CMAKE_CURRENT_SOURCE_DIR}/data/QPramp/example3/xref.csv ${CMAKE_CURRENT_SOURCE_DIR}/data/QPramp/example3/xout.csv 5e-5)

add_test(NAME TestLtiMpc_OsqpModel4 COMMAND TestLtiMpc ${CMAKE_CURRENT_SOURCE_DIR}/data/OSQP/model4 ${CMAKE_CURRENT_SOURCE_DIR}/data/OSQP/model4 10 100)
add_test(NAME TestLtiMpc_OsqpModel4_Verify_Uout COMMAND python ${CMAKE_CURRENT_SOURCE_DIR}/scripts/compare_csvs.py ${CMAKE_CURRENT_SOURCE_DIR}/data/OSQP/model4/uref.csv ${CMAKE_CURRENT_SOURCE_DIR}/data/OSQP/model4/uout.csv 6e-4)
add_test(NAME TestLtiMpc_OsqpModel4_Verify_Xout COMMAND python ${CMAKE_CURRENT_SOURCE_DIR}/scripts/compare_csvs.py ${CMAKE_CURRENT_SOURCE_DIR}/data/OSQP/model4/xref.csv ${CMAKE_CURRENT_SOURCE_DIR}/data/OSQP/model4/xout.csv 2e-3)

add_test(NAME TestLtiMpc_OsqpModel5 COMMAND TestLtiMpc ${CMAKE_CURRENT_SOURCE_DIR}/data/OSQP/model5 ${CMAKE_CURRENT_SOURCE_DIR}/data/OSQP/model5 10 100)
add_test(NAME TestLtiMpc_OsqpModel5_Verify_Uout COMMAND python ${CMAKE_CURRENT_SOURCE_DIR}/scripts/compare_csvs.py ${CMAKE_CURRENT_SOURCE_DIR}/data/OSQP/model5/uref.csv ${CMAKE_CURRENT_SOURCE_DIR}/data/OSQP/model5/uout.csv 8e-4)
add_test(NAME TestLtiMpc_OsqpModel5_Verify_Xout COMMAND python ${CMAKE_CURRENT_SOURCE_DIR}/scripts/compare_csvs.py ${CMAKE_CURRENT_SOURCE_DIR}/data/OSQP/model5/xref.csv ${CMAKE_CURRENT_SOURCE_DIR}/data/OSQP/model5/xout.csv 3e-3)

add_test(NAME TestLtiMpc_OsqpModel6 COMMAND TestLtiMpc ${CMAKE_CURRENT_SOURCE_DIR}/data/OSQP/model6 ${CMAKE_CURRENT_SOURCE_DIR}/data/OSQP/model6 10 100)
add_test(NAME TestLtiMpc_OsqpModel6_Verify_Uout COMMAND python ${CMAKE_CURRENT_SOURCE_DIR}/scripts/compare_csvs.py ${CMAKE_CURRENT_SOURCE_DIR}/data/OSQP/model6/uref.csv ${CMAKE_CURRENT_SOURCE_DIR}/data/OSQP/model6/uout.csv 2e-2)
add_test(NAME TestLtiMpc_OsqpModel6_Verify_Xout COMMAND python ${CMAKE_CURRENT_SOURCE_DIR}/scripts/compare_csvs.py ${CMAKE_CURRENT_SOURCE_DIR}/data/OSQP/model6/xref.csv ${CMAKE_CURRENT_SOURCE_DIR}/data/OSQP/model6/xout.csv 1e-1)


add_test(NAME TestLtvMpc_QPrampExample1 COMMAND TestLtvMpc ${CMAKE_CURRENT_SOURCE_DIR}/data/QPramp/example1 ${CMAKE_CURRENT_SOURCE_DIR}/data/QPramp/example1 10 100)
add_test(NAME TestLtvMpc_QPrampExample1_Verify_Uout COMMAND python ${CMAKE_CURRENT_SOURCE_DIR}/scripts/compare_csvs.py ${CMAKE_CURRENT_SOURCE_DIR}/data/QPramp/example1/uref.csv ${CMAKE_CURRENT_SOURCE_DIR}/data/QPramp/example1/uout.csv 4e-5)
add_test(NAME TestLtvMpc_QPrampExample1_Verify_Xout COMMAND python ${CMAKE_CURRENT_SOURCE_DIR}/scripts/compare_csvs.py ${CMAKE_CURRENT_SOURCE_DIR}/data/QPramp/example1/xref.csv ${CMAKE_CURRENT_SOURCE_DIR}/data/QPramp/example1/xout.csv 5e-5)

add_test(NAME TestLtvMpc_QPrampExample2 COMMAND TestLtvMpc ${CMAKE_CURRENT_SOURCE_DIR}/data/QPramp/example2 ${CMAKE_CURRENT_SOURCE_DIR}/data/QPramp/example2 10 100)
add_test(NAME TestLtvMpc_QPrampExample2_Verify_Uout COMMAND python ${CMAKE_CURRENT_SOURCE_DIR}/scripts/compare_csvs.py ${CMAKE_CURRENT_SOURCE_DIR}/data/QPramp/example2/uref.csv ${CMAKE_CURRENT_SOURCE_DIR}/data/QPramp/example2/uout.csv 3e-2)
add_test(NAME TestLtvMpc_QPrampExample2_Verify_Xout COMMAND python ${CMAKE_CURRENT_SOURCE_DIR}/scripts/compare_csvs.py ${CMAKE_CURRENT_SOURCE_DIR}/data/QPramp/example2/xref.csv ${CMAKE_CURRENT_SOURCE_DIR}/data/QPramp/example2/xout.csv 5e-3)

add_test(NAME TestLtvMpc_QPrampExample3 COMMAND TestLtvMpc ${CMAKE_CURRENT_SOURCE_DIR}/data/QPramp/example3 ${CMAKE_CURRENT_SOURCE_DIR}/data/QPramp/example3 10 100)
add_test(NAME TestLtvMpc_QPrampExample3_Verify_Uout COMMAND python ${CMAKE_CURRENT_SOURCE_DIR}/scripts/compare_csvs.py ${CMAKE_CURRENT_SOURCE_DIR}/data/QPramp/example3/uref.csv ${CMAKE_CURRENT_SOURCE_DIR}/data/QPramp/example3/uout.csv 6e-6)
add_test(NAME TestLtvMpc_QPrampExample3_Verify_Xout COMMAND python ${CMAKE_CURRENT_SOURCE_DIR}/scripts/compare_csvs.py ${CMAKE_CURRENT_SOURCE_DIR}/data/QPramp/example3/xref.csv ${CMAKE_CURRENT_SOURCE_DIR}/data/QPramp/example3/xout.csv 5e-5)

add_test(NAME TestLtvMpc_OsqpModel4 COMMAND TestLtvMpc ${CMAKE_CURRENT_SOURCE_DIR}/data/OSQP/model4 ${CMAKE_CURRENT_SOURCE_DIR}/data/OSQP/model4 10 100)
add_test(NAME TestLtvMpc_OsqpModel4_Verify_Uout COMMAND python ${CMAKE_CURRENT_SOURCE_DIR}/scripts/compare_csvs.py ${CMAKE_CURRENT_SOURCE_DIR}/data/OSQP/model4/uref.csv ${CMAKE_CURRENT_SOURCE_DIR}/data/OSQP/model4/uout.csv 7e-4)
add_test(NAME TestLtvMpc_OsqpModel4_Verify_Xout COMMAND python ${CMAKE_CURRENT_SOURCE_DIR}/scripts/compare_csvs.py ${CMAKE_CURRENT_SOURCE_DIR}/data/OSQP/model4/xref.csv ${CMAKE_CURRENT_SOURCE_DIR}/data/OSQP/model4/xout.csv 2e-3)

add_test(NAME TestLtvMpc_OsqpModel5 COMMAND TestLtvMpc ${CMAKE_CURRENT_SOURCE_DIR}/data/OSQP/model5 ${CMAKE_CURRENT_SOURCE_DIR}/data/OSQP/model5 10 100)
add_test(NAME TestLtvMpc_OsqpModel5_Verify_Uout COMMAND python ${CMAKE_CURRENT_SOURCE_DIR}/scripts/compare_csvs.py ${CMAKE_CURRENT_SOURCE_DIR}/data/OSQP/model5/uref.csv ${CMAKE_CURRENT_SOURCE_DIR}/data/OSQP/model5/uout.csv 7e-4)
add_test(NAME TestLtvMpc_OsqpModel5_Verify_Xout COMMAND python ${CMAKE_CURRENT_SOURCE_DIR}/scripts/compare_csvs.py ${CMAKE_CURRENT_SOURCE_DIR}/data/OSQP/model5/xref.csv ${CMAKE_CURRENT_SOURCE_DIR}/data/OSQP/model5/xout.csv 2e-3)

add_test(NAME TestLtvMpc_OsqpModel6 COMMAND TestLtvMpc ${CMAKE_CURRENT_SOURCE_DIR}/data/OSQP/model6 ${CMAKE_CURRENT_SOURCE_DIR}/data/OSQP/model6 10 100)
add_test(NAME TestLtvMpc_OsqpModel6_Verify_Uout COMMAND python ${CMAKE_CURRENT_SOURCE_DIR}/scripts/compare_csvs.py ${CMAKE_CURRENT_SOURCE_DIR}/data/OSQP/model6/uref.csv ${CMAKE_CURRENT_SOURCE_DIR}/data/OSQP/model6/uout.csv 3e-2)
add_test(NAME TestLtvMpc_OsqpModel6_Verify_Xout COMMAND python ${CMAKE_CURRENT_SOURCE_DIR}/scripts/compare_csvs.py ${CMAKE_CURRENT_SOURCE_DIR}/data/OSQP/model6/xref.csv ${CMAKE_CURRENT_SOURCE_DIR}/data/OSQP/model6/xout.csv 2e-1)


add_test(NAME TestMmcMpc_Horizon1 COMMAND TestMmcMpc ${CMAKE_CURRENT_SOURCE_DIR}/data/TestMmcMpc 1 10000)
add_test(NAME TestMmcMpc_Horizon2 COMMAND TestMmcMpc ${CMAKE_CURRENT_SOURCE_DIR}/data/TestMmcMpc 2 10000)
add_test(NAME TestMmcMpc_Horizon5 COMMAND TestMmcMpc ${CMAKE_CURRENT_SOURCE_DIR}/data/TestMmcMpc 5 10000)
add_test(NAME TestMmcMpc_Horizon10 COMMAND TestMmcMpc ${CMAKE_CURRENT_SOURCE_DIR}/data/TestMmcMpc 10 10000)

add_subdirectory(mex)
