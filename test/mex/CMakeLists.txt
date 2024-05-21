if (CMAKE_HOST_WIN32)
    message(STATUS "TestMexMmcMpc configured for Windows")
    set(TEST_MEX_MATLAB_FLAGS -nodesktop -nosplash -wait)
else()
    set(TEST_MEX_MATLAB_FLAGS -nodisplay)
endif()

add_test(NAME TestMexMmcMpc_Horizons_1_2_5_10 COMMAND matlab ${TEST_MEX_MATLAB_FLAGS} -r 
    "addpath('${CMAKE_BINARY_DIR}/mex'), \
     addpath('${CMAKE_CURRENT_SOURCE_DIR}'), \
     TestMexMmcMpc('${CMAKE_SOURCE_DIR}/test/data/TestMexMmcMpc', 1, 10000), \
     clear mex, \
     TestMexMmcMpc('${CMAKE_SOURCE_DIR}/test/data/TestMexMmcMpc', 2, 10000), \
     clear mex, \
     TestMexMmcMpc('${CMAKE_SOURCE_DIR}/test/data/TestMexMmcMpc', 5, 10000), \
     clear mex, \
     TestMexMmcMpc('${CMAKE_SOURCE_DIR}/test/data/TestMexMmcMpc', 10, 10000), \
     exit")

set(TEST_MEX_SIMULATE_MMC_COMPARISON_LIMIT 8e-3)

add_test(NAME TestMexMmcMpc_Compare_xoutN1_With_TestMmcMpc COMMAND python ${CMAKE_SOURCE_DIR}/test/scripts/compare_csvs.py 
    ${CMAKE_SOURCE_DIR}/test/data/TestMexMmcMpc/xoutN1.csv ${CMAKE_SOURCE_DIR}/test/data/TestMmcMpc/xoutN1.csv ${TEST_MEX_SIMULATE_MMC_COMPARISON_LIMIT})
add_test(NAME TestMexMmcMpc_Compare_uoutN1_With_TestMmcMpc COMMAND python ${CMAKE_SOURCE_DIR}/test/scripts/compare_csvs.py 
    ${CMAKE_SOURCE_DIR}/test/data/TestMexMmcMpc/uoutN1.csv ${CMAKE_SOURCE_DIR}/test/data/TestMmcMpc/uoutN1.csv ${TEST_MEX_SIMULATE_MMC_COMPARISON_LIMIT})

add_test(NAME TestMexMmcMpc_Compare_xoutN2_With_TestMmcMpc COMMAND python ${CMAKE_SOURCE_DIR}/test/scripts/compare_csvs.py 
    ${CMAKE_SOURCE_DIR}/test/data/TestMexMmcMpc/xoutN2.csv ${CMAKE_SOURCE_DIR}/test/data/TestMmcMpc/xoutN2.csv ${TEST_MEX_SIMULATE_MMC_COMPARISON_LIMIT})
add_test(NAME TestMexMmcMpc_Compare_uoutN2_With_TestMmcMpc COMMAND python ${CMAKE_SOURCE_DIR}/test/scripts/compare_csvs.py 
    ${CMAKE_SOURCE_DIR}/test/data/TestMexMmcMpc/uoutN2.csv ${CMAKE_SOURCE_DIR}/test/data/TestMmcMpc/uoutN2.csv ${TEST_MEX_SIMULATE_MMC_COMPARISON_LIMIT})

add_test(NAME TestMexMmcMpc_Compare_xoutN5_With_TestMmcMpc COMMAND python ${CMAKE_SOURCE_DIR}/test/scripts/compare_csvs.py 
    ${CMAKE_SOURCE_DIR}/test/data/TestMexMmcMpc/xoutN5.csv ${CMAKE_SOURCE_DIR}/test/data/TestMmcMpc/xoutN5.csv ${TEST_MEX_SIMULATE_MMC_COMPARISON_LIMIT})
add_test(NAME TestMexMmcMpc_Compare_uoutN5_With_TestMmcMpc COMMAND python ${CMAKE_SOURCE_DIR}/test/scripts/compare_csvs.py 
    ${CMAKE_SOURCE_DIR}/test/data/TestMexMmcMpc/uoutN5.csv ${CMAKE_SOURCE_DIR}/test/data/TestMmcMpc/uoutN5.csv ${TEST_MEX_SIMULATE_MMC_COMPARISON_LIMIT})

add_test(NAME TestMexMmcMpc_Compare_xoutN10_With_TestMmcMpc COMMAND python ${CMAKE_SOURCE_DIR}/test/scripts/compare_csvs.py 
    ${CMAKE_SOURCE_DIR}/test/data/TestMexMmcMpc/xoutN10.csv ${CMAKE_SOURCE_DIR}/test/data/TestMmcMpc/xoutN10.csv ${TEST_MEX_SIMULATE_MMC_COMPARISON_LIMIT})
add_test(NAME TestMexMmcMpc_Compare_uoutN10_With_TestMmcMpc COMMAND python ${CMAKE_SOURCE_DIR}/test/scripts/compare_csvs.py 
    ${CMAKE_SOURCE_DIR}/test/data/TestMexMmcMpc/uoutN10.csv ${CMAKE_SOURCE_DIR}/test/data/TestMmcMpc/uoutN10.csv ${TEST_MEX_SIMULATE_MMC_COMPARISON_LIMIT})


add_test(NAME TestMexLtiMpc_AllModels COMMAND matlab ${TEST_MEX_MATLAB_FLAGS} -r 
    "addpath('${CMAKE_BINARY_DIR}/mex'), \
     addpath('${CMAKE_CURRENT_SOURCE_DIR}'), \
     TestMexLtiMpc('${CMAKE_SOURCE_DIR}/test/data/QPramp/example1', '${CMAKE_SOURCE_DIR}/test/data/QPramp/example1', 10, 100), \
     clear mex, \
     TestMexLtiMpc('${CMAKE_SOURCE_DIR}/test/data/QPramp/example2', '${CMAKE_SOURCE_DIR}/test/data/QPramp/example2', 10, 100), \
     clear mex, \
     TestMexLtiMpc('${CMAKE_SOURCE_DIR}/test/data/QPramp/example3', '${CMAKE_SOURCE_DIR}/test/data/QPramp/example3', 10, 100), \
     clear mex, \
     TestMexLtiMpc('${CMAKE_SOURCE_DIR}/test/data/OSQP/model4', '${CMAKE_SOURCE_DIR}/test/data/OSQP/model4', 10, 100), \
     clear mex, \
     TestMexLtiMpc('${CMAKE_SOURCE_DIR}/test/data/OSQP/model5', '${CMAKE_SOURCE_DIR}/test/data/OSQP/model5', 10, 100), \
     clear mex, \
     TestMexLtiMpc('${CMAKE_SOURCE_DIR}/test/data/OSQP/model6', '${CMAKE_SOURCE_DIR}/test/data/OSQP/model6', 10, 100), \
     clear mex, \
     exit")

add_test(NAME TestMexLtiMpc_QPrampExample1_Verify_Uout COMMAND python ${CMAKE_SOURCE_DIR}/test/scripts/compare_csvs.py ${CMAKE_SOURCE_DIR}/test/data/QPramp/example1/uref.csv ${CMAKE_SOURCE_DIR}/test/data/QPramp/example1/uout.csv 4e-5)
add_test(NAME TestMexLtiMpc_QPrampExample1_Verify_Xout COMMAND python ${CMAKE_SOURCE_DIR}/test/scripts/compare_csvs.py ${CMAKE_SOURCE_DIR}/test/data/QPramp/example1/xref.csv ${CMAKE_SOURCE_DIR}/test/data/QPramp/example1/xout.csv 5e-5)

add_test(NAME TestMexLtiMpc_QPrampExample2_Verify_Uout COMMAND python ${CMAKE_SOURCE_DIR}/test/scripts/compare_csvs.py ${CMAKE_SOURCE_DIR}/test/data/QPramp/example2/uref.csv ${CMAKE_SOURCE_DIR}/test/data/QPramp/example2/uout.csv 3e-2)
add_test(NAME TestMexLtiMpc_QPrampExample2_Verify_Xout COMMAND python ${CMAKE_SOURCE_DIR}/test/scripts/compare_csvs.py ${CMAKE_SOURCE_DIR}/test/data/QPramp/example2/xref.csv ${CMAKE_SOURCE_DIR}/test/data/QPramp/example2/xout.csv 5e-3)

add_test(NAME TestMexLtiMpc_QPrampExample3_Verify_Uout COMMAND python ${CMAKE_SOURCE_DIR}/test/scripts/compare_csvs.py ${CMAKE_SOURCE_DIR}/test/data/QPramp/example3/uref.csv ${CMAKE_SOURCE_DIR}/test/data/QPramp/example3/uout.csv 6e-6)
add_test(NAME TestMexLtiMpc_QPrampExample3_Verify_Xout COMMAND python ${CMAKE_SOURCE_DIR}/test/scripts/compare_csvs.py ${CMAKE_SOURCE_DIR}/test/data/QPramp/example3/xref.csv ${CMAKE_SOURCE_DIR}/test/data/QPramp/example3/xout.csv 5e-5)

add_test(NAME TestMexLtiMpc_OsqpModel4_Verify_Uout COMMAND python ${CMAKE_SOURCE_DIR}/test/scripts/compare_csvs.py ${CMAKE_SOURCE_DIR}/test/data/OSQP/model4/uref.csv ${CMAKE_SOURCE_DIR}/test/data/OSQP/model4/uout.csv 6e-4)
add_test(NAME TestMexLtiMpc_OsqpModel4_Verify_Xout COMMAND python ${CMAKE_SOURCE_DIR}/test/scripts/compare_csvs.py ${CMAKE_SOURCE_DIR}/test/data/OSQP/model4/xref.csv ${CMAKE_SOURCE_DIR}/test/data/OSQP/model4/xout.csv 2e-3)

add_test(NAME TestMexLtiMpc_OsqpModel5_Verify_Uout COMMAND python ${CMAKE_SOURCE_DIR}/test/scripts/compare_csvs.py ${CMAKE_SOURCE_DIR}/test/data/OSQP/model5/uref.csv ${CMAKE_SOURCE_DIR}/test/data/OSQP/model5/uout.csv 8e-4)
add_test(NAME TestMexLtiMpc_OsqpModel5_Verify_Xout COMMAND python ${CMAKE_SOURCE_DIR}/test/scripts/compare_csvs.py ${CMAKE_SOURCE_DIR}/test/data/OSQP/model5/xref.csv ${CMAKE_SOURCE_DIR}/test/data/OSQP/model5/xout.csv 3e-3)

add_test(NAME TestMexLtiMpc_OsqpModel6_Verify_Uout COMMAND python ${CMAKE_SOURCE_DIR}/test/scripts/compare_csvs.py ${CMAKE_SOURCE_DIR}/test/data/OSQP/model6/uref.csv ${CMAKE_SOURCE_DIR}/test/data/OSQP/model6/uout.csv 2e-2)
add_test(NAME TestMexLtiMpc_OsqpModel6_Verify_Xout COMMAND python ${CMAKE_SOURCE_DIR}/test/scripts/compare_csvs.py ${CMAKE_SOURCE_DIR}/test/data/OSQP/model6/xref.csv ${CMAKE_SOURCE_DIR}/test/data/OSQP/model6/xout.csv 1e-1)
