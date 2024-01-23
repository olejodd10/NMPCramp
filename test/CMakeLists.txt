add_subdirectory(src)


add_test(NAME TestNMPCramp_Compiles COMMAND make)


add_test(NAME TestLmpcConstant_BadPath_ShouldFail COMMAND TestLmpcConstant /bad/path /bad/path2 10 100)
set_tests_properties(TestLmpcConstant_BadPath_ShouldFail PROPERTIES WILL_FAIL TRUE)

add_test(NAME TestLmpcConstant_MissingArguments_ShouldFail COMMAND TestLmpcConstant /bad/path /bad/path2 10)
set_tests_properties(TestLmpcConstant_MissingArguments_ShouldFail PROPERTIES WILL_FAIL TRUE)
