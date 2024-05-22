add_executable(TestLtiMpc TestLtiMpc.c)
target_link_libraries(TestLtiMpc PRIVATE LtiMpc Simulate Csv Timer)

add_executable(TestLtvMpc TestLtvMpc.c)
target_link_libraries(TestLtvMpc PRIVATE LtvMpc Simulate Csv Timer)

add_executable(TestMmcMpc TestMmcMpc.c)
target_link_libraries(TestMmcMpc PRIVATE MmcMpc Simulate Csv MmcModel MmcMpcUtils Timer m)
