add_executable(TestLtiMpc TestLtiMpc.c)
target_link_libraries(TestLtiMpc PRIVATE LtiMpc Simulate Timer)

add_executable(TestLtvMpc TestLtvMpc.c)
target_link_libraries(TestLtvMpc PRIVATE LtvMpc Simulate Timer)

add_executable(TestMmcMpc TestMmcMpc.c)
target_link_libraries(TestMmcMpc PRIVATE MmcMpc Simulate MmcModel MmcMpcUtils Timer m)
