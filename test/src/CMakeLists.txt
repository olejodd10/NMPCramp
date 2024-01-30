add_executable(TestLmpcConstant TestLmpcConstant.c)
target_link_libraries(TestLmpcConstant PRIVATE LmpcConstant)

add_executable(TestLmpcVarying TestLmpcVarying.c ${CMAKE_SOURCE_DIR}/src/Simulate.c ${CMAKE_SOURCE_DIR}/src/Csv.c)
target_link_libraries(TestLmpcVarying PRIVATE SdqpLmpcVarying)
