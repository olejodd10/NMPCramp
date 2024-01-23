add_library(SdqpLmpcConstant STATIC Ramp.c SdqpLmpcConstant.c LinAlg.c IndexedVectors.c IterableSet.c)
target_include_directories(SdqpLmpcConstant PUBLIC ${CMAKE_SOURCE_DIR}/include)

add_library(LmpcConstant STATIC LmpcConstant.c Simulate.c Csv.c)
target_link_libraries(LmpcConstant PUBLIC SdqpLmpcConstant)
