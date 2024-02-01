add_library(SdqpLmpcConstant STATIC Ramp.c SdqpLmpcConstant.c LinAlg.c IndexedVectors.c IterableSet.c)
target_include_directories(SdqpLmpcConstant PUBLIC ${CMAKE_SOURCE_DIR}/include)

add_library(LmpcConstant STATIC LmpcConstant.c Simulate.c Csv.c)
target_link_libraries(LmpcConstant PUBLIC SdqpLmpcConstant)

add_library(SdqpLmpcVarying STATIC Ramp.c SdqpLmpcVarying.c LinAlg.c IndexedVectors.c IterableSet.c)
target_include_directories(SdqpLmpcVarying PUBLIC ${CMAKE_SOURCE_DIR}/include)

add_library(LmpcVarying OBJECT Simulate.c Csv.c)
target_link_libraries(LmpcVarying PUBLIC SdqpLmpcVarying)
