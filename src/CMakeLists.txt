add_library(Ramp OBJECT Ramp.c LinAlg.c IndexedVectors.c IterableSet.c)
target_include_directories(Ramp PUBLIC ${CMAKE_SOURCE_DIR}/include)
target_link_libraries(Ramp PUBLIC m)

add_library(Simulate OBJECT Simulate.c Csv.c)
target_include_directories(Simulate PUBLIC ${CMAKE_SOURCE_DIR}/include)

add_library(Timer OBJECT Timer.c)
target_include_directories(Timer PUBLIC ${CMAKE_SOURCE_DIR}/include)


add_library(SdqpLmpcConstant STATIC SdqpLmpcConstant.c)
target_link_libraries(SdqpLmpcConstant PUBLIC Ramp)

add_library(SdqpLmpcVarying STATIC SdqpLmpcVarying.c)
target_link_libraries(SdqpLmpcVarying PUBLIC Ramp)

add_library(SdqpLmpcMmc STATIC SdqpLmpcMmc.c)
target_link_libraries(SdqpLmpcMmc PUBLIC Ramp)


add_library(MmcModel OBJECT MmcModel.c LinAlg.c)
target_include_directories(MmcModel PUBLIC ${CMAKE_SOURCE_DIR}/include)

add_library(SimulateMmcUtils OBJECT Extrapolate.c MmcTrajectory.c)
target_include_directories(SimulateMmcUtils PUBLIC ${CMAKE_SOURCE_DIR}/include)
target_link_libraries(SimulateMmcUtils PUBLIC m)
