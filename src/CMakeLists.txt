add_library(Ramp OBJECT Ramp.c LinAlg.c IndexedVectors.c IterableSet.c)
target_include_directories(Ramp PUBLIC ${CMAKE_SOURCE_DIR}/include)
target_link_libraries(Ramp PUBLIC m)

add_library(Simulate OBJECT Simulate.c Csv.c)
target_include_directories(Simulate PUBLIC ${CMAKE_SOURCE_DIR}/include)

add_library(Timer OBJECT Timer.c)
target_include_directories(Timer PUBLIC ${CMAKE_SOURCE_DIR}/include)


add_library(LtiMpc STATIC LtiMpc.c)
target_link_libraries(LtiMpc PUBLIC Ramp)

add_library(LtvMpc STATIC LtvMpc.c Csr.c)
target_link_libraries(LtvMpc PUBLIC Ramp)

add_library(MmcMpc STATIC MmcMpc.c)
target_link_libraries(MmcMpc PUBLIC Ramp)


add_library(MmcModel OBJECT MmcModel.c LinAlg.c)
target_include_directories(MmcModel PUBLIC ${CMAKE_SOURCE_DIR}/include)

add_library(MmcMpcUtils OBJECT Extrapolate.c MmcTrajectory.c)
target_include_directories(MmcMpcUtils PUBLIC ${CMAKE_SOURCE_DIR}/include)
target_link_libraries(MmcMpcUtils PUBLIC m)
