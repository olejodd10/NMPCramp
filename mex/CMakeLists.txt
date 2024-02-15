find_package(Matlab)

matlab_add_mex(NAME MexMmcModel SRC MexMmcModel.c MexUtils.c LINK_TO m MmcModel R2018a)
matlab_add_mex(NAME MexSdqpLmpcMmc SRC MexSdqpLmpcMmc.c MexUtils.c LINK_TO m SdqpLmpcMmc R2018a)
