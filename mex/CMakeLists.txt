find_package(Matlab 9.7 EXACT)

matlab_add_mex(NAME MexMmcModel SRC MexMmcModel.c MexUtils.c LINK_TO m MmcModel R2018a)
matlab_add_mex(NAME MexMmcMpc SRC MexMmcMpc.c MexUtils.c LINK_TO m MmcMpc R2018a)
