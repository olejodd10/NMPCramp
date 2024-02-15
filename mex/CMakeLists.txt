find_package(Matlab)

matlab_add_mex(NAME MexMmcModel SRC MexMmcModel.c MexUtils.c LINK_TO m MmcModel R2018a)
