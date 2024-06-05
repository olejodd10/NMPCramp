project_path = '..';
defs = [];

def = legacy_code('initialize');
def.SFunctionName = 'MmcSFunctionA';
% specification
def.StartFcnSpec = 'void mmc_s_function_start(int16 p1, double p2, double p3, double p4, double p5, double p6, double p7, double p8, double p9, double p10, double p11, double p12[4], double p13[4], double p14, double p15[2], double p16[2])';
def.OutputFcnSpec = 'int16 y3 = mmc_s_function(int16 p1, double u1, double u2, double u3, double u4, double u5, double u6, double u7[4], double y1[p1][4], double y2[p1][2])';
def.TerminateFcnSpec = 'void mmc_s_function_terminate(void)';
def.IncPaths      = {strcat(project_path, '/include'), strcat(project_path, '/mex')};
def.SrcPaths      = {strcat(project_path, '/src'), strcat(project_path, '/mex')};
def.HeaderFiles   = {'MmcSFunction.h'};
def.SourceFiles   = {'MmcTrajectory.c', 'MmcSFunction.c', 'MmcModel.c' 'SdqpLmpcMmc.c', 'LinAlg.c', 'Ramp.c', 'IndexedVectors.c', 'IterableSet.c'};
def.Options.convertNDArrayToRowMajor = true;
defs = [defs; def];

def.SFunctionName = 'MmcSFunctionB';
defs = [defs; def];

def.SFunctionName = 'MmcSFunctionC';
defs = [defs; def];

legacy_code('generate_for_sim', defs);
legacy_code('slblock_generate', defs);
