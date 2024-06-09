# RaMPC

This project develops a ramp function-based solver for QP problems having a positive semi-definite Hessian matrix. Complete solvers $LtiMpc$, $LtvMpc$ and $MmcMpc$ are implemented for a common linear time-variant model predictive control (MPC) formulation, a common linear time-invariant MPC formulation and an MPC formulation for modular multilevel converters (MMC), respectively.

The project is a continuation of my specialization project, which implemented a ramp function-based solver for a condensed formulation of QP problems with a positive definite Hessian. That solver is available [here](https://github.com/olejodd10/QPramp).

Developers interested in implementation details are encouraged to read my master's thesis.

## Building on Linux/WSL

```
cd /path/to/RaMPC
mkdir build
cd build
cmake ..
cmake --build .
```

## Building in PowerShell

```
cd /path/to/RaMPC
mkdir build
cd build
cmake .. -G "MinGW Makefiles"
cmake --build .
```

## Testing

```
cd /path/to/RaMPC/build
mkdir ../test/data/TestMmcMpc 
mkdir ../test/data/TestMexMmcMpc
ctest
```

## Using the C Example Program

```
/path/to/RaMPC/build/test/src/TestLtiMpc "/path/to/input/folder" \
"/path/to/output/folder" "10" "100"
```

## Using the MATLAB Example Program

```
addpath(’/path/to/RaMPC/build/mex’); % The MEX file MexLtiMpc is here
addpath(’/path/to/RaMPC/test/mex’); % The MATLAB function TestMexLtiMpc is here
TestMexLtiMpc(’/path/to/input/folder’, ’/path/to/output/folder’, 10, 100);
```
