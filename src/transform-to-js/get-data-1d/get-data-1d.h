#ifndef GET_DATA_1D
#define GET_DATA_1D

#include <math.h>
#include <napi.h>

#include <iostream>
#include <string>
#include <vector>

#include "../../fdtd/1d-pml/fdtd-pml-1d.h"

// Fdtd method in 1D case.
Napi::Value GetData1D(const Napi::CallbackInfo &info);

#endif