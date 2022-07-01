#ifndef GET_DATA_2D
#define GET_DATA_2D


#include <math.h>
#include <napi.h>
#include <iostream>
#include <string>
#include <vector>

#include "../../fdtd/2d-pml/fdtd-pml-2d.h"

// Fdtd method in 2D case.
Napi::Value GetData2D(const Napi::CallbackInfo &info);

#endif