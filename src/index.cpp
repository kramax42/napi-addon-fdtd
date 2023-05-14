#include <napi.h>

#include <string>

#include "./transform-to-js/fdtd-1d/fdtd-1d.h"
#include "./transform-to-js/fdtd-2d/fdtd-2d.h"
#include "./transform-to-js/fdtd-2d-tf-sf/fdtd-2d-tf-sf.h"

// Callback method when module is registered with Node.js.
Napi::Object Init(Napi::Env env, Napi::Object exports) {
    Napi::Object new_exports = exports;
    new_exports = Fdtd1D::Init(env, new_exports);
    new_exports = Fdtd2D::Init(env, new_exports);
    new_exports = Fdtd2DTFSF::Init(env, new_exports);
    
    return Fdtd1D::Init(env, new_exports);
}

// Register `FDTD` module which calls `Init` method.
NODE_API_MODULE(FDTD, Init)