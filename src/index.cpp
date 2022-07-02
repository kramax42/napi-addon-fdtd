#include <napi.h>
#include <string>

#include "./transform-to-js/get-data-1d/get-data-1d.h"
#include "./transform-to-js/get-data-2d/get-data-2d.h"

// Callback method when module is registered with Node.js.
Napi::Object Init(Napi::Env env, Napi::Object exports)
{
  exports.Set(Napi::String::New(env, "getData1D"),
              Napi::Function::New(env, GetData1D));

  exports.Set(Napi::String::New(env, "getData2D"),
              Napi::Function::New(env, GetData2D));

  return exports;
}

// Register `FDTD` module which calls `Init` method.
NODE_API_MODULE(FDTD, Init)