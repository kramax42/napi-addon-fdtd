// tutorial - how to make node addons
// https://medium.com/jspoint/a-simple-guide-to-load-c-c-code-into-node-js-javascript-applications-3fcccf54fd32
// https://github.com/course-one/native-addon-starter-kit/blob/master/src/index.cpp
// https://github.com/nodejs/node-addon-examples

// NAPI DESTINATION
// C:/Users/BOSS/AppData/Local/node-gyp/Cache/15.11.0/include/node
// /usr/include/node
// ${workspaceFolder}/**
// https://stackoverflow.com/questions/47616834/visual-studio-code-cannot-open-source-file-iostream

#include <math.h>
#include <napi.h>
#include <iostream>
#include <string>
#include <vector>

#include "./transform-to-js/get-data-1d/get-data-1d.h"
#include "./transform-to-js/get-data-2d/get-data-2d.h"


// Callback method when module is registered with Node.js.
Napi::Object Init(Napi::Env env, Napi::Object exports)
{
  exports.Set(Napi::String::New(env, "getData1D"),
              Napi::Function::New(env, GetData1D));

  exports.Set(Napi::String::New(env, "getData2D"),
              Napi::Function::New(env, GetData2D));

  // Return `exports` object (always).
  return exports;
}

// Register `FDTD` module which calls `Init` method.
NODE_API_MODULE(FDTD, Init)