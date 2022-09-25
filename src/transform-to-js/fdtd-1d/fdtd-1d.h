#ifndef FDTD_1D_H
#define FDTD_1D_H

#include <napi.h>

#include "../../fdtd/1d-pml/fdtd-pml-1d.h"

class Fdtd1D : public Napi::ObjectWrap<Fdtd1D> {
   public:
    static Napi::Object Init(Napi::Env env, Napi::Object exports);
    Fdtd1D(const Napi::CallbackInfo& info);

   private:
    Napi::Value GetNextTimeLayer(const Napi::CallbackInfo& info);
    FdtdPml1D fdtd;
};

#endif