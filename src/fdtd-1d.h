#ifndef FDTD_1D_H
#define FDTD_1D_H

#include <napi.h>

#include "./fdtd/1d-pml/fdtd-pml-1d.h"

class Fdtd1D : public Napi::ObjectWrap<Fdtd1D> {
   public:
    static Napi::Object Init(Napi::Env env, Napi::Object exports);
    Fdtd1D(const Napi::CallbackInfo& info);

   private:
    // Napi::Value GetValue(const Napi::CallbackInfo& info);
    // Napi::Value PlusOne(const Napi::CallbackInfo& info);
    // Napi::Value Multiply(const Napi::CallbackInfo& info);
    Napi::Value GetNextTimeLayer(const Napi::CallbackInfo &info);
    // Napi::Value CreateNewItem(const Napi::CallbackInfo& info);

    // double value;
    FdtdPml1D fdtd;
};

#endif