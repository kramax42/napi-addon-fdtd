#ifndef FDTD_2D_TF_SF_H
#define FDTD_2D_TF_SF_H

#include <napi.h>

#include "../../fdtd/2d-upml-tf-sf/fdtd-2d-upml-tf-sf.h"

class Fdtd2DTFSF : public Napi::ObjectWrap<Fdtd2DTFSF> {
   public:
    static Napi::Object Init(Napi::Env env, Napi::Object exports);
    Fdtd2DTFSF(const Napi::CallbackInfo& info);

   private:
    Napi::Value GetNextTimeLayer(const Napi::CallbackInfo& info);
    TFSF fdtd;
    // Data return type('Ez' = 0 | 'Hy' = 1 |'Hx' = 2 |'Energy' = 3)
    // int data_return_type = 0;
};

#endif