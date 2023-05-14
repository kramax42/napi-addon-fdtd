// #ifndef FDTD_2D_H
// #define FDTD_2D_H

// #include <napi.h>

// #include "../../fdtd/2d-pml/fdtd-pml-2d.h"

// class Fdtd2D : public Napi::ObjectWrap<Fdtd2D> {
//    public:
//     static Napi::Object Init(Napi::Env env, Napi::Object exports);
//     Fdtd2D(const Napi::CallbackInfo& info);

//    private:
//     Napi::Value GetNextTimeLayer(const Napi::CallbackInfo& info);
//     FdtdPml2D fdtd;
//     // Data return type('Ez' = 0 | 'Hy' = 1 |'Hx' = 2 |'Energy' = 3)
//     int data_return_type = 0;
// };

// #endif