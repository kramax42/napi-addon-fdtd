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

#include "./FDTD/2D/FDTD_2D.h"
#include "./FDTD/3D/FDTD_3D.h"
#include "./FDTD/3D_DIFRACTION/FDTD_3D_DIFRACTION.h"
// #include "./FDTD/3D_DIFRACTION/get-data-difraction-3d.cpp"
#include "./FDTD/3D_INTERFERENCE/FDTD_3D_INTERFERENCE.h"
#include "globals.h"  // Nx, Ny

// https://stackoverflow.com/questions/12573816/what-is-an-undefined-reference-unresolved-external-symbol-error-and-how-do-i-fix/12574420#12574420
// https://www.it-swarm.com.ru/ru/c%2B%2B/chto-takoe-neopredelennaya-ssylka-nerazreshennaya-vneshnyaya-oshibka-simvola-i-kak-ee-ispravit/1069256308/
extern "C" {
#include "./FDTD/pure-c/FDTD-2D.h"
}

// Difraction FDTD-3D.
Napi::Value GetDataDifraction3D(const Napi::CallbackInfo &info) {

  Napi::Env env = info.Env();

  // 0 - lambda.
  // 1 - reload checker.
  // 2 - refractive index matrix(flatten).
  // 3 - refractive index matrix(flatten) size (for 2x2 is 2).
  // 4 - data return type('Ez' = 0 | 'Hy' = 1 |'Hx' = 2 |'Energy' = 3)
  const Napi::Array input_array_condition = info[0].As<Napi::Array>();
  
  // Reload params checker.
  bool reload_check = static_cast<bool>(info[1].As<Napi::Boolean>());

  // Refraction index matrix transformation JS -> C++.
  const Napi::Array refr_index_matrix_js = info[2].As<Napi::Array>();

  // Must be even.
  int refr_index_matrix_size = static_cast<int>(info[3].As<Napi::Number>());

  // Temporary matrix.
  std::vector<std::vector<double>> temp_matrix;

  // Data return type('Ez' = 0 | 'Hy' = 1 |'Hx' = 2 |'Energy' = 3)
  int data_return_type = static_cast<int>(info[4].As<Napi::Number>());

  // Params transformation JS -> C++.
  double lambda = (double)input_array_condition[(uint32_t)0].As<Napi::Number>();
  double beamsize = (double)input_array_condition[1].As<Napi::Number>();
  double n1 = (double)input_array_condition[2].As<Napi::Number>();
  double n2 = (double)input_array_condition[3].As<Napi::Number>();

  // Transform input flatten matrix into 2-dimensional matrix.
  for (int i = 0; i < refr_index_matrix_size; i++) {
    temp_matrix.push_back(std::vector<double>());
    for (int j = 0; j < refr_index_matrix_size; j++) {
      temp_matrix[i].push_back(
          (float)refr_index_matrix_js[i * refr_index_matrix_size + j]
              .As<Napi::Number>());
    }
  }

  // Matrix size  coefficient.
  size_t coeff = Nx / refr_index_matrix_size;

  // Initialization refractive index matrix.
  std::vector<std::vector<double>> refr_index_matrix;
  for (int i = 0; i < Nx; i++) {
    refr_index_matrix.push_back(std::vector<double>());
    for (int j = 0; j < Ny; j++) {
      refr_index_matrix[i].push_back(n1);
    }
  }

  // Filling refractive index matrix.
  for (int i = 0; i < refr_index_matrix_size; i++) {
    for (int j = 0; j < refr_index_matrix_size; j++) {
      for (int k = 0; k < coeff; k++) {
        for (int f = 0; f < coeff; f++) {
          refr_index_matrix[i * coeff + k][j * coeff + f] = temp_matrix[i][j];
        }
      }
    }
  }

  //temp

// std::cout << "start" << std::endl;
// for (int i = 0; i < Nx; i++) {
//     std::cout << std::endl;
//     for (int j = 0; j < Ny; j++) {
//       std::cout << refr_index_matrix[i][j] << " ";
//     }
//   }
// std::cout << "sstop" << std::endl;
  // end temp



  static FDTD_3D_DIFRACTION fdtd_3D =
      FDTD_3D_DIFRACTION(lambda, beamsize, n1, n2, refr_index_matrix);
  if ((fdtd_3D.getLambda() != lambda) || (fdtd_3D.getBeamsize() != beamsize) ||
      (fdtd_3D.getN1() != n1) || reload_check) {
    fdtd_3D.setLambda(lambda);
    fdtd_3D.setBeamsize(beamsize);
    fdtd_3D.setN1(n1);
    fdtd_3D.setN2(n2);
    fdtd_3D.setParams(refr_index_matrix);
  }

  std::vector<double> vect_X = {};
  std::vector<double> vect_Y = {};
  std::vector<double> vect_Ez = {};
  std::vector<double> vect_Hy = {};
  std::vector<double> vect_Hx = {};
  std::vector<double> vect_Energy = {};

  fdtd_3D.calcNextLayer(vect_X, vect_Y, vect_Ez, vect_Hy, vect_Hx, vect_Energy);

  // Matrix sizes.
  size_t Nx = fdtd_3D.getNx() / fdtd_3D.getStep();
  size_t Ny = fdtd_3D.getNy() / fdtd_3D.getStep();

  // Creating JS arrays to store C++ arrays.
  Napi::Array js_data_X = Napi::Array::New(env, Nx * Ny);
  Napi::Array js_data_Y = Napi::Array::New(env, Nx * Ny);
  Napi::Array js_data_Ez = Napi::Array::New(env, Nx * Ny);
  Napi::Array js_data_Hy = Napi::Array::New(env, Nx * Ny);
  Napi::Array js_data_Hx = Napi::Array::New(env, Nx * Ny);
  Napi::Array js_data_Energy = Napi::Array::New(env, Nx * Ny);

  // Filling JS arrays with C++ arrays data.
  for (size_t i = 0; i < Nx * Ny; i++) {
    js_data_X[i] = Napi::Number::New(env, vect_X[i]);
    js_data_Y[i] = Napi::Number::New(env, vect_Y[i]);
    js_data_Ez[i] = Napi::Number::New(env, vect_Ez[i]);
    js_data_Hy[i] = Napi::Number::New(env, vect_Hy[i]);
    js_data_Hx[i] = Napi::Number::New(env, vect_Hx[i]);
    js_data_Energy[i] = Napi::Number::New(env, vect_Energy[i]);
  }

  // Creating JS object to return.
  Napi::Object data = Napi::Array::New(env);
  data.Set("dataX", js_data_X);
  data.Set("dataY", js_data_Y);
  data.Set("row", Nx);
  data.Set("col", Ny);
  data.Set("currentTick", fdtd_3D.getCurrentTick());

  switch (data_return_type) {
    case 0:
      data.Set("dataEz", js_data_Ez);
      break;
    case 1:
      data.Set("dataHy", js_data_Hy);
      break;
    case 2:
      data.Set("dataHx", js_data_Hy);
      break;
    case 3:
      data.Set("dataEnergy", js_data_Energy);
      break;

    default:
      break;
  }

  return data;
}

// Interference FDTD_3D.
Napi::Value getFDTD_3D_INTERFERENCE(const Napi::CallbackInfo &info) {
  Napi::Env env = info.Env();

  const Napi::Array inputArrayCondition = info[0].As<Napi::Array>();
  // 0 - lambda.
  // 1 - beamsize.
  // 2 - n1.

  // Reload current params?
  bool reload = static_cast<bool>(info[1].As<Napi::Boolean>());

  int first = 0;                                                          //????
  double lambda = (double)inputArrayCondition[first].As<Napi::Number>();  //????
  double beamsize = (double)(inputArrayCondition[1].As<Napi::Number>());
  double n1 = (double)(inputArrayCondition[2].As<Napi::Number>());

  static FDTD_3D_INTERFERENCE fdtd_3D =
      FDTD_3D_INTERFERENCE(lambda, beamsize, n1);
  if ((fdtd_3D.getLambda() != lambda) || (fdtd_3D.getBeamsize() != beamsize) ||
      (fdtd_3D.getN1() != n1) || reload) {
    fdtd_3D.setLambda(lambda);
    fdtd_3D.setBeamsize(beamsize);
    fdtd_3D.setN1(n1);
    fdtd_3D.setParams();
  }

  vector<double> vect_X = {};
  vector<double> vect_Y = {};
  vector<double> vect_Ez = {};
  vector<double> vect_Hy = {};
  vector<double> vect_Hx = {};
  vector<double> vect_Energy = {};

  vect_X.clear();
  vect_Y.clear();
  vect_Ez.clear();
  vect_Hy.clear();
  vect_Hx.clear();
  vect_Energy.clear();

  fdtd_3D.calcNextLayer(vect_X, vect_Y, vect_Ez, vect_Hy, vect_Hx, vect_Energy);

  size_t Nx = fdtd_3D.getNx() / fdtd_3D.getStep();
  size_t Ny = fdtd_3D.getNy() / fdtd_3D.getStep();

  // Creating arrays.
  Napi::Array js_data_X = Napi::Array::New(env, Nx * Ny);
  Napi::Array js_data_Y = Napi::Array::New(env, Nx * Ny);
  Napi::Array js_data_Ez = Napi::Array::New(env, Nx * Ny);
  Napi::Array js_data_Hy = Napi::Array::New(env, Nx * Ny);
  Napi::Array js_data_Hx = Napi::Array::New(env, Nx * Ny);
  Napi::Array js_data_Energy = Napi::Array::New(env, Nx * Ny);

  // Temporary variables.
  Napi::Number elem;

  for (size_t i = 0; i < Nx * Ny; i++) {
    elem = Napi::Number::New(env, vect_X[i]);
    js_data_X[i] = elem;

    elem = Napi::Number::New(env, vect_Y[i]);
    js_data_Y[i] = elem;

    elem = Napi::Number::New(env, vect_Ez[i]);
    js_data_Ez[i] = elem;

    elem = Napi::Number::New(env, vect_Hy[i]);
    js_data_Hy[i] = elem;

    elem = Napi::Number::New(env, vect_Hx[i]);
    js_data_Hx[i] = elem;

    elem = Napi::Number::New(env, vect_Energy[i]);
    js_data_Energy[i] = elem;
  }

  Napi::Object data = Napi::Array::New(env);
  data.Set("dataX", js_data_X);
  data.Set("dataY", js_data_Y);
  data.Set("dataEz", js_data_Ez);
  data.Set("dataHy", js_data_Hy);
  data.Set("dataHx", js_data_Hx);
  data.Set("dataEnergy", js_data_Energy);
  data.Set("row", Nx);
  data.Set("col", Ny);
  data.Set("currentTick", fdtd_3D.getCurrentTick());

  return data;
}

// Plain FDTD_3D.
Napi::Value getFDTD_3D(const Napi::CallbackInfo &info) {
  Napi::Env env = info.Env();

  const Napi::Array inputArrayCondition = info[0].As<Napi::Array>();
  // 0 - lambda.
  // 1 - beamsize.
  // 2 - n1.

  // Reload current params?
  bool reload = static_cast<bool>(info[1].As<Napi::Boolean>());

  int first = 0;                                                          //????
  double lambda = (double)inputArrayCondition[first].As<Napi::Number>();  //????
  double beamsize = (double)(inputArrayCondition[1].As<Napi::Number>());
  double n1 = (double)(inputArrayCondition[2].As<Napi::Number>());

  static FDTD_3D fdtd_3D = FDTD_3D(lambda, beamsize, n1);
  if ((fdtd_3D.getLambda() != lambda) || (fdtd_3D.getBeamsize() != beamsize) ||
      (fdtd_3D.getN1() != n1) || reload) {
    fdtd_3D.setLambda(lambda);
    fdtd_3D.setBeamsize(beamsize);
    fdtd_3D.setN1(n1);
    fdtd_3D.setParams();
  }

  vector<double> vect_X = {};
  vector<double> vect_Y = {};
  vector<double> vect_Ez = {};
  vector<double> vect_Hy = {};
  vector<double> vect_Hx = {};
  vector<double> vect_Energy = {};

  vect_X.clear();
  vect_Y.clear();
  vect_Ez.clear();
  vect_Hy.clear();
  vect_Hx.clear();
  vect_Energy.clear();

  fdtd_3D.calcNextLayer(vect_X, vect_Y, vect_Ez, vect_Hy, vect_Hx, vect_Energy);

  size_t Nx = fdtd_3D.getNx() / fdtd_3D.getStep();
  size_t Ny = fdtd_3D.getNy() / fdtd_3D.getStep();

  // Creating arrays.
  Napi::Array js_data_X = Napi::Array::New(env, Nx * Ny);
  Napi::Array js_data_Y = Napi::Array::New(env, Nx * Ny);
  Napi::Array js_data_Ez = Napi::Array::New(env, Nx * Ny);
  Napi::Array js_data_Hy = Napi::Array::New(env, Nx * Ny);
  Napi::Array js_data_Hx = Napi::Array::New(env, Nx * Ny);
  Napi::Array js_data_Energy = Napi::Array::New(env, Nx * Ny);

  // Temporary variables.
  Napi::Number elem;

  for (size_t i = 0; i < Nx * Ny; i++) {
    elem = Napi::Number::New(env, vect_X[i]);
    js_data_X[i] = elem;

    elem = Napi::Number::New(env, vect_Y[i]);
    js_data_Y[i] = elem;

    elem = Napi::Number::New(env, vect_Ez[i]);
    js_data_Ez[i] = elem;

    elem = Napi::Number::New(env, vect_Hy[i]);
    js_data_Hy[i] = elem;

    elem = Napi::Number::New(env, vect_Hx[i]);
    js_data_Hx[i] = elem;

    elem = Napi::Number::New(env, vect_Energy[i]);
    js_data_Energy[i] = elem;
  }

  Napi::Object data = Napi::Array::New(env);
  data.Set("dataX", js_data_X);
  data.Set("dataY", js_data_Y);
  data.Set("dataEz", js_data_Ez);
  data.Set("dataHy", js_data_Hy);
  data.Set("dataHx", js_data_Hx);
  data.Set("dataEnergy", js_data_Energy);
  data.Set("row", Nx);
  data.Set("col", Ny);
  data.Set("currentTick", fdtd_3D.getCurrentTick());

  return data;
}

// FDTD method in 2D case.
Napi::Value GetData2D(const Napi::CallbackInfo &info) {
  Napi::Env env = info.Env();

  // Params - ([lambda, tau, refractive_index], reload)
  const Napi::Array input_array_condition = info[0].As<Napi::Array>();

  // Reload params checker.
  bool reload_check = static_cast<bool>(info[1].As<Napi::Boolean>());

  int nil = 0;  //!!!! MUST BE REFACTORED !!!!!!!
  float lambda = (float)input_array_condition[nil].As<Napi::Number>();
  float tau = (float)(input_array_condition[1].As<Napi::Number>());
  float refractive_index = (float)(input_array_condition[2].As<Napi::Number>());

  // Containers to storage coordinates.
  vector<double> vect_X = {};
  vector<double> vect_Y = {};

  // Using static to save save data for different function call.
  static FDTD_2D fdtd = FDTD_2D(lambda, tau, refractive_index);

  // Сhecking for changes in preconditions.
  if ((fdtd.GetLambda() != lambda) || (fdtd.GetTau() != tau) ||
      (fdtd.GetRefractiveIndex() != refractive_index) || reload_check) {
    fdtd.SetLambda(lambda);
    fdtd.SetTau(tau);
    fdtd.SetRefractiveIndex(refractive_index);
    fdtd.SetParams();
  }

  fdtd.CalcNextLayer(vect_X, vect_Y);

  // Creating JS data for response.
  Napi::Array js_data_X = Napi::Array::New(env, Nx);
  Napi::Array js_data_Y = Napi::Array::New(env, Nx);

  for (size_t i = 0; i < Nx; i++) {
    js_data_X[i] = Napi::Number::New(env, vect_X[i]);
    js_data_Y[i] = Napi::Number::New(env, vect_Y[i]);
  }

  Napi::Object data = Napi::Array::New(env);
  data.Set("dataX", js_data_X);
  data.Set("dataY", js_data_Y);
  data.Set("col", Nx);
  data.Set("currentTick", fdtd.GetCurrentTick());

  return data;
}

// For 2D pure c code.-----------------------------------------
size_t ticks_checker = 0;
DATA_STRUCT *data_struct;

void set_params(DATA_STRUCT *data, float lambda, float tau, float n1);

void init_2D_pure_c() {
  data_struct = (DATA_STRUCT *)malloc(sizeof(DATA_STRUCT));
}

void free_2D_pure_c(DATA_STRUCT *data) {
  free(data->eps);
  free(data->H1);
  free(data->H2);
  free(data->E1);
  free(data->E2);
  free(data);
  std::cout << "FDTD 2d data CLEARED!!";
}

void clear_FDTD_2D_data_pure_c(const Napi::CallbackInfo &info) {
  free_2D_pure_c(data_struct);
}

Napi::Value getFDTD_2D_pure_c(const Napi::CallbackInfo &info) {
  Napi::Env env = info.Env();

  const Napi::Array inputArrayCondition = info[0].As<Napi::Array>();
  // 0 - lambda
  // 1 - tau
  // 2 - n1

  // Reload current params?
  bool reload = static_cast<bool>(info[1].As<Napi::Boolean>());

  int first = 0;                                                        //????
  float lambda = (float)inputArrayCondition[first].As<Napi::Number>();  //????
  float tau = (float)(inputArrayCondition[1].As<Napi::Number>());
  float n1 = (float)(inputArrayCondition[2].As<Napi::Number>());

  // pure c implementation.

  if (ticks_checker == 0) {
    init_2D_pure_c();
    set_params(data_struct, lambda, tau, n1);
  }

  double *vect_X = (double *)malloc(data_struct->Nx * sizeof(double));
  double *vect_Y = (double *)malloc(data_struct->Nx * sizeof(double));

  //   static FDTD_2D fdtd = FDTD_2D(lambda, tau, n1);
  if ((data_struct->lambda != lambda) || (data_struct->tau != tau) ||
      (data_struct->n1 != n1) || reload) {
    ticks_checker = 0;
    set_params(data_struct, lambda, tau, n1);
  }

  // for (int i = 0; i < 400; ++i) {
  calculate_next_time_layer(data_struct, vect_X, vect_Y);
  // }
  ticks_checker = data_struct->ticks;

  size_t Nx = data_struct->Nx;

  // Creating arrays.
  Napi::Array js_data_X = Napi::Array::New(env, Nx);
  Napi::Array js_data_Y = Napi::Array::New(env, Nx);

  // Temporary variables.
  Napi::Number elem;

  for (size_t j = 0; j < Nx; j++) {
    elem = Napi::Number::New(env, vect_X[j]);
    js_data_X[j] = elem;

    elem = Napi::Number::New(env, vect_Y[j]);
    js_data_Y[j] = elem;
  }
  std::cout << ticks_checker << "\n";
  Napi::Object data = Napi::Array::New(env);
  data.Set("dataX", js_data_X);
  data.Set("dataY", js_data_Y);
  data.Set("col", Nx);
  data.Set("currentTick", data_struct->ticks);

  /* Free memory */
  free(vect_X);
  free(vect_Y);

  return data;
}

// Callback method when module is registered with Node.js.
Napi::Object Init(Napi::Env env, Napi::Object exports) {
  exports.Set(Napi::String::New(env, "getFdtd2D"),
              Napi::Function::New(env, GetData2D));

  exports.Set(Napi::String::New(env, "clearFDTD_2D_data"),
              Napi::Function::New(env, clear_FDTD_2D_data_pure_c));

  exports.Set(Napi::String::New(env, "getFDTD_3D"),
              Napi::Function::New(env, getFDTD_3D));

  exports.Set(Napi::String::New(env, "getFDTD_3D_DIFRACTION"),
              Napi::Function::New(env, GetDataDifraction3D));

  exports.Set(Napi::String::New(env, "getFDTD_3D_INTERFERENCE"),
              Napi::Function::New(env, getFDTD_3D_INTERFERENCE));

  // Return `exports` object (always).
  return exports;
}

// Register `FDTD` module which calls `Init` method.
NODE_API_MODULE(FDTD, Init)