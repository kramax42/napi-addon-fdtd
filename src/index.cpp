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
#include "./FDTD/3D_INTERFERENCE/FDTD_3D_INTERFERENCE.h"
#include "globals.h" // Nx, Ny

// https://stackoverflow.com/questions/12573816/what-is-an-undefined-reference-unresolved-external-symbol-error-and-how-do-i-fix/12574420#12574420
// https://www.it-swarm.com.ru/ru/c%2B%2B/chto-takoe-neopredelennaya-ssylka-nerazreshennaya-vneshnyaya-oshibka-simvola-i-kak-ee-ispravit/1069256308/
extern "C" {
#include "./FDTD/pure-c/FDTD-2D.h"
}

// Difraction FDTD_3D.
Napi::Value getFDTD_3D_DIFRACTION(const Napi::CallbackInfo &info) {
  Napi::Env env = info.Env();

  const Napi::Array inputArrayCondition = info[0].As<Napi::Array>();
  // 0 - lambda.
  // 1 - beamsize.
  // 2 - refractive index matrix(flatten).
  // 3 - refractive index matrix(flatten) size for 2x2 - 2.
  // 4 - data return type('Ez' = 0 | 'Hy' = 1 |'Hx' = 2 |'Energy' = 3)

  // Reload current params?
  bool reload = static_cast<bool>(info[1].As<Napi::Boolean>());

  // Refraction index matrix transformation JS -> C++.
  const Napi::Array refrIndexMatrixJS = info[2].As<Napi::Array>();

  // Должени быть кратным 400 !!!!!
  int refrIndexMatrixSize = static_cast<int>(info[3].As<Napi::Number>());
  std::vector<vector<double>> tempMatrix;

  // Data return type('Ez' = 0 | 'Hy' = 1 |'Hx' = 2 |'Energy' = 3)
  int dataReturnType = static_cast<int>(info[4].As<Napi::Number>());

  // Transform input flatten matrix to 2 dimensional.
  for (int i = 0; i < refrIndexMatrixSize; i++) {
    tempMatrix.push_back(std::vector<double>());
    for (int j = 0; j < refrIndexMatrixSize; j++) {
      tempMatrix[i].push_back(
          (double)refrIndexMatrixJS[i * refrIndexMatrixSize + j]
              .As<Napi::Number>());
    }
  }

  // // Output input matrix as 2-dimesioned.
  // for (int i = 0; i < refrIndexMatrixSize; i++)
  // {
  //     std::cout << std::endl;
  //     for (int j = 0; j < refrIndexMatrixSize; j++)
  //     {
  //         std::cout << tempMatrix[i][j] << "\t";
  //     }
  // }

  std::vector<std::vector<double>> refrIndexMatrix;
  size_t koeff = Nx / refrIndexMatrixSize;

  // Initialization refractive index matrix.
  for (int i = 0; i < Nx; i++) {
    refrIndexMatrix.push_back(std::vector<double>());
    for (int j = 0; j < Ny; j++) {
      refrIndexMatrix[i].push_back(1);
    }
  }

  for (int i = 0; i < refrIndexMatrixSize; i++) {
    for (int j = 0; j < refrIndexMatrixSize; j++) {
      for (int k = 0; k < koeff; k++) {
        for (int f = 0; f < koeff; f++) {
          refrIndexMatrix[i * koeff + k][j * koeff + f] = tempMatrix[i][j];
        }
      }
    }
  }

  // int q2 = 0;
  // int q1 = 0;
  // for (int i = 0; i < Nx; i++)
  // {
  //     //std::cout << std::endl;
  //     for (int j = 0; j < Ny; j++)
  //     {
  //         //  std::cout << refrIndexMatrix[i][j] << " ";
  //         if (refrIndexMatrix[i][j] == 2)
  //             q2++;
  //         else
  //             q1++;
  //     }
  // }
  // std::cout << q1 << " ";
  // std::cout << q2 << " ";

  // for (int i = 0; i < Nx; i++)
  // {
  //     std::cout << std::endl;
  //     for (int j = 0; j < Ny; j++)
  //     {
  //         std::cout << refrIndexMatrix[i][j] << " ";
  //     }
  // }
  // std::cout << std::endl;

  // Params transformation JS -> C++.
  int first = 0;                                                          //????
  double lambda = (double)inputArrayCondition[first].As<Napi::Number>();  //????
  double beamsize = (double)inputArrayCondition[1].As<Napi::Number>();
  double n1 = (double)inputArrayCondition[2].As<Napi::Number>();
  double n2 = (double)inputArrayCondition[3].As<Napi::Number>();

  static FDTD_3D_DIFRACTION fdtd_3D =
      FDTD_3D_DIFRACTION(lambda, beamsize, n1, n2, refrIndexMatrix);
  if ((fdtd_3D.getLambda() != lambda) || (fdtd_3D.getBeamsize() != beamsize) ||
      (fdtd_3D.getN1() != n1) || reload) {
    // std::cout << "Works!! " << reload << std::endl;
    fdtd_3D.setLambda(lambda);
    fdtd_3D.setBeamsize(beamsize);
    fdtd_3D.setN1(n1);
    fdtd_3D.setN2(n2);
    fdtd_3D.setParams(refrIndexMatrix);
  }

  vector<double> vectX = {};
  vector<double> vectY = {};
  vector<double> vectEz = {};
  vector<double> vectHy = {};
  vector<double> vectHx = {};
  vector<double> vectEnergy = {};

  vectX.clear();
  vectY.clear();
  vectEz.clear();
  vectHy.clear();
  vectHx.clear();
  vectEnergy.clear();

  fdtd_3D.calcNextLayer(vectX, vectY, vectEz, vectHy, vectHx, vectEnergy);

  size_t Nx = fdtd_3D.getNx() / fdtd_3D.getStep();
  size_t Ny = fdtd_3D.getNy() / fdtd_3D.getStep();
  // size_t Nx = vectX.size();
  // size_t Ny = vectY.size();

  // Creating arrays.
  Napi::Array jsDataX = Napi::Array::New(env, Nx * Ny);
  Napi::Array jsDataY = Napi::Array::New(env, Nx * Ny);
  Napi::Array jsDataEz = Napi::Array::New(env, Nx * Ny);
  Napi::Array jsDataHy = Napi::Array::New(env, Nx * Ny);
  Napi::Array jsDataHx = Napi::Array::New(env, Nx * Ny);
  Napi::Array jsDataEnergy = Napi::Array::New(env, Nx * Ny);

  // Temporary variables.
  Napi::Number elem;

  for (size_t i = 0; i < Nx * Ny; i++) {
    elem = Napi::Number::New(env, vectX[i]);
    jsDataX[i] = elem;

    elem = Napi::Number::New(env, vectY[i]);
    jsDataY[i] = elem;

    elem = Napi::Number::New(env, vectEz[i]);
    jsDataEz[i] = elem;

    elem = Napi::Number::New(env, vectHy[i]);
    jsDataHy[i] = elem;

    elem = Napi::Number::New(env, vectHx[i]);
    jsDataHx[i] = elem;

    elem = Napi::Number::New(env, vectEnergy[i]);
    jsDataEnergy[i] = elem;
  }

  Napi::Object data = Napi::Array::New(env);
  data.Set("dataX", jsDataX);
  data.Set("dataY", jsDataY);
  // data.Set("dataEz", jsDataEz);
  // data.Set("dataHy", jsDataHy);
  // data.Set("dataHx", jsDataHx);
  // data.Set("dataEnergy", jsDataEnergy);
  data.Set("row", Nx);
  data.Set("col", Ny);
  data.Set("currentTick", fdtd_3D.getCurrentTick());
  std::cout << "from C++: " << dataReturnType << std::endl;
  switch (dataReturnType) {
    case 0:
      data.Set("dataEz", jsDataEz);
      break;
    case 1:
      data.Set("dataHy", jsDataHy);
      break;
    case 2:
      data.Set("dataHx", jsDataHy);
      break;
    case 3:
      data.Set("dataEnergy", jsDataEnergy);
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
    // std::cout << "Works!! " << reload << std::endl;
    fdtd_3D.setLambda(lambda);
    fdtd_3D.setBeamsize(beamsize);
    fdtd_3D.setN1(n1);
    fdtd_3D.setParams();
  }

  vector<double> vectX = {};
  vector<double> vectY = {};
  vector<double> vectEz = {};
  vector<double> vectHy = {};
  vector<double> vectHx = {};
  vector<double> vectEnergy = {};

  vectX.clear();
  vectY.clear();
  vectEz.clear();
  vectHy.clear();
  vectHx.clear();
  vectEnergy.clear();

  fdtd_3D.calcNextLayer(vectX, vectY, vectEz, vectHy, vectHx, vectEnergy);

  size_t Nx = fdtd_3D.getNx() / fdtd_3D.getStep();
  size_t Ny = fdtd_3D.getNy() / fdtd_3D.getStep();
  // size_t Nx = vectX.size();
  // size_t Ny = vectY.size();

  // Creating arrays.
  Napi::Array jsDataX = Napi::Array::New(env, Nx * Ny);
  Napi::Array jsDataY = Napi::Array::New(env, Nx * Ny);
  Napi::Array jsDataEz = Napi::Array::New(env, Nx * Ny);
  Napi::Array jsDataHy = Napi::Array::New(env, Nx * Ny);
  Napi::Array jsDataHx = Napi::Array::New(env, Nx * Ny);
  Napi::Array jsDataEnergy = Napi::Array::New(env, Nx * Ny);

  // Temporary variables.
  Napi::Number elem;

  for (size_t i = 0; i < Nx * Ny; i++) {
    elem = Napi::Number::New(env, vectX[i]);
    jsDataX[i] = elem;

    elem = Napi::Number::New(env, vectY[i]);
    jsDataY[i] = elem;

    elem = Napi::Number::New(env, vectEz[i]);
    jsDataEz[i] = elem;

    elem = Napi::Number::New(env, vectHy[i]);
    jsDataHy[i] = elem;

    elem = Napi::Number::New(env, vectHx[i]);
    jsDataHx[i] = elem;

    elem = Napi::Number::New(env, vectEnergy[i]);
    jsDataEnergy[i] = elem;
  }

  Napi::Object data = Napi::Array::New(env);
  data.Set("dataX", jsDataX);
  data.Set("dataY", jsDataY);
  data.Set("dataEz", jsDataEz);
  data.Set("dataHy", jsDataHy);
  data.Set("dataHx", jsDataHx);
  data.Set("dataEnergy", jsDataEnergy);
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

  // static FDTD_3D_DIFRACTION fdtd_3D = FDTD_3D_DIFRACTION(lambda, beamsize,
  // n1, 1.5);
  static FDTD_3D fdtd_3D = FDTD_3D(lambda, beamsize, n1);
  if ((fdtd_3D.getLambda() != lambda) || (fdtd_3D.getBeamsize() != beamsize) ||
      (fdtd_3D.getN1() != n1) || reload) {
    // std::cout << "Works!! " << reload << std::endl;
    fdtd_3D.setLambda(lambda);
    fdtd_3D.setBeamsize(beamsize);
    fdtd_3D.setN1(n1);
    fdtd_3D.setParams();
  }

  vector<double> vectX = {};
  vector<double> vectY = {};
  vector<double> vectEz = {};
  vector<double> vectHy = {};
  vector<double> vectHx = {};
  vector<double> vectEnergy = {};

  vectX.clear();
  vectY.clear();
  vectEz.clear();
  vectHy.clear();
  vectHx.clear();
  vectEnergy.clear();

  fdtd_3D.calcNextLayer(vectX, vectY, vectEz, vectHy, vectHx, vectEnergy);

  size_t Nx = fdtd_3D.getNx() / fdtd_3D.getStep();
  size_t Ny = fdtd_3D.getNy() / fdtd_3D.getStep();
  // size_t Nx = vectX.size();
  // size_t Ny = vectY.size();

  // Creating arrays.
  Napi::Array jsDataX = Napi::Array::New(env, Nx * Ny);
  Napi::Array jsDataY = Napi::Array::New(env, Nx * Ny);
  Napi::Array jsDataEz = Napi::Array::New(env, Nx * Ny);
  Napi::Array jsDataHy = Napi::Array::New(env, Nx * Ny);
  Napi::Array jsDataHx = Napi::Array::New(env, Nx * Ny);
  Napi::Array jsDataEnergy = Napi::Array::New(env, Nx * Ny);

  // Temporary variables.
  Napi::Number elem;

  for (size_t i = 0; i < Nx * Ny; i++) {
    elem = Napi::Number::New(env, vectX[i]);
    jsDataX[i] = elem;

    elem = Napi::Number::New(env, vectY[i]);
    jsDataY[i] = elem;

    elem = Napi::Number::New(env, vectEz[i]);
    jsDataEz[i] = elem;

    elem = Napi::Number::New(env, vectHy[i]);
    jsDataHy[i] = elem;

    elem = Napi::Number::New(env, vectHx[i]);
    jsDataHx[i] = elem;

    elem = Napi::Number::New(env, vectEnergy[i]);
    jsDataEnergy[i] = elem;
  }

  Napi::Object data = Napi::Array::New(env);
  data.Set("dataX", jsDataX);
  data.Set("dataY", jsDataY);
  data.Set("dataEz", jsDataEz);
  data.Set("dataHy", jsDataHy);
  data.Set("dataHx", jsDataHx);
  data.Set("dataEnergy", jsDataEnergy);
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
  vector<double> vect_x = {};
  vector<double> vect_y = {};

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

  fdtd.CalcNextLayer(vect_x, vect_y);

  // Creating JS data for response.
  Napi::Array js_data_x = Napi::Array::New(env, Nx);
  Napi::Array js_data_y = Napi::Array::New(env, Nx);

  for (size_t i = 0; i < Nx; i++) {
    js_data_x[i] = Napi::Number::New(env, vect_x[i]);
    js_data_y[i] = Napi::Number::New(env, vect_y[i]);
  }

  Napi::Object data = Napi::Array::New(env);
  data.Set("dataX", js_data_x);
  data.Set("dataY", js_data_y);
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

  double *vect_x = (double *)malloc(data_struct->Nx * sizeof(double));
  double *vect_y = (double *)malloc(data_struct->Nx * sizeof(double));

  //   static FDTD_2D fdtd = FDTD_2D(lambda, tau, n1);
  if ((data_struct->lambda != lambda) || (data_struct->tau != tau) ||
      (data_struct->n1 != n1) || reload) {
    ticks_checker = 0;
    set_params(data_struct, lambda, tau, n1);
  }

  // for (int i = 0; i < 400; ++i) {
  calculate_next_time_layer(data_struct, vect_x, vect_y);
  // }
  ticks_checker = data_struct->ticks;

  size_t Nx = data_struct->Nx;

  // Creating arrays.
  Napi::Array jsDataX = Napi::Array::New(env, Nx);
  Napi::Array jsDataY = Napi::Array::New(env, Nx);

  // Temporary variables.
  Napi::Number elem;

  for (size_t j = 0; j < Nx; j++) {
    elem = Napi::Number::New(env, vect_x[j]);
    jsDataX[j] = elem;

    elem = Napi::Number::New(env, vect_y[j]);
    jsDataY[j] = elem;
  }
  std::cout << ticks_checker << "\n";
  Napi::Object data = Napi::Array::New(env);
  data.Set("dataX", jsDataX);
  data.Set("dataY", jsDataY);
  data.Set("col", Nx);
  data.Set("currentTick", data_struct->ticks);

  /* Free memory */
  free(vect_x);
  free(vect_y);

  return data;
}

// Callback method when module is registered with Node.js.
Napi::Object Init(Napi::Env env, Napi::Object exports) {
  exports.Set(Napi::String::New(env, "getFDTD_2D"),
              Napi::Function::New(env, GetData2D));
  //_pure_c

  exports.Set(Napi::String::New(env, "clearFDTD_2D_data"),
              Napi::Function::New(env, clear_FDTD_2D_data_pure_c));

  exports.Set(Napi::String::New(env, "getFDTD_3D"),
              Napi::Function::New(env, getFDTD_3D));

  exports.Set(Napi::String::New(env, "getFDTD_3D_DIFRACTION"),
              Napi::Function::New(env, getFDTD_3D_DIFRACTION));

  exports.Set(Napi::String::New(env, "getFDTD_3D_INTERFERENCE"),
              Napi::Function::New(env, getFDTD_3D_INTERFERENCE));

  // Return `exports` object (always).
  return exports;
}

// Register `FDTD` module which calls `Init` method.
NODE_API_MODULE(FDTD, Init)