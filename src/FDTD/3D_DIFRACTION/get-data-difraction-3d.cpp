// #include <napi.h>
// #include <string>
// #include <vector>
// #include <iostream>
// #include "../../globals.h" // Nx, Ny
// #include "./FDTD_3D_DIFRACTION.h"

// // Difraction FDTD-3D.
// Napi::Value GetDataDifraction3D(const Napi::CallbackInfo &info) {
//   Napi::Env env = info.Env();

//   const Napi::Array input_array_condition = info[0].As<Napi::Array>();
//   // 0 - lambda.
//   // 1 - reload checker.
//   // 2 - refractive index matrix(flatten).
//   // 3 - refractive index matrix(flatten) size for 2x2 - 2.
//   // 4 - data return type('Ez' = 0 | 'Hy' = 1 |'Hx' = 2 |'Energy' = 3)

//   // Reload params checker.
//   bool reload_check = static_cast<bool>(info[1].As<Napi::Boolean>());

//   // Refraction index matrix transformation JS -> C++.
//   const Napi::Array refr_index_matrix_js = info[2].As<Napi::Array>();

//   // Must be odd.
//   int refr_index_matrix_size = static_cast<int>(info[3].As<Napi::Number>());
//   std::cout << refr_index_matrix_size << "!!!" << endl;
//   std::vector<std::vector<double>> temp_matrix;

//   // Data return type('Ez' = 0 | 'Hy' = 1 |'Hx' = 2 |'Energy' = 3)
//   int data_return_type = static_cast<int>(info[4].As<Napi::Number>());

//   // Transform input flatten matrix into 2-dimensional matrix.
//   for (int i = 0; i < refr_index_matrix_size; i++) {
//     temp_matrix.push_back(std::vector<double>());
//     for (int j = 0; j < refr_index_matrix_size; j++) {
//       temp_matrix[i].push_back(
//           (float)refr_index_matrix_js[i * refr_index_matrix_size + j]
//               .As<Napi::Number>());
//     }
//   }

//   // // Output input matrix as 2-dimesioned.
//   // for (int i = 0; i < refr_index_matrix_size; i++)
//   // {
//   //     std::cout << std::endl;
//   //     for (int j = 0; j < refr_index_matrix_size; j++)
//   //     {
//   //         std::cout << temp_matrix[i][j] << "\t";
//   //     }
//   // }

//   // Matrix size koefficient.
//   size_t koeff = Nx / refr_index_matrix_size;

//   // Initialization refractive index matrix.
//   std::vector<std::vector<double>> refr_index_matrix;
//   for (int i = 0; i < Nx; i++) {
//     refr_index_matrix.push_back(std::vector<double>());
//     for (int j = 0; j < Ny; j++) {
//       refr_index_matrix[i].push_back(1);
//     }
//   }

//   // Filling refractive index matrix.
//   for (int i = 0; i < refr_index_matrix_size; i++) {
//     for (int j = 0; j < refr_index_matrix_size; j++) {
//       for (int k = 0; k < koeff; k++) {
//         for (int f = 0; f < koeff; f++) {
//           refr_index_matrix[i * koeff + k][j * koeff + f] = temp_matrix[i][j];
//         }
//       }
//     }
//   }

//   // int q2 = 0;
//   // int q1 = 0;
//   // for (int i = 0; i < Nx; i++)
//   // {
//   //     //std::cout << std::endl;
//   //     for (int j = 0; j < Ny; j++)
//   //     {
//   //         //  std::cout << refr_index_matrix[i][j] << " ";
//   //         if (refr_index_matrix[i][j] == 2)
//   //             q2++;
//   //         else
//   //             q1++;
//   //     }
//   // }
//   // std::cout << q1 << " ";
//   // std::cout << q2 << " ";

//   for (int i = 0; i < Nx; i++) {
//     std::cout << std::endl;
//     for (int j = 0; j < Ny; j++) {
//       std::cout << refr_index_matrix[i][j] << " ";
//     }
//   }
//   std::cout << std::endl;

//   // Params transformation JS -> C++.
//   int first = 0;  //????
//   double lambda =
//       (double)input_array_condition[first].As<Napi::Number>();  //????
//   double beamsize = (double)input_array_condition[1].As<Napi::Number>();
//   double n1 = (double)input_array_condition[2].As<Napi::Number>();
//   double n2 = (double)input_array_condition[3].As<Napi::Number>();

//   static FDTD_3D_DIFRACTION fdtd_3D =
//       FDTD_3D_DIFRACTION(lambda, beamsize, n1, n2, refr_index_matrix);
//   if ((fdtd_3D.getLambda() != lambda) || (fdtd_3D.getBeamsize() != beamsize) ||
//       (fdtd_3D.getN1() != n1) || reload_check) {
//     // std::cout << "Works!! " << reload << std::endl;
//     fdtd_3D.setLambda(lambda);
//     fdtd_3D.setBeamsize(beamsize);
//     fdtd_3D.setN1(n1);
//     fdtd_3D.setN2(n2);
//     fdtd_3D.setParams(refr_index_matrix);
//   }

//   std::vector<double> vectX = {};
//   std::vector<double> vectY = {};
//   std::vector<double> vectEz = {};
//   std::vector<double> vectHy = {};
//   std::vector<double> vectHx = {};
//   std::vector<double> vectEnergy = {};

//   vectX.clear();
//   vectY.clear();
//   vectEz.clear();
//   vectHy.clear();
//   vectHx.clear();
//   vectEnergy.clear();

//   fdtd_3D.calcNextLayer(vectX, vectY, vectEz, vectHy, vectHx, vectEnergy);

//   size_t Nx = fdtd_3D.getNx() / fdtd_3D.getStep();
//   size_t Ny = fdtd_3D.getNy() / fdtd_3D.getStep();
//   // size_t Nx = vectX.size();
//   // size_t Ny = vectY.size();

//   // Creating arrays.
//   Napi::Array jsDataX = Napi::Array::New(env, Nx * Ny);
//   Napi::Array jsDataY = Napi::Array::New(env, Nx * Ny);
//   Napi::Array jsDataEz = Napi::Array::New(env, Nx * Ny);
//   Napi::Array jsDataHy = Napi::Array::New(env, Nx * Ny);
//   Napi::Array jsDataHx = Napi::Array::New(env, Nx * Ny);
//   Napi::Array jsDataEnergy = Napi::Array::New(env, Nx * Ny);

//   // Temporary variables.
//   Napi::Number elem;

//   for (size_t i = 0; i < Nx * Ny; i++) {
//     elem = Napi::Number::New(env, vectX[i]);
//     jsDataX[i] = elem;

//     elem = Napi::Number::New(env, vectY[i]);
//     jsDataY[i] = elem;

//     elem = Napi::Number::New(env, vectEz[i]);
//     jsDataEz[i] = elem;

//     elem = Napi::Number::New(env, vectHy[i]);
//     jsDataHy[i] = elem;

//     elem = Napi::Number::New(env, vectHx[i]);
//     jsDataHx[i] = elem;

//     elem = Napi::Number::New(env, vectEnergy[i]);
//     jsDataEnergy[i] = elem;
//   }

//   Napi::Object data = Napi::Array::New(env);
//   data.Set("dataX", jsDataX);
//   data.Set("dataY", jsDataY);

//   data.Set("row", Nx);
//   data.Set("col", Ny);
//   data.Set("currentTick", fdtd_3D.getCurrentTick());

//   std::cout << "from C++: " << data_return_type << std::endl;

//   switch (data_return_type) {
//     case 0:
//       data.Set("dataEz", jsDataEz);
//       break;
//     case 1:
//       data.Set("dataHy", jsDataHy);
//       break;
//     case 2:
//       data.Set("dataHx", jsDataHy);
//       break;
//     case 3:
//       data.Set("dataEnergy", jsDataEnergy);
//       break;

//     default:
//       break;
//   }

//   return data;
// }

