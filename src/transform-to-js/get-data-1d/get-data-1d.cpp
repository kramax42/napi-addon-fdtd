#include "./get-data-1d.h"

// Fdtd method in 1D case.
Napi::Value GetData1D(const Napi::CallbackInfo &info)
{
  Napi::Env env = info.Env();

  // 0 - conditions - [omega, tau]
  // 1 - reload checker.
  // 2 - material vector.
  // 3 - cols (material vector size). rows x rows
  // 4 - epsilon array.
  // 5 - mu array.
  // 6 - sigma array.
  // 7 - relative source position array.

  const Napi::Array input_array_condition = info[0].As<Napi::Array>();
  double omega = (float)input_array_condition[(uint32_t)0].As<Napi::Number>();
  double tau = (float)(input_array_condition[1].As<Napi::Number>());

  // Reload params checker.
  bool reload_check = static_cast<bool>(info[1].As<Napi::Boolean>());

  // Material------ matrix transformation JS -> C++.
  const Napi::Array material_matrix_js = info[2].As<Napi::Array>();
  const Napi::Array eps_js = info[4].As<Napi::Array>();
  const Napi::Array mu_js = info[5].As<Napi::Array>();
  const Napi::Array sigma_js = info[6].As<Napi::Array>();



  // Must be even.
  int material_vector_size = static_cast<int>(info[3].As<Napi::Number>());

  // Temporary vector.
  std::vector<int> temp_vector;

  // Transform input vector to C++ vector.

  for (int j = 0; j < material_vector_size; j++)
  {
    temp_vector.push_back(
        (int)material_matrix_js[j]
            .As<Napi::Number>());
  }

  const size_t grid_size = FdtdPml1D::GetGridSize();

  const Napi::Array src_position_array = info[7].As<Napi::Array>();
  double relative_src_position = (float)src_position_array[(uint32_t)0].As<Napi::Number>();
  // double relative_src_position = static_cast<double>(info[7].As<Napi::Number>());
  size_t src_position = static_cast<size_t>(relative_src_position * grid_size);

  // Matrix size  coefficient.
  size_t coeff = grid_size / material_vector_size;

  // Initialization eps, mu, sigma vectors.
  std::vector<double> eps_vector;
  std::vector<double> mu_vector;
  std::vector<double> sigma_vector;

  // Filling eps, mu, sigma vectors.
  for (int i = 0; i < material_vector_size; i++)
  {
    for (int j = 0; j < coeff; j++)
    {
      int index = temp_vector[i];
      eps_vector.push_back(static_cast<double>(eps_js[index].As<Napi::Number>()));
      mu_vector.push_back(static_cast<double>(mu_js[index].As<Napi::Number>()));
      sigma_vector.push_back(static_cast<double>(sigma_js[index].As<Napi::Number>()));
    }
  }

  

  // Using static variable to save save data for different function call.
  static FdtdPml1D fdtd = FdtdPml1D(tau,
                  omega,
                  eps_vector,
                  mu_vector,
                  sigma_vector,
                  src_position);

  // (fdtd.GetLambda() != omega) ||
  // if (
  //     (fdtd.GetTau() != tau)
  //     // || (fdtd.GetSourcePosition() != src)
  //     || reload_check)
  // {
  //   // fdtd.setLambda(omega);

  //   fdtd.setSourcePosition(src);
  //   // fdtd.setParams(eps_vector, sigma_vector, source_position_vector);
  //   fdtd.setParams(tau, lambda, epsilon_vector, mu_ve);
  // }
  if (reload_check)
  {
    fdtd.SetParams(tau,
                  omega,
                  eps_vector,
                  mu_vector,
                  sigma_vector,
                  src_position);
  }



  // Containers to storage coordinates.
  std::vector<double> vect_x = {};
  std::vector<double> vect_ex = {};
  std::vector<double> vect_hy = {};

  fdtd.Calculation(vect_x, vect_ex, vect_hy);

  // Creating JS data for response.
  Napi::Array js_data_x = Napi::Array::New(env, grid_size);
  Napi::Array js_data_ex = Napi::Array::New(env, grid_size);
  Napi::Array js_data_hy = Napi::Array::New(env, grid_size);

  for (size_t i = 0; i < grid_size; i++)
  {
    js_data_x[i] = Napi::Number::New(env, vect_x[i]);
    js_data_ex[i] = Napi::Number::New(env, vect_ex[i]);
    js_data_hy[i] = Napi::Number::New(env, vect_hy[i]);
  }

  double maxEx = *std::max_element(std::begin(vect_ex), std::end(vect_ex));
  double minEx = *std::min_element(std::begin(vect_ex), std::end(vect_ex));

  double maxHy = *std::max_element(std::begin(vect_hy), std::end(vect_hy));
  double minHy = *std::min_element(std::begin(vect_hy), std::end(vect_hy));

  Napi::Object data = Napi::Array::New(env);
  data.Set("dataX", js_data_x);
  data.Set("dataEx", js_data_ex);
  data.Set("dataHy", js_data_hy);
  data.Set("col", grid_size);
  data.Set("currentTick", fdtd.GetCurrentTimeStep());

  data.Set("maxEx", maxEx);
  data.Set("minEx", minEx);
  data.Set("maxHy", maxHy);
  data.Set("minHy", minHy);

  return data;
}