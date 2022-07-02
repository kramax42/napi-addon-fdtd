#include "./get-data-2d.h"

Napi::Value GetData2D(const Napi::CallbackInfo &info)
{
  Napi::Env env = info.Env();

  // 0 - conditions - [lambda, beamsize]
  // 1 - reload checker.
  // 2 - material matrix(flatten).
  // 3 - rows (material matrix size). rows x rows
  // 4 - epsilon array.
  // 5 - mu array.
  // 6 - sigma array.
  // 7 - data return type(number)   ('Ez' = 0 | 'Hy' = 1 |'Hx' = 2 |'Energy' = 3)
  // 8 - relative source position array.
  const Napi::Array input_array_condition = info[0].As<Napi::Array>();

  // Reload params checker.
  bool reload_check = static_cast<bool>(info[1].As<Napi::Boolean>());

  // Material matrix transformation JS -> C++.
  const Napi::Array material_matrix_js = info[2].As<Napi::Array>();
  const Napi::Array eps_js = info[4].As<Napi::Array>();
  const Napi::Array mu_js = info[5].As<Napi::Array>();
  const Napi::Array sigma_js = info[6].As<Napi::Array>();

  // Must be even.
  int material_matrix_size = static_cast<int>(info[3].As<Napi::Number>());

  // Temporary matrix.
  std::vector<std::vector<int>> temp_matrix;

  // Data return type('Ez' = 0 | 'Hy' = 1 |'Hx' = 2 |'Energy' = 3)
  int data_return_type = static_cast<int>(info[7].As<Napi::Number>());

  // Params transformation JS -> C++.
  double lambda = (double)input_array_condition[(uint32_t)0].As<Napi::Number>();
  double beamsize = (double)input_array_condition[1].As<Napi::Number>();

  // Transform input flatten matrix into 2-dimensional matrix.
  for (int i = 0; i < material_matrix_size; i++)
  {
    temp_matrix.push_back(std::vector<int>());
    for (int j = 0; j < material_matrix_size; j++)
    {
      temp_matrix[i].push_back(
          (int)material_matrix_js[i * material_matrix_size + j]
              .As<Napi::Number>());
    }
  }

  const size_t rows = FdtdPml2D::GetRows();
  const size_t cols = FdtdPml2D::GetCols();

  // Matrix size  coefficient.
  size_t coeff = rows / material_matrix_size;

  // Initialization eps, mu, sigma matrixes.
  std::vector<std::vector<double>> eps_matrix;
  std::vector<std::vector<double>> mu_matrix;
  std::vector<std::vector<double>> sigma_matrix;
  for (int i = 0; i < rows; i++)
  {
    eps_matrix.push_back(std::vector<double>());
    mu_matrix.push_back(std::vector<double>());
    sigma_matrix.push_back(std::vector<double>());
    for (int j = 0; j < cols; j++)
    {
      eps_matrix[i].push_back(0);
      mu_matrix[i].push_back(0);
      sigma_matrix[i].push_back(0);
    }
  }

  // Filling eps, mu, sigma matrixes.
  for (int i = 0; i < material_matrix_size; i++)
  {
    for (int j = 0; j < material_matrix_size; j++)
    {
      for (int k = 0; k < coeff; k++)
      {
        for (int f = 0; f < coeff; f++)
        {
          int index = temp_matrix[i][j];

          // Rotate matrix on 90 degree for correctness in numerical method.
          // eps_matrix[j * coeff + f][i * coeff + k] = static_cast<double>(eps_js[index].As<Napi::Number>());
          // mu_matrix[j * coeff + f][i * coeff + k] = static_cast<double>(mu_js[index].As<Napi::Number>());
          // sigma_matrix[j * coeff + f][i * coeff + k] = static_cast<double>(sigma_js[index].As<Napi::Number>());

          // Without rotate.
          eps_matrix[i * coeff + k][j * coeff + f] = static_cast<double>(eps_js[index].As<Napi::Number>());
          mu_matrix[i * coeff + k][j * coeff + f] = static_cast<double>(mu_js[index].As<Napi::Number>());
          sigma_matrix[i * coeff + k][j * coeff + f] = static_cast<double>(sigma_js[index].As<Napi::Number>());
        }
      }
    }
  }

  const Napi::Array src_position_array = info[8].As<Napi::Array>();
  double relative_src_position_x = (float)src_position_array[(uint32_t)0].As<Napi::Number>();
  double relative_src_position_y = (float)src_position_array[(uint32_t)1].As<Napi::Number>();
  // double relative_src_position = static_cast<double>(info[7].As<Napi::Number>());
  size_t src_position_row = static_cast<size_t>(relative_src_position_y * rows);
  size_t src_position_col = static_cast<size_t>(relative_src_position_x * cols);

  static FdtdPml2D fdtd = FdtdPml2D(eps_matrix, mu_matrix, sigma_matrix, src_position_row, src_position_col);

  // (fdtd_3D.getLambda() != lambda) || (fdtd_3D.getBeamsize() != beamsize) ||
  if (reload_check)
  {
    fdtd.SetParams(eps_matrix, mu_matrix, sigma_matrix, src_position_row, src_position_col);
  }

  std::vector<double> vect_X = {};
  std::vector<double> vect_Y = {};
  std::vector<double> vect_Ez = {};
  std::vector<double> vect_Hy = {};
  std::vector<double> vect_Hx = {};
  std::vector<double> vect_Energy = {};

  double max = 0.001;
  double min = -0.001;

  fdtd.CalcNextLayer(vect_X, vect_Y, vect_Ez, vect_Hy, vect_Hx, vect_Energy, max, min);

  // Matrix sizes.
  size_t client_rows = rows / fdtd.GetStep();
  size_t client_cols = cols / fdtd.GetStep();

  const size_t js_arrays_size = client_rows * client_cols;

  // Creating JS arrays to store C++ arrays.
  Napi::Array js_data_X = Napi::Array::New(env, js_arrays_size);
  Napi::Array js_data_Y = Napi::Array::New(env, js_arrays_size);
  Napi::Array js_data_Ez = Napi::Array::New(env, js_arrays_size);
  Napi::Array js_data_Hy = Napi::Array::New(env, js_arrays_size);
  Napi::Array js_data_Hx = Napi::Array::New(env, js_arrays_size);
  Napi::Array js_data_Energy = Napi::Array::New(env, js_arrays_size);

  // Filling JS arrays with C++ arrays data.
  for (size_t i = 0; i < js_arrays_size; i++)
  {
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
  data.Set("rows", client_rows);
  data.Set("cols", client_cols);
  data.Set("timeStep", fdtd.GetCurrentTimeStep());

  switch (data_return_type)
  {
  case 0:
    data.Set("dataEz", js_data_Ez);
    // max = *std::max_element(std::begin(vect_Ez), std::end(vect_Ez));
    // min = *std::min_element(std::begin(vect_Ez), std::end(vect_Ez));
    break;
  case 1:
    data.Set("dataHy", js_data_Hy);
    // max = *std::max_element(std::begin(vect_Hy), std::end(vect_Hy));
    // min = *std::min_element(std::begin(vect_Hy), std::end(vect_Hy));
    break;
  case 2:
    data.Set("dataHx", js_data_Hx);
    // max = *std::max_element(std::begin(vect_Hx), std::end(vect_Hx));
    // min = *std::min_element(std::begin(vect_Hx), std::end(vect_Hx));
    break;
  case 3:
    data.Set("dataEnergy", js_data_Energy);
    // max = *std::max_element(std::begin(vect_Energy), std::end(vect_Energy));
    // min = *std::min_element(std::begin(vect_Energy), std::end(vect_Energy));
    break;

  default:
    // max = *std::max_element(std::begin(vect_Ez), std::end(vect_Ez));
    // min = *std::min_element(std::begin(vect_Ez), std::end(vect_Ez));
    break;
  }
  // Fill max and min values.
  data.Set("max", max);
  data.Set("min", min);

  return data;
}