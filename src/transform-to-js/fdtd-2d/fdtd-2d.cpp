#include <iostream> 
#include <cmath>
#include <array>
#include "fdtd-2d.h"

Napi::Object Fdtd2D::Init(Napi::Env env, Napi::Object exports) {
    // This method is used to hook the accessor and method callbacks.
    Napi::Function func =
        DefineClass(env,
                    "Fdtd2D",
                    {
                        InstanceMethod("getNextTimeLayer", &Fdtd2D::GetNextTimeLayer),
                    });

    Napi::FunctionReference* constructor = new Napi::FunctionReference();

    // Create a persistent reference to the class constructor. This will allow
    // a function called on a class prototype and a function
    // called on instance of a class to be distinguished from each other.
    *constructor = Napi::Persistent(func);
    exports.Set("Fdtd2D", func);

    // Store the constructor as the add-on instance data. This will allow this
    // add-on to support multiple instances of itself running on multiple worker
    // threads, as well as multiple instances of itself running in different
    // contexts on the same thread.
    //
    // By default, the value set on the environment here will be destroyed when
    // the add-on is unloaded using the `delete` operator, but it is also
    // possible to supply a custom deleter.
    env.SetInstanceData(constructor);

    return exports;
}

Fdtd2D::Fdtd2D(const Napi::CallbackInfo& info)
    : Napi::ObjectWrap<Fdtd2D>(info), fdtd() {
    Napi::Env env = info.Env();

    // Object argument:
    // {
    //   lambda,
    //   beamsize,
    //   isReload,
    //   materialVector, // material matrix(flatten)
    //   eps,
    //   mu,
    //   sigma,
    //   dataReturnType, // ('Ez' = 0 | 'Hy' = 1 |'Hx' = 2 |'Energy' = 3)
    //   srcPosition
    // }

    if (info.Length() != 1) {
        Napi::TypeError::New(env, "Wrong number of arguments")
            .ThrowAsJavaScriptException();
        return;
    }

    if (!info[0].IsObject()) {
        Napi::TypeError::New(env, "Wrong argument! Should be an object.")
            .ThrowAsJavaScriptException();
        return;
    }

    Napi::Object argObj = info[0].As<Napi::Object>().ToObject();

    // lambda
    if (!argObj.Has("lambda") || !argObj.Get("lambda").IsNumber()) {
        Napi::TypeError::New(
            env,
            "Object should contain 'lambda' property and it should be number.")
            .ThrowAsJavaScriptException();
        return;
    }
    double lambda = (double)argObj.Get("lambda").As<Napi::Number>();

    // beamsize
    if (!argObj.Has("beamsize") || !argObj.Get("beamsize").IsNumber()) {
        Napi::TypeError::New(
            env,
            "Object should contain 'beamsize' property and it should be number.")
            .ThrowAsJavaScriptException();
        return;
    }
    double beamsize = (double)argObj.Get("beamsize").As<Napi::Number>();

    // src_type
    if (!argObj.Has("srcType") || !argObj.Get("srcType").IsNumber()) {
        Napi::TypeError::New(
            env,
            "Object should contain 'srcType' property and it should be number.")
            .ThrowAsJavaScriptException();
        return;
    }
    int src_type = (int)argObj.Get("srcType").As<Napi::Number>();

    // Reload params checker.
    if (!argObj.Has("isReload") || !argObj.Get("isReload").IsBoolean()) {
        Napi::TypeError::New(
            env,
            "Object should contain 'isReload' property and it should be boolean.")
            .ThrowAsJavaScriptException();
        return;
    }
    bool reload_check = static_cast<bool>(argObj.Get("isReload").As<Napi::Boolean>());

    // Material matrix(flatten).
    if (!argObj.Has("materialVector") || !argObj.Get("materialVector").IsArray()) {
        Napi::TypeError::New(
            env,
            "Object should contain 'materialVector' property and it should be array.")
            .ThrowAsJavaScriptException();
        return;
    }

    // Material matrix transformation JS -> C++.
    const Napi::Array material_matrix_js = argObj.Get("materialVector").As<Napi::Array>();
    int material_matrix_size = std::sqrt(material_matrix_js.Length());

    // Eps vector.
    if (!argObj.Has("eps") || !argObj.Get("eps").IsArray()) {
        Napi::TypeError::New(
            env,
            "Object should contain 'eps' property and it should be array.")
            .ThrowAsJavaScriptException();
        return;
    }

    // eps transformation JS -> C++.
    const Napi::Array eps_js = argObj.Get("eps").As<Napi::Array>();

    // mu vector.
    if (!argObj.Has("mu") || !argObj.Get("mu").IsArray()) {
        Napi::TypeError::New(
            env,
            "Object should contain 'mu' property and it should be array.")
            .ThrowAsJavaScriptException();
        return;
    }

    // mu transformation JS -> C++.
    const Napi::Array mu_js = argObj.Get("mu").As<Napi::Array>();

    // sigma vector.
    if (!argObj.Has("sigma") || !argObj.Get("sigma").IsArray()) {
        Napi::TypeError::New(
            env,
            "Object should contain 'sigma' property and it should be array.")
            .ThrowAsJavaScriptException();
        return;
    }

    // sigma transformation JS -> C++.
    const Napi::Array sigma_js = argObj.Get("sigma").As<Napi::Array>();

    // Data return type('Ez' = 0 | 'Hy' = 1 |'Hx' = 2 |'Energy' = 3) old
    // Data return type('Hz' = 0 | 'Ex' = 1 |'Ey' = 2 |)
    if (!argObj.Has("dataReturnType") || !argObj.Get("dataReturnType").IsNumber()) {
        Napi::TypeError::New(
            env,
            "Object should contain 'dataReturnType' property and it should be number(0, 1, 2).")
            .ThrowAsJavaScriptException();
        return;
    }

    int data_return_type = static_cast<int>(argObj.Get("dataReturnType").As<Napi::Number>());
    this->data_return_type = data_return_type;

    // std::cout << "before ";//  << material_matrix_size;

    // Temporary matrix.
    
    // const size_t rows = FdtdPml2D::GetRows();
    // const size_t cols = FdtdPml2D::GetCols();
    //  const size_t nx = 220;
    // static const size_t ny = 220;
    static const size_t rows = 220;
    static const size_t cols = 220;

    // std::array<std::array<size_t, cols>, rows> temp_matrix;
    std::vector<std::vector<int>> temp_matrix;
    std::vector<std::vector<int>> temp_matrix_2;

    // std::cout << "--" << material_matrix_size << "--";

    // Transform input flatten matrix into 2-dimensional matrix.
    for (int i = 0; i < material_matrix_size; i++) {
        temp_matrix.push_back(std::vector<int>());
        
        for (int j = 0; j < material_matrix_size; j++) {
            temp_matrix[i].push_back(
                (int)material_matrix_js[i * material_matrix_size + j]
                    .As<Napi::Number>());
            // temp_matrix[i][j] = 
            //     (int)material_matrix_js[i * material_matrix_size + j]
            //         .As<Napi::Number>();
        }
    }


    // Matrix size  coefficient.
    size_t coeff = rows / material_matrix_size;

    
    std::vector<double> eps;
    std::vector<double> mu;
    std::vector<double> sigma;

    for( size_t i = 0; i < eps_js.Length(); ++i) {
        eps.push_back(static_cast<double>(eps_js[i].As<Napi::Number>()));
        mu.push_back(static_cast<double>(mu_js[i].As<Napi::Number>()));
        sigma.push_back(static_cast<double>(sigma_js[i].As<Napi::Number>()));
    }

    // Initialization eps, mu, sigma matrixes.
    std::vector<std::vector<double>> eps_matrix;
    std::vector<std::vector<double>> mu_matrix;
    std::vector<std::vector<double>> sigma_matrix;
    for (int i = 0; i < rows; i++) {
        eps_matrix.push_back(std::vector<double>());
        mu_matrix.push_back(std::vector<double>());
        sigma_matrix.push_back(std::vector<double>());

        temp_matrix_2.push_back(std::vector<int>());

        for (int j = 0; j < cols; j++) {
            eps_matrix[i].push_back(0);
            mu_matrix[i].push_back(0);
            sigma_matrix[i].push_back(0);

            temp_matrix_2[i].push_back(0);
        }
    }

    // Filling eps, mu, sigma matrixes.
    for (int i = 0; i < material_matrix_size; i++) {
        for (int j = 0; j < material_matrix_size; j++) {
            for (int k = 0; k < coeff; k++) {
                for (int f = 0; f < coeff; f++) {
                    int index = temp_matrix[i][j];

                    // temp_matrix_2[i * coeff + k][j * coeff + f] = index;
                    temp_matrix_2[j * coeff + f][i * coeff + k] = index;

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

    // srcPositionArray.
    if (!argObj.Has("srcPosition") || !argObj.Get("srcPosition").IsArray()) {
        Napi::TypeError::New(
            env,
            "Object should contain 'srcPosition' property and it should be array.")
            .ThrowAsJavaScriptException();
        return;
    }

    // srcPosition transformation JS -> C++.
    const Napi::Array src_position_array = argObj.Get("srcPosition").As<Napi::Array>();
    double relative_src_position_x = (float)src_position_array[(uint32_t)0].As<Napi::Number>();
    double relative_src_position_y = (float)src_position_array[(uint32_t)1].As<Napi::Number>();

    size_t src_position_row = static_cast<size_t>(relative_src_position_y * rows);
    size_t src_position_col = static_cast<size_t>(relative_src_position_x * cols);

    // fdtd.SetParams(eps_matrix, mu_matrix, sigma_matrix, src_position_row, src_position_col);
    // fdtd.InitializeFdtd();
    fdtd.InitializeFdtd(temp_matrix_2, eps, mu, sigma, src_position_row, src_position_col, src_type);
}

// Fdtd method in 1D case.
Napi::Value Fdtd2D::GetNextTimeLayer(const Napi::CallbackInfo& info) {
    Napi::Env env = info.Env();

    
    if (info.Length() != 1) {
        Napi::TypeError::New(env, "Wrong number of arguments")
            .ThrowAsJavaScriptException();
        // return;
    }

    if (!info[0].IsObject()) {
        Napi::TypeError::New(env, "Wrong argument! Should be an object.")
            .ThrowAsJavaScriptException();
        // return;
    }

    Napi::Object argObj = info[0].As<Napi::Object>().ToObject();

    // Data return type('Hz' = 0 | 'Ex' = 1 |'Ey' = 2 |)
    if (!argObj.Has("dataReturnType") || !argObj.Get("dataReturnType").IsNumber()) {
        Napi::TypeError::New(
            env,
            "Object should contain 'dataReturnType' property and it should be number(0, 1, 2).")
            .ThrowAsJavaScriptException();
        // return;
    }

    int data_return_type = static_cast<int>(argObj.Get("dataReturnType").As<Napi::Number>());

    std::vector<double> vect_X = {};
    std::vector<double> vect_Y = {};
    std::vector<double> vect_Hz = {};
    // std::vector<double> vect_Ez = {};
    // std::vector<double> vect_Hy = {};
    // std::vector<double> vect_Hx = {};
    // std::vector<double> vect_Energy = {};

    double max = 0.001;
    double min = -0.001;

    // this->fdtd.CalcNextLayer(vect_X, vect_Y, vect_Ez, vect_Hy, vect_Hx, vect_Energy, max, min);
    this->fdtd.CalcNextLayer();

    //hz - 1, ex - 2, ey - 3
    struct FdtdPml2D::Output fdtd_output = this->fdtd.GetValues(data_return_type);
    size_t rows = fdtd_output.rows;
    size_t cols = fdtd_output.cols;

    // const size_t rows = FdtdPml2D::GetRows();
    // const size_t cols = FdtdPml2D::GetCols();
    // const size_t rows = 100;
    // const size_t cols = 100;

    // Matrix sizes.
    // size_t client_rows = rows / fdtd.GetStep();
    // size_t client_cols = cols / fdtd.GetStep();
    size_t client_rows = rows;
    size_t client_cols = cols;

    const size_t js_arrays_size = client_rows * client_cols;

    // Creating JS arrays to store C++ arrays.
    Napi::Array js_data_X = Napi::Array::New(env, js_arrays_size);
    Napi::Array js_data_Y = Napi::Array::New(env, js_arrays_size);
    Napi::Array js_data_Hz = Napi::Array::New(env, js_arrays_size);
    Napi::Array js_data_field = Napi::Array::New(env, js_arrays_size);
    // Napi::Array js_data_Ez = Napi::Array::New(env, js_arrays_size);
    // Napi::Array js_data_Hy = Napi::Array::New(env, js_arrays_size);
    // Napi::Array js_data_Hx = Napi::Array::New(env, js_arrays_size);
    // Napi::Array js_data_Energy = Napi::Array::New(env, js_arrays_size);

    // Filling JS arrays with C++ arrays data.
    for (size_t i = 0; i < js_arrays_size; i++) {
        js_data_X[i] = Napi::Number::New(env, fdtd_output.X[i]);
        js_data_Y[i] = Napi::Number::New(env, fdtd_output.Y[i]);
        js_data_Hz[i] = Napi::Number::New(env, fdtd_output.Hz[i]);
        js_data_field[i] = Napi::Number::New(env, fdtd_output.field[i]);
        // js_data_Ez[i] = Napi::Number::New(env, vect_Ez[i]);
        // js_data_Hy[i] = Napi::Number::New(env, vect_Hy[i]);
        // js_data_Hx[i] = Napi::Number::New(env, vect_Hx[i]);
        // js_data_Energy[i] = Napi::Number::New(env, vect_Energy[i]);
    }

    // Creating JS object to return.
    Napi::Object data = Napi::Array::New(env);
    data.Set("dataX", js_data_X);
    data.Set("dataY", js_data_Y);
    data.Set("rows", client_rows);
    data.Set("cols", client_cols);
    data.Set("dataEz", js_data_Hz);
    data.Set("field", js_data_field);
    data.Set("max", fdtd_output.maxHz);
    data.Set("min", fdtd_output.minHz);
    data.Set("timeStep", fdtd.GetCurrentTimeStep());

    // switch (this->data_return_type) {
    //     case 0:
    //         data.Set("dataEz", js_data_Ez);
    //         // max = *std::max_element(std::begin(vect_Ez), std::end(vect_Ez));
    //         // min = *std::min_element(std::begin(vect_Ez), std::end(vect_Ez));
    //         break;
    //     case 1:
    //         data.Set("dataHy", js_data_Hy);
    //         // max = *std::max_element(std::begin(vect_Hy), std::end(vect_Hy));
    //         // min = *std::min_element(std::begin(vect_Hy), std::end(vect_Hy));
    //         break;
    //     case 2:
    //         data.Set("dataHx", js_data_Hx);
    //         // max = *std::max_element(std::begin(vect_Hx), std::end(vect_Hx));
    //         // min = *std::min_element(std::begin(vect_Hx), std::end(vect_Hx));
    //         break;
    //     case 3:
    //         data.Set("dataEnergy", js_data_Energy);
    //         // max = *std::max_element(std::begin(vect_Energy), std::end(vect_Energy));
    //         // min = *std::min_element(std::begin(vect_Energy), std::end(vect_Energy));
    //         break;

    //     default:
    //         // max = *std::max_element(std::begin(vect_Ez), std::end(vect_Ez));
    //         // min = *std::min_element(std::begin(vect_Ez), std::end(vect_Ez));
    //         break;
    // }
    // // Fill max and min values.
    // data.Set("max", max);
    // data.Set("min", min);
    
    return data;
}