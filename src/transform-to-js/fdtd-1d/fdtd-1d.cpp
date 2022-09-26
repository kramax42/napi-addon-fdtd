#include "fdtd-1d.h"

Napi::Object Fdtd1D::Init(Napi::Env env, Napi::Object exports) {
    // This method is used to hook the accessor and method callbacks.
    Napi::Function func =
        DefineClass(env,
                    "Fdtd1D",
                    {
                        InstanceMethod("getNextTimeLayer", &Fdtd1D::GetNextTimeLayer),
                    });

    Napi::FunctionReference* constructor = new Napi::FunctionReference();

    // Create a persistent reference to the class constructor. This will allow
    // a function called on a class prototype and a function
    // called on instance of a class to be distinguished from each other.
    *constructor = Napi::Persistent(func);
    exports.Set("Fdtd1D", func);

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

Fdtd1D::Fdtd1D(const Napi::CallbackInfo& info)
    : Napi::ObjectWrap<Fdtd1D>(info), fdtd() {
    Napi::Env env = info.Env();

    // Object argument:
    // {
    //   omega,
    //   tau,
    //   isReload,
    //   materialVector,
    //   eps,
    //   mu,
    //   sigma,
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

    // Omega
    if (!argObj.Has("omega") || !argObj.Get("omega").IsNumber()) {
        Napi::TypeError::New(
            env,
            "Object should contain 'omega' property and it should be number.")
            .ThrowAsJavaScriptException();
        return;
    }
    double omega = (double)argObj.Get("omega").As<Napi::Number>();

    // Tau
    if (!argObj.Has("tau") || !argObj.Get("tau").IsNumber()) {
        Napi::TypeError::New(
            env,
            "Object should contain 'tau' property and it should be number.")
            .ThrowAsJavaScriptException();
        return;
    }
    double tau = (double)argObj.Get("tau").As<Napi::Number>();

    // Reload params checker.
    if (!argObj.Has("isReload") || !argObj.Get("isReload").IsBoolean()) {
        Napi::TypeError::New(
            env,
            "Object should contain 'isReload' property and it should be boolean.")
            .ThrowAsJavaScriptException();
        return;
    }
    bool reload_check = static_cast<bool>(argObj.Get("isReload").As<Napi::Boolean>());

    // Material vector.
    if (!argObj.Has("materialVector") || !argObj.Get("materialVector").IsArray()) {
        Napi::TypeError::New(
            env,
            "Object should contain 'materialVector' property and it should be array.")
            .ThrowAsJavaScriptException();
        return;
    }

    // Material matrix transformation JS -> C++.
    const Napi::Array material_matrix_js = argObj.Get("materialVector").As<Napi::Array>();
    int material_vector_size = material_matrix_js.Length();

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

    // Temporary vector.
    std::vector<int> temp_vector;

    // Transform input JS vector to C++ vector.
    for (int j = 0; j < material_vector_size; j++) {
        temp_vector.push_back(
            (int)material_matrix_js[j]
                .As<Napi::Number>());
    }

    const size_t grid_size = FdtdPml1D::GetGridSize();

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
    double relative_src_position = (float)src_position_array[(uint32_t)0].As<Napi::Number>();
    size_t src_position = static_cast<size_t>(relative_src_position * grid_size);

    // Matrix size  coefficient.
    size_t coeff = grid_size / material_vector_size;

    // Initialization eps, mu, sigma vectors.
    std::vector<double> eps_vector;
    std::vector<double> mu_vector;
    std::vector<double> sigma_vector;

    // Filling eps, mu, sigma vectors.
    for (int i = 0; i < material_vector_size; i++) {
        for (int j = 0; j < coeff; j++) {
            int index = temp_vector[i];
            eps_vector.push_back(static_cast<double>(eps_js[index].As<Napi::Number>()));
            mu_vector.push_back(static_cast<double>(mu_js[index].As<Napi::Number>()));
            sigma_vector.push_back(static_cast<double>(sigma_js[index].As<Napi::Number>()));
        }
    }

    fdtd.SetParams(tau,
                   omega,
                   eps_vector,
                   mu_vector,
                   sigma_vector,
                   src_position);
}

// Fdtd method in 1D case.
Napi::Value Fdtd1D::GetNextTimeLayer(const Napi::CallbackInfo& info) {
    Napi::Env env = info.Env();

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
    // if (reload_check) {
    //     fdtd.SetParams(tau,
    //                    omega,
    //                    eps_vector,
    //                    mu_vector,
    //                    sigma_vector,
    //                    src_position);
    // }

    // Containers to storage coordinates.
    std::vector<double> vect_x = {};
    std::vector<double> vect_ex = {};
    std::vector<double> vect_hy = {};

    this->fdtd.Calculation(vect_x, vect_ex, vect_hy);

    const size_t grid_size = FdtdPml1D::GetGridSize();

    // Creating JS data for response.
    Napi::Array js_data_x = Napi::Array::New(env, grid_size);
    Napi::Array js_data_ex = Napi::Array::New(env, grid_size);
    Napi::Array js_data_hy = Napi::Array::New(env, grid_size);

    for (size_t i = 0; i < grid_size; i++) {
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