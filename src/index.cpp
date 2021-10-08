//tutorial - how to make node addons
//https://medium.com/jspoint/a-simple-guide-to-load-c-c-code-into-node-js-javascript-applications-3fcccf54fd32
//https://github.com/course-one/native-addon-starter-kit/blob/master/src/index.cpp
//https://github.com/nodejs/node-addon-examples


// NAPI DESTINATION
//C:/Users/BOSS/AppData/Local/node-gyp/Cache/15.11.0/include/node

#include <napi.h>
#include <string>
#include <vector>
#include <math.h>
#include <iostream>


#include "./FDTD/2D/FDTD_2D.h"
#include "./FDTD/3D/FDTD_3D.h"
#include "./FDTD/3D_DIFRACTION/FDTD_3D_DIFRACTION.h"
#include "./FDTD/3D_INTERFERENCE/FDTD_3D_INTERFERENCE.h"

// Difraction FDTD_3D.
Napi::Value getFDTD_3D_DIFRACTION(const Napi::CallbackInfo &info)
{
    Napi::Env env = info.Env();

    const Napi::Array inputArrayCondition = info[0].As<Napi::Array>();
    // 0 - lambda.
    // 1 - beamsize.
    // 2 - n1.
    // 3 - n2.

    // Reload current params?
    bool reload = static_cast<bool>(info[1].As<Napi::Boolean>());

    // Refraction index matrix transformation JS -> C++.
    const Napi::Array refrIndexMatrixJS = info[2].As<Napi::Array>();

    // Должени быть кратным 400 !!!!!
    int refrIndexMatrixSize = static_cast<int>(info[3].As<Napi::Number>());
    std::vector<vector<double>> tempMatrix;

    for (int i = 0; i < refrIndexMatrixSize; i++)
    {
        tempMatrix.push_back(std::vector<double>());
        for (int j = 0; j < refrIndexMatrixSize; j++)
        {
            tempMatrix[i].push_back((double)refrIndexMatrixJS[i * refrIndexMatrixSize + j].As<Napi::Number>());
        }
    }

    static const size_t Nxx = 400;
    static const size_t Nyx = 400;
    std::vector<std::vector<double>> refrIndexMatrix;
    size_t koeff = Nxx / refrIndexMatrixSize;

    for (int i = 0; i < Nxx; i++)
    {
        refrIndexMatrix.push_back(std::vector<double>());
        for (int j = 0; j < Nyx; j++)
        {
            refrIndexMatrix[i].push_back(1);
        }
    }

    // for (int i = 0; i < Nx; i++)
    // {
    //     for (int j = 0; j < Ny; j++)
    //     {
    //         refrIndexMatrix[i][j] = 1;
    //     }
    // }

    for (int i = 0; i < refrIndexMatrixSize; i++)
    {
        for (int j = 0; j < refrIndexMatrixSize; j++)
        {
            for (int k = 0; k < koeff; k++)
            {
                for (int f = 0; f < koeff; f++)
                {
                    refrIndexMatrix[i * koeff + k][j*koeff + f] = tempMatrix[i][j];
                }
            }
        }
    }

    // for (int i = 0; i < Nxx; i++)
    // {
    //     std::cout << std::endl;
    //     for (int j = 0; j < Nyx; j++)
    //     {
    //         std::cout << refrIndexMatrix[i][j] << " ";
    //     }
    // }
    // std::cout << std::endl;

    // Params transformation JS -> C++.
    int first = 0;                                                         //????
    double lambda = (double)inputArrayCondition[first].As<Napi::Number>(); //????
    double beamsize = (double)(inputArrayCondition[1].As<Napi::Number>());
    double n1 = (double)(inputArrayCondition[2].As<Napi::Number>());
    double n2 = (double)(inputArrayCondition[3].As<Napi::Number>());

    static FDTD_3D_DIFRACTION fdtd_3D = FDTD_3D_DIFRACTION(lambda, beamsize, n1, n2, refrIndexMatrix);
    if ((fdtd_3D.getLambda() != lambda) || (fdtd_3D.getBeamsize() != beamsize) || (fdtd_3D.getN1() != n1) || reload)
    {
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
            //size_t Nx = vectX.size();
            //size_t Ny = vectY.size();

            // Creating arrays.
            Napi::Array jsDataX = Napi::Array::New(env, Nx * Ny);
            Napi::Array jsDataY = Napi::Array::New(env, Nx * Ny);
            Napi::Array jsDataEz = Napi::Array::New(env, Nx * Ny);
            Napi::Array jsDataHy = Napi::Array::New(env, Nx * Ny);
            Napi::Array jsDataHx = Napi::Array::New(env, Nx * Ny);
            Napi::Array jsDataEnergy = Napi::Array::New(env, Nx * Ny);

            // Temporary variables.
            Napi::Number elem;

            for (size_t i = 0; i < Nx * Ny; i++)
            {
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

        // Interference FDTD_3D.
        Napi::Value getFDTD_3D_INTERFERENCE(const Napi::CallbackInfo &info)
        {
            Napi::Env env = info.Env();

            const Napi::Array inputArrayCondition = info[0].As<Napi::Array>();
            // 0 - lambda.
            // 1 - beamsize.
            // 2 - n1.

            // Reload current params?
            bool reload = static_cast<bool>(info[1].As<Napi::Boolean>());

            int first = 0;                                                         //????
            double lambda = (double)inputArrayCondition[first].As<Napi::Number>(); //????
            double beamsize = (double)(inputArrayCondition[1].As<Napi::Number>());
            double n1 = (double)(inputArrayCondition[2].As<Napi::Number>());

            static FDTD_3D_INTERFERENCE fdtd_3D = FDTD_3D_INTERFERENCE(lambda, beamsize, n1);
            if ((fdtd_3D.getLambda() != lambda) || (fdtd_3D.getBeamsize() != beamsize) || (fdtd_3D.getN1() != n1) || reload)
            {
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
            //size_t Nx = vectX.size();
            //size_t Ny = vectY.size();

            // Creating arrays.
            Napi::Array jsDataX = Napi::Array::New(env, Nx * Ny);
            Napi::Array jsDataY = Napi::Array::New(env, Nx * Ny);
            Napi::Array jsDataEz = Napi::Array::New(env, Nx * Ny);
            Napi::Array jsDataHy = Napi::Array::New(env, Nx * Ny);
            Napi::Array jsDataHx = Napi::Array::New(env, Nx * Ny);
            Napi::Array jsDataEnergy = Napi::Array::New(env, Nx * Ny);

            // Temporary variables.
            Napi::Number elem;

            for (size_t i = 0; i < Nx * Ny; i++)
            {
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
        Napi::Value getFDTD_3D(const Napi::CallbackInfo &info)
        {
            Napi::Env env = info.Env();

            const Napi::Array inputArrayCondition = info[0].As<Napi::Array>();
            // 0 - lambda.
            // 1 - beamsize.
            // 2 - n1.

            // Reload current params?
            bool reload = static_cast<bool>(info[1].As<Napi::Boolean>());

            int first = 0;                                                         //????
            double lambda = (double)inputArrayCondition[first].As<Napi::Number>(); //????
            double beamsize = (double)(inputArrayCondition[1].As<Napi::Number>());
            double n1 = (double)(inputArrayCondition[2].As<Napi::Number>());

            //static FDTD_3D_DIFRACTION fdtd_3D = FDTD_3D_DIFRACTION(lambda, beamsize, n1, 1.5);
            static FDTD_3D fdtd_3D = FDTD_3D(lambda, beamsize, n1);
            if ((fdtd_3D.getLambda() != lambda) || (fdtd_3D.getBeamsize() != beamsize) || (fdtd_3D.getN1() != n1) || reload)
            {
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
            //size_t Nx = vectX.size();
            //size_t Ny = vectY.size();

            // Creating arrays.
            Napi::Array jsDataX = Napi::Array::New(env, Nx * Ny);
            Napi::Array jsDataY = Napi::Array::New(env, Nx * Ny);
            Napi::Array jsDataEz = Napi::Array::New(env, Nx * Ny);
            Napi::Array jsDataHy = Napi::Array::New(env, Nx * Ny);
            Napi::Array jsDataHx = Napi::Array::New(env, Nx * Ny);
            Napi::Array jsDataEnergy = Napi::Array::New(env, Nx * Ny);

            // Temporary variables.
            Napi::Number elem;

            for (size_t i = 0; i < Nx * Ny; i++)
            {
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

        Napi::Value getFDTD_2D(const Napi::CallbackInfo &info)
        {
            Napi::Env env = info.Env();

            const Napi::Array inputArrayCondition = info[0].As<Napi::Array>();
            // 0 - lambda
            // 1 - tau
            // 2 - n1

            // Reload current params?
            bool reload = static_cast<bool>(info[1].As<Napi::Boolean>());

            int first = 0;                                                         //????
            double lambda = (double)inputArrayCondition[first].As<Napi::Number>(); //????
            double tau = (double)(inputArrayCondition[1].As<Napi::Number>());
            double n1 = (double)(inputArrayCondition[2].As<Napi::Number>());

            vector<double> vectX = {};
            vector<double> vectY = {};

            static FDTD_2D fdtd = FDTD_2D(lambda, tau, n1);
            if ((fdtd.getLambda() != lambda) || (fdtd.getTau() != tau) || (fdtd.getN1() != n1) || reload)
            {
                //    std::cout << "Works!! " << reload << std::endl;
                fdtd.setLambda(lambda);
                fdtd.setTau(tau);
                fdtd.setN1(n1);
                fdtd.setParams();
            }

            vectX.clear();
            vectY.clear();
            fdtd.calcNextLayer(vectX, vectY);

            size_t Nx = vectX.size();

            // Creating arrays.
            Napi::Array jsDataX = Napi::Array::New(env, Nx);
            Napi::Array jsDataY = Napi::Array::New(env, Nx);

            // Temporary variables.
            Napi::Number elem;

            for (size_t j = 0; j < Nx; j++)
            {
                elem = Napi::Number::New(env, vectX[j]);
                jsDataX[j] = elem;

                elem = Napi::Number::New(env, vectY[j]);
                jsDataY[j] = elem;
            }

            Napi::Object data = Napi::Array::New(env);
            data.Set("dataX", jsDataX);
            data.Set("dataY", jsDataY);
            data.Set("col", Nx);
            data.Set("currentTick", fdtd.getCurrentTick());

            return data;
        }

        // Callback method when module is registered with Node.js.
        Napi::Object Init(Napi::Env env, Napi::Object exports)
        {
            exports.Set(Napi::String::New(env, "getFDTD_2D"),
                        Napi::Function::New(env, getFDTD_2D));

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