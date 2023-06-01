// 2D FDTD with UPML and TF/SF.

// Based on code from book:
// k_Fz_1_new
#ifndef TF_SF_H
#define TF_SF_H

#include <math.h>
#include <napi.h>

#include <algorithm>  // std::fill_n
#include <vector>
#include <array>
#include <iostream>
#include <fstream>


class TFSF {
    //  Matlab array index start
    const size_t mmm = 0;

    //  Physical constants
    //  Absolute vacuum permittivity
    const double epsilon_0 = 8.854*1e-12;

    const double pi = 3.1415926535897932385;

    // Absolute vacuum permeability
    const double mu_0 = 4 * pi * 1e-7;


    // Light speed in vacuum
    const double light_speed = 2.99792458e8;

    //  Main parameters
    //  Calculation area length per x and y axes
    static constexpr double area_width = 2.5;
    static constexpr double area_height = 2.5;

    // Uniform grid points for x and y axes
    // static const size_t nx = 500;
    // static const size_t ny = 500;
    // static const size_t nx = 30;
    // static const size_t ny = 30;
    static const size_t nx = 220;
    static const size_t ny = 220;

    // Perfect match layer (upml_width) thickness in uniform grid cells
    static const size_t upml_width = 10;
    // static const size_t upml_width = 3;

    //  End time [s]
    double t_end = 15e-9;

    // Excitation source amplitude and frequency [Hz]
    double E0 = 1.0;     
    double frequency = 2.0e+9;   // 2 GHz


    //Width of alinea between total field area and calculation area border
    //(scattered field interface width) [m]
    static constexpr double len_tf_dx = 0.5;
    static constexpr double len_tf_dy = 0.5;

    // Grid calculation
    static constexpr double dx = area_width / nx; //C
    static constexpr double dy = area_width / ny; //C

    // % C++ ignore
    // X[nx];
    // Y[ny];

    // % C++ ignore
    // for i = 1:nx
    //     X(i) = (i-1)*dx - area_width/2;
    // end

    // % C++ ignore
    // for i = 1:ny
    //     Y(i) = (i-1)*dy - area_height/2;
    // end


    // Time step from CFL condition
    const double CFL_FACTOR = 0.99;
    // const double dt = ( 1/light_speed/sqrt( 1/(dx^2) + 1/(dy^2) ) )*CFL_FACTOR;
    const double dt = CFL_FACTOR / (light_speed * sqrt(std::pow(1 / dx, 2) + std::pow(1 / dy, 2))); // C


    // number_of_iterations = round(t_end/dt);

    // Geometry matrix 
    // Fill with 0 in code below
    std::array<std::array<double, ny>, nx> Index; // C
    std::array<std::array<double, ny>, nx> IndexX; // C
    std::array<std::array<double, ny>, nx> IndexY; // C

    // Calculate size of total field area in TF/SF formalism
    // static constexpr size_t nx_a = round(len_tf_dx/dx);
    // static constexpr size_t nx_b = round((area_width-len_tf_dx)/dx);
    // static constexpr size_t ny_a = round(len_tf_dy/dy);
    // static constexpr size_t ny_b = round((area_height-len_tf_dy)/dy);
    // static constexpr size_t nx_a = 100;
    // static constexpr size_t nx_b = 400;
    // static constexpr size_t ny_a = 100;
    // static constexpr size_t ny_b = 400;

    // static constexpr size_t nx_a = 6;
    // static constexpr size_t nx_b = 24;
    // static constexpr size_t ny_a = 6;
    // static constexpr size_t ny_b = 24;

    static constexpr size_t nx_a = 60;
    static constexpr size_t nx_b = 160;
    static constexpr size_t ny_a = 60;
    static constexpr size_t ny_b = 160;

    // static constexpr size_t nx_a = 15;
    // static constexpr size_t nx_b = 205;
    // static constexpr size_t ny_a = 15;
    // static constexpr size_t ny_b = 205;

    // static constexpr size_t nx_a = 10;
    // static constexpr size_t nx_b = 209;
    // static constexpr size_t ny_a = 10;
    // static constexpr size_t ny_b = 209;

    // Pre-allocate 1D fields for TF/SF interface 
    // TM components
    std::array<double, nx_b+1> Ez_1D; // C
    std::array<double, nx_b+1> Fz_1D; // C
    std::array<double, nx_b> Hy_1D; // C

    std::array<double, nx_b+1> k_Fz_a; // C
    std::array<double, nx_b+1> k_Fz_b; // C
    std::array<double, nx_b> k_Hy_a; // C
    std::array<double, nx_b> k_Hy_b; // C

    //  TE components
    std::array<double, nx_b+1> Hz_1D; // C
    std::array<double, nx_b> Ey_1D; // C


    std::array<double, nx_b+1> k_Hz_a; // C half
    std::array<double, nx_b+1> k_Hz_b; // C half
    std::array<double, nx_b> k_Ey_a; // C
    std::array<double, nx_b> k_Ey_b; // C

    const double ka_max = 1;
    const double m = 4;
    const double R_err = 1e-16;


    //  auxiliary
    std::array<std::array<double, ny-2>, nx-2> k_Fz_1_new; // C
    std::array<std::array<double, ny-2>, nx-2> k_Fz_2_new; // C
    std::array<std::array<double, ny-2>, nx-2> k_Ez_1_new; // C
    std::array<std::array<double, ny-2>, nx-2> k_Ez_2_new; // C
    std::array<std::array<double, ny-1>, nx> k_Gx_1_new; // C
    std::array<std::array<double, ny-1>, nx> k_Gx_2_new; // C
    std::array<std::array<double, ny-1>, nx> k_Hx_1_new; // C
    std::array<std::array<double, ny-1>, nx> k_Hx_2_new; // C
    std::array<std::array<double, ny-1>, nx> k_Ex_1_new; // C
    std::array<std::array<double, ny-1>, nx> k_Ex_2_new; // C
    
    std::array<std::array<double, ny>, nx-1> k_Gy_1_new; // C
    std::array<std::array<double, ny>, nx-1> k_Gy_2_new; // C half  

    std::array<std::array<double, ny>, nx-1> k_Ey_1_new; // C   
    std::array<std::array<double, ny>, nx-1> k_Ey_2_new; // C half  
    std::array<std::array<double, ny>, nx-1> k_Hy_1_new; // C half
    std::array<std::array<double, ny>, nx-1> k_Hy_2_new; // C half  

    std::array<std::array<double, ny-2>, nx-2> k_Hz_1_new; // C
    std::array<std::array<double, ny-2>, nx-2> k_Hz_2_new; // C half

    std::array<std::array<double, ny-2>, nx-2> M0; // c
    std::array<std::array<double, ny-2>, nx-2> M1; // c 
    std::array<std::array<double, ny-1>, nx> M2; // c
    std::array<std::array<double, ny>, nx-1> M3; // c
    

    std::array<double, nx_b-1> Fz_1D_r; // C
    std::array<std::array<double, ny-2>, nx-2> Fz_r1;

    std::array<std::array<double, ny-2>, nx-2> Fz_r; //
    std::array<std::array<double, ny-2>, nx-2> Tz_r;

    std::array<std::array<double, ny>, nx-1> Gy_r;
    std::array<std::array<double, ny-1>, nx> Gx_r;
    std::array<std::array<double, ny-1>, nx> Gx_r1;
    std::array<std::array<double, ny>, nx-1> Gy_r1;

    std::array<std::array<double, ny-2>, nx-2> K_a_new;
    std::array<std::array<double, ny-2>, nx-2> K_b_new;

    // Allocate 2D arrays
    // TM physical and auxiliary fields
    std::array<std::array<double, ny>, nx> Fz;
    std::array<std::array<double, ny>, nx> Tz;
    std::array<std::array<double, ny-1>, nx> Gx;
    std::array<std::array<double, ny>, nx-1> Gy;    
    std::array<std::array<double, ny>, nx> Ez;    
    std::array<std::array<double, ny-1>, nx> Hx;    
    std::array<std::array<double, ny>, nx-1> Hy;    

    // TE physical and auxiliary fields
    std::array<std::array<double, ny>, nx> Wz;
    std::array<std::array<double, ny-1>, nx> Mx;
    std::array<std::array<double, ny>, nx-1> My;
    std::array<std::array<double, ny>, nx> Hz;
    std::array<std::array<double, ny-1>, nx> Ex;
    std::array<std::array<double, ny>, nx-1> Ey;

    // Allocate UPML FDTD 1D coefficient arrays
    // TM coefficients
    std::array<double, ny> k_Fz_1; // C
    std::array<double, ny> k_Fz_2; // C
    std::array<double, nx> k_Ez_1; // C
    std::array<double, nx> k_Ez_2; // C
    std::array<double, ny-1> k_Gx_1; // C
    std::array<double, ny-1> k_Gx_2; // C
    std::array<double, nx> k_Hx_1; // C
    std::array<double, nx> k_Hx_2; // C

    std::array<double, ny-1> k_Gy_1; // C
    std::array<double, ny-1> k_Gy_2; // C half

    std::array<double, ny> k_Hy_1; // C half
    std::array<double, ny> k_Hy_2; // C half

    // TM coefficients
    std::array<double, nx> k_Hz_1; // C half
    std::array<double, nx> k_Hz_2; // C half
    std::array<double, nx> k_Ex_1; // C
    std::array<double, nx> k_Ex_2; // C
    std::array<double, ny> k_Ey_1; // C
    std::array<double, ny> k_Ey_2; // C half


    // Number of materials (include vacuum backrgound)
    static const size_t number_of_materials = 2;
    //Materials matrix 
    // Material = zeros(number_of_materials,3);
    
    // % eps, mu, sigma
    // % Background relative permittivity, relative permeability and absolute
    // % conductivity
    std::array<std::array<double, 3>, number_of_materials> Material = {{
        {1.0, 1.0, 0.0}, //vacuum
        {3, 3, 9.0e+70}
    }};

    std::array<double, number_of_materials> K_a; // C
    std::array<double, number_of_materials> K_a1; // C
    std::array<double, number_of_materials> K_b; // C
    std::array<double, number_of_materials> K_b1; // C

    // tmp here start
    size_t mat_count = 3;

// mat_conf = [
//     80, 15, 100, 20;
//     100, 10, 100, 20;
//     120, 10, 100, 20;
//     140, 10, 100, 20;
//     ];
    
    // std::array<std::array<size_t, 4>, 4> mat_conf = {{
    //     {80, 15, 100, 20},
    //     {100, 10, 100, 20},
    //     {120, 10, 100, 20},
    //     {140, 10, 100, 20}

    // }};

    std::array<std::array<size_t, 4>, 3> mat_conf = {{
        {100, 25, 80, 50},
        {100, 25, 80, 50},
        {100, 25, 80, 50},
        // {120, 10, 100, 20},
        // {140, 10, 100, 20}

    }};

    // tmp here end

    double k_Ez_a;
    double k_Ez_b;

    double eta;
    double sigma_max;

    size_t T = 0;

public:
    TFSF();
    void SetParams();
    void CalcNextLayer();
    size_t GetCurrentTimeStep();

    struct Output {
	    // std::array<std::array<double, ny>, nx> Ez;
        // std::array<std::array<double, ny>, nx> Hz;
        std::array<double, nx * ny> Ez;
        std::array<double, nx * ny> Hz;
        std::array<size_t, nx * ny> X;
        std::array<size_t, nx * ny> Y;
        size_t rows;
        size_t cols;
        double maxEz;
        double minEz;
        double maxHz;
        double minHz;
    };

    struct Output GetValues();
};

#endif