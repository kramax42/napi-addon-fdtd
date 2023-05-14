#include "fdtd-2d-upml-tf-sf.h"

TFSF::TFSF() {
    SetParams();
}

size_t TFSF::GetCurrentTimeStep() {
    return T;
}

struct TFSF::Output TFSF::GetValues() {
    struct Output output; 

    // output.Ez = Ez;
    // output.Hz = Hz;  
    output.rows = nx;
    output.cols = ny;

    // size_t flatten_array_size = nx * ny;

    for(size_t i = 0; i < nx; i += 1) {
        for(size_t j = 0; j < ny; j += 1) {
            output.Ez[i*ny + j] = Ez[i][j];
            output.Hz[i*ny + j] = Hz[i][j];
            output.X[i*ny + j] = i;
            output.Y[i*ny + j] = j;
        }
    }

    output.maxEz = *std::max_element(std::begin(output.Ez), std::end(output.Ez));
    output.minEz = *std::min_element(std::begin(output.Ez), std::end(output.Ez));

    output.maxHz = *std::max_element(std::begin(output.Hz), std::end(output.Hz));
    output.minHz = *std::min_element(std::begin(output.Hz), std::end(output.Hz));
	
	return output;
}


void TFSF::CalcNextLayer() {
    // Calculate incident 1D plain waves for TF/SF implementation
    // TE mode (Hz)

    // H
    for(size_t i=0; i < nx_b; i += 1) {
        Ey_1D[i] = k_Ey_a[i]*Ey_1D[i] - k_Ey_b[i]*(Hz_1D[i+1]-Hz_1D[i]); // C
    }
    
    // std::cout << "Ey_1d\n";

    Hz_1D[0] = E0*std::sin(2 * pi * frequency * T * dt);
    
    // H
    for(size_t i=1; i < nx_b; i += 1) {
        Hz_1D[i] = k_Hz_a[i] * Hz_1D[i] - k_Hz_b[i]*(Ey_1D[i]-Ey_1D[i-1]); // C
    }


    //  TM mode (Ez)  
    for(size_t i=0; i < nx_b; i += 1) {
        Hy_1D[i] = k_Hy_a[i]*Hy_1D[i] + k_Hy_b[i]*(Ez_1D[i+1]-Ez_1D[i]); // C
    }


    Ez_1D[0] = E0*std::sin(2*pi*frequency*T*dt);

    for(size_t i=0; i < nx_b-1; i += 1) {
        Fz_1D_r[i] = Fz_1D[i+1]; // C
    }
    
    for(size_t i=1; i < nx_b; i += 1) {
        Fz_1D[i] = k_Fz_a[i]*Fz_1D[i] + k_Fz_b[i]*(Hy_1D[i] - Hy_1D[i-1]); // C
    }

    for(size_t i=1; i < nx_b; i += 1) {
        Ez_1D[i] = k_Ez_a*Ez_1D[i] + k_Ez_b*( Fz_1D[i] - Fz_1D_r[i-1]); // C
    }


    // Calculate Ez (TM mode) and Hz (TE mode) total fields 
    // TE: Wz -> Hz    

    // h
    // C
    for(size_t i=1; i < nx-1; i += 1) {
        for(size_t j=1; j < ny-1; j += 1) {
            Fz_r1[i-1][j-1] = Wz[i][j];
        }
    }

    // h   
        
    for(size_t i=1; i < nx-1; i += 1) {
        for(size_t j=1; j < ny-1; j += 1) {

            
             // C
            Wz[i][j] = k_Fz_1_new[i-1][j-1] * Wz[i][j] + 
                        k_Fz_2_new[i-1][j-1] * ((Ex[i][j]-Ex[i][j-1])/dy -
                        (Ey[i][j]-Ey[i-1][j])/dx );

                        
            // C
            Hz[i][j] = k_Hz_1_new[i-1][j-1]*Hz[i][j] +
                        k_Hz_2_new[i-1][j-1] * ( Wz[i][j]-Fz_r1[i-1][j-1] )/M1[i-1][j-1];
        }
    }
    


    // TM: Fz -> Tz -> Ez    
    // h      
    for(size_t i=1; i < nx-1; i += 1) {
        for(size_t j=1; j < ny-1; j += 1) {
            Fz_r[i-1][j-1] = Fz[i][j];
            Tz_r[i-1][j-1] = Tz[i][j];
        }
    }


    // h
    for(size_t i=1; i < nx-1; i += 1) {
        for(size_t j=1; j < ny-1; j += 1) {
            Fz[i][j] = k_Fz_1_new[i-1][j-1] * Fz[i][j] + k_Fz_2_new[i-1][j-1]*( (Hy[i][j] -
                        Hy[i-1][j])/dx - (Hx[i][j] -
                        Hx[i][j-1])/dy);
        }
    }
    
    // h
    for(size_t i=0; i < nx-2; i += 1) {
        for(size_t j=0; j < ny-2; j += 1) {
            Tz[i+1][j+1] = K_a_new[i][j] * Tz[i+1][j+1] + 
                        K_b_new[i][j]*( Fz[i+1][j+1] - Fz_r[i][j]);
        }
    }
    
    // h
    for(size_t i=1; i < nx-1; i += 1) {
        for(size_t j=1; j < ny-1; j += 1) {
            Ez[i][j] = k_Ez_1_new[i-1][j-1]*Ez[i][j] + 
                        k_Ez_2_new[i-1][j-1]*( Tz[i][j] - Tz_r[i-1][j-1] )/ M0[i-1][j-1];                   
        }
    }

    // Calculate scattered field Ez and Hz in TF/SF
    // TE mode             
    size_t i = nx_a;

    //h   
    for(size_t j=ny_a; j <= ny_b; j += 1) {     
        Hz[i][j] = Hz[i][j] + dt / (mu_0*Material[Index[i][j]][1]*dx) * Ey_1D[0];
    }

    //h
    i = nx_b;
    for(size_t j=ny_a; j <= ny_b; j += 1) {     
        Hz[i][j] = Hz[i][j] -  dt/(mu_0*Material[Index[i][j]][1]*dx)
        *Ey_1D[nx_b-nx_a+1]; // C
    }

    // TM mode 
    // h
    i = nx_a;
    for(size_t j=ny_a; j <= ny_b; j += 1) {                      
        Ez[i][j] = Ez[i][j] - dt / (epsilon_0*Material[Index[i][j]][0]*dx)*Hy_1D[0];
    }

    // h
    i = nx_b;
    for(size_t j=ny_a; j <= ny_b; j += 1) {                      
        Ez[i][j] = Ez[i][j] + dt/(epsilon_0*Material[Index[i][j]][0]*dx)*Hy_1D[nx_b-nx_a+1];
    }

    // Calculate Hx and Ex total fields 
    // TE mode
    // h
    for(size_t i=0; i < nx; i += 1) {
        for(size_t j=0; j < ny-1; j += 1) {
            Gx_r1[i][j] = Mx[i][j];
        }
    }
    
    // h
    // C
    for(size_t i=0; i < nx; i += 1) {
        for(size_t j=0; j < ny-1; j += 1) {
            Mx[i][j] = k_Gx_1_new[i][j] * Mx[i][j] + 
                      k_Gx_2_new[i][j]*( Hz[i][j+1]-Hz[i][j] );
        }
    }

    // h
    // 
    for(size_t i=0; i < nx; i += 1) {
        for(size_t j=0; j < ny-1; j += 1) {
            Ex[i][j] = K_a1[IndexX[i][j]] * (K_b1[IndexX[i][j]] *
                      Ex[i][j] + k_Ex_1_new[i][j]*Mx[i][j]-k_Ex_2_new[i][j]*Gx_r1[i][j]);

            // if(T == 9 && i == 1 && j == 1 ) {
                // std::cout << "a: " << k_Fz_1_new[i-1][j-1] * Wz[i][j] << "\n";
                // std::cout << "b: " << ((Ex[i][j]-Ex[i][j-1])/dy -
                //         (Ey[i][j]-Ey[i-1][j])/dx ) << "\n";
                // std::cout << "ex: " << (K_b1[IndexX[i][j]] *
                    //   Ex[i][j] + k_Ex_1_new[i][j]*Mx[i][j]-k_Ex_2_new[i][j]*Gx_r1[i][j]) << "\n";
            // }                        
        }
    }

    // TM mode
    // h
    for(size_t i=0; i < nx; i += 1) {
        for(size_t j=0; j < ny-1; j += 1) {
            Gx_r[i][j] = Gx[i][j];
        }
    }

    // // h
    for(size_t i=0; i < nx; i += 1) {
        for(size_t j=0; j < ny-1; j += 1) {
            Gx[i][j] = k_Gx_1_new[i][j] * Gx[i][j] - 
                      k_Gx_2_new[i][j] * (Ez[i][j+1] - Ez[i][j]);
        }
    }

    // // h
    for(size_t i=0; i < nx; i += 1) {
        for(size_t j=0; j < ny-1; j += 1) {
            Hx[i][j] = Hx[i][j] + (k_Hx_1_new[i][j] * Gx[i][j] - k_Hx_2_new[i][j]*Gx_r[i][j]) / M2[i][j];
        }
    }


    // // Calculate Hy and Ey total fields 
    // TE mode
    // h
    for(size_t i=0; i < nx-1; i += 1) {
        for(size_t j=0; j < ny; j += 1) {
            Gy_r1[i][j] = My[i][j];
        }
    }
    
    // h
    // C
    for(size_t i=0; i < nx-1; i += 1) {
        for(size_t j=0; j < ny; j += 1) {
            My[i][j] = k_Gy_1_new[i][j]*My[i][j] - 
                      k_Gy_2_new[i][j] * (Hz[i+1][j]-Hz[i][j]);
        }
    }

    
    // h
    // C
    for(size_t i=0; i < nx-1; i += 1) {
        for(size_t j=0; j < ny; j += 1) {
            Ey[i][j] = K_a1[IndexY[i][j]] * (K_b1[IndexY[i][j]] *
                      Ey[i][j] + k_Ey_1_new[i][j] * My[i][j] - k_Ey_2_new[i][j] * Gy_r1[i][j]);
        }
    }



    
    // ???
    // TM mode
    // h
    // C
    for(size_t i=0; i < nx-1; i += 1) {
        for(size_t j=0; j < ny; j += 1) {
            Gy_r[i][j] = Gy[i][j];
        }
    }

    // h
    // C
    for(size_t i=0; i < nx-1; i += 1) {
        for(size_t j=0; j < ny; j += 1) {
            Gy[i][j] = k_Gy_1_new[i][j] * Gy[i][j] + 
                      k_Gy_2_new[i][j] * (Ez[i+1][j]-Ez[i][j]);
        }
    }

    // h
    // C
    for(size_t i=0; i < nx-1; i += 1) {
        for(size_t j=0; j < ny; j += 1) {
            Hy[i][j] = Hy[i][j] + 
                    (k_Hy_1_new[i][j] * Gy[i][j] - k_Hy_2_new[i][j] * Gy_r[i][j]) / M3[i][j];
        }
    }


    // // Calculate scattered field Hx and Ex in TF/SF
    // // TE mode
    // h
    // C
    size_t j = ny_a-1;
    for(size_t i=nx_a; i <= nx_b; i += 1) {
        Ex[i][j] = Ex[i][j] -
            2 * dt / dy * K_a1[Index[i][j]] * Hz_1D[i-nx_a+1];

    }

    // h
    j = ny_b;
    for(size_t i=nx_a; i <= nx_b; i += 1) {
        Ex[i][j] = Ex[i][j] + 
            2 * dt / dy * K_a1[Index[i][j]] * Hz_1D[i-nx_a+1];

        
    }

    // TM mode  
    //h          
    j = ny_a-1;
    for(size_t i=nx_a; i <= nx_b; i += 1) {
        Hx[i][j] = Hx[i][j] + 
            dt / (mu_0 * dy * Material[Index[i][j]][1]) * Ez_1D[i-nx_a+1];
    }

    // h
    j = ny_b;
    for(size_t i=nx_a; i <= nx_b; i += 1) {
        Hx[i][j] = Hx[i][j] -
            dt / (mu_0 * dy * Material[Index[i][j]][1]) * Ez_1D[i-nx_a+1];
    }

    // Calculate scattered field Hy and Ey in TF/SF
    // TE mode
    // h
    i = nx_a-1;
    for(size_t j=ny_a; j <= ny_b; j += 1) {    
        Ey[i][j] = Ey[i][j] + 
        2 * dt / dx * K_a1[Index[i][j]] * Hz_1D[1];
    } // K_a1????

    // h
    i = nx_b;
    for(size_t j=ny_a; j <= ny_b; j += 1) {    
        Ey[i][j] = Ey[i][j] - 
            2 * dt / dx * K_a1[Index[i][j]] * Hz_1D[nx_b-nx_a+1];

      
    }// K_a1????

    // TM mode   
    // h                  
    i = nx_a-1;
    for(size_t j=ny_a; j <= ny_b; j += 1) {    
        Hy[i][j] = Hy[i][j] -
            dt / (mu_0 * dx * Material[Index[i][j]][1]) * Ez_1D[1];
    }


    // h
    // C
    i = nx_b;
    for(size_t j=ny_a; j <= ny_b; j += 1) {  
        Hy[i][j] = Hy[i][j] + 
            dt / (mu_0 * dx * Material[Index[i][j]][1]) * Ez_1D[nx_b-nx_a+1];

        //   if(T == 9  ) {
        //     std::cout << j << ": " << dt / (mu_0 * dx * Material[Index[i][j]][1]) * Ez_1D[nx_b-nx_a+1] << " ";
    
        //     }   
        //     std::cout << "\n";
    }


        if(T == 9) {
        // for(size_t i=0; i < Fz_r.size(); i += 1) {
        //     for(size_t j=0; j < Fz_r[0].size(); j += 1) {
        //         std::cout << Fz_r[i][j] << " "; 
        //     }
        //     std::cout << "\n";
        // }
        // << "Fz_r1[0][j]" << "\n"; 
        std::ofstream myfile;
 myfile.open("Wz.txt");
  if (myfile.is_open())
  {
    //  myfile << "This is a line.\n";
    // myfile << "This is another line.\n";

    for(size_t i=0; i < Hz.size(); i += 1) {
        for(size_t j=0; j < Hz[0].size(); j += 1) {
   
            myfile << Hz[i][j] << " ";
        }
        myfile << '\n';
    }
    
// std::cout << "filfe" << '\n';
    myfile.close();
  }
  else std::cout << "Unable to open file";
    }


    T++;
}


void TFSF::SetParams() {
    
    // Free space FDTD coefficients
    k_Fz_1.fill(1.0);
    k_Fz_2.fill(dt);
    
    k_Ez_1.fill(1.0);
    k_Ez_2.fill(1.0/epsilon_0);
    k_Gx_1.fill(1.0);
    k_Gx_2.fill(dt/dy);
    k_Hx_1.fill(1.0/mu_0);
    k_Hx_2.fill(1.0/mu_0);
    k_Gy_1.fill(1.0);
    k_Gy_2.fill(dt/dx);
    k_Hy_1.fill(1.0/mu_0);
    k_Hy_2.fill(1.0/mu_0);
    k_Hz_1.fill(1.0);
    k_Hz_2.fill(1.0/mu_0);
    k_Ex_1.fill(2.0);
    k_Ex_2.fill(2.0);
    k_Ey_1.fill(2.0);
    k_Ey_2.fill(2.0);
 

    
    k_Ez_a = (2.0*epsilon_0 * Material[0][0] - Material[0][2]*dt)/(2.0*epsilon_0*Material[0][0] + Material[0][2]*dt);
    k_Ez_b =  2.0/(2.0*epsilon_0*Material[0][0] + Material[0][2]*dt);

    
    // // Free space plane wave coefficients
    k_Hy_a.fill(1.0);
    k_Hy_b.fill(dt/(mu_0*Material[0][1]*dx));
    k_Fz_a.fill(1.0);
    k_Fz_b.fill(dt/dx);

    k_Hz_a.fill(1.0);
    k_Hz_b.fill(dt/(mu_0*Material[0][1]*dx));
    k_Ey_a.fill(1.0);
    k_Ey_b.fill(dt/(epsilon_0*Material[0][0]*dx));



    eta = std::sqrt(mu_0*Material[0][1]/epsilon_0/Material[0][0]);

    //  Field transformation coefficients in upml_width areas for TM ans TE modes
    //  Along x-axis
    sigma_max = -(m+1)*std::log(R_err)/(2*eta*upml_width*dx);

    for(size_t i=0; i < upml_width; i += 1) {

        size_t end;

        double sigma_x = sigma_max * std::pow((double)(upml_width-i)/upml_width, m);
        double ka_x = 1+(ka_max-1)*std::pow((double)(upml_width-i)/upml_width, m);

        k_Ez_1[i] = (2*epsilon_0*ka_x-sigma_x*dt)/(2*epsilon_0*ka_x+sigma_x*dt);
        end = k_Ez_1.size() - 1;
        k_Ez_1[end-i] = k_Ez_1[i];

        k_Ez_2[i] = 2.0/(2*epsilon_0*ka_x + sigma_x*dt);
        end = k_Ez_2.size() - 1;
        k_Ez_2[end-i] = k_Ez_2[i];

        k_Hx_1[i] = (2*epsilon_0*ka_x+sigma_x*dt)/(2*epsilon_0*mu_0);
        end = k_Hx_1.size() - 1;
        k_Hx_1[end-i] = k_Hx_1[i];

        k_Hx_2[i] = (2*epsilon_0*ka_x-sigma_x*dt)/(2*epsilon_0*mu_0);
        end = k_Hx_2.size() - 1;
        k_Hx_2[end-i] = k_Hx_2[i];

        k_Hz_1[i] = (2*epsilon_0*ka_x-sigma_x*dt)/(2*epsilon_0*ka_x+sigma_x*dt);
        end = k_Hz_1.size() - 1;
        k_Hz_1[end-i] = k_Hz_1[i];

        k_Hz_2[i] = 2*epsilon_0/(2*epsilon_0*ka_x+sigma_x*dt)/mu_0;
        end = k_Hz_2.size() - 1;
        k_Hz_2[end-i] = k_Hz_2[i];

        k_Ex_1[i] = (2*epsilon_0*ka_x+sigma_x*dt)/epsilon_0;
        end = k_Ex_1.size() - 1;
        k_Ex_1[end-i] = k_Ex_1[i];

        k_Ex_2[i] = (2*epsilon_0*ka_x-sigma_x*dt)/epsilon_0;
        end = k_Ex_2.size() - 1;
        k_Ex_2[end-i] = k_Ex_2[i];


       sigma_x = sigma_max * std::pow((double)(upml_width-i-0.5)/upml_width, m);
       ka_x = 1+(ka_max-1)*std::pow((double)(upml_width-i-0.5)/upml_width, m);

       k_Gy_1[i] = (2*epsilon_0*ka_x-sigma_x*dt)/(2*epsilon_0*ka_x+sigma_x*dt);
       end = k_Gy_1.size() - 1;
       k_Gy_1[end-i] = k_Gy_1[i];

       k_Gy_2[i] = 2*epsilon_0*dt/(2*epsilon_0*ka_x+sigma_x*dt)/dx;
       end = k_Gy_2.size() - 1;
       k_Gy_2[end-i] = k_Gy_2[i];
    }




// % Along y-axis

    sigma_max = -(m+1)*std::log(R_err)/(2*eta*upml_width*dy);

    for(size_t i=0; i < upml_width; i += 1) {

        size_t end;

        double sigma_y = sigma_max * std::pow((double)(upml_width-i)/upml_width, m);
        double ka_y = 1+(ka_max-1)*std::pow((double)(upml_width-i)/upml_width, m);

        k_Fz_1[i] = (2*epsilon_0*ka_y-sigma_y*dt)/(2*epsilon_0*ka_y+sigma_y*dt);
        end = k_Fz_1.size() - 1;
        k_Fz_1[end-i] = k_Fz_1[i];
    
        k_Fz_2[i] = 2*epsilon_0*dt/(2*epsilon_0*ka_y+sigma_y*dt);
        end = k_Fz_2.size() - 1;
        k_Fz_2[end-i] = k_Fz_2[i];

        k_Hy_1[i] = (2*epsilon_0*ka_y+sigma_y*dt)/(2*epsilon_0*mu_0);
        end = k_Hy_1.size() - 1;
        k_Hy_1[end-i] = k_Hy_1[i];

        k_Hy_2[i] = (2*epsilon_0*ka_y-sigma_y*dt)/(2*epsilon_0*mu_0);
        end = k_Hy_2.size() - 1;
        k_Hy_2[end-i] = k_Hy_2[i];

        k_Ey_1[i] = (2*epsilon_0*ka_y+sigma_y*dt)/epsilon_0;
        end = k_Ey_1.size() - 1;
        k_Ey_1[end-i] = k_Ey_1[i];

        k_Ey_2[i] = (2*epsilon_0*ka_y-sigma_y*dt)/epsilon_0;
        end = k_Ey_2.size() - 1;
        k_Ey_2[end-i] = k_Ey_2[i];


        sigma_y = sigma_max * std::pow((double)(upml_width-i-0.5)/upml_width, m);
        ka_y = 1+(ka_max-1)*std::pow((double)(upml_width-i-0.5)/upml_width, m);

         k_Gx_1[i] = (2*epsilon_0*ka_y-sigma_y*dt)/(2*epsilon_0*ka_y+sigma_y*dt);
        end = k_Gx_1.size() - 1;
         k_Gx_1[end-i] = k_Gx_1[i];

         k_Gx_2[i] = 2*epsilon_0*dt/(2*epsilon_0*ka_y+sigma_y*dt)/dy;
         end = k_Gx_2.size() - 1;
         k_Gx_2[end-i] = k_Gx_2[i];    
    }


//    for (auto & row : arr) {
//         for (auto & col : row) {
//             std::cout << col << ' ';
//         }
//         std::cout << std::endl;
//     }


// % Vectorize transformation coefficients
size_t rows;
size_t cols;

//nx-2 ny-2
rows = k_Fz_1_new.size();
cols = k_Fz_1_new[0].size();
std::fill_n(&k_Fz_1_new[0][0], rows * cols, k_Fz_1[upml_width]); // C
std::fill_n(&k_Fz_2_new[0][0], rows * cols, k_Fz_2[upml_width]);

// left right border
for(size_t i=0; i < rows; i += 1) {
    for(size_t j=0; j < upml_width-1; j += 1) {
        // left
        k_Fz_1_new[i][j] = k_Fz_1[j+1]; // C
        k_Fz_2_new[i][j] = k_Fz_2[j+1];

        // right
        k_Fz_1_new[i][cols-j-1] = k_Fz_1[j+1]; // C
        k_Fz_2_new[i][cols-j-1] = k_Fz_2[j+1];
    }
}



rows = k_Ez_1_new.size();
cols = k_Ez_1_new[0].size();
std::fill_n(&k_Ez_1_new[0][0], rows * cols, k_Ez_1[upml_width]);
std::fill_n(&k_Ez_2_new[0][0], rows * cols, k_Ez_2[upml_width]);

// up down border
for(size_t i=0; i < upml_width-1; i += 1) {
    for(size_t j=0; j < cols; j += 1) {

        // up
        k_Ez_1_new[i][j] = k_Ez_1[i+1];
        k_Ez_2_new[i][j] = k_Ez_2[i+1];

        //down
        k_Ez_1_new[rows-i-1][j] = k_Ez_1[i+1];
        k_Ez_2_new[rows-i-1][j] = k_Ez_2[i+1];
    }
}


rows = k_Gx_1_new.size();
cols = k_Gx_1_new[0].size();
std::fill_n(&k_Gx_1_new[0][0], rows * cols, k_Gx_1[upml_width]);
std::fill_n(&k_Gx_2_new[0][0], rows * cols, k_Gx_2[upml_width]);

// left right border
for(size_t i=0; i < rows; i += 1) {
    for(size_t j=0; j < upml_width; j += 1) {
        // left
        k_Gx_1_new[i][j] = k_Gx_1[j];
        k_Gx_2_new[i][j] = k_Gx_2[j];

        // right
        k_Gx_1_new[i][cols-j-1] = k_Gx_1[j];
        k_Gx_2_new[i][cols-j-1] = k_Gx_2[j];
    }
}


rows = k_Hx_1_new.size();
cols = k_Hx_1_new[0].size();
std::fill_n(&k_Hx_1_new[0][0], rows * cols, k_Hx_1[upml_width]);
std::fill_n(&k_Hx_2_new[0][0], rows * cols, k_Hx_2[upml_width]);
std::fill_n(&k_Ex_1_new[0][0], rows * cols, k_Ex_1[upml_width]);
std::fill_n(&k_Ex_2_new[0][0], rows * cols, k_Ex_2[upml_width]);

// up down border
for(size_t i=0; i < upml_width; i += 1) {
    for(size_t j=0; j < cols; j += 1) {
        // up
        k_Hx_1_new[i][j] = k_Hx_1[i];
        k_Hx_2_new[i][j] = k_Hx_2[i];

        k_Ex_1_new[i][j] = k_Ex_1[i];
        k_Ex_2_new[i][j] = k_Ex_2[i];

        //down
        k_Hx_1_new[rows-i-1][j] = k_Hx_1[i];
        k_Hx_2_new[rows-i-1][j] = k_Hx_2[i];
        k_Ex_1_new[rows-i-1][j] = k_Ex_1[i];
        k_Ex_2_new[rows-i-1][j] = k_Ex_2[i];
    }
}




rows = k_Gy_1_new.size();
cols = k_Gy_1_new[0].size();
std::fill_n(&k_Gy_1_new[0][0], rows * cols, k_Gy_1[upml_width]);
std::fill_n(&k_Gy_2_new[0][0], rows * cols, k_Gy_2[upml_width]);

// up down border
for(size_t i=0; i < upml_width; i += 1) {
    for(size_t j=0; j < cols; j += 1) {

        // up
        k_Gy_1_new[i][j] = k_Gy_1[i];
        k_Gy_2_new[i][j] = k_Gy_2[i];

        //down
        k_Gy_1_new[rows-i-1][j] = k_Gy_1[i];
        k_Gy_2_new[rows-i-1][j] = k_Gy_2[i];
    }
}



rows = k_Hy_1_new.size();
cols = k_Hy_1_new[0].size();
std::fill_n(&k_Hy_1_new[0][0], rows * cols, k_Hy_1[upml_width]);
std::fill_n(&k_Hy_2_new[0][0], rows * cols, k_Hy_2[upml_width]);
std::fill_n(&k_Ey_1_new[0][0], rows * cols, k_Ey_1[upml_width]);
std::fill_n(&k_Ey_2_new[0][0], rows * cols, k_Ey_2[upml_width]);

// left right border
for(size_t i=0; i < rows; i += 1) {
    for(size_t j=0; j < upml_width; j += 1) {
        // left
        k_Hy_1_new[i][j] = k_Hy_1[j];
        k_Hy_2_new[i][j] = k_Hy_2[j];

        k_Ey_1_new[i][j] = k_Ey_1[j];
        k_Ey_2_new[i][j] = k_Ey_2[j];

        // right
        k_Hy_1_new[i][cols-j-1] = k_Hy_1[j];
        k_Hy_2_new[i][cols-j-1] = k_Hy_2[j];

        k_Ey_1_new[i][cols-j-1] = k_Ey_1[j];
        k_Ey_2_new[i][cols-j-1] = k_Ey_2[j];

    }
}


rows = k_Hz_1_new.size();
cols = k_Hz_1_new[0].size();
std::fill_n(&k_Hz_1_new[0][0], rows * cols, k_Hz_1[upml_width]);
std::fill_n(&k_Hz_2_new[0][0], rows * cols, k_Hz_2[upml_width]);

// up down border
for(size_t i=0; i < upml_width-1; i += 1) {
    for(size_t j=0; j < cols; j += 1) {

        // up
        k_Hz_1_new[i][j] = k_Hz_1[i+1];
        k_Hz_2_new[i][j] = k_Hz_2[i+1];

        //down
        k_Hz_1_new[rows-i-1][j] = k_Hz_1[i+1];
        k_Hz_2_new[rows-i-1][j] = k_Hz_2[i+1];
    }
}



std::fill_n(&Index[0][0], nx * ny, 0.0); // C
std::fill_n(&IndexX[0][0], nx * ny, 0.0); // C
std::fill_n(&IndexY[0][0], nx * ny, 0.0); // C
// % User material config coordinates in cells
// % [i, cells_i, j, cells_j]
// mat_conf = [
//     [160, 9, 240, 20],
//     [180, 15, 240, 20],
//     [200, 10, 240, 20],
//     [220, 10, 240, 20],
//     [240, 10, 240, 20],
//     [260, 10, 240, 20],
//     [280, 10, 240, 20],
//     [300, 10, 240, 20],
//     [320, 10, 240, 20]
//     ];

// mat_count = size(mat_conf, 1);
// % /////////////////////////////////////////////////////


// Fill geometry matrix for 1/2 cylinder (no vectorization here), but for TE
// mode we need three different Index for Ex, Ey and Hz due to leapfrog
for(size_t i=0; i < nx; i += 1) {
    for(size_t j=0; j < ny; j += 1) {
        for(size_t mat_idx=0; mat_idx < mat_count; mat_idx += 1) {
        
            //conf = mat_conf(mat_idx);
            size_t x = mat_conf[mat_idx][0];
            size_t x_cells = mat_conf[mat_idx][1];
            size_t y = mat_conf[mat_idx][2];
            size_t y_cells = mat_conf[mat_idx][3];

            if ((i > x) &&  (i < (x + x_cells))
                && (j > y) &&  (j < (y + y_cells))) {
                Index[i][j] = 1;
            }
            if ((i+0.5 > x) &&  (i+0.5 < (x + x_cells))
                && (j > y) &&  (j < (y + y_cells))) {
                IndexX[i][j] = 1;
            }
            if ((i > x) &&  (i < (x + x_cells))
                && (j+0.5 > y) &&  (j+0.5 < (y + y_cells))) {
                IndexY[i][j] = 1;
            }
        }
    }
}


for(size_t i=0; i < upml_width; i += 1) {

        size_t end;

        double sigma_x = sigma_max * std::pow((double)(upml_width-i)/upml_width, m);
        double ka_x = 1+(ka_max-1)*std::pow((double)(upml_width-i)/upml_width, m);

        end = k_Fz_a.size() - 1;
        k_Fz_a[end-i] = (2*epsilon_0*ka_x-sigma_x*dt)/(2*epsilon_0*ka_x+sigma_x*dt);

        end = k_Fz_b.size() - 1;
        k_Fz_b[end-i] = 2*epsilon_0*dt/(2*epsilon_0*ka_x+sigma_x*dt)/dx;

        end = k_Hz_a.size() - 1;
        k_Hz_a[end-i] = (2*epsilon_0*ka_x - sigma_x*dt)/(2*epsilon_0*ka_x + sigma_x*dt);

        end = k_Hz_b.size() - 1;
        k_Hz_b[end-i] = 2*epsilon_0*dt/(2*epsilon_0*ka_x+sigma_x*dt)/(mu_0*Material[0][1]*dx);

        sigma_x = sigma_max * std::pow((double)(upml_width-i-0.5)/upml_width, m);
        ka_x = 1+(ka_max-1)*std::pow((double)(upml_width-i-0.5)/upml_width, m);

        end = k_Hy_a.size() - 1;
        k_Hy_a[end-i] = (2*epsilon_0*ka_x - sigma_x*dt)/(2*epsilon_0*ka_x + sigma_x*dt);

        end = k_Hy_b.size() - 1;
        k_Hy_b[end-i] = 2*epsilon_0*dt/(2*epsilon_0*ka_x + sigma_x*dt)/(mu_0*Material[0][1]*dx);

        end = k_Ey_a.size() - 1;
        k_Ey_a[end-i] = (2*epsilon_0*ka_x - sigma_x*dt)/(2*epsilon_0*ka_x + sigma_x*dt);

        end = k_Ey_b.size() - 1;
        k_Ey_b[end-i] = 2*dt/(2*epsilon_0*ka_x + sigma_x*dt)/(Material[0][0]*dx);
}


           
// end

// %% Another transformation coefficients

for(size_t i=0; i < number_of_materials; i += 1) {
    K_a[i] = (2*epsilon_0*Material[i][0] - 
        Material[i][2]*dt) / (2*epsilon_0*Material[i][0] +
        Material[i][2]*dt);

    K_b[i] = 2*epsilon_0*Material[i][0] /
    (2*epsilon_0*Material[i][0] + Material[i][2]*dt);

    K_a1[i] = 1/(2*epsilon_0*Material[i][0] + Material[i][2]*dt);

    K_b1[i] = 2*epsilon_0*Material[i][0] - Material[i][2]*dt;
}


rows = K_a_new.size();
cols = K_a_new[0].size();
std::fill_n(&K_a_new[0][0], rows * cols, 1);
std::fill_n(&K_b_new[0][0], rows * cols, 1);

 
// % Replace [1,2] -> [1,-1]
for(size_t i=0; i < nx-2; i += 1) {
    for(size_t j=0; j < ny-2; j += 1) {
        
      if(K_a_new[i][j] == 2) {
        K_a_new[i][j] = K_a[1];
      }

      if(K_b_new[i][j] == 2) {
        K_b_new[i][j] = K_b[1];
      }
    }
}


for(size_t i=0; i < nx-2; i += 1) {
    for(size_t j=0; j < ny-2; j += 1) {
        M0[i][j] = Material[Index[i+1][j+1]][0];
        M1[i][j] = Material[Index[i+1][j+1]][1];
    }
}

for(size_t i=0; i < nx; i += 1) {
    for(size_t j=0; j < ny-1; j += 1) {
        M2[i][j] = Material[Index[i][j]][1];
    }
}

for(size_t i=0; i < nx-1; i += 1) {
    for(size_t j=0; j < ny; j += 1) {
        M3[i][j] = Material[Index[i][j]][1];
    }
}



for(size_t j=0; j < K_b.size(); j += 1) {
        // std::cout << K_b[j] << " "; 
 }
//  std::cout << "\n" << "K_b[j]" << "\n"; 

std::ofstream myfile;
 myfile.open("k_Fz_1_new3.txt");
  if (myfile.is_open())
  {
     myfile << "This is a line.\n";
    myfile << "This is another line.\n";

    // for(size_t i=0; i < nx; i += 1) {
    //     for(size_t j=0; j < ny; j += 1) {
   
    //         std::cout << Index[i][j] << " ";
    //     }
    //     std::cout << '\n';
    // }
    
// std::cout << "filfe" << '\n';
    myfile.close();
  }
  else std::cout << "Unable to open file";
for(size_t i=0; i < M3.size(); i += 1) {
        for(size_t j=0; j < M3[0].size(); j += 1) {
   
            // std::cout << M3[i][j] << " ";
        }
        // std::cout << '\n';
}



    // std::cout << "ka_max" << ka_max << std::endl;
    }
