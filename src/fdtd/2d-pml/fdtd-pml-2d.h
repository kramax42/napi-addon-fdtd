// 2D FDTD with PML.

// Based on code from book:
// Electromagnetic simulation using the fdtd method with python.
// Chapter 3.

#ifndef FDTD_PML_2D
#define FDTD_PML_2D

#include <math.h>
#include <napi.h>

#include <algorithm>  // std::fill_n
#include <vector>
#include <array>


#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#define  MEDIACONSTANT    (2)
#define  NUMBEROFITERATIONCONSTANT    (1000)          // was 300
#define  BIGLINESIZE    (8192) 


class FdtdPml2D {


///////////////////////////////////////////////////////
       // char  ch;
    double  eps[MEDIACONSTANT] = {1.0, 1.0};       // index=0 is for vacuum, index=1 is for the metallic cylinder
    double  sig[MEDIACONSTANT] = {0.0, 1.0e+7};
    double  mur[MEDIACONSTANT] = {1.0, 1.0};
    double  sim[MEDIACONSTANT] = {0.0, 0.0};
    double  ca[MEDIACONSTANT];
    double  cb[MEDIACONSTANT];
    double  da[MEDIACONSTANT];
    double  db[MEDIACONSTANT];
    double  source[NUMBEROFITERATIONCONSTANT];
    int  n = 0;
    int i,j, icenter, jcenter;
    double  eaf,haf, diam,rad,dist2;
    double  temporaryi,temporaryj,temporary;
    double  cc,muz,epsz,pi,freq,lambda,omega,dx,dt;
    double  rmax, orderbc;
    static const size_t ie = 220;   //number of grid cells in x-direction
    static const size_t je = 220; //number of grid cells in y-direction
    int ib, jb, is, js, nmax, iebc, jebc, ibbc, jbbc, iefbc, jefbc, ibfbc, jbfbc;
    int  media;
    double  rtau, tau, delay;
    double  **ex;        //fields in main grid 
    double  **ey;      
    double  **hz;      
    double  **exbcf;     //fields in front PML region
    double  **eybcf;   
    double  **hzxbcf;  
    double  **hzybcf;  
    double  **exbcb;     //fields in back PML region
    double  **eybcb;   
    double  **hzxbcb;  
    double  **hzybcb;  
    double  **exbcl;     //fields in left PML region
    double  **eybcl;   
    double  **hzxbcl;  
    double  **hzybcl;  
    double  **exbcr;     //fields in right PML region
    double  **eybcr;   
    double  **hzxbcr;  
    double  **hzybcr;  
    double  **caex;     // main grid coefficents
    double  **cbex;
    double  **caey;
    double  **cbey;
    double  **dahz;
    double  **dbhz;
    double  **caexbcf;   // pml coefficients
    double  **cbexbcf;
    double  **caexbcb;
    double  **cbexbcb;
    double  **caexbcl;
    double  **cbexbcl;
    double  **caexbcr;
    double  **cbexbcr;
    double  **caeybcf;
    double  **cbeybcf;
    double  **caeybcb;
    double  **cbeybcb;
    double  **caeybcl;
    double  **cbeybcl;
    double  **caeybcr;
    double  **cbeybcr;
    double  **dahzxbcf;
    double  **dbhzxbcf;
    double  **dahzxbcb;
    double  **dbhzxbcb;
    double  **dahzxbcl;
    double  **dbhzxbcl;
    double  **dahzxbcr;
    double  **dbhzxbcr;
    double  **dahzybcf;
    double  **dbhzybcf;
    double  **dahzybcb;
    double  **dbhzybcb;
    double  **dahzybcl;
    double  **dbhzybcl;
    double  **dahzybcr;
    double  **dbhzybcr;
    double  delbc,sigmam,bcfactor;
    double  ca1,cb1,y1,y2, sigmay,sigmays,da1,db1;
    double  x1,x2,sigmax,sigmaxs;
    double  minimumValue, maximumValue;
    char  filename[BIGLINESIZE];
    FILE *filePointer;
    double  scaleValue;
    int  iValue,plottingInterval,centery,centerx;     


    int src_type;


   public:

    FdtdPml2D() {}

     size_t GetCurrentTimeStep() {
        return n;
     }
    
    struct Output {
	    // std::array<std::array<double, ny>, nx> Ez;
        // std::array<std::array<double, ny>, nx> Hz;
        // std::array<double, nx * ny> Ez;
        std::array<double, ie * je> Hz;
        std::array<double, ie * je> field;
        std::array<size_t, ie * je> X;
        std::array<size_t, ie * je> Y;
        size_t rows;
        size_t cols;
        double max;
        double min;
        double maxHz;
        double minHz;
    };

    //dataReturnType: 1 - Hz, 2-Ex, 3-Ey
    struct Output GetValues(int dataReturnType) {
        struct Output output; 

        // output.Ez = Ez;
        // output.Hz = Hz;  
        output.rows = ie;
        output.cols = je;

        // size_t flatten_array_size = ie * nyje
        double min = -0.000001;
        double max = 0.000001;
        for(size_t i = 0; i < ie; i += 1) {
            for(size_t j = 0; j < je; j += 1) {
                output.Hz[i*je + j] = hz[i][j];
                if(hz[i][j] > max) {
                    max = hz[i][j];
                }
                if(hz[i][j] < min) {
                    min = hz[i][j];
                }
                output.X[i*je + j] = i;
                output.Y[i*je + j] = j;


                double val;
                if(dataReturnType == 1) {
                    val =  hz[i][j];
                } else if(dataReturnType == 2) {
                    val =  ex[i][j];
                } else if(dataReturnType == 3) {
                    val =  ey[i][j];
                }

                output.field[i*je + j] = val;
                if(val > max) {
                    max = val;
                }
                if(val < min) {
                    min = val;
                }
                
            }
        }

        // output.maxEz = *std::max_element(std::begin(output.Ez), std::end(output.Ez));
        // output.minEz = *std::min_element(std::begin(output.Ez), std::end(output.Ez));

        // output.maxHz = *std::max_element(std::begin(output.Hz), std::end(output.Hz));
        // output.minHz = *std::min_element(std::begin(output.Hz), std::end(output.Hz));
        output.maxHz = max;
        output.minHz = min;
        output.max = max;
        output.min = min;
        
        return output;
    }



    void CalcNextLayer() {
        

        //***********************************************************************
        //     Update electric fields (EX and EY) in main grid
        //***********************************************************************

        // Note the 4 edges, ey left (i=0), ey right (i=ie), ex bottom (j=0) and ex top (j=je) are evaluated in the pml section

        for (i = 0; i < ie; i++) { 
            for (j = 1; j < je; j++) {  // dont do ex at j=0 or j=je, it will be done in the PML section
                ex[i][j] = caex[i][j] * ex[i][j] + cbex[i][j] * ( hz[i][j] - hz[i][j-1] );
            } /* jForLoop */     
        } /* iForLoop */     

        for (i = 1; i < ie; i++) {        // dont do ey at i=0 or i=ie,  it will be done in the PML section 
            for (j = 0; j < je; j++) {   
                ey[i][j] = caey[i][j] * ey[i][j] + cbey[i][j] * ( hz[i-1][j] - hz[i][j] );
            } /* jForLoop */     
        } /* iForLoop */     
          
         
        //***********************************************************************
        //     Update EX in PML regions
        //***********************************************************************


        //     FRONT
         
        for (i = 0; i < iefbc; i++) { 
            for (j = 1; j < jebc; j++) {   // don't Evaluate exbcf at j=0, as it is the PEC
                // note: sign change in the second term from main grid!! ... (due to the Exponential time stepping algorithm?)
                exbcf[i][j] = caexbcf[i][j] * exbcf[i][j] - cbexbcf[i][j] * (hzxbcf[i][j-1] + hzybcf[i][j-1] - hzxbcf[i][j] - hzybcf[i][j]);
            } /* jForLoop */     
        } /* iForLoop */     

        for (j=0, i = 0; i < ie; i++) {    // fill in the edge for ex at j=0 main grid  (ties in the pml with the main grid)
            ex[i][j] = caex[i][j] * ex[i][j] - cbex[i][j] * (hzxbcf[iebc+i][jebc-1] + hzybcf[iebc+i][jebc-1] - hz[i][j]);
        } /* iForLoop */     


        //     BACK
         
        for (i = 0; i < iefbc; i++) { 
            for (j = 1; j < jebc; j++) {   // don't Evaluate exbcb at j=jebc, as it is the PEC, also dont eval at j=0 as this point is the same as j=je on the main grid
                exbcb[i][j] = caexbcb[i][j] * exbcb[i][j] - cbexbcb[i][j] * (hzxbcb[i][j-1] + hzybcb[i][j-1] - hzxbcb[i][j] - hzybcb[i][j]);
            } /* jForLoop */     
        } /* iForLoop */     

        for (j=je, i = 0; i < ie; i++) {    // fill in the edge for ex at j=je main grid  (ties in the pml with the main grid)
            ex[i][j] = caex[i][j] * ex[i][j] - cbex[i][j] * (hz[i][j-1] - hzxbcb[iebc+i][0] - hzybcb[iebc+i][0]);   
        } /* iForLoop */     


        //     LEFT

        for (i = 0; i < iebc; i++) {  
            // don't Evaluate exbcl at j=0, j=0 is a special case, it needs data from the "front grid" see below
            // likewise, don't Evaluate exbcl at j=je, j=je is a special case, it needs data from the "back grid" see below
            for (j = 1; j < je; j++) {
                exbcl[i][j] = caexbcl[i][j] * exbcl[i][j] - cbexbcl[i][j] * (hzxbcl[i][j-1] + hzybcl[i][j-1] - hzxbcl[i][j] - hzybcl[i][j]);
            } /* jForLoop */     
        } /* iForLoop */     
         
        for (j=0, i = 0; i < iebc; i++) { // exbcl at j=0 case, uses data from the "front grid"
            exbcl[i][j] = caexbcl[i][j] * exbcl[i][j] - cbexbcl[i][j] * (hzxbcf[i][jebc-1] + hzybcf[i][jebc-1] - hzxbcl[i][j] - hzybcl[i][j]);
        } /* iForLoop */     

        for (j=je, i = 0; i < iebc; i++) { // exbcl at j=je case, uses data from the "back grid"
            exbcl[i][j] = caexbcl[i][j] * exbcl[i][j] - cbexbcl[i][j] * (hzxbcl[i][j-1] + hzybcl[i][j-1] - hzxbcb[i][0] - hzybcb[i][0]);
        } /* iForLoop */     


        //     RIGHT

        for (i = 0; i < iebc; i++) {  
            // don't Evaluate exbcr at j=0, j=0 is a special case, it needs data from the "front grid" see below
            // likewise, don't Evaluate exbcr at j=je, j=je is a special case, it needs data from the "back grid" see below
            for (j = 1; j < je; j++) {
                exbcr[i][j] = caexbcr[i][j] * exbcr[i][j] - cbexbcr[i][j] * (hzxbcr[i][j-1] + hzybcr[i][j-1] - hzxbcr[i][j] - hzybcr[i][j]);
            } /* jForLoop */     
        } /* iForLoop */     
         
        for (j=0, i = 0; i < iebc; i++) { // exbcr at j=0 case, uses data from the "front grid" (on the right side)
            exbcr[i][j] = caexbcr[i][j] * exbcr[i][j] - cbexbcr[i][j] * (hzxbcf[iebc+ie + i][jebc-1] + hzybcf[iebc+ie + i][jebc-1] -  hzxbcr[i][j] - hzybcr[i][j]);
        } /* iForLoop */     

        for (j=je, i = 0; i < iebc; i++) { // exbcr at j=je case, uses data from the "back grid"
            exbcr[i][j] = caexbcr[i][j] * exbcr[i][j] - cbexbcr[i][j] * (hzxbcr[i][j-1] + hzybcr[i][j-1] - hzxbcb[iebc+ie + i][0] - hzybcb[iebc+ie + i][0]);
        } /* iForLoop */     


        //***********************************************************************
        //     Update EY in PML regions
        //***********************************************************************
         
        //     FRONT

        for (i = 1; i < iefbc; i++) {     // don't Evaluate eybcf at i=0 or iefbc, as it is the PEC 
            for (j = 0; j < jebc; j++) {  
                // note: sign change in the second term from main grid!!
                eybcf[i][j] = caeybcf[i][j] * eybcf[i][j] - cbeybcf[i][j] * (hzxbcf[i][j] + hzybcf[i][j] - hzxbcf[i-1][j] - hzybcf[i-1][j]);
            } /* jForLoop */     
        } /* iForLoop */     
         
         
        //     BACK
         
        for (i = 1; i < iefbc; i++) {     // don't Evaluate eybcb at i=0 or iefbc, as it is the PEC 
            for (j = 0; j < jebc; j++) {  
                eybcb[i][j] = caeybcb[i][j] * eybcb[i][j] - cbeybcb[i][j] * (hzxbcb[i][j] + hzybcb[i][j] - hzxbcb[i-1][j] - hzybcb[i-1][j]);
            } /* jForLoop */     
        } /* iForLoop */     

                        
        //     LEFT

        for (i = 1; i < iebc; i++) {    // don't Evaluate eybcb at i=0, as it is the PEC  
            for (j = 0; j < je; j++) {
                eybcl[i][j] = caeybcl[i][j] * eybcl[i][j] - cbeybcl[i][j] * (hzxbcl[i][j] + hzybcl[i][j] - hzxbcl[i-1][j] - hzybcl[i-1][j] );
            } /* jForLoop */     
        } /* iForLoop */     
         
        for (i=0, j = 0; j < je; j++) {  // fill in the edge for ey at i=0 main grid  (ties in the pml with the main grid) 
            ey[i][j] = caey[i][j] * ey[i][j] - cbey[i][j] * (hz[i][j] - hzxbcl[iebc-1][j] - hzybcl[iebc-1][j]);
        } /* jForLoop */     

                        
        //     RIGHT

        for (i = 1; i < iebc; i++) {  // don't Evaluate eybcr at i=iebc, as it is the PEC, also dont eval at i=0 as this point is the same as i=ie on the main grid
            for (j = 0; j < je; j++) {
                eybcr[i][j] = caeybcr[i][j] * eybcr[i][j] - cbeybcr[i][j] * (hzxbcr[i][j] + hzybcr[i][j] - hzxbcr[i-1][j] - hzybcr[i-1][j] );
            } /* jForLoop */     
        } /* iForLoop */     
                         
        for (i=ie, j = 0; j < je; j++) {  // fill in the edge for ey at i=ie main grid  (ties in the pml with the main grid) 
            ey[i][j] = caey[i][j] * ey[i][j] - cbey[i][j] * (hzxbcr[0][j] + hzybcr[0][j] - hz[i-1][j] );
        } /* jForLoop */     
                              

         
        //***********************************************************************
        //     Update magnetic fields (HZ) in main grid
        //***********************************************************************
         

        for (i = 0; i < ie; i++) { 
            for (j = 0; j < je; j++) {
                hz[i][j] = dahz[i][j] * hz[i][j] + dbhz[i][j] * ( ex[i][j+1] - ex[i][j] + ey[i][j] - ey[i+1][j] );
            } /* jForLoop */     
        } /* iForLoop */     


        //  rtau = 160.0e-12;
    // tau = rtau / dt;
    // delay = 3 * tau;
    // for (i = 0; i < nmax; i++) {
    //     source[i] = 0.0;
    // } /* iForLoop */     
    
    // int nn;
    // for (nn = 0; nn < (int  )(7.0 * tau); nn++) {
    //     temporary = (double  )nn - delay;
    //     source[nn] = sin( omega * (temporary) * dt) * exp(-( (temporary * temporary)/(tau * tau) ) );
    // } /* forLoop */  
         

        const double t0 = 20;
        // Beam width.
        const double tau = 30;
        double src = 1*std::sin(2*3.14*freq*dt*n) * std::exp(-1.0 * std::pow((t0 - n) / tau, 2));
        // double source = 20.0 * ((n - t0) / tau) * std::exp(-1.0 * std::pow((n - t0) / tau, 2));
        
        // double src = -20.0 * ((n - t0) / tau) * std::exp(-1.0 * std::pow((n - t0) / tau, 2));

        // std::sin(2*3.14*freq*dt*n)

        // hz[is][js] = src;
        // hz[is][js] = source;
        // hz[is][js] = src;

        if(src_type == 1) {
            hz[is][js] += 10*std::sin(2*3.14*freq*dt*n);// sin
        } else {
            hz[is][js] = source[n]; //gaussin
        }
        // 
        
                       
         
        //***********************************************************************
        //     Update HZX in PML regions
        //***********************************************************************
         
        //     FRONT
         
        for (i = 0; i < iefbc; i++) { 
            for (j = 0; j < jebc; j++) {
                // note: sign change in the second term from main grid!!
                hzxbcf[i][j] = dahzxbcf[i][j] * hzxbcf[i][j] - dbhzxbcf[i][j] * (eybcf[i+1][j] - eybcf[i][j] );
            } /* jForLoop */     
        } /* iForLoop */     
        
        //     BACK
         
        for (i = 0; i < iefbc; i++) { 
            for (j = 0; j < jebc; j++) {
                hzxbcb[i][j] = dahzxbcb[i][j] * hzxbcb[i][j] - dbhzxbcb[i][j] * (eybcb[i+1][j] - eybcb[i][j] );
            } /* jForLoop */     
        } /* iForLoop */     

        //     LEFT

        for (i = 0; i < (iebc-1); i++) {   // don't evaluate hzxbcl at i=iebc-1 as it needs ey from main grid, see below
            for (j = 0; j < je; j++) {
                hzxbcl[i][j] = dahzxbcl[i][j] * hzxbcl[i][j] - dbhzxbcl[i][j] * (eybcl[i+1][j] - eybcl[i][j]);
            } /* jForLoop */     
        } /* iForLoop */     

        for (i=(iebc-1), j = 0; j < je; j++) {   // fix-up hzxbcl at i=iebc-1
            hzxbcl[i][j] = dahzxbcl[i][j] * hzxbcl[i][j] - dbhzxbcl[i][j] * (ey[0][j] - eybcl[i][j] );
        } /* jForLoop */     

        //     RIGHT

        for (i = 1; i < iebc; i++) {   // don't evaluate hzxbcl at i=0 as it needs ey from main grid, see below
            for (j = 0; j < je; j++) {
                hzxbcr[i][j] = dahzxbcr[i][j] * hzxbcr[i][j] - dbhzxbcr[i][j] * (eybcr[i+1][j] - eybcr[i][j] );
            } /* jForLoop */     
        } /* iForLoop */     

        for (i=0, j = 0; j < je; j++) {   // fix-up hzxbcl at i=0
            hzxbcr[i][j] = dahzxbcr[i][j] * hzxbcr[i][j] - dbhzxbcr[i][j] * (eybcr[i+1][j] - ey[ie][j] );
        } /* jForLoop */     
                       
         
         
        //***********************************************************************
        //     Update HZY in PML regions
        //***********************************************************************

        //     FRONT
         
        for (i = 0; i < iefbc; i++) { 
            for (j = 0; j < (jebc-1); j++) {  // don't evaluate hzxbcf at j=jebc-1 as it needs data from main,left,right grids, see below 
                // note: sign change in the second term from main grid!!
                hzybcf[i][j] = dahzybcf[i][j] * hzybcf[i][j] - dbhzybcf[i][j] * (exbcf[i][j] - exbcf[i][j+1] );
            } /* jForLoop */     
        } /* iForLoop */     

        for (j = (jebc-1), i = 0; i < iebc; i++) {  // fix-up hzybcf at j=jebc-1, with left grid
            hzybcf[i][j] = dahzybcf[i][j] * hzybcf[i][j] - dbhzybcf[i][j] * (exbcf[i][j] - exbcl[i][0] );
        } /* iForLoop */     
         
        for (j = (jebc-1), i = 0; i < ie; i++) {  // fix-up hzybcf at j=jebc-1, with main grid
            hzybcf[iebc+i][j] = dahzybcf[iebc+i][j] * hzybcf[iebc+i][j] - dbhzybcf[iebc+i][j] * (exbcf[iebc+i][j] - ex[i][0]);
        } /* iForLoop */     
         
        for (j = (jebc-1), i = 0; i < iebc; i++) {  // fix-up hzybcf at j=jebc-1, with right grid
            hzybcf[iebc+ie + i][j] = dahzybcf[iebc+ie + i][j] * hzybcf[iebc+ie + i][j] - dbhzybcf[iebc+ie + i][j] * (exbcf[iebc+ie + i][j] -  exbcr[i][0]);
        } /* iForLoop */     

        //     BACK

        for (i = 0; i < iefbc; i++) { 
            for (j = 1; j < jebc; j++) {  // don't evaluate hzxbcb at j=0 as it needs data from main,left,right grids, see below 
                hzybcb[i][j] = dahzybcb[i][j] * hzybcb[i][j] - dbhzybcb[i][j] * (exbcb[i][j] - exbcb[i][j+1] );
            } /* jForLoop */     
        } /* iForLoop */     

        for (j = 0, i = 0; i < iebc; i++) {  // fix-up hzybcb at j=0, with left grid
            hzybcb[i][j] = dahzybcb[i][j] * hzybcb[i][j] - dbhzybcb[i][j] * (exbcl[i][je] - exbcb[i][j+1] );
        } /* iForLoop */     
               
        for (j = 0, i = 0; i < ie; i++) {  // fix-up hzybcb at j=0, with main grid
            hzybcb[iebc+i][j] = dahzybcb[iebc+i][j] * hzybcb[iebc+i][j] - dbhzybcb[iebc+i][j] * (ex[i][je] - exbcb[iebc+i][j+1]);
        } /* iForLoop */     

        for (j = 0, i = 0; i < iebc; i++) {  // fix-up hzybcb at j=0, with right grid
            hzybcb[iebc+ie + i][j] = dahzybcb[iebc+ie + i][j] * hzybcb[iebc+ie + i][j] - dbhzybcb[iebc+ie + i][j] * (exbcr[i][je] - exbcb[iebc+ie + i][j+1] );
        } /* iForLoop */     

        //     LEFT
         
        for (i = 0; i < iebc; i++) {
            for (j = 0; j < je; j++) {
                hzybcl[i][j] = dahzybcl[i][j] * hzybcl[i][j] - dbhzybcl[i][j] * (exbcl[i][j] - exbcl[i][j+1] );
            } /* jForLoop */     
        } /* iForLoop */     

        //     RIGHT
         
        for (i = 0; i < iebc; i++) {
            for (j = 0; j < je; j++) {
                hzybcr[i][j] = dahzybcr[i][j] * hzybcr[i][j] - dbhzybcr[i][j] * (exbcr[i][j] - exbcr[i][j+1] );
            } /* jForLoop */     
        } /* iForLoop */     

        n++;
    }

// standard C memory allocation for 2-D array
double  **AllocateMemory (int  imax, int  jmax, double  initialValue)
{
    int  i,j;                 
    double  **pointer;
    pointer = (double **)malloc(imax * sizeof(double *));
    for (i = 0; i < imax; i++) {
        pointer[i] = (double *)malloc(jmax * sizeof(double));
        for (j = 0; j < jmax; j++) {
            pointer[i][j] = initialValue;
        } /* jForLoop */     
    } /* iForLoop */     
    return(pointer);
}    

    // template <size_t ie, size_t je>
    // template <typename TwoD>
    // void  InitializeFdtd (TwoD&  material_matrix)
    void  InitializeFdtd (std::vector<std::vector<int>> &material_matrix, std::vector<double> &eps,
                   std::vector<double> &mur,
                   std::vector<double> &sig, int src_position_row, int src_position_col, int src_type)
                   // src_type: 1 - gaussian, 2 - sin
    // void  InitializeFdtd (size_t (&material_matrix)[ie][je])
    
    // void  InitializeFdtd ()
    {
   
        this->src_type = src_type;

    //***********************************************************************
    //     Printing/Plotting variables
    //***********************************************************************
    minimumValue = -0.001;
    maximumValue =  0.001;   
    plottingInterval = 0;
    centery = 25;
    centerx = 15;

    //***********************************************************************
    //     Fundamental constants
    //***********************************************************************
    pi  = (acos(-1.0));
    cc = 2.99792458e8;                  //speed of light in free space (meters/second)
    muz = 4.0 * pi * 1.0e-7;            //permeability of free space
    epsz = 1.0 / (cc * cc * muz);       //permittivity of free space

    freq = 6.0e+12;                  //center frequency of source excitation (Hz)
    lambda = cc / freq;             //center wavelength of source excitation
    omega = 2.0 * pi * freq;        //center frequency in radians  


    double eta1 = sqrt(mur[0]*muz/(eps[0]*epsz));
    double eta2 = sqrt(mur[1]*muz/(eps[1]*epsz));
    double reflection_coef = abs((eta2-eta1) / (eta2 + eta1));
    double transmittance_coef = 1-reflection_coef;

    std::cout << "eta1: " << eta1 << std::endl;
    std::cout << "eta2: " << eta2 << std::endl;
    std::cout << "reflection_coef: " << reflection_coef << std::endl;
    std::cout << "transmittance_coef: " << 1-abs(reflection_coef) << std::endl;


    //***********************************************************************
    //     Grid parameters
    //***********************************************************************
     
    // ie = 100;               //number of grid cells in x-direction
    // je = 50;                //number of grid cells in y-direction
     
    ib = ie + 1;            // one extra is needed for fields on the boundaries (ie Ex on top boundary, Ey on right boundary)
    jb = je + 1;            // ditto
     
    // is = 15;                //location of z-directed hard source
    // js = je / 2;            //location of z-directed hard source

    // swap later
    js = src_position_row;                //location of z-directed hard source
    is = src_position_col;            //location of z-directed hard source
     
    // grid size = dx * nx = 3.0e-6m * 220 = 660-6m = 660 mkm
    dx = 3.0e-6;            //space increment of square lattice  (meters)
    dt = dx / (2.0 * cc);   //time step,  seconds, courant limit, Taflove1995 page 177
     
    nmax = NUMBEROFITERATIONCONSTANT;             //total number of time steps
     
    iebc = 10;               //thickness of left and right PML region
    jebc = 10;               //thickness of front and back PML region
    rmax = 0.00001;         // R(0) reflection coefficient (in %)  Nikolova part4 p.25
    orderbc = 2;            // m, grading order, optimal values: 2 <= m <= 6,  Nikolova part4 p.29
    ibbc = iebc + 1;
    jbbc = jebc + 1;
    iefbc = ie + 2 * iebc;  // for front and bottom (width of region)
    jefbc = je + 2 * jebc;  // not used
    ibfbc = iefbc + 1;      // one extra for Ey on right boundary
    jbfbc = jefbc + 1;      // not used
   

    //***********************************************************************
    //     Material parameters
    //***********************************************************************

    media = MEDIACONSTANT;        // number of different medias, ie 2: vacuum, metallicCylinder


    //***********************************************************************
    //     Wave excitation
    //***********************************************************************

    rtau = 160.0e-15;
    // rtau = 80.0e-15;
    tau = rtau / dt;
    delay = 3 * tau;
    // delay = 1 * tau;
    for (i = 0; i < nmax; i++) {
        source[i] = 0.0;
    } /* iForLoop */     
    
    int nn;
    // for (nn = 0; nn < (int  )(7.0 * tau); nn++) {
    for (nn = 0; nn < nmax; nn++) {
        temporary = (double  )nn - delay;
        source[nn] = 10*sin( omega * (temporary) * dt) * exp(-( (temporary * temporary)/(tau * tau) ) );
        // source[nn] = 10*sin( omega * (temporary) * dt);
    } /* forLoop */     


    //***********************************************************************
    //     Field arrays
    //***********************************************************************

    ex = AllocateMemory(ie,jb, 0.0 );           //fields in main grid 
    ey = AllocateMemory(ib,je, 0.0 );
    hz = AllocateMemory(ie,je, 0.0 );
     
    exbcf = AllocateMemory(iefbc,jebc, 0.0 );   //fields in front PML region
    eybcf = AllocateMemory(ibfbc,jebc, 0.0 );
    hzxbcf = AllocateMemory(iefbc,jebc, 0.0 );
    hzybcf = AllocateMemory(iefbc,jebc, 0.0 );
     
    exbcb = AllocateMemory(iefbc,jbbc, 0.0 );   //fields in back PML region
    eybcb = AllocateMemory(ibfbc,jebc, 0.0 );
    hzxbcb = AllocateMemory(iefbc,jebc, 0.0 );
    hzybcb = AllocateMemory(iefbc,jebc, 0.0 );
     
    exbcl = AllocateMemory(iebc,jb, 0.0 );      //fields in left PML region
    eybcl = AllocateMemory(iebc,je, 0.0 );
    hzxbcl = AllocateMemory(iebc,je, 0.0 );
    hzybcl = AllocateMemory(iebc,je, 0.0 );
     
    exbcr = AllocateMemory(iebc,jb, 0.0 );      //fields in right PML region
    eybcr = AllocateMemory(ibbc,je, 0.0 );
    hzxbcr = AllocateMemory(iebc,je, 0.0 );
    hzybcr = AllocateMemory(iebc,je, 0.0 );


    //***********************************************************************
    //     Updating coefficients
    //***********************************************************************

    for (i = 0; i < media; i++) {
        eaf   = dt * sig[i] / (2.0 * epsz * eps[i] );    // Taflove1995 p.67
        ca[i] = (1.0 - eaf) / (1.0 + eaf);               // ditto
        cb[i] = dt / epsz / eps[i] / dx / (1.0 + eaf);   // ditto
        haf   = dt * sim[i] / (2.0 * muz * mur[i]);      // ditto
        da[i] = (1.0 - haf) / (1.0 + haf);               // ditto
        db[i] = dt / muz / mur[i] / dx / (1.0 + haf);    // ditto
    } /* iForLoop */     
    

    //***********************************************************************
    //     Geometry specification (main grid)
    //***********************************************************************

    //     Initialize entire main grid to free space

    caex = AllocateMemory(ie,jb, ca[0]);     
    cbex = AllocateMemory(ie,jb, cb[0]);
                                   
    caey = AllocateMemory(ib,je, ca[0] );
    cbey = AllocateMemory(ib,je, cb[0] );

    dahz = AllocateMemory(ie,je, da[0] );
    dbhz = AllocateMemory(ie,je, db[0] );

    //     Add metal cylinder

    diam = 20;          // diameter of cylinder: 6 cm
    rad = diam / 2.0;     // radius of cylinder: 3 cm
    icenter = (4 * ie) / 5;   // i-coordinate of cylinder's center
    jcenter = je / 2;     // j-coordinate of cylinder's center
    for (i = 0; i < ie; i++) {
        for (j = 0; j < je; j++) {
            temporaryi = (double  )(i - icenter);
            temporaryj = (double  )(j - jcenter);
            dist2 = (temporaryi + 0.5) * (temporaryi + 0.5) + (temporaryj) * (temporaryj);

            caex[i][j] = ca[material_matrix[i][j]];
            cbex[i][j] = cb[material_matrix[i][j]];

            caey[i][j] = ca[material_matrix[i][j]];
            cbey[i][j] = cb[material_matrix[i][j]];

            // std::cout << i*je + j << ": " << ca[material_matrix[i][j]] << std::endl;


            // if (dist2 <= (rad * rad)) {
            //     caex[i][j] = ca[1];
            //     cbex[i][j] = cb[1];
            // } /* if */     

            // // This looks tricky! Why can't caey/cbey use the same 'if' statement as caex/cbex above ?? 
            // dist2 = (temporaryj + 0.5) * (temporaryj + 0.5) + (temporaryi) * (temporaryi);
            // if (dist2 <= (rad * rad)) {
            //     caey[i][j] = ca[1];
            //     cbey[i][j] = cb[1];
            // } /* if */     
        } /* jForLoop */     
    } /* iForLoop */     


    //***********************************************************************
    //     Fill the PML regions
    //***********************************************************************

    caexbcf = AllocateMemory(iefbc,jebc, 0.0 );
    cbexbcf = AllocateMemory(iefbc,jebc, 0.0 );
    caexbcb = AllocateMemory(iefbc,jbbc, 0.0 );
    cbexbcb = AllocateMemory(iefbc,jbbc, 0.0 );
    caexbcl = AllocateMemory(iebc,jb, 0.0 );
    cbexbcl = AllocateMemory(iebc,jb, 0.0 );
    caexbcr = AllocateMemory(iebc,jb, 0.0 ); 
    cbexbcr = AllocateMemory(iebc,jb, 0.0 ); 
    caeybcf = AllocateMemory(ibfbc,jebc, 0.0 );
    cbeybcf = AllocateMemory(ibfbc,jebc, 0.0 );
    caeybcb = AllocateMemory(ibfbc,jebc, 0.0 );
    cbeybcb = AllocateMemory(ibfbc,jebc, 0.0 );
    caeybcl = AllocateMemory(iebc,je, 0.0 );
    cbeybcl = AllocateMemory(iebc,je, 0.0 );
    caeybcr = AllocateMemory(ibbc,je, 0.0 );
    cbeybcr = AllocateMemory(ibbc,je, 0.0 );
    dahzxbcf = AllocateMemory(iefbc,jebc, 0.0 );
    dbhzxbcf = AllocateMemory(iefbc,jebc, 0.0 );
    dahzxbcb = AllocateMemory(iefbc,jebc, 0.0 );
    dbhzxbcb = AllocateMemory(iefbc,jebc, 0.0 );
    dahzxbcl = AllocateMemory(iebc,je, 0.0 );
    dbhzxbcl = AllocateMemory(iebc,je, 0.0 );
    dahzxbcr = AllocateMemory(iebc,je, 0.0 );
    dbhzxbcr = AllocateMemory(iebc,je, 0.0 );
    dahzybcf = AllocateMemory(iefbc,jebc, 0.0 );
    dbhzybcf = AllocateMemory(iefbc,jebc, 0.0 );
    dahzybcb = AllocateMemory(iefbc,jebc, 0.0 );
    dbhzybcb = AllocateMemory(iefbc,jebc, 0.0 );
    dahzybcl = AllocateMemory(iebc,je, 0.0 );
    dbhzybcl = AllocateMemory(iebc,je, 0.0 );
    dahzybcr = AllocateMemory(iebc,je, 0.0 );
    dbhzybcr = AllocateMemory(iebc,je, 0.0 );

    delbc = (double  )iebc * dx;    // width of PML region (in mm)
                                                                                                      
    // SigmaMaximum, using polynomial grading (Nikolova part 4, p.30), rmax=reflectionMax in percent
    sigmam =-log(rmax/100.0) * epsz * cc * (orderbc + 1) / (2 * delbc);

    // bcfactor comes from the polynomial grading equation: sigma_x = sigmaxMaximum * (x/d)^m, where d=width of PML, m=gradingOrder, (Nikolova part4, p.28)
    //  IMPORTANT: The conductivity (sigma) must use the "average" value at each mesh point as follows:
    //  sigma_x = sigma_Maximum/dx * Integral_from_x0_to_x1 of (x/d)^m dx,  where x0=currentx-0.5, x1=currentx+0.5   (Nikolova part 4, p.32)
    //  integrating gives: sigma_x = (sigmaMaximum / (dx * d^m * m+1)) * ( x1^(m+1) - x0^(m+1) )     (Nikolova part 4, p.32)
    //  the first part is "bcfactor", so, sigma_x = bcfactor * ( x1^(m+1) - x0^(m+1) )   (Nikolova part 4, p.32)
    // note: it's not exactly clear what the term eps[0] is for. It's probably to cover the case in which eps[0] is not equal to one (ie the main grid area next to the pml boundary is not vacuum)
    bcfactor = eps[0] * sigmam / ( dx * (pow(delbc,orderbc)) * (orderbc + 1));

    //     FRONT region 

    for (i = 0; i < iefbc; i++) {   // IS THIS EVER USED? -- coef for the pec ex at j=0 set so that ex_t+1 = ex_t, but ex is never evaluated at j=0.
        caexbcf[i][0] = 1.0;
        cbexbcf[i][0] = 0.0;
    } /* iForLoop */     

    for (j = 1; j < jebc; j++) {     // calculate the coefs for the PML layer (except for the boundary at the main grid, which is a special case (see below))
        y1 = ((double  )(jebc - j) + 0.5) * dx;       // upper bounds for point j      (re-adujsted for C indexes!)
        y2 = ((double  )(jebc - j) - 0.5) * dx;       // lower bounds for point j
        sigmay = bcfactor * (pow(y1,(orderbc+1)) - pow(y2,(orderbc+1)) );   //   polynomial grading
        ca1 = exp (-sigmay * dt / (epsz * eps[0]));   // exponential time step, Taflove1995 p.77,78
        cb1 = (1.0 - ca1) / (sigmay * dx);            // ditto, but note sign change from Taflove1995
        for (i = 0; i < iefbc; i++) {
            caexbcf[i][j] = ca1;
            cbexbcf[i][j] = cb1;
        } /* iForLoop */     
    } /* jForLoop */     
    
    sigmay = bcfactor * pow( (0.5 * dx), (orderbc+1) );  // calculate for the front edge of the pml at j=0 in the main grid  (half vacuum (sigma=0) and half pml)
    ca1 = exp( -sigmay * dt / (epsz * eps[0]) );
    cb1 = ( 1 - ca1) / (sigmay * dx);
    for (i = 0; i < ie; i++) { 
        caex[i][0] = ca1;
        cbex[i][0] = cb1;
    } /* iForLoop */     

    for (i = 0; i < iebc; i++) {   // this continues the front edge into the left and right grids
        caexbcl[i][0] = ca1;
        cbexbcl[i][0] = cb1;
        caexbcr[i][0] = ca1;
        cbexbcr[i][0] = cb1;
    } /* iForLoop */     
    
    for (j = 0; j < jebc; j++) {    // for ey and hz  (which are offset spacially 1/2 dx from ex)
        y1 = ((double  )(jebc - j) + 0.0) * dx;       // upper bounds for point j  
        y2 = ((double  )(jebc - j) - 1.0) * dx;       // lower bounds for point j
        sigmay = bcfactor * (pow(y1,(orderbc+1)) - pow(y2,(orderbc+1)) );   //   polynomial grading
        sigmays = sigmay * (muz / (epsz*eps[0]) );    // Taflove1995 p.182  (for no reflection: sigmaM = sigmaE * mu0/eps0)
        da1 = exp (-sigmays * dt / muz);              // exponential time step, Taflove1995 p.77,78
        db1 = (1.0 - da1) / (sigmays * dx);
        for (i = 0; i < iefbc; i++) {
            dahzybcf[i][j] = da1;
            dbhzybcf[i][j] = db1;
            dahzxbcf[i][j] = da[0];   // important note: hzx is Perpendicular to the front pml and so is not attenuated (sigma=0) (looks like vacuum) 
            dbhzxbcf[i][j] = db[0];   // ditto
        } /* iForLoop */     
        for (i = 0; i < ibfbc; i++) {
            caeybcf[i][j]  = ca[0];   // important note: ey is Perpendicular to the front pml and so is not attenuated (sigma=0) (looks like vacuum)
            cbeybcf[i][j]  = cb[0];   // ditto
        } /* iForLoop */     
    } /* jForLoop */     


    //     BACK region 

    for (j=jebc, i = 0; i < iefbc; i++) {   // IS THIS EVER USED? -- coef for the pec ex at j=jebc set to ex_t+1 = ex_t
        caexbcb[i][j] = 1.0;
        cbexbcb[i][j] = 0.0;
    } /* iForLoop */     

    for (j = 1; j < jebc; j++) {     // calculate the coefs for the PML layer (except for the boundary at the main grid, which is a special case (see below))
        y1 = ((double  )(j) + 0.5) * dx;       // upper bounds for point j         (re-adujsted for C indexes!)
        y2 = ((double  )(j) - 0.5) * dx;       // lower bounds for point j
        sigmay = bcfactor * (pow(y1,(orderbc+1)) - pow(y2,(orderbc+1)) );   //   polynomial grading
        ca1 = exp (-sigmay * dt / (epsz * eps[0]));   // exponential time step
        cb1 = (1.0 - ca1) / (sigmay * dx);            // ditto, but note sign change from Taflove
        for (i = 0; i < iefbc; i++) {
            caexbcb[i][j] = ca1;
            cbexbcb[i][j] = cb1;
        } /* iForLoop */     
    } /* jForLoop */     

    sigmay = bcfactor * pow( (0.5 * dx), (orderbc+1) );  // calculate for the front edge of the pml at j=0 in the main grid  (half vacuum (sigma=0) and half pml)
    ca1 = exp( -sigmay * dt / (epsz * eps[0]) );
    cb1 = ( 1 - ca1) / (sigmay * dx);
    for (j=je, i = 0; i < ie; i++) { 
        caex[i][j] = ca1;
        cbex[i][j] = cb1;
    } /* iForLoop */     

    for (j=je, i = 0; i < iebc; i++) {  // this continues the back edge into the left and right grids 
        caexbcl[i][j] = ca1;
        cbexbcl[i][j] = cb1;
        caexbcr[i][j] = ca1;
        cbexbcr[i][j] = cb1;
    } /* iForLoop */     

    for (j = 0; j < jebc; j++) {    // for ey and hz  (which are offset spacially 1/2 dx from ex)
        y1 = ((double  )(j) + 1.0) * dx;       // upper bounds for point j  
        y2 = ((double  )(j) + 0.0) * dx;       // lower bounds for point j  
        sigmay = bcfactor * (pow(y1,(orderbc+1)) - pow(y2,(orderbc+1)) );   //   polynomial grading
        sigmays = sigmay * (muz / (epsz*eps[0]) );
        da1 = exp (-sigmays * dt / muz);
        db1 = (1.0 - da1) / (sigmays * dx);
        for (i = 0; i < iefbc; i++) {
            dahzybcb[i][j] = da1;
            dbhzybcb[i][j] = db1;
            dahzxbcb[i][j] = da[0];   // important note: hzx is Perpendicular to the back pml and so is not attenuated (sigma=0) (looks like vacuum) 
            dbhzxbcb[i][j] = db[0];   // ditto
        } /* iForLoop */     
        for (i = 0; i < ibfbc; i++) {
            caeybcb[i][j]  = ca[0];   // important note: ey is Perpendicular to the back pml and so is not attenuated (sigma=0) (looks like vacuum)
            cbeybcb[i][j]  = cb[0];   // ditto
        } /* iForLoop */     
    } /* jForLoop */     

    //     LEFT region 


    for (i=0, j = 0; j < je; j++) {    // IS THIS EVER USED? -- coef for the pec ey at i=0 set to ey_t+1 = ey_t 
        caeybcl[i][j] = 1.0;
        cbeybcl[i][j] = 0.0;
    } /* jForLoop */     
    
    for (i = 1; i < iebc; i++) {    // calculate the coefs for the PML layer (except for the boundary at the main grid, which is a special case (see below)) 
        x1 = ((double  )(iebc - i) + 0.5) * dx;     // upper bounds for point i    (re-adujsted for C indexes!)
        x2 = ((double  )(iebc - i) - 0.5) * dx;     // lower bounds for point i   
        sigmax = bcfactor * ( pow(x1,(orderbc+1)) - pow(x2,(orderbc+1)) );
        ca1 = exp(-sigmax * dt / (epsz * eps[0]) );
        cb1 = (1.0 - ca1) / (sigmax * dx);
        for (j = 0; j < je; j++) {
            caeybcl[i][j] = ca1;
            cbeybcl[i][j] = cb1;
        } /* jForLoop */     
        for (j = 0; j < jebc; j++) {   // fill in the front left and back left corners for ey
            caeybcf[i][j] = ca1;    
            cbeybcf[i][j] = cb1;
            caeybcb[i][j] = ca1;
            cbeybcb[i][j] = cb1;
        } /* jForLoop */     
    } /* iForLoop */     
    
    sigmax = bcfactor * pow( (0.5 * dx), (orderbc+1) );  // calculate for the left edge of the pml at x=0 in the main grid  (half vacuum (sigma=0) and half pml)
    ca1 = exp(-sigmax * dt / (epsz * eps[0]) );
    cb1 = (1.0 - ca1) / (sigmax * dx);
    for (i=0, j = 0; j < je; j++) {
        caey[i][j] = ca1;
        cbey[i][j] = cb1;
    } /* jForLoop */     
    for (i=iebc, j = 0; j < jebc; j++) {   // continue the left edge into the front and back grids
        caeybcf[i][j] = ca1;
        cbeybcf[i][j] = cb1;
        caeybcb[i][j] = ca1;
        cbeybcb[i][j] = cb1;
    } /* jForLoop */     
    
        
    for (i = 0; i < iebc; i++) {          // for ex and hz  (which are offset spacially 1/2 dx from ey) 
        x1 = ((double  )(iebc - i) + 0.0) * dx;     // upper bounds for point i    (re-adujsted for C indexes!)
        x2 = ((double  )(iebc - i) - 1.0) * dx;     // lower bounds for point i   
        sigmax = bcfactor * ( pow(x1,(orderbc+1)) - pow(x2,(orderbc+1)) );
        sigmaxs = sigmax * (muz / (epsz * eps[0]) );
        da1 = exp(-sigmaxs * dt / muz);
        db1 = (1 - da1) / (sigmaxs * dx);
        for (j = 0; j < je; j++) {
            dahzxbcl[i][j] = da1;
            dbhzxbcl[i][j] = db1;
            dahzybcl[i][j] = da[0];   // important note: hzy is Perpendicular to the left pml and so is not attenuated (sigma=0) (looks like vacuum) 
            dbhzybcl[i][j] = db[0];
        } /* jForLoop */     
        for (j = 0; j < jebc; j++) {    // fill in the front left and back left corners for hzx 
            dahzxbcf[i][j] = da1;
            dbhzxbcf[i][j] = db1;
            dahzxbcb[i][j] = da1;
            dbhzxbcb[i][j] = db1;
        } /* jForLoop */     
        for (j = 1; j < je; j++) {    // important note: ex is Perpendicular to the left pml and so is not attenuated (sigma=0) (looks like vacuum) 
            caexbcl[i][j] =ca[0];
            cbexbcl[i][j] =cb[0];
        } /* jForLoop */     
    } /* iForLoop */     
        


    //     RIGHT region 

    for (i=iebc, j = 0; j < je; j++) {    // IS THIS EVER USED? -- coef for the pec ey at i=iebc set to ey_t+1 = ey_t 
        caeybcr[i][j] = 1.0;
        cbeybcr[i][j] = 0.0;
    } /* jForLoop */     

    for (i = 1; i < iebc; i++) {    // calculate the coefs for the PML layer (except for the boundary at the main grid, which is a special case (see below)) 
        x1 = ((double  )(i) + 0.5) * dx;     // upper bounds for point i        (re-adujsted for C indexes!)
        x2 = ((double  )(i) - 0.5) * dx;     // lower bounds for point i   
        sigmax = bcfactor * ( pow(x1,(orderbc+1)) - pow(x2,(orderbc+1)) );
        ca1 = exp(-sigmax * dt / (epsz * eps[0]) );
        cb1 = (1.0 - ca1) / (sigmax * dx);
        for (j = 0; j < je; j++) {
            caeybcr[i][j] = ca1;
            cbeybcr[i][j] = cb1;
        } /* jForLoop */     
        for (j = 0; j < jebc; j++) {   // fill in the front right and back right corners for ey
            caeybcf[i+iebc+ie][j] = ca1;    
            cbeybcf[i+iebc+ie][j] = cb1;
            caeybcb[i+iebc+ie][j] = ca1;
            cbeybcb[i+iebc+ie][j] = cb1;
        } /* jForLoop */     
    } /* iForLoop */     


    sigmax = bcfactor * pow( (0.5 * dx), (orderbc+1) );  // calculate for the right edge of the pml at x=ic in the main grid  (half vacuum (sigma=0) and half pml)
    ca1 = exp(-sigmax * dt / (epsz * eps[0]) );
    cb1 = (1.0 - ca1) / (sigmax * dx);
    for (i=ie, j = 0; j < je; j++) {
        caey[i][j] = ca1;
        cbey[i][j] = cb1;
    } /* jForLoop */     
    for (i=iebc+ie, j = 0; j < jebc; j++) {   // continue the right edge into the front and back grids
        caeybcf[i][j] = ca1;
        cbeybcf[i][j] = cb1;
        caeybcb[i][j] = ca1;
        cbeybcb[i][j] = cb1;
    } /* jForLoop */     


    for (i = 0; i < iebc; i++) {          // for ex and hz  (which are offset spacially 1/2 dx from ey) 
        x1 = ((double  )(i) + 1.0) * dx;     // upper bounds for point i         (re-adujsted for C indexes!)
        x2 = ((double  )(i) + 0.0) * dx;     // lower bounds for point i   
        sigmax = bcfactor * ( pow(x1,(orderbc+1)) - pow(x2,(orderbc+1)) );
        sigmaxs = sigmax * (muz / (epsz * eps[0]) );
        da1 = exp(-sigmaxs * dt / muz);
        db1 = (1 - da1) / (sigmaxs * dx);
        for (j = 0; j < je; j++) {
            dahzxbcr[i][j] = da1;
            dbhzxbcr[i][j] = db1;
            dahzybcr[i][j] = da[0];   // important note: hzy is Perpendicular to the right pml and so is not attenuated (sigma=0) (looks like vacuum) 
            dbhzybcr[i][j] = db[0];
        } /* jForLoop */     
        for (j = 0; j < jebc; j++) {    // fill in the front right and back right corners for hzx 
            dahzxbcf[i+ie+iebc][j] = da1;
            dbhzxbcf[i+ie+iebc][j] = db1;
            dahzxbcb[i+ie+iebc][j] = da1;
            dbhzxbcb[i+ie+iebc][j] = db1;
        } /* jForLoop */     
        for (j = 1; j < je; j++) {    // important note: ex is Perpendicular to the right pml and so is not attenuated (sigma=0) (looks like vacuum) 
            caexbcr[i][j] =ca[0];
            cbexbcr[i][j] =cb[0];
        } /* jForLoop */     
    } /* iForLoop */     
        

    }
                       
};

#endif