// dllmain.cpp : Defines the entry point for the DLL application.
// #include "pch.h"
#include "Eigen/Dense"
#include "json.hpp"
#include "ASRE_Timo_functions.h"
#include <cstring>

// #ifdef _MSVC //Define the macro and DLLMain needed on Windows
//     #include <windows.h>
//     #define DLLEXPORT __declspec(dllexport)
//     BOOL APIENTRY DllMain( HMODULE hModule,
//                         DWORD  ul_reason_for_call,
//                         LPVOID lpReserved
//                         )
//     {
//         switch (ul_reason_for_call)
//         {
//         case DLL_PROCESS_ATTACH:
//         case DLL_THREAD_ATTACH:
//         case DLL_THREAD_DETACH:
//         case DLL_PROCESS_DETACH:
//             break;
//         }
//         return TRUE;
//     }
// #else
//     #define DLLEXPORT __declspec(dllexport)
// #endif

extern "C" {

    //__declspec(dllexport) int run() {
    DLLEXPORT int run(int nnode, double* meshX, double* meshY, double* meshZ, double* dispV, double* dispL, double* dispT,
        double Eb, double EoverG, double EsNominal, double nis, double dfoot, double bfoot, double ni_foot, double mu_int, double qz_foot,
        double d_NA,
        const char* solver,
        const char* output,
        double* result_array,
        int result_size
        ) {

        #ifdef PRINT_INT_RESULTS
        std::cout << "ASRE run starts" << std::endl;
        std::cout << "meshX" << std::endl;
        std::cout << meshX << std::endl;
        std::cout << "output type" << std::endl;
        std::cout << output << std::endl;
        #endif

        VectorXd x_global_foot = VectorXd::Zero(nnode);
        VectorXd y_global_foot = VectorXd::Zero(nnode);
        VectorXd z_global_foot = VectorXd::Zero(nnode);
        VectorXd ucat = VectorXd::Zero(nnode * 6);
        for (int i = 0; i < nnode; i++) {
            x_global_foot(i) = meshX[i];
            y_global_foot(i) = meshY[i];
            z_global_foot(i) = meshZ[i];
            ucat(i * 6) = dispL[i];
            ucat(i * 6 + 1) = dispT[i];
            ucat(i * 6 + 2) = dispV[i];
        }
        double h_el_foot = std::sqrt((x_global_foot(0) - x_global_foot(nnode - 1)) * (x_global_foot(0) - x_global_foot(nnode - 1)) +
            (y_global_foot(0) - y_global_foot(nnode - 1)) * (y_global_foot(0) - y_global_foot(nnode - 1))) / (nnode - 1);
        
        #ifdef PRINT_INT_RESULTS
        //std::cout << "x_global_foot: " << std::endl;
        //std::cout << x_global_foot << std::endl;
        //std::cout << "ucat: " << std::endl;
        //std::cout << ucat << std::endl;
        //std::cout << "h_el_foot: " << std::endl;
        //std::cout << h_el_foot << std::endl;
        #endif

        int nelement = nnode - 1;
        int dimK = nnode * 6;
        MatrixXd KKfoot = MatrixXd::Zero(dimK, dimK);
        for (int i = 0; i < nelement; i++) {
            double dx = std::sqrt((meshX[i] - meshX[i + 1]) * (meshX[i] - meshX[i + 1]) + (meshY[i] - meshY[i + 1]) * (meshY[i] - meshY[i + 1]));
            MatrixXd KBern3Delt = KBern3D_foot_TIM_dNA(Eb, dfoot, bfoot, dx, EoverG, ni_foot, d_NA);
            KKfoot(seq((i) * 6, (i) * 6 + 11, 1),
                seq((i) * 6, (i) * 6 + 11, 1)) =
                KKfoot(seq((i) * 6, (i) * 6 + 11, 1),
                    seq((i) * 6, (i) * 6 + 11, 1)) + KBern3Delt;
        }

        #ifdef PRINT_INT_RESULTS
        std::cout << "KKfoot: " << std::endl;
        std::cout << KKfoot << std::endl;
        #endif
        
        double dx = std::sqrt((meshX[0] - meshX[0 + 1]) * (meshX[0] - meshX[0 + 1]) + (meshY[0] - meshY[0 + 1]) * (meshY[0] - meshY[0 + 1]));
        
        MatrixXd FLEX_3DOF = _calFlexVaziri(z_global_foot, x_global_foot, y_global_foot, EsNominal, nnode, h_el_foot, nis, bfoot);

        double kh_gazetas;
        double kv_gazetas;

        soilsprings_static_foot_gazetas(&kh_gazetas, &kv_gazetas, EsNominal, nis, bfoot, h_el_foot);
        MatrixXd Ks_3DOF = FLEX_3DOF.inverse();
        MatrixXd Ks = MatrixXd::Zero(nnode * 6, nnode * 6);
        MatrixXd FLEX = MatrixXd::Zero(nnode * 6, nnode * 6);
        for (int i = 0; i < nnode; i++) {
            for (int j = 0; j < nnode; j++) {
                Ks(seq(i * 6, 2 + i * 6, 1), seq(j * 6, 2 + j * 6, 1)) =
                    Ks_3DOF(seq(i * 3, 2 + i * 3, 1), seq(j * 3, 2 + j * 3, 1));
                FLEX(seq(i * 6, 2 + i * 6, 1), seq(j * 6, 2 + j * 6, 1)) =
                    FLEX_3DOF(seq(i * 3, 2 + i * 3, 1), seq(j * 3, 2 + j * 3, 1));
            }
            // ADDITIONAL STIFFNESS IN THE DOF 4 TO REMOVE UNCOSTRAINED ROTATION ABOUT THE BEAM AXIS
            Ks(i * 6 + 3, i * 6 + 3) = kv_gazetas * 1000;
            FLEX(i * 6 + 3, i * 6 + 3) = 1.0/(kv_gazetas * 1000.0);
        }
        Ks_3DOF.resize(0, 0);

        #ifdef PRINT_INT_RESULTS
        //FLEX = Ks.inverse();
        //std::cout << "FLEX after" << std::endl;
        //std::cout << FLEX << std::endl;
        #endif
        VectorXd P_el = VectorXd::Zero(6 * nnode);
        calc_P_el(&P_el, nnode, qz_foot, h_el_foot);

        #ifdef PRINT_INT_RESULTS
        //std::cout << "P_el: " << std::endl;
        //std::cout << P_el << std::endl;
        #endif
        // Because the deadload is unlikely to cause plastic deformation, use elastic solution to estimate initial displacement
        MatrixXd Amatrix = Ks + KKfoot;
        PartialPivLU<MatrixXd> lu = PartialPivLU<MatrixXd>(Amatrix);
        VectorXd u_P_el = lu.solve(P_el);
        Amatrix.resize(0, 0);
        #ifdef PRINT_INT_RESULTS
        //std::cout << "u_P_el: " << std::endl;
        //std::cout << u_P_el << std::endl;
        #endif



        VectorXd uinc;

        if (strcmp(solver, "elastic") == 0){
            //Elastic solution
            #ifdef PRINT_INT_RESULTS
            // std::cout << "elastic" << std::endl;
            //std::cout << "KKfoot: " << std::endl;
            //std::cout << KKfoot << std::endl;
            //std::cout << "Ks: " << std::endl;
            //std::cout << Ks << std::endl;
            #endif
            VectorXd u_cat_el = lu.solve(P_el + Ks * ucat);

            #ifdef PRINT_INT_RESULTS
            //std::cout << "u_cat_el: " << std::endl;
            //std::cout << u_cat_el << std::endl;
            #endif
            uinc = u_cat_el;
        } else if (strcmp(solver, "elasto-plastic") == 0) {
            //Prepare for elastioPlaticIteration to solve ucat induced displacement
            #ifdef PRINT_INT_RESULTS
            std::cout << "elasto-plastic" << std::endl;
            #endif
            VectorXd Kstar = VectorXd::Zero(nnode * 6);
            MatrixXd Lstar = MatrixXd::Zero(nnode * 6, nnode * 6);
            Lstar = FLEX;
            MatrixXd Stiffness = KKfoot;
            for (int i = 0; i < nnode; i++) {
                Kstar(i * 6 + 0) = 1.0 / FLEX(i * 6+0, i * 6 + 0);
                Kstar(i * 6 + 1) = 1.0 / FLEX(i * 6 + 1, i * 6 + 1);
                Kstar(i * 6 + 2) = 1.0 / FLEX(i * 6 + 2, i * 6 + 2);
                // Kstar(i * 6 + 3) = 1.0 / FLEX(i * 6 + 3, i * 6 + 3);
                Lstar(i * 6 + 0, i * 6 + 0) = 0.0;
                Lstar(i * 6 + 1, i * 6 + 1) = 0.0;
                Lstar(i * 6 + 2, i * 6 + 2) = 0.0;
                // Lstar(i * 6 + 3, i * 6 + 3) = 0.0;
                Stiffness(i * 6 + 0, i * 6 + 0) += Kstar(i * 6 + 0);
                Stiffness(i * 6 + 1, i * 6 + 1) += Kstar(i * 6 + 1);
                Stiffness(i * 6 + 2, i * 6 + 2) += Kstar(i * 6 + 2);
                // Stiffness(i * 6 + 3, i * 6 + 3) += Kstar(i * 6 + 3);
            }

                #ifdef PRINT_INT_RESULTS
                /*std::cout << "Kstar: " << std::endl;
                std::cout << Kstar << std::endl;
                std::cout << "Lstar: " << std::endl;
                std::cout << Lstar << std::endl;
                std::cout << "Stiffness: " << std::endl;
                std::cout << Stiffness << std::endl;*/
                #endif

            uinc = u_P_el;
            VectorXd uip = VectorXd::Zero(uinc.size());
            double lim_t_int = 0;
            double lim_c_int = INFINITY;

            #ifdef PRINT_INT_RESULTS
            /*std::cout << "ucat: " << std::endl;
            std::cout << ucat << std::endl;*/
            #endif
            elastoPlasticIterationLDLT(Stiffness, KKfoot, Kstar,
                FLEX, Lstar, uinc, ucat, P_el,
                mu_int, lim_t_int, lim_c_int,
                h_el_foot * bfoot,
                uip);
            
            #ifdef PRINT_INT_RESULTS
            /*std::cout << "uinc-u_P_el" << std::endl;
            std::cout << uinc-u_P_el << std::endl;*/
            #endif
        } else {
            std::cout << "Invalid solver type" << std::endl;
            return -1;
        }
        
        if (strcmp(output, "disp") == 0) {
            VectorXd result = uinc-u_P_el;
            // for (int i=0; i < result.size(); i++) {
            //     // std::cout << result(i) <<"," <<result.data()[i] << std::endl;
            //     result_array[i] = result(i);
            // }
            std::copy(result.data(), result.data() + result.size(), result_array);
            return 0;
        } else if (strcmp(output, "strain") == 0) {
            VectorXd F_M_deltaT_el_M, F_N_deltaT_el_M, F_S_deltaT_el_M;
            calInternalForces(&F_M_deltaT_el_M, &F_N_deltaT_el_M, &F_S_deltaT_el_M,
                    u_P_el, uinc, Eb, EoverG, h_el_foot, dfoot, bfoot, ni_foot, nnode);
            VectorXd epsilon_vector = calculateStrain(&F_S_deltaT_el_M, &F_M_deltaT_el_M,
                    &F_N_deltaT_el_M, Eb, EoverG, bfoot, dfoot, ni_foot, d_NA);
            VectorXd result = epsilon_vector;
            std::copy(result.data(), result.data() + result.size(), result_array);
            return 0;
        } else if (strcmp(output, "strain+disp") == 0) {
            VectorXd F_M_deltaT_el_M, F_N_deltaT_el_M, F_S_deltaT_el_M;
            calInternalForces(&F_M_deltaT_el_M, &F_N_deltaT_el_M, &F_S_deltaT_el_M,
                    u_P_el, uinc, Eb, EoverG, h_el_foot, dfoot, bfoot, ni_foot, nnode);
            VectorXd epsilon_vector = calculateStrain(&F_S_deltaT_el_M, &F_M_deltaT_el_M,
                    &F_N_deltaT_el_M, Eb, EoverG, bfoot, dfoot, ni_foot, d_NA);
            VectorXd strain = epsilon_vector;
            VectorXd disp = uinc - u_P_el;
            VectorXd result(strain.size() + disp.size());
            result << strain, disp;
            std::copy(result.data(), result.data() + result.size(), result_array);
            return 0;    
        } else if (strcmp(output, "strain+disp+force") == 0) {
            VectorXd F_M_deltaT_el_M, F_N_deltaT_el_M, F_S_deltaT_el_M;
            calInternalForces(&F_M_deltaT_el_M, &F_N_deltaT_el_M, &F_S_deltaT_el_M,
                    u_P_el, uinc, Eb, EoverG, h_el_foot, dfoot, bfoot, ni_foot, nnode);
            VectorXd epsilon_vector = calculateStrain(&F_S_deltaT_el_M, &F_M_deltaT_el_M,
                    &F_N_deltaT_el_M, Eb, EoverG, bfoot, dfoot, ni_foot, d_NA);
            VectorXd strain = epsilon_vector;
            VectorXd disp = uinc - u_P_el;
            VectorXd result(F_M_deltaT_el_M.size() + F_N_deltaT_el_M.size() + F_S_deltaT_el_M.size()+
                    strain.size() + disp.size());
            result << strain, disp, F_M_deltaT_el_M, F_N_deltaT_el_M, F_S_deltaT_el_M;
            std::copy(result.data(), result.data() + result.size(), result_array);
            return 0;
        } else {
            return -1;
        }
        
        




        /*MatrixXd KKextra;
        MatrixXd KKpg = calKKpg(&KKsoil_rid, kv_gazetas, KKfoot, &KKextra, nnode);

        MatrixXd SS = KKfoot + KKextra;
        MatrixXd lamdasts, lamdastd, Kst;
        calLamdastsLamdastdKst(&lamdasts, &lamdastd, &Kst, FLEX, nnode);*/
    }
}

