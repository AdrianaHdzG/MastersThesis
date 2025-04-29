// dllmain.cpp : Defines the entry point for the DLL application.
// #include "pch.h"
#include "Eigen/Dense"
#include "json.hpp"
#include "ASRE_Frame_functions.h"
#include <cstring>

extern "C" {

    //__declspec(dllexport) int run() {
    DLLEXPORT int run(double* frame_stiffness, int frame_system_size, int* footing_nodes_ind,
        int footing_nodes_size,
        double* footing_coord_x, double* footing_coord_y, double* footing_coord_z,
        double* footing_ele_length, double* footing_ele_width, 
        double EsNominal, double nis, double mu_int,
        double* dispV, double* dispL, double* dispT,
        double* external_load,
        const char* solver,
        double* result_array,
        int result_size
        ) {

        #ifdef PRINT_INT_RESULTS
        std::cout << "ASRE run starts" << std::endl;
        std::cout << "solver type" << std::endl;
        std::cout << solver << std::endl;
        #endif

        // Get the 
        VectorXd x_global_foot = VectorXd::Zero(footing_nodes_size);
        VectorXd y_global_foot = VectorXd::Zero(footing_nodes_size);
        VectorXd z_global_foot = VectorXd::Zero(footing_nodes_size);

        VectorXd footing_ele_length_eigen = VectorXd::Map(footing_ele_length, footing_nodes_size);
        VectorXd footing_ele_width_eigen = VectorXd::Map(footing_ele_width, footing_nodes_size);
        VectorXi footing_nodes_ind_eigen = VectorXi::Map(footing_nodes_ind, footing_nodes_size);
        
        for (int i = 0; i < footing_nodes_size; i++) {
            x_global_foot(i) = footing_coord_x[i];
            y_global_foot(i) = footing_coord_y[i];
            z_global_foot(i) = footing_coord_z[i];
        }
        MatrixXd FLEX_3DOF = _calFlexVaziri(z_global_foot, x_global_foot, 
            y_global_foot, EsNominal, footing_nodes_size, footing_ele_length, nis, footing_ele_width);
        
        #ifdef PRINT_INT_RESULTS
        std::cout << "FLEX_3DOF: " << std::endl;
        std::cout << FLEX_3DOF << std::endl;
        #endif
        
        int dimK = frame_system_size;
        MatrixXd KKfoot = MatrixXd::Map(frame_stiffness, dimK, dimK);

        #ifdef PRINT_INT_RESULTS
        std::cout << "KKfoot: " << std::endl;
        std::cout << KKfoot << std::endl;
        #endif
        
        double kh_gazetas;
        double kv_gazetas;

        MatrixXd Ks_3DOF = FLEX_3DOF.inverse();
        MatrixXd Ks = MatrixXd::Zero(frame_system_size, frame_system_size);
        std::vector<int> groundNodeDOF1, groundNodeDOF2, groundNodeDOF3, groundNodeDOF;
        for (int i = 0; i < footing_nodes_size; i++) {
            groundNodeDOF.push_back(footing_nodes_ind[i] * 6 + 0);
            groundNodeDOF.push_back(footing_nodes_ind[i] * 6 + 1);
            groundNodeDOF.push_back(footing_nodes_ind[i] * 6 + 2);
            // groundNodeDOF.push_back(footing_nodes_ind[i] * 6 + 3);
            // groundNodeDOF.push_back(footing_nodes_ind[i] * 6 + 4);
            // groundNodeDOF.push_back(footing_nodes_ind[i] * 6 + 5);
            groundNodeDOF1.push_back(footing_nodes_ind[i] * 6 + 0);
            groundNodeDOF2.push_back(footing_nodes_ind[i] * 6 + 1);
            groundNodeDOF3.push_back(footing_nodes_ind[i] * 6 + 2);
        }

        Ks(groundNodeDOF, groundNodeDOF) = Ks_3DOF;

        // MatrixXd FLEX = MatrixXd::Zero(frame_system_size, frame_system_size);
        // for (int i = 0; i < footing_nodes_size; i++) {
        //     int node_ind_i = footing_nodes_ind[i];
        //     for (int j = 0; j < footing_nodes_size; j++) {
        //         int node_ind_j = footing_nodes_ind[j];
        //         Ks(seq(node_ind_i * 6, 2 + node_ind_i * 6, 1), seq(node_ind_j * 6, 2 + node_ind_j * 6, 1)) =
        //             Ks_3DOF(seq(i * 3, 2 + i * 3, 1), seq(j * 3, 2 + j * 3, 1));
        //         // FLEX(seq(node_ind_i * 6, 2 + node_ind_i * 6, 1), seq(node_ind_j * 6, 2 + node_ind_j * 6, 1)) =
        //         //     FLEX_3DOF(seq(i * 3, 2 + i * 3, 1), seq(j * 3, 2 + j * 3, 1));
        //     }
        //     // ADDITIONAL STIFFNESS IN THE DOF 4 TO REMOVE UNCOSTRAINED ROTATION ABOUT THE BEAM AXIS
        //     soilsprings_static_foot_gazetas(&kh_gazetas, &kv_gazetas, EsNominal, nis, footing_ele_width[i], footing_ele_length[i]);
        //     Ks(node_ind_i * 6 + 3, node_ind_i* 6 + 3) = kv_gazetas * 1000;
        //     // FLEX(node_ind_i * 6 + 3, node_ind_i * 6 + 3) = 1.0/(kv_gazetas * 1000.0);
        // }
        Ks_3DOF.resize(0, 0);
        
        #ifdef PRINT_INT_RESULTS
        //FLEX = Ks.inverse();
        //std::cout << "FLEX after" << std::endl;
        //std::cout << FLEX << std::endl;
        #endif
        VectorXd P_el = VectorXd::Map(external_load, frame_system_size);

        #ifdef PRINT_INT_RESULTS
        std::cout << "P_el: " << std::endl;
        std::cout << P_el << std::endl;
        #endif
        
        // Because the deadload is unlikely to cause plastic deformation, use elastic solution to estimate initial displacement
        MatrixXd Amatrix = Ks + KKfoot;
        #ifdef PRINT_INT_RESULTS
        std::cout << "Ks: " << std::endl;
        std::cout << Ks << std::endl;
        std::cout << "Amatrix: " << std::endl;
        std::cout << Amatrix << std::endl;
        #endif

        PartialPivLU<MatrixXd> lu = PartialPivLU<MatrixXd>(Amatrix);
        VectorXd u_P_el = lu.solve(P_el);
        Amatrix.resize(0, 0);
        #ifdef PRINT_INT_RESULTS
        std::cout << "u_P_el: " << std::endl;
        std::cout << u_P_el << std::endl;
        #endif


        
        VectorXd uinc;

        if (strcmp(solver, "elastic") == 0){
            //Elastic solution
            #ifdef PRINT_INT_RESULTS
            std::cout << "elastic" << std::endl;
            std::cout << "frame_system_size: " << frame_system_size << std::endl;
            //std::cout << "KKfoot: " << std::endl;
            //std::cout << KKfoot << std::endl;
            //std::cout << "Ks: " << std::endl;
            //std::cout << Ks << std::endl;
            #endif
            VectorXd ucat = VectorXd::Zero(frame_system_size);
            for (int i = 0; i < footing_nodes_size; i++) {
                int node_ind = footing_nodes_ind[i];
                ucat(node_ind * 6) = dispL[i];
                ucat(node_ind * 6 + 1) = dispT[i];
                ucat(node_ind * 6 + 2) = dispV[i];
            }
            VectorXd u_cat_el = lu.solve(P_el + Ks * ucat);

            #ifdef PRINT_INT_RESULTS
            std::cout << "u_cat_el: " << std::endl;
            std::cout << u_cat_el << std::endl;
            #endif
            uinc = u_cat_el;
        } else if (strcmp(solver, "elasto-plastic") == 0) {
            //Prepare for elastioPlaticIteration to solve ucat induced displacement
            #ifdef PRINT_INT_RESULTS
            std::cout << "elasto-plastic" << std::endl;
            #endif
            VectorXd Kstar = VectorXd::Zero(frame_system_size);
            MatrixXd Lstar = MatrixXd::Zero(frame_system_size, frame_system_size);
            Lstar(groundNodeDOF, groundNodeDOF) = FLEX_3DOF;
            MatrixXd Stiffness = KKfoot;
            for (int i = 0; i < groundNodeDOF.size(); i++) {
                Kstar(groundNodeDOF[i]) = 1.0/FLEX_3DOF(i,i); 
                Lstar(groundNodeDOF[i],groundNodeDOF[i]) = 0.0;
                Stiffness(groundNodeDOF[i],groundNodeDOF[i]) += Kstar(groundNodeDOF[i]);
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

            VectorXd patchArea = footing_ele_length_eigen.cwiseProduct(footing_ele_width_eigen);
            #ifdef PRINT_INT_RESULTS
            std::cout << "patchArea: " << std::endl;
            std::cout << patchArea << std::endl;
            #endif

            VectorXd ucat = VectorXd::Zero(frame_system_size);
            for (int i = 0; i < footing_nodes_size; i++) {
                int node_ind = footing_nodes_ind[i];
                ucat(node_ind * 6) = dispL[i];
                ucat(node_ind * 6 + 1) = dispT[i];
                ucat(node_ind * 6 + 2) = dispV[i];
            }

            elastoPlasticIterationLDLT(Stiffness, KKfoot, Kstar,
                FLEX_3DOF, Lstar, uinc, ucat, P_el,
                mu_int, lim_t_int, lim_c_int,
                patchArea,
                uip, groundNodeDOF1, groundNodeDOF2, groundNodeDOF3, groundNodeDOF);  
            #ifdef PRINT_INT_RESULTS
            /*std::cout << "uinc-u_P_el" << std::endl;
            std::cout << uinc-u_P_el << std::endl;*/
            #endif
        } else {
            std::cout << "Invalid solver type" << std::endl;
            return -1;
        }
        
        VectorXd result = uinc-u_P_el;
        std::copy(result.data(), result.data() + result.size(), result_array);
        return 0;
        // if (strcmp(output, "disp") == 0) {
        //     VectorXd result = uinc-u_P_el;
        //     // for (int i=0; i < result.size(); i++) {
        //     //     // std::cout << result(i) <<"," <<result.data()[i] << std::endl;
        //     //     result_array[i] = result(i);
        //     // }
        //     std::copy(result.data(), result.data() + result.size(), result_array);
        //     return 0;
        // } else if (strcmp(output, "strain") == 0) {
        //     VectorXd F_M_deltaT_el_M, F_N_deltaT_el_M, F_S_deltaT_el_M;
        //     calInternalForces(&F_M_deltaT_el_M, &F_N_deltaT_el_M, &F_S_deltaT_el_M,
        //             u_P_el, uinc, Eb, EoverG, h_el_foot, dfoot, bfoot, ni_foot, nnode);
        //     VectorXd epsilon_vector = calculateStrain(&F_S_deltaT_el_M, &F_M_deltaT_el_M,
        //             &F_N_deltaT_el_M, Eb, EoverG, bfoot, dfoot, ni_foot);
        //     VectorXd result = epsilon_vector;
        //     std::copy(result.data(), result.data() + result.size(), result_array);
        //     return 0;
        // } else if (strcmp(output, "strain+disp") == 0) {
        //     VectorXd F_M_deltaT_el_M, F_N_deltaT_el_M, F_S_deltaT_el_M;
        //     calInternalForces(&F_M_deltaT_el_M, &F_N_deltaT_el_M, &F_S_deltaT_el_M,
        //             u_P_el, uinc, Eb, EoverG, h_el_foot, dfoot, bfoot, ni_foot, nnode);
        //     VectorXd epsilon_vector = calculateStrain(&F_S_deltaT_el_M, &F_M_deltaT_el_M,
        //             &F_N_deltaT_el_M, Eb, EoverG, bfoot, dfoot, ni_foot);
        //     VectorXd strain = epsilon_vector;
        //     VectorXd disp = uinc - u_P_el;
        //     VectorXd result(strain.size() + disp.size());
        //     result << strain, disp;
        //     std::copy(result.data(), result.data() + result.size(), result_array);
        //     return 0;
        // } else {
        //     return -1;
        // }
        
        




        /*MatrixXd KKextra;
        MatrixXd KKpg = calKKpg(&KKsoil_rid, kv_gazetas, KKfoot, &KKextra, nnode);

        MatrixXd SS = KKfoot + KKextra;
        MatrixXd lamdasts, lamdastd, Kst;
        calLamdastsLamdastdKst(&lamdasts, &lamdastd, &Kst, FLEX, nnode);*/
    }
}

