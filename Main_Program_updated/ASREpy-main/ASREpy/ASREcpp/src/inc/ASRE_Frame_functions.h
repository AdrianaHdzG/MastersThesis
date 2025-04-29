#pragma once
#ifndef ASRE_Frame_functions
#define ASRE_Frame_functions
#define _USE_MATH_DEFINES

#include "Eigen/Dense"
#include "Eigen/LU"
#include "json.hpp"
#include <math.h>
#include <fstream>
#include <string>
#include <iostream>
#include <chrono>
#include "Eigen/SparseCore"
using namespace Eigen;
// using json = nlohmann::json;
using namespace std::chrono;

// #if defined _MSVC //Define the macro and DLLMain needed on Windows
#if defined(_MSVC)
    // #include <windows.h>
    #define DLLEXPORT __declspec(dllexport)
    // BOOL APIENTRY DllMain( HMODULE hModule,
    //                     DWORD  ul_reason_for_call,
    //                     LPVOID lpReserved
    //                     )
    // {
    //     switch (ul_reason_for_call)
    //     {
    //     case DLL_PROCESS_ATTACH:
    //     case DLL_THREAD_ATTACH:
    //     case DLL_THREAD_DETACH:
    //     case DLL_PROCESS_DETACH:
    //         break;
    //     }
    //     return TRUE;
    // }
#elif defined(_GNC)
    //  GCC
    #define DLLEXPORT __attribute__((visibility("default")))
    #define DLLIMPORT
#else
    #define DLLEXPORT
    #define DLLIMPORT
    #pragma warning Unknown dynamic link import/export semantics.
#endif


MatrixXd KBern3D_foot_TIM(double E, double d_foot, double b_foot, double dx, double EGratio,
        double ni_str) {
        double A = b_foot * d_foot;
        double a = std::max(b_foot, d_foot) / 2;
        double b = std::min(b_foot, d_foot) / 2;
        // double I11 = a * (b * b * b) * ((double)16 / 3 - 3.36 * b / a * (1 - (b * b * b * b) / 12 / (a * a * a * a)));
        double I22 = b_foot * pow(d_foot, 3) / 12;
        double I33 = d_foot * pow(b_foot, 3) / 12; 
        double I11 = I22 + I33;
        Vector3d Xi(0, 0, 0);
        Vector3d Xf(dx, 0, 0);
        double L = (Xi - Xf).norm();
        double k = 10 * (1 + ni_str) / (12 + 11 * ni_str);
        double G = E / EGratio;
        double phi2 = 12 * E * I22 / k / G / A / (L * L);
        double phi3 = 12 * E * I33 / k / G / A / (L * L);
        double phi2Bar = 1/(1+phi2);
        double phi3Bar = 1/(1+phi3);
        Vector3d Z(0, 0, 1);
        Vector3d x1 = (Xf - Xi) / L;
        Vector3d x2 = Z.cross(x1) / (Z.cross(x1).norm());
        Vector3d x3 = x1.cross(x2) / (x1.cross(x2).norm());
        MatrixXd Rot(3, 3);
        Rot << x1(0), x1(1), x1(2),
            x2(0), x2(1), x2(2),
            x3(0), x3(1), x3(2);
        MatrixXd MatRot = MatrixXd::Zero(12, 12);
        MatRot.block(0, 0, 3, 3) = Rot;
        MatRot.block(3, 3, 3, 3) = Rot;
        MatRot.block(6, 6, 3, 3) = Rot;
        MatRot.block(9, 9, 3, 3) = Rot;

        MatrixXd K = MatrixXd::Zero(12, 12);
        K(0, 0) = E * A / L;
        K(0, 6) = (-1) * E * A / L;
        K(1, 1) = 12 * E * I33 * phi3Bar/ pow(L, 3);
        K(1, 5) = 6 * E * I33 *phi3Bar/ pow(L, 2);
        K(1, 7) = -12 * E * I33 *phi3Bar/ pow(L, 3);
        K(1, 11) = 6 * E * I33 *phi3Bar/ pow(L, 2);
        K(2, 2) = 12 * E * I22 / pow(L, 3) * phi2Bar;
        K(2, 4) = -6 * E * I22 / pow(L, 2) * phi2Bar;
        K(2, 8) = -12 * E * I22 / pow(L, 3) * phi2Bar;
        K(2, 10) = -6 * E * I22 / pow(L, 2) * phi2Bar;
        K(3, 3) = G * I11 / L;
        K(3, 9) = (-1) * G * I11 / L;
        K(4, 2) = -6 * E * I22 / pow(L, 2) * phi2Bar;
        K(4, 4) = (4 + phi2) * E * I22 / L * phi2Bar;
        K(4, 8) = 6 * E * I22 / pow(L, 2) * phi2Bar;
        K(4, 10) = (2 - phi2) * E * I22 / L * phi2Bar;
        K(5, 1) = 6 * E * I33 *phi3Bar/ pow(L, 2);
        K(5, 5) = (4 + phi3) * phi3Bar * E * I33 / L;
        K(5, 7) = -6 * phi3Bar *E * I33 / pow(L, 2);
        K(5, 11) = (2 - phi3) * phi3Bar * E * I33 / L;
        K(6, 0) = -1 * E * A / L;
        K(6, 6) = E * A / L;
        K(7, 1) = -12 * E * I33 * phi3Bar/ pow(L, 3);
        K(7, 5) = -6 * E * I33 * phi3Bar/ pow(L, 2);
        K(7, 7) = 12 * E * I33 * phi3Bar/ pow(L, 3);
        K(7, 11) = -6 * E * I33 *phi3Bar/ pow(L, 2);
        K(8, 2) = -12 * E * I22 * phi2Bar/ pow(L, 3);
        K(8, 4) = 6 * E * I22 / pow(L, 2) *phi2Bar;
        K(8, 8) = 12 * E * I22 / pow(L, 3) *phi2Bar;
        K(8, 10) = 6 * E * I22 / pow(L, 2) * phi2Bar;

        K(9, 3) = -1 * G * I11 / L;
        K(9, 9) = G * I11 / L;

        K(10, 2) = -6 * E * I22 / pow(L, 2) * phi2Bar;
        K(10, 4) = (2 - phi2) * E * I22 / L * phi2Bar;
        K(10, 8) = 6 * E * I22 / pow(L, 2) * phi2Bar;
        K(10, 10) = (4 + phi2) * E * I22 / L * phi2Bar;

        K(11, 1) = 6 * E *phi3Bar * I33 / pow(L, 2);
        K(11, 5) = (2 - phi3) * phi3Bar * E * I33 / L;
        K(11, 7) = -6 * phi3Bar * E * I33 / pow(L, 2);
        K(11, 11) = (4 + phi3) * E * phi3Bar * I33 / L;

        MatrixXd Kelem = MatRot.transpose() * K * MatRot;
        return Kelem;
    }

void f_int_mind_dz_p(const VectorXd& u, VectorXd v, VectorXd w, VectorXd x, VectorXd y,
    const VectorXd& z, VectorXd* fdz, VectorXd* fdx, int nnodes_foot, double nis) {
    VectorXd X = x - u;
    VectorXd Y = y - v;
    VectorXd Z1 = z - w;
    VectorXd Z2 = z + w;
    VectorXd R1 = (X.cwiseProduct(X) + Y.cwiseProduct(Y) + Z1.cwiseProduct(Z1)).cwiseSqrt();
    VectorXd R2 = (X.cwiseProduct(X) + Y.cwiseProduct(Y) + Z2.cwiseProduct(Z2)).cwiseSqrt();
    //VectorXd Tx1 = (Y + Z1 + R1).cwiseQuotient(X).array().atan();
    VectorXd Tx2 = (Y + Z2 + R2).cwiseQuotient(X).array().atan();
    //VectorXd Ty1 = (X + Z1 + R1).cwiseQuotient(Y).array().atan();
    //VectorXd Ty2 = (X + Z2 + R2).cwiseQuotient(Y).array().atan();
    VectorXd Tz1 = (X + Y + R1).cwiseQuotient(Z1).array().atan();
    VectorXd Tz2 = (X + Y + R2).cwiseQuotient(Z2).array().atan();

    //VectorXd Txy = Y.cwiseQuotient(X).array().atan();
    //VectorXd Tyx = X.cwiseQuotient(Y).array().atan();

    VectorXd fdz_result = (3 - 4 * nis) * (Y.cwiseProduct((VectorXd)(X + R1).array().log()) +
        X.cwiseProduct((VectorXd)(Y + R1).array().log()) + 2 * Z1.cwiseProduct(Tz1)
        - 2 * Z2.cwiseProduct(Tz2)) - 2 * Z1.cwiseProduct(Tz1) +
        (2 * X.cwiseProduct(Y).cwiseProduct(z).cwiseProduct(w)).cwiseQuotient(R2).cwiseProduct(
            VectorXd::Ones(nnodes_foot).cwiseQuotient(Y.cwiseProduct(Y) + Z2.cwiseProduct(Z2)) +
            VectorXd::Ones(nnodes_foot).cwiseQuotient(X.cwiseProduct(X) + Z2.cwiseProduct(Z2))
        ) + (8 * (1 - nis) * (1 - nis) - 3 + 4 * nis) *
        (Y.cwiseProduct((VectorXd)(X + R2).array().log()) +
            X.cwiseProduct((VectorXd)(Y + R2).array().log()) + 2 * Z2.cwiseProduct(Tz2));
    *fdz = fdz_result;

    VectorXd fdx_result = (4 * nis - 3) * Z1.cwiseProduct((VectorXd)(Y + R2).array().log()) -
        Z1.cwiseProduct((VectorXd)(Y + R1).array().log()) - (2 * w.cwiseProduct(z).cwiseProduct(Y).
            cwiseProduct(Z2)).cwiseQuotient(R2.cwiseProduct(X.cwiseProduct(X) + Z2.cwiseProduct(Z2))) -
        4 * (1 - nis) * (1 - 2 * nis) * (Y.cwiseProduct((VectorXd)(Z2 + R2).array().log()) +
            Z2.cwiseProduct((VectorXd)(Y + R2).array().log()) + 2 * X.cwiseProduct(Tx2));
    *fdx = fdx_result;

}

void f_int_mind_dz_q(const VectorXd& u, VectorXd v, VectorXd w, VectorXd x, VectorXd y,
    const VectorXd& z, VectorXd* fdz, VectorXd* fdx, double nis) {
    VectorXd X = x - u;
    VectorXd Y = y - v;
    VectorXd Z1 = z - w;
    VectorXd Z2 = z + w;
    VectorXd R1 = (X.cwiseProduct(X) + Y.cwiseProduct(Y) + Z1.cwiseProduct(Z1)).cwiseSqrt();
    VectorXd R2 = (X.cwiseProduct(X) + Y.cwiseProduct(Y) + Z2.cwiseProduct(Z2)).cwiseSqrt();
    //VectorXd Tx1 = (Y + Z1 + R1).cwiseQuotient(X).array().atan();
    VectorXd Tx2 = (Y + Z2 + R2).cwiseQuotient(X).array().atan();
    //VectorXd Ty1 = (X + Z1 + R1).cwiseQuotient(Y).array().atan();
    //VectorXd Ty2 = (X + Z2 + R2).cwiseQuotient(Y).array().atan();
    VectorXd Tz1 = (X + Y + R1).cwiseQuotient(Z1).array().atan();
    VectorXd Tz2 = (X + Y + R2).cwiseQuotient(Z2).array().atan();

    //VectorXd Txy = Y.cwiseQuotient(X).array().atan();
    //VectorXd Tyx = X.cwiseQuotient(Y).array().atan();

    VectorXd fdz_result = (4 * nis - 3) * Z1.cwiseProduct((VectorXd)(Y + R2).array().log()) -
        Z1.cwiseProduct((VectorXd)(Y + R1).array().log()) +
        (2 * w.cwiseProduct(z).cwiseProduct(Z2).cwiseProduct(Y))
        .cwiseQuotient((X.cwiseProduct(X) + Z2.cwiseProduct(Z2)).cwiseProduct(R2)) +
        4 * (1 - nis) * (1 - 2 * nis) * (Y.cwiseProduct((VectorXd)(Z2 + R2).array().log()) +
            Z2.cwiseProduct((VectorXd)(Y + R2).array().log()) + 2 * X.cwiseProduct(Tx2));

    *fdz = fdz_result;

    VectorXd fdx_result = (3 - 4 * nis) * (Y.cwiseProduct((VectorXd)(X + R1).array().log()) +
        X.cwiseProduct((VectorXd)(Y + R1).array().log()) + 2 * Z1.cwiseProduct(Tz1) +
        Y.cwiseProduct((VectorXd)(X + R2).array().log()) + 2 * Z2.cwiseProduct(Tz2)) +
        Y.cwiseProduct((VectorXd)(X + R2).array().log()) +
        X.cwiseProduct((VectorXd)(Y + R2).array().log()) +
        2 * Z2.cwiseProduct(Tz2) + Y.cwiseProduct((VectorXd)(X + R1).array().log()) +
        2 * Z1.cwiseProduct(Tz1) + (2 * w.cwiseProduct(z).cwiseProduct(X).cwiseProduct(Y)).cwiseQuotient(
            R2.cwiseProduct((X.cwiseProduct(X) + Z2.cwiseProduct(Z2)))
        ) + 4 * (1 - nis) * (1 - 2 * nis) * (
            X.cwiseProduct((VectorXd)(Y + R2).array().log()) - 2 * Z2.cwiseProduct(Tx2)
            );
    *fdx = fdx_result;

}

// Disp at point i when uniform force applied around point j
MatrixXd int_factor_mindlin_j_1_vect_Vasiri(const VectorXd& Xi, VectorXd Yi,
    VectorXd Zi, VectorXd Xj, VectorXd Yj, VectorXd Zj, double Gs, double Es,
    int nnodes_foot, double nis, double h_el_foot, double bfoot) {
    MatrixXd result(nnodes_foot, 3);
    Yi = -1 * Yi;
    Yj = -1 * Yj;
    double beta_soil = 8 * M_PI * Es * (1 - nis) / (1 + nis);
    VectorXd u1 = Xj - VectorXd::Ones(nnodes_foot) * h_el_foot / 2;
    VectorXd u2 = Xj + VectorXd::Ones(nnodes_foot) * h_el_foot / 2;
    VectorXd v1 = Yj - VectorXd::Ones(nnodes_foot) * bfoot / 2;
    VectorXd v2 = Yj + VectorXd::Ones(nnodes_foot) * bfoot / 2;

    VectorXd q = VectorXd::Ones(nnodes_foot).cwiseQuotient((u2 - u1).cwiseProduct(v2 - v1));

    VectorXd fdz_q_u1_v1 = VectorXd::Zero(nnodes_foot);
    VectorXd fdx_q_u1_v1 = VectorXd::Zero(nnodes_foot);
    VectorXd fdz_q_u1_v2 = VectorXd::Zero(nnodes_foot);
    VectorXd fdx_q_u1_v2 = VectorXd::Zero(nnodes_foot);
    VectorXd fdz_q_u2_v1 = VectorXd::Zero(nnodes_foot);
    VectorXd fdx_q_u2_v1 = VectorXd::Zero(nnodes_foot);
    VectorXd fdz_q_u2_v2 = VectorXd::Zero(nnodes_foot);
    VectorXd fdx_q_u2_v2 = VectorXd::Zero(nnodes_foot);

    f_int_mind_dz_q(u1, v1, Zi, Xi, Yi, Zi, &fdz_q_u1_v1, &fdx_q_u1_v1, nis);
    f_int_mind_dz_q(u1, v2, Zi, Xi, Yi, Zi, &fdz_q_u1_v2, &fdx_q_u1_v2, nis);
    f_int_mind_dz_q(u2, v1, Zi, Xi, Yi, Zi, &fdz_q_u2_v1, &fdx_q_u2_v1, nis);
    f_int_mind_dz_q(u2, v2, Zi, Xi, Yi, Zi, &fdz_q_u2_v2, &fdx_q_u2_v2, nis);


    VectorXd deltaz_q = (q / beta_soil).cwiseProduct((fdz_q_u2_v2 - fdz_q_u1_v2) - (fdz_q_u2_v1 - fdz_q_u1_v1));
    VectorXd deltax_q = (q / beta_soil).cwiseProduct((fdx_q_u2_v2 - fdx_q_u1_v2) - (fdx_q_u2_v1 - fdx_q_u1_v1));

    //cout<<(fdz_q_u2_v2-fdz_q_u1_v2)-(fdz_q_u2_v1-fdz_q_u1_v1)<<endl;
    //cout<<deltaz_q<<deltax_q<<endl;
    result.col(0) = deltax_q;
    result.col(1) = VectorXd::Zero(nnodes_foot);
    result.col(2) = deltaz_q;
    return result;
}

MatrixXd int_factor_mindlin_j_2_vect_Vasiri(VectorXd Xi, VectorXd Yi,
    VectorXd Zi, VectorXd Xj, VectorXd Yj, VectorXd Zj, double Gs, double Es,
    int nnodes_foot, double nis, double h_el_foot, double bfoot) {
    MatrixXd result(nnodes_foot, 3);
    double beta_soil = 8 * M_PI * Es * (1 - nis) / (1 + nis);
    VectorXd v1 = Xj - VectorXd::Ones(nnodes_foot) * h_el_foot / 2;
    VectorXd v2 = Xj + VectorXd::Ones(nnodes_foot) * h_el_foot / 2;
    VectorXd u1 = Yj - VectorXd::Ones(nnodes_foot) * bfoot / 2;
    VectorXd u2 = Yj + VectorXd::Ones(nnodes_foot) * bfoot / 2;

    VectorXd q = VectorXd::Ones(nnodes_foot).cwiseQuotient((u2 - u1).cwiseProduct(v2 - v1));

    VectorXd fdz_q_u1_v1 = VectorXd::Zero(nnodes_foot);
    VectorXd fdx_q_u1_v1 = VectorXd::Zero(nnodes_foot);
    VectorXd fdz_q_u1_v2 = VectorXd::Zero(nnodes_foot);
    VectorXd fdx_q_u1_v2 = VectorXd::Zero(nnodes_foot);
    VectorXd fdz_q_u2_v1 = VectorXd::Zero(nnodes_foot);
    VectorXd fdx_q_u2_v1 = VectorXd::Zero(nnodes_foot);
    VectorXd fdz_q_u2_v2 = VectorXd::Zero(nnodes_foot);
    VectorXd fdx_q_u2_v2 = VectorXd::Zero(nnodes_foot);

    f_int_mind_dz_q(u1, v1, Zi, Yi, Xi, Zi, &fdz_q_u1_v1, &fdx_q_u1_v1, nis);
    f_int_mind_dz_q(u1, v2, Zi, Yi, Xi, Zi, &fdz_q_u1_v2, &fdx_q_u1_v2, nis);
    f_int_mind_dz_q(u2, v1, Zi, Yi, Xi, Zi, &fdz_q_u2_v1, &fdx_q_u2_v1, nis);
    f_int_mind_dz_q(u2, v2, Zi, Yi, Xi, Zi, &fdz_q_u2_v2, &fdx_q_u2_v2, nis);

    VectorXd deltaz_q = (q / beta_soil).cwiseProduct((fdz_q_u2_v2 - fdz_q_u1_v2) - (fdz_q_u2_v1 - fdz_q_u1_v1));
    VectorXd deltax_q = (q / beta_soil).cwiseProduct((fdx_q_u2_v2 - fdx_q_u1_v2) - (fdx_q_u2_v1 - fdx_q_u1_v1));

    result.col(1) = deltax_q;
    result.col(0) = VectorXd::Zero(nnodes_foot);
    result.col(2) = deltaz_q;
    return result;
}

MatrixXd int_factor_mindlin_j_3_vect_Vasiri(VectorXd Xi, VectorXd Yi,
    VectorXd Zi, VectorXd Xj, VectorXd Yj, VectorXd Zj, double Gs, double Es,
    int nnodes_foot, double nis, double h_el_foot, double bfoot) {
    MatrixXd result(nnodes_foot, 3);
    Yi = -Yi;
    Yj = -Yj;
    double beta_soil = 8 * M_PI * Es * (1 - nis) / (1 + nis);
    VectorXd u1 = Xj - VectorXd::Ones(nnodes_foot) * h_el_foot / 2;
    VectorXd u2 = Xj + VectorXd::Ones(nnodes_foot) * h_el_foot / 2;
    VectorXd v1 = Yj - VectorXd::Ones(nnodes_foot) * bfoot / 2;
    VectorXd v2 = Yj + VectorXd::Ones(nnodes_foot) * bfoot / 2;

    VectorXd p = VectorXd::Ones(nnodes_foot).cwiseQuotient((u2 - u1).cwiseProduct(v2 - v1));

    VectorXd fdz_p_u1_v1 = VectorXd::Zero(nnodes_foot);
    VectorXd fdx_p_u1_v1 = VectorXd::Zero(nnodes_foot);
    VectorXd fdz_p_u1_v2 = VectorXd::Zero(nnodes_foot);
    VectorXd fdx_p_u1_v2 = VectorXd::Zero(nnodes_foot);
    VectorXd fdz_p_u2_v1 = VectorXd::Zero(nnodes_foot);
    VectorXd fdx_p_u2_v1 = VectorXd::Zero(nnodes_foot);
    VectorXd fdz_p_u2_v2 = VectorXd::Zero(nnodes_foot);
    VectorXd fdx_p_u2_v2 = VectorXd::Zero(nnodes_foot);

    f_int_mind_dz_p(u1, v1, Zi, Xi, Yi, Zi, &fdz_p_u1_v1, &fdx_p_u1_v1, nnodes_foot, nis);
    f_int_mind_dz_p(u1, v2, Zi, Xi, Yi, Zi, &fdz_p_u1_v2, &fdx_p_u1_v2, nnodes_foot, nis);
    f_int_mind_dz_p(u2, v1, Zi, Xi, Yi, Zi, &fdz_p_u2_v1, &fdx_p_u2_v1, nnodes_foot, nis);
    f_int_mind_dz_p(u2, v2, Zi, Xi, Yi, Zi, &fdz_p_u2_v2, &fdx_p_u2_v2, nnodes_foot, nis);

    VectorXd deltaz_p = (p / beta_soil).cwiseProduct((fdz_p_u2_v2 - fdz_p_u1_v2) - (fdz_p_u2_v1 - fdz_p_u1_v1));
    VectorXd deltax_p = (p / beta_soil).cwiseProduct((fdx_p_u2_v2 - fdx_p_u1_v2) - (fdx_p_u2_v1 - fdx_p_u1_v1));

    result.col(0) = deltax_p;
    result.col(1) = VectorXd::Zero(nnodes_foot);
    result.col(2) = deltaz_p;
    return result;
}

MatrixXd _calFlexVaziri(VectorXd Zi_nod, VectorXd Xi_nod, VectorXd Yi_nod, double Es,
    int nnodes, double* h_el_foot, double nis, double *bfoot) {
    MatrixXd FLEX = MatrixXd::Zero(nnodes * 3, nnodes * 3);
    double Gs = Es / 2 / (1 + nis);
    for (int i = 0; i < nnodes; i++) {
        VectorXd Zj_nod = VectorXd::Ones(nnodes) * Zi_nod(i);
        VectorXd Xj_nod = VectorXd::Ones(nnodes) * Xi_nod(i);
        VectorXd Yj_nod = VectorXd::Ones(nnodes) * Yi_nod(i);

        MatrixXd uij_1 = int_factor_mindlin_j_1_vect_Vasiri(Xi_nod, Yi_nod, Zi_nod,
            Xj_nod, Yj_nod, Zj_nod, Gs, Es, nnodes, nis, h_el_foot[i], bfoot[i]);
        MatrixXd uij_2 = int_factor_mindlin_j_2_vect_Vasiri(Xi_nod, Yi_nod, Zi_nod,
            Xj_nod, Yj_nod, Zj_nod, Gs, Es, nnodes, nis, h_el_foot[i], bfoot[i]);
        MatrixXd uij_3 = int_factor_mindlin_j_3_vect_Vasiri(Xi_nod, Yi_nod, Zi_nod,
            Xj_nod, Yj_nod, Zj_nod, Gs, Es, nnodes, nis, h_el_foot[i], bfoot[i]);
        //cout<<uij_3<<endl;
        FLEX(seq(0, nnodes * 3 - 1, 3), i * 3) = uij_1.col(0);
        FLEX(seq(1, nnodes * 3 - 1, 3), i * 3) = uij_1.col(1);
        FLEX(seq(2, nnodes * 3 - 1, 3), i * 3) = uij_1.col(2);

        FLEX(seq(0, nnodes * 3 - 1, 3), i * 3 + 1) = uij_2.col(0);
        FLEX(seq(1, nnodes * 3 - 1, 3), i * 3 + 1) = uij_2.col(1);
        FLEX(seq(2, nnodes * 3 - 1, 3), i * 3 + 1) = uij_2.col(2);

        FLEX(seq(0, nnodes * 3 - 1, 3), i * 3 + 2) = uij_3.col(0);
        FLEX(seq(1, nnodes * 3 - 1, 3), i * 3 + 2) = uij_3.col(1);
        FLEX(seq(2, nnodes * 3 - 1, 3), i * 3 + 2) = uij_3.col(2);
    }
    double r1 = -0.3136 * exp(-0.531 * (h_el_foot[0] / bfoot[0])) + 0.2174 * pow(h_el_foot[0] / bfoot[0], 0.1755) + 1.327;
    double r2 = -0.3136 * exp(-0.531 * (bfoot[0] / h_el_foot[0])) + 0.2174 * pow(bfoot[0] / h_el_foot[0], 0.1755) + 1.327;
    double r3 = r1;

    FLEX(0, 0) *= r1;
    FLEX(1, 1) *= r2;
    FLEX(2, 2) *= r3;

    r1 = -0.3136 * exp(-0.531 * (h_el_foot[nnodes - 1] / bfoot[nnodes - 1])) + 0.2174 * pow(h_el_foot[nnodes - 1] / bfoot[nnodes - 1], 0.1755) + 1.327;
    r2 = -0.3136 * exp(-0.531 * (bfoot[nnodes - 1] / h_el_foot[nnodes - 1])) + 0.2174 * pow(bfoot[nnodes - 1] / h_el_foot[nnodes - 1], 0.1755) + 1.327;
    r3 = r1;

    FLEX((nnodes - 1) * 3, (nnodes - 1) * 3) *= r1;
    FLEX(1 + (nnodes - 1) * 3, 1 + (nnodes - 1) * 3) *= r2;
    FLEX(2 + (nnodes - 1) * 3, 2 + (nnodes - 1) * 3) *= r3;
    return FLEX;
}

extern "C" {
    // Export calFlex
    DLLEXPORT double* calFlexVaziri(double* Zi, double* Xi, double* Yi, double Es,
        int nnodes, double* h_el_foot, double nis, double* bfoot) {

        // Convert array to vector
        VectorXd Zi_nod = Eigen::Map<VectorXd>(Zi, nnodes);
        VectorXd Xi_nod = Eigen::Map<VectorXd>(Xi, nnodes);
        VectorXd Yi_nod = Eigen::Map<VectorXd>(Yi, nnodes);

        MatrixXd FLEX = MatrixXd::Zero(nnodes * 3, nnodes * 3);
        double Gs = Es / 2 / (1 + nis);
        for (int i = 0; i < nnodes; i++) {
            VectorXd Zj_nod = VectorXd::Ones(nnodes) * Zi_nod(i);
            VectorXd Xj_nod = VectorXd::Ones(nnodes) * Xi_nod(i);
            VectorXd Yj_nod = VectorXd::Ones(nnodes) * Yi_nod(i);

            MatrixXd uij_1 = int_factor_mindlin_j_1_vect_Vasiri(Xi_nod, Yi_nod, Zi_nod,
                Xj_nod, Yj_nod, Zj_nod, Gs, Es, nnodes, nis, h_el_foot[i], bfoot[i]);
            MatrixXd uij_2 = int_factor_mindlin_j_2_vect_Vasiri(Xi_nod, Yi_nod, Zi_nod,
                Xj_nod, Yj_nod, Zj_nod, Gs, Es, nnodes, nis, h_el_foot[i], bfoot[i]);
            MatrixXd uij_3 = int_factor_mindlin_j_3_vect_Vasiri(Xi_nod, Yi_nod, Zi_nod,
                Xj_nod, Yj_nod, Zj_nod, Gs, Es, nnodes, nis, h_el_foot[i], bfoot[i]);
            //cout<<uij_3<<endl;
            FLEX(seq(0, nnodes * 3 - 1, 3), i * 3) = uij_1.col(0);
            FLEX(seq(1, nnodes * 3 - 1, 3), i * 3) = uij_1.col(1);
            FLEX(seq(2, nnodes * 3 - 1, 3), i * 3) = uij_1.col(2);

            FLEX(seq(0, nnodes * 3 - 1, 3), i * 3 + 1) = uij_2.col(0);
            FLEX(seq(1, nnodes * 3 - 1, 3), i * 3 + 1) = uij_2.col(1);
            FLEX(seq(2, nnodes * 3 - 1, 3), i * 3 + 1) = uij_2.col(2);

            FLEX(seq(0, nnodes * 3 - 1, 3), i * 3 + 2) = uij_3.col(0);
            FLEX(seq(1, nnodes * 3 - 1, 3), i * 3 + 2) = uij_3.col(1);
            FLEX(seq(2, nnodes * 3 - 1, 3), i * 3 + 2) = uij_3.col(2);
        }
        // double r1 = -0.3136 * exp(-0.531 * (h_el_foot / bfoot)) + 0.2174 * pow(h_el_foot / bfoot, 0.1755) + 1.327;
        // double r2 = -0.3136 * exp(-0.531 * (bfoot / h_el_foot)) + 0.2174 * pow(bfoot / h_el_foot, 0.1755) + 1.327;
        // double r3 = r1;

        // FLEX(0, 0) *= r1;
        // FLEX(1, 1) *= r2;
        // FLEX(2, 2) *= r3;
        // FLEX((nnodes - 1) * 3, (nnodes - 1) * 3) *= r1;
        // FLEX(1 + (nnodes - 1) * 3, 1 + (nnodes - 1) * 3) *= r2;
        // FLEX(2 + (nnodes - 1) * 3, 2 + (nnodes - 1) * 3) *= r3;
        std::cout<<FLEX<<std::endl;
        return FLEX.data();
    }
}    

void soilsprings_static_foot_gazetas(double* ptr_kh_gazetas, double* ptr_kv_gazetas, double Es, double nis, double bfoot, double h_el_foot) {
    double Gs = Es / 2 / (1 + nis);
    double ll = std::max(bfoot * 0.5, h_el_foot * 0.5);
    double bb = std::min(bfoot * 0.5, h_el_foot * 0.5);
    *ptr_kv_gazetas = 2 * Gs * ll / (1 - nis) * (0.73 + 1.54 * pow(ll * bb / ll / ll, 0.75));
    *ptr_kh_gazetas = 2 * Gs * ll / (2 - nis) * (2 + 2.5 * pow(ll * bb / ll / ll, 0.85));
}

MatrixXd calKKpg(MatrixXd* ptr_KKsoil_rid, double kv, MatrixXd Kfoot, MatrixXd* ptr_KKextra, int nnodes_foot) {
    MatrixXd KKsoil = MatrixXd::Zero(nnodes_foot * 6, nnodes_foot * 6);
    for (int i = 0; i < nnodes_foot; i++) {
        for (int j = 0; j < nnodes_foot; j++) {
            KKsoil(seq(i * 6, 2 + i * 6, 1), seq(j * 6, 2 + j * 6, 1)) =
                (*ptr_KKsoil_rid)(seq(i * 3, 2 + i * 3, 1), seq(j * 3, 2 + j * 3, 1));
        }
    }

    *ptr_KKextra = MatrixXd::Zero(nnodes_foot * 6, nnodes_foot * 6);
    for (int i = 1; i < 6 * nnodes_foot; i += 6) {
        (*ptr_KKextra)(i, i) = kv * 1000;
        (*ptr_KKextra)(i + 2, i + 2) = kv * 1000;
    }

    return Kfoot + KKsoil + (*ptr_KKextra);
}

void calc_P_el(VectorXd* ptr_P_el_pass, int nnodes_foot, double qz_foot, double h_el_foot) {
    // Point load
    /*if (nelementi_foot % 2 == 0) {
        (*ptr_P_el_pass)(nelementi_foot / 2 * 6) = Fx_foot_centr;
        (*ptr_P_el_pass)(2 + nelementi_foot / 2 * 6) = Fz_foot_centr;
        (*ptr_P_el_pass)(4 + nelementi_foot / 2 * 6) = My_foot_centr;
    }
    else {
        (*ptr_P_el_pass)(nelementi_foot / 2 * 6) = Fx_foot_centr * 0.5;
        (*ptr_P_el_pass)(2 + nelementi_foot / 2 * 6) = Fz_foot_centr * 0.5;
        (*ptr_P_el_pass)(4 + nelementi_foot / 2 * 6) = My_foot_centr * 0.5;
        (*ptr_P_el_pass)(nelementi_foot / 2 * 6 + 6) = Fx_foot_centr * 0.5;
        (*ptr_P_el_pass)(2 + nelementi_foot / 2 * 6 + 6) = Fz_foot_centr * 0.5;
        (*ptr_P_el_pass)(4 + nelementi_foot / 2 * 6 + 6) = My_foot_centr * 0.5;
    }*/
    // Uniformly distributed load
    int nelementi_foot = nnodes_foot - 1;
    double qx_foot=0; // Assume only vertical dead load
    VectorXd q_ind = VectorXd::Ones(nnodes_foot * 6);
    q_ind(0) = 0.5;
    q_ind(2) = 0.5;
    q_ind(nelementi_foot * 6) = 0.5;
    q_ind(2 + nelementi_foot * 6) = 0.5;
    // Reduce the load at each node is q_foot * h_el_foot, the load at two ends are half;
    (*ptr_P_el_pass)(seq(0, ptr_P_el_pass->size() - 1, 6)) = -1 * (
        q_ind(seq(0, ptr_P_el_pass->size() - 1, 6)) * qx_foot * h_el_foot
        + (*ptr_P_el_pass)(seq(0, ptr_P_el_pass->size() - 1, 6)));
    (*ptr_P_el_pass)(seq(2, ptr_P_el_pass->size() - 1, 6)) = -1 * (
        q_ind(seq(2, ptr_P_el_pass->size() - 1, 6)) * qz_foot * h_el_foot
        + (*ptr_P_el_pass)(seq(2, ptr_P_el_pass->size() - 1, 6)));
}

void calLamdastsLamdastdKst(MatrixXd* ptr_lamdasts, MatrixXd* ptr_lamdastd, MatrixXd* ptr_Kst, MatrixXd F, int nnodes_foot) {
    MatrixXd lamdasts_rid = F;
    lamdasts_rid.diagonal() -= lamdasts_rid.diagonal();
    MatrixXd lamdastd_rid = MatrixXd::Zero(F.rows(), F.cols());
    lamdastd_rid.diagonal() = F.diagonal();

    MatrixXd kst_rid = lamdastd_rid.inverse();
    (*ptr_Kst) = MatrixXd::Zero(6 * nnodes_foot, 6 * nnodes_foot);
    (*ptr_lamdasts) = MatrixXd::Zero(6 * nnodes_foot, 6 * nnodes_foot);
    (*ptr_lamdastd) = MatrixXd::Zero(6 * nnodes_foot, 6 * nnodes_foot);
    for (int i = 0; i < nnodes_foot; i++) {
        for (int j = 0; j < nnodes_foot; j++) {
            (*ptr_Kst)(seq(6 * i, 2 + 6 * i, 1), seq(j * 6, 2 + j * 6, 1)) =
                kst_rid(seq(i * 3, 2 + i * 3, 1), seq(j * 3, 2 + j * 3, 1));
            (*ptr_lamdasts)(seq(6 * i, 2 + 6 * i, 1), seq(j * 6, 2 + j * 6, 1)) =
                lamdasts_rid(seq(i * 3, 2 + i * 3, 1), seq(j * 3, 2 + j * 3, 1));
            (*ptr_lamdastd)(seq(6 * i, 2 + 6 * i, 1), seq(j * 6, 2 + j * 6, 1)) =
                lamdastd_rid(seq(i * 3, 2 + i * 3, 1), seq(j * 3, 2 + j * 3, 1));
        }
    }
}

void script_ep_solutions_freeFooting_exter(VectorXd* ptr_P_el,
    MatrixXd* ptr_SS, MatrixXd* ptr_Kst, MatrixXd* ptr_lamdasts, MatrixXd* ptr_lamdastd,
    MatrixXd* ptr_vect_UINC, MatrixXd* ptr_vect_dUIP, double mu, int nnodes_foot) {
    // solution for External load and self weight
    // cout<<"ep T1"<<ptr_P_el->rows()<<endl;
    std::string log_file = "C:\\Users\\Jinyan\\Desktop\\NGI_GIBV\\ASRETimoshenko\\ASRETimoshenko\\log_file.txt";
    double flim2 = 0;
    double flim1 = INFINITY;
    double eps_err = 0.05;
    int nelementi_foot = nnodes_foot - 1;
    int nfoot = 1;
    int num_incr_P = 10;
    VectorXd UINC_iprevious(6 * nfoot * (nelementi_foot + 1));
    //VectorXd UINCIP_iprevious(6 * nfoot * (nelementi_foot + 1));
    VectorXd P_elINC_iprevious(6 * nfoot * (nelementi_foot + 1));
    //MatrixXd dUCAT_vect = diff( *ptr_dfft_gruppo_vect);
    //MatrixXd dUCAT_vect_memory = dUCAT_vect;
    if (!(*ptr_P_el).array().isZero()) {
        int i_inc = 0;
        UINC_iprevious = VectorXd::Zero(6 * nfoot * nnodes_foot);
        P_elINC_iprevious = VectorXd::Zero(6 * nfoot * nnodes_foot);
        //dUCAT_vect = VectorXd::Zero(dUCAT_vect.size());
        VectorXd dP_el = (*ptr_P_el) / num_incr_P;
        // cout<<"ep T2"<<dP_el.rows()<<endl;
        VectorXd dU = VectorXd::Zero(ptr_SS->rows());
        MatrixXd A = *ptr_SS + *ptr_Kst + (*ptr_Kst) * (*ptr_lamdasts) * (*ptr_SS);
        PartialPivLU<MatrixXd> lu = PartialPivLU<MatrixXd>(A);
        for (int i = 0; i < num_incr_P; i++) {
            VectorXd dUIP = VectorXd::Zero(6 * nfoot * nnodes_foot);
            VectorXd dUCAP = VectorXd::Zero(6 * nfoot * nnodes_foot);
            //VectorXd dUCAT = VectorXd::Zero(6 * nfoot * nnodes_foot);
            // cout<<"ep T3"<<dUIP.rows() <<endl;
            int n_while = 0;
            while (true) {
                n_while++;
                VectorXd b = dP_el + (*ptr_Kst) * (*ptr_lamdasts) * (dP_el)+*ptr_Kst * dUIP;
                // cout<<"ep T4"<<dU.rows() <<endl;
                dU = lu.solve(b);
                // cout<<"ep T5"<<dU.rows() <<endl;
                VectorXd UINC = UINC_iprevious + dU;
                VectorXd P_elINC = P_elINC_iprevious + dP_el;
                VectorXd dF = (-1) * (*ptr_SS) * dU + dP_el;
                //cout<<"dF" << dF <<endl;
                VectorXd react_iprevious = (-1) * (*ptr_SS) * UINC_iprevious + P_elINC_iprevious;
                VectorXd react_inc = (-1) * (*ptr_SS) * UINC + P_elINC;
                VectorXd react_iprevious_3dof = react_iprevious(seq(2, react_iprevious.size() - 1, 6));
                VectorXd react_inc_3dof = react_inc(seq(2, react_inc.size() - 1, 6));
                VectorXd react_iprevious_1dof = react_iprevious(seq(0, react_iprevious.size() - 1, 6));
                VectorXd react_inc_1dof = react_inc(seq(0, react_inc.size() - 1, 6));

                VectorXd flim2vect = flim2 * VectorXd::Ones(nfoot * nnodes_foot);
                VectorXd flim1vect = flim1 * VectorXd::Ones(nfoot * nnodes_foot);
                VectorXd dF_3dof = dF(seq(2, dF.size() - 1, 6));
                VectorXd dF_1dof = dF(seq(0, dF.size() - 1, 6));


                VectorXd flim2Mreact_ipre = flim2vect - react_iprevious_3dof;
                VectorXd flim1Mreact_ipre = flim1vect - react_iprevious_3dof;

                dF_3dof = (react_inc_3dof.array() < flim2).select(flim2Mreact_ipre, dF_3dof);
                dF_3dof = (react_inc_3dof.array() > flim1).select(flim1Mreact_ipre, dF_3dof);
                dF(seq(2, dF.size() - 1, 6)) = dF_3dof;
                if (i == 0) {
                    VectorXd fmax_1dof = mu * (react_inc_3dof.array() > 0).select(react_inc_3dof, 0);
                    //auto ind_ip_1_up = ((react_inc_1dof - fmax_1dof).array() > 0 );
                    //auto ind_ip_1_dw = ((react_inc_1dof  + fmax_1dof).array() < 0 );
                    dF_1dof = ((react_inc_1dof - fmax_1dof).array() > 0).select(0, dF_1dof);
                    dF_1dof = ((react_inc_1dof + fmax_1dof).array() < 0).select(0, dF_1dof);
                    dF(seq(0, dF.size() - 1, 6)) = dF_1dof;
                    //cout<<"dF_1dof"<<dF_1dof<<endl;

                }
                else {
                    //double fmax_1dof = (mu * (react_inc_3dof.array() > 0).select(react_inc_3dof, 0)).maxCoeff();
                    VectorXd fmax_1dof_vect = mu * (react_inc_3dof.array() > 0).select(react_inc_3dof, 0);
                    VectorXd fmaxMinusReact = fmax_1dof_vect - react_iprevious_1dof;
                    VectorXd minusFmaxMinusReact = (-1) * fmax_1dof_vect - react_iprevious_1dof;
                    dF_1dof = (react_inc_1dof.array() > fmax_1dof_vect.array()).select(fmaxMinusReact, dF_1dof);
                    dF_1dof = (react_inc_1dof.array() < (-1) * fmax_1dof_vect.array()).select(minusFmaxMinusReact, dF_1dof);
                    dF(seq(0, dF.size() - 1, 6)) = dF_1dof;
                }
                //cout<<dF<<endl;
                /*VectorXd CI = (((*ptr_lamdasts) + (*ptr_lamdastd)) * dF  + dUIP).cwiseQuotient(dU);
                CI = (CI.array().isNaN()).select(1, CI);
                for (int i = 0; i < seq(1, CI.size()-1, 6).size(); i++) {
                    CI(seq(1, CI.size()-1, 6)[i]) = 1;
                    CI(seq(1, CI.size()-1, 6)[i] + 2) = 1;
                    CI(seq(1, CI.size()-1, 6)[i] + 3) = 1;
                    CI(seq(1, CI.size()-1, 6)[i] + 4) = 1;
                }*/
                VectorXd dUCAP = (*ptr_lamdasts) * dF;
                VectorXd dUIP_old = dUIP;
                VectorXd dUd = (*ptr_lamdastd) * dF;
                //cout<<dUIP<<endl;
                dUIP(seq(2, dUIP.size() - 1, 6)) = dU(seq(2, dUIP.size() - 1, 6)) - dUd(seq(2, dUIP.size() - 1, 6))
                    - dUCAP(seq(2, dUIP.size() - 1, 6));
                dUIP(seq(0, dUIP.size() - 1, 6)) = dU(seq(0, dUIP.size() - 1, 6)) - dUd(seq(0, dUIP.size() - 1, 6))
                    - dUCAP(seq(0, dUIP.size() - 1, 6));

                //dUCAP = ((*ptr_lamdasts) * dF + beta * dUCAP)/(1 + beta);
                //cout<<dUIP<<endl;
                VectorXd perr = (dUIP - dUIP_old).cwiseQuotient(dUIP);
                perr = perr.array().isNaN().select(0, perr);

                //cout<<"" <<endl<<(CI - VectorXd::Ones(CI.size())).cwiseAbs().maxCoeff()<<endl;
                //cout<<"" <<endl<<perr.cwiseAbs().maxCoeff()<<endl;
                if (n_while > 2000) {
                    std::cout << "external while loops more than 2000 times" << std::endl;
                    break;
                }
                if (perr.cwiseAbs().maxCoeff() < eps_err) {
                    i_inc += 1;
                    (*ptr_vect_UINC).col(i_inc - 1) = UINC;
                    (*ptr_vect_dUIP).col(i_inc - 1) = dUIP;
                    UINC_iprevious = UINC;
                    P_elINC_iprevious = P_elINC;
                    //cout<<"i_inc is " << i_inc << " loop finished" <<endl;
//                    cout<<n_while<<endl;
                    break;
                }
            }
        }
    }
    else {
        *ptr_vect_UINC = (*ptr_vect_UINC).block(0, 0, ptr_vect_UINC->rows(), 1); //use vect_UINC as vect_Utotf_P and vect_Utotf_free_ep_P
        *ptr_vect_dUIP = (*ptr_vect_dUIP).block(0, 0, ptr_vect_dUIP->rows(), 1);
    }
}


void script_ep_solutions_freeFooting_tunnel(MatrixXd* ptr_dfft_gruppo_vect,
    VectorXd* ptr_P_el, MatrixXd* ptr_vect_UINC, MatrixXd* ptr_vect_dUIP,
    MatrixXd* ptr_SS, MatrixXd* ptr_Kst, MatrixXd* ptr_lamdasts, MatrixXd* ptr_lamdastd
    , MatrixXd* vect_Utotf_P, double dvm, double mu, int nnodes_foot) {
    //cout<<(*ptr_vect_UINC)<<endl;
    double flim2 = 0;
    double flim1 = INFINITY;
    double eps_err = 0.05;
    int nelementi_foot = nnodes_foot - 1;
    int nfoot = 1;
    int num_incr_TUNNEL = 20;
    if (dvm != 0) {
        int i_inc = 0;
        VectorXd UINC_iprevious = (*vect_Utotf_P).col(vect_Utotf_P->cols() - 1);
        //cout<<UINC_iprevious<<endl;
        VectorXd P_elINC_iprevious = (*ptr_P_el);
        MatrixXd dUCAT_vect = (*ptr_dfft_gruppo_vect)/num_incr_TUNNEL;
        VectorXd dP_el = VectorXd::Zero(6 * nfoot * nnodes_foot);
        VectorXd dU = VectorXd::Zero(ptr_SS->rows());
        MatrixXd A = *ptr_SS + *ptr_Kst + (*ptr_Kst) * (*ptr_lamdasts) * (*ptr_SS);
        PartialPivLU<MatrixXd> lu = PartialPivLU<MatrixXd>(A);
        for (int ii = 0; ii < num_incr_TUNNEL; ii++) {
            VectorXd dUIP = VectorXd::Zero(6 * nfoot * nnodes_foot);
            VectorXd dUCAP = VectorXd::Zero(6 * nfoot * nnodes_foot);
            VectorXd dUCAT = dUCAT_vect.col(ii);
            int n_while = 0;
            while (true) {
                n_while++;
                //                VectorXd dU = (*ptr_SS + *ptr_Kst + (*ptr_Kst)*(*ptr_lamdasts)*(*ptr_SS)).
                //                        partialPivLu().solve(
                //                        dP_el + *ptr_Kst*((*ptr_lamdasts)*dP_el + dUCAT + dUIP)
                //                );
                
                VectorXd b = dP_el + *ptr_Kst * ((*ptr_lamdasts) * dP_el + dUCAT + dUIP);
                dU = lu.solve(b);
                VectorXd UINC = UINC_iprevious + dU;
                VectorXd P_elINC = P_elINC_iprevious + dP_el;
                VectorXd dF = (-1) * (*ptr_SS) * dU + dP_el;
                VectorXd react_iprevious = (-1) * (*ptr_SS) * UINC_iprevious + P_elINC_iprevious;
                VectorXd react_inc = (-1) * (*ptr_SS) * UINC + P_elINC;
                VectorXd react_iprevious_3dof = react_iprevious(seq(2, react_iprevious.size() - 1, 6));
                VectorXd react_inc_3dof = react_inc(seq(2, react_inc.size() - 1, 6));
                VectorXd react_iprevious_1dof = react_iprevious(seq(0, react_iprevious.size() - 1, 6));
                VectorXd react_inc_1dof = react_inc(seq(0, react_inc.size() - 1, 6));
                VectorXd flim2vect = flim2 * VectorXd::Ones(nfoot * nnodes_foot);
                VectorXd flim1vect = flim1 * VectorXd::Ones(nfoot * nnodes_foot);
                VectorXd dF_3dof = dF(seq(2, dF.size() - 1, 6));
                VectorXd dF_1dof = dF(seq(0, dF.size() - 1, 6));
                VectorXd flim2Mreact_ipre = flim2vect - react_iprevious_3dof;
                VectorXd flim1Mreact_ipre = flim1vect - react_iprevious_3dof;
                dF_3dof = (react_inc_3dof.array() < flim2).select(flim2Mreact_ipre, dF_3dof);
                dF_3dof = (react_inc_3dof.array() > flim1).select(flim1Mreact_ipre, dF_3dof);
                dF(seq(2, dF.size() - 1, 6)) = dF_3dof;

                if (ii == 0) {
                    VectorXd fmax_1dof = mu * (react_inc_3dof.array() > 0).select(react_inc_3dof, 0);
                    //auto ind_ip_1_up = ((react_inc_1dof - fmax_1dof).array() > 0 );
                    //auto ind_ip_1_dw = ((react_inc_1dof  + fmax_1dof).array() < 0 );
                    dF_1dof = ((react_inc_1dof - fmax_1dof).array() > 0).select(0, dF_1dof);
                    dF_1dof = ((react_inc_1dof + fmax_1dof).array() < 0).select(0, dF_1dof);
                    dF(seq(0, dF.size() - 1, 6)) = dF_1dof;
                    //cout<<"dF_1dof"<<dF_1dof<<endl;

                }
                else {
                    VectorXd fmax_1dof_vect = mu * (react_inc_3dof.array() > 0).select(react_inc_3dof, 0);
                    VectorXd fmaxMinusReact = fmax_1dof_vect - react_iprevious_1dof;
                    VectorXd minusFmaxMinusReact = (-1) * fmax_1dof_vect - react_iprevious_1dof;
                    dF_1dof = (react_inc_1dof.array() > fmax_1dof_vect.array()).select(fmaxMinusReact, dF_1dof);
                    dF_1dof = (react_inc_1dof.array() < (-1) * fmax_1dof_vect.array()).select(minusFmaxMinusReact, dF_1dof);
                    dF(seq(0, dF.size() - 1, 6)) = dF_1dof;
                }

                /*VectorXd CI = (((*ptr_lamdasts) + (*ptr_lamdastd)) * dF + dUCAT + dUIP).cwiseQuotient(dU);
                CI = (CI.array().isNaN()).select(1, CI);
                for (int i = 0; i < seq(1, CI.size()-1, 6).size(); i++) {
                    CI(seq(1, CI.size()-1, 6)[i]) = 1;
                    CI(seq(1, CI.size()-1, 6)[i] + 2) = 1;
                    CI(seq(1, CI.size()-1, 6)[i] + 3) = 1;
                    CI(seq(1, CI.size()-1, 6)[i] + 4) = 1;
                }*/


                VectorXd dUIP_old = dUIP;
                VectorXd dUCAP = (*ptr_lamdasts) * dF;
                VectorXd dUd = (*ptr_lamdastd) * dF;
                dUIP(seq(2, dUIP.size() - 1, 6)) = dU(seq(2, dUIP.size() - 1, 6)) - dUd(seq(2, dUIP.size() - 1, 6))
                    - dUCAP(seq(2, dUIP.size() - 1, 6)) - dUCAT(seq(2, dUIP.size() - 1, 6));
                dUIP(seq(0, dUIP.size() - 1, 6)) = dU(seq(0, dUIP.size() - 1, 6)) - dUd(seq(0, dUIP.size() - 1, 6))
                    - dUCAP(seq(0, dUIP.size() - 1, 6)) - dUCAT(seq(0, dUIP.size() - 1, 6));
                //dUCAP = ((*ptr_lamdasts) * dF + beta * dUCAP)/(1 + beta);
                VectorXd perr = (dUIP - dUIP_old).cwiseQuotient(dUIP);
                perr = perr.array().isNaN().select(0, perr);
                if (n_while > 1500) {
                    std::cout << "tunnel while loops more than 1500 times" << std::endl;
                    break;
                }
                if (perr.cwiseAbs().maxCoeff() < eps_err) {
                    i_inc += 1;
                    (*ptr_vect_UINC).col(i_inc - 1) = UINC;
                    (*ptr_vect_dUIP).col(i_inc - 1) = dUIP;
                    UINC_iprevious = UINC;
                    P_elINC_iprevious = P_elINC;
                    //cout<<"i_inc is " << i_inc << " loop finished" <<endl;
//                    cout<<n_while<<endl;
                    break;
                }
            }
        }
    }
}

void elastoPlasticIterationLDLT(MatrixXd& Stiffness, MatrixXd& openKKfoot,VectorXd& Kstar,
                            MatrixXd& LL, MatrixXd& Lstar, VectorXd& uinc, VectorXd& ucat,VectorXd& P_el,
                            double mu, double lim_t, double lim_c,
                            VectorXd& patchArea,
                            VectorXd& uip, std::vector<int>& groundNodeDOF1, std::vector<int>& groundNodeDOF2, 
                            std::vector<int>& groundNodeDOF3, std::vector<int>& groundNodeDOF){
    
    #ifdef PRINT_INT_RESULTS
    VectorXi eig_vec;
    eig_vec = VectorXi::Map(groundNodeDOF.data(), groundNodeDOF.size());
    std::cout<<"groundNodeDOF.size(): "<<std::endl;
    std::cout<<groundNodeDOF.size()<<std::endl;
    std::cout<<"groundNodeDOF: "<<std::endl;
    std::cout<<eig_vec<<std::endl;
    eig_vec.resize(0);
    #endif
     auto start = high_resolution_clock::now();
    //  int load_step = 20;
     VectorXd ducat = ucat/20.0;
     int iter; VectorXd residual1; VectorXi residual3;
     VectorXd du; VectorXd rhs;
     VectorXd ducap = VectorXd::Zero(ucat.size());
     VectorXd duip = VectorXd::Zero(ucat.size());
     VectorXd f_prev, f_prev_v, f_prev_x, f_prev_y, f_curr, f_curr_v;
     VectorXd df, df_v, f_h_x, f_h_y, f_h_curr,f_h_x_slide, f_h_y_slide; 
     VectorXd df_t_limit, df_c_limit, f_h_max, df_h_x, df_h_y;
     VectorXd df_h_x_slide, df_h_y_slide, du_prim;
     VectorXd uinc_prev = uinc;
     VectorXd f_v_capacity = -lim_c*patchArea;
     VectorXd CI, duip_old,duip_prim;
     VectorXd uip_prev = VectorXd::Zero(ucat.size());
     VectorXd ucap;
     VectorXd ucap_prev = Lstar*(P_el - openKKfoot*uinc);
     VectorXd react_ep, react_ep_3dof;
     double beta = 4.0;
    //  SimplicialLDLT<SparseMatrix<double>, Eigen::Upper> iterSolver;
    //  iterSolver.compute(Stiffness);
    LDLT<MatrixXd> solver;
    solver.compute(Stiffness);
    auto stop = high_resolution_clock::now();
    auto duration = duration_cast<microseconds>(stop - start);
    //  std::cout << "eigen sparse ldlt compute time: "<<duration.count()/1000000.0 << std::endl;
     for (int i=0; i < 20; i++){
         iter = 0;
         int monitor_ind = 0;
         int n_iter =1000;
         residual1 = VectorXd::Zero(n_iter);
         residual3 = VectorXi::Zero(n_iter);
         du = VectorXd::Zero(ucat.size());
         double smallResidualNo = 0;
         double residual_prev1 = 0;
         while (iter < n_iter){
            rhs = Kstar.cwiseProduct(ducap + ducat + duip);
            du = solver.solve(rhs);
            
            f_prev = P_el - openKKfoot*uinc_prev;
            f_prev_v = f_prev(groundNodeDOF3);
            f_prev_x = f_prev(groundNodeDOF1);
            f_prev_y = f_prev(groundNodeDOF2);
            uinc = uinc_prev + du;
            f_curr = P_el - openKKfoot*uinc;
            f_curr_v = f_curr(groundNodeDOF3);

            df = f_curr-f_prev;
            df_v = df(groundNodeDOF3);
            df_t_limit = VectorXd::Zero(f_prev_v.size()) - f_prev_v;
            df_v = (f_curr_v.array()>lim_t).select(df_t_limit, df_v);
            df_c_limit = f_v_capacity - f_prev_v;
            df_v = (f_curr_v.array() < f_v_capacity.array()).select(df_c_limit, df_v);
            
            f_curr_v = df_v + f_prev_v;
            f_h_max = f_curr_v*(-mu);
            f_h_max = (f_h_max.array() > 0).select(f_h_max, VectorXd::Zero(f_h_max.size()));
            f_h_x = f_curr(groundNodeDOF1);
            f_h_y = f_curr(groundNodeDOF2);
            f_h_curr = (f_h_x.cwiseProduct(f_h_x) + f_h_y.cwiseProduct(f_h_y)).cwiseSqrt();
            f_h_x_slide = f_h_x.cwiseProduct(f_h_max.cwiseQuotient(f_h_curr).cwiseAbs());
            f_h_y_slide = f_h_y.cwiseProduct(f_h_max.cwiseQuotient(f_h_curr).cwiseAbs());
            f_h_x = (f_h_curr.cwiseAbs().array()>f_h_max.cwiseAbs().array()).select(f_h_x_slide, f_h_x);
            f_h_y = (f_h_curr.cwiseAbs().array()>f_h_max.cwiseAbs().array()).select(f_h_y_slide, f_h_y);
            df_h_x = df(groundNodeDOF1);
            df_h_y = df(groundNodeDOF2);
            df_h_x_slide = f_h_x - f_prev_x;
            df_h_y_slide = f_h_y - f_prev_y;
            df_h_x = (f_h_curr.cwiseAbs().array()>f_h_max.cwiseAbs().array()).select(df_h_x_slide, df_h_x);
            df_h_y = (f_h_curr.cwiseAbs().array()>f_h_max.cwiseAbs().array()).select(df_h_y_slide, df_h_y);
            df(groundNodeDOF1) = df_h_x;
            df(groundNodeDOF2) = df_h_y;
            df(groundNodeDOF3) = df_v;
            du = (du.cwiseAbs().array()<1e-8).select(0, du);
            VectorXd du_soil_el = VectorXd::Zero(df.size());
            du_soil_el(groundNodeDOF) = LL*(df(groundNodeDOF));
            du_prim = ducat + duip + du_soil_el;
            du_prim = (du_prim.cwiseAbs().array()<1e-8).select(0, du_prim);
            CI = (du.array()==0).select(1.0, du_prim.cwiseQuotient(du));
            // CI = CI(seq(0, CI.size()-1, 3));
            double CImax = (CI - VectorXd::Ones(CI.size())).cwiseAbs().maxCoeff();
            duip_old = duip;
            VectorXd diagFLEX = LL.diagonal();
            VectorXd du_soil_el_nonLocal = VectorXd::Zero(df.size());
            du_soil_el_nonLocal(groundNodeDOF) = diagFLEX.cwiseProduct(df(groundNodeDOF));
            duip = du - ducap - ducat-du_soil_el_nonLocal;
            // duip(groundNodeDOF) = duip(groundNodeDOF) - diagFLEX.cwiseProduct(df(groundNodeDOF));

            ducap = (Lstar*df + beta*ducap)/(1.0+beta);
            duip = (duip.cwiseAbs().array()<1e-8).select(0, duip);
            double perr = (duip - duip_old).norm();
            uip = uip_prev + duip;
            ucap = ucap_prev + ducap;
            react_ep = P_el - openKKfoot*uinc;
            react_ep_3dof = react_ep(groundNodeDOF3);

            residual1(iter) = (Stiffness*uinc - (P_el + Kstar.cwiseProduct(ucap+ducat*((double)i+1.0)+uip))).norm();
            residual3(iter) = (react_ep_3dof.array()>0.5).count();
            if (CImax<0.05 && perr<0.0001 && residual1(iter)<0.0001 && residual3(iter)==0){
                uinc_prev = uinc;
                uip_prev = uip;
                ucap_prev = ucap;
                // std::cout<<uinc(groundNodeDOF).head(50)<<std::endl;
                break;
            }
            if (residual1(iter)<1e-6) {
                uinc_prev = uinc;
                uip_prev = uip;
                ucap_prev = ucap;
                break;
            }
            if (iter ==n_iter-1){
                uinc_prev = uinc;
                uip_prev = uip;
                ucap_prev = ucap;
                // std::cout<<"iter reaches 500 at load step: "<< i<<std::endl;
                // std::cout<<"residual1"<<residual1(iter)<<std::endl;
                // std::cout<<"residual3"<<residual3(iter)<<std::endl;
            }
            if (std::abs(residual1(iter)-residual_prev1) < 1e-7){
                smallResidualNo++;
                residual_prev1 = residual1(iter);
            } else {
                smallResidualNo = 0;
                residual_prev1 = residual1(iter);
            }
            if (smallResidualNo>=10 && CImax<1e-4){
                uinc_prev = uinc;
                uip_prev = uip;
                ucap_prev = ucap;
                break;
            }
            if (residual1(iter)<0.00015 && CImax<1e-7){
                uinc_prev = uinc;
                uip_prev = uip;
                ucap_prev = ucap;
                break;
            }
//--------------------------monitor for debug----------------------------------------------
            // if ((iter%20) ==0){
            //     // uinc_prev = uinc;
            //     // uip_prev = uip;
            //     // ucap_prev = ucap;
            //     std::cout<< i  <<" step iter reaches 20 at monitor_ind: "<< monitor_ind<<std::endl;
            //     auto stop = high_resolution_clock::now();
            //     auto duration = duration_cast<microseconds>(stop - start);
            //     std::cout << duration.count()/1000000.0 << std::endl;
            //     std::cout << CI << std::endl;
            //     // std::cout<<"monitor: "<<monitor_ind*20<<" "<<monitor_ind*20+20<<std::endl;
            //     // std::cout<<"residual1"<<residual1.segment(monitor_ind*20, monitor_ind*20+20)<<std::endl;
            //     // std::cout<<"residual3"<<residual3.segment(monitor_ind*20, monitor_ind*20+20)<<std::endl;
            //     monitor_ind++;
            // }
            #ifdef PRINT_INT_RESULTS
            std::cout<<"load step: "<< i <<" iter: "<< iter <<" residual1: "
                     <<residual1(iter)<< " residual3" << residual3(iter) 
                     <<" CImax: "<<CImax<<" perr: "<<perr<<std::endl;
            #endif
            
//-----------------------------------------------------------------------------------------
            iter = iter + 1;

         }
            // std::cout<< i  <<" step iter finishes at iter: "<< iter<<std::endl;
            // auto stop = high_resolution_clock::now();
            // auto duration = duration_cast<microseconds>(stop - start);
            // std::cout <<"time from iter starts: "<< duration.count()/1000000.0 << std::endl;
            // std::cout<<"uinc: "<<std::endl;
            // std::cout<<uinc.head(50)<<std::endl;
            // std::cout<<"uinc(groundNodeDOF): "<<std::endl;
            // std::cout<<uinc(groundNodeDOF).head(50)<<std::endl;
     }

}       

void calInternalForces(VectorXd* F_M_deltaT_el_M_ptr, VectorXd* F_N_deltaT_el_M_ptr,
    VectorXd* F_S_deltaT_el_M_ptr, VectorXd& vect_Utotf_P,
    VectorXd& vect_Utotf_T, double Efoot, double EGratio,
    double dx, double dfoot, double bfoot, double ni_foot, int nnode) {
    int nelementi_foot = nnode - 1;
    MatrixXd KBern3Delt = KBern3D_foot_TIM(Efoot, dfoot, bfoot, dx, EGratio, ni_foot);
    // std::cout<<"KBern3Delt"<<std::endl;
    // std::cout<<KBern3Delt<<std::endl;
    MatrixXd Fin_P_1 = MatrixXd::Zero(6, nelementi_foot),
        Fin_P_2 = MatrixXd::Zero(6, nelementi_foot);
    MatrixXd Fin_T_1 = MatrixXd::Zero(6, nelementi_foot),
        Fin_T_2 = MatrixXd::Zero(6, nelementi_foot);

    //VectorXd Fin_P = KBern3Delt * (VectorXd)vect_Utotf_P(seq(6 * 0, 11 + 6 * 0, 1));
    //Fin_P_1(seqN(0, 6), 0) = Fin_P(seqN(0, 6));
    //cout<< (VectorXd)vect_Utotf_T(seq(6 * 0, 11 + 6 * 0, 1))<<endl;
    for (int j = 0; j < nelementi_foot; j++) {
        VectorXd Fin_P = KBern3Delt * (VectorXd)vect_Utotf_P(seq(6 * j, 11 + 6 * j, 1));
        Fin_P_1(seqN(0, 6), j) = Fin_P(seqN(0, 6));
        Fin_P_2(seqN(0, 6), j) = Fin_P(seqN(6, 6));
        VectorXd Fin_T = KBern3Delt * (VectorXd)vect_Utotf_T(seq(6 * j, 11 + 6 * j, 1));
        Fin_T_1(seqN(0, 6), j) = Fin_T(seqN(0, 6));
        Fin_T_2(seqN(0, 6), j) = Fin_T(seqN(6, 6));
    }

    MatrixXd Fin_deltaT_1 = Fin_T_1 - Fin_P_1;
    MatrixXd Fin_deltaT_2 = Fin_T_2 - Fin_P_2;
    VectorXd F_M_deltaT_el(2 * nelementi_foot), F_N_deltaT_el(2 * nelementi_foot),
        F_S_deltaT_el(2 * nelementi_foot);

    for (int i = 0; i < 2 * nelementi_foot; i += 2) {
        F_M_deltaT_el(i) = -Fin_deltaT_1(4, i / 2);
        F_M_deltaT_el(i + 1) = Fin_deltaT_2(4, i / 2);
        F_N_deltaT_el(i) = -Fin_deltaT_1(0, i / 2);
        F_N_deltaT_el(i + 1) = Fin_deltaT_2(0, i / 2);
        F_S_deltaT_el(i) = -Fin_deltaT_1(2, i / 2);
        F_S_deltaT_el(i + 1) = Fin_deltaT_2(2, i / 2);
    }
    *F_M_deltaT_el_M_ptr = F_M_deltaT_el;
    *F_N_deltaT_el_M_ptr = F_N_deltaT_el;
    *F_S_deltaT_el_M_ptr = F_S_deltaT_el;
}

VectorXd calculateStrain(VectorXd* F_S_deltaT_el_ptr, VectorXd* F_M_deltaT_el_ptr,
    VectorXd* F_N_deltaT_el_ptr, double Efoot, double EGratio, double bfoot, double dfoot,
    double ni_foot) {
    double I_fot = bfoot * pow(dfoot, 3) / 12.0;
    double c_shear = 1.5;
    double s_shear = 0;
    double d_na = dfoot / 2.0;
    VectorXd epsilon_dmax = c_shear * (*F_S_deltaT_el_ptr).cwiseAbs() /
        ((10 + ni_foot * 10) / (12 + 11 * ni_foot) * bfoot * dfoot) / (2 * Efoot / EGratio);
    VectorXd epsilon_bmax_top = (*F_M_deltaT_el_ptr) * (-1 * (dfoot - d_na)) / (Efoot * I_fot);
    VectorXd epsilon_bmax_bot = (*F_M_deltaT_el_ptr) * d_na / (Efoot * I_fot);
    VectorXd epsilon_bmax_pri = (*F_M_deltaT_el_ptr).cwiseAbs() * s_shear / (Efoot * I_fot);
    VectorXd epsilon_h = (*F_N_deltaT_el_ptr) / (Efoot * (bfoot * dfoot));
    VectorXd epsilon_br_top = epsilon_bmax_top + epsilon_h;
    VectorXd epsilon_br_bot = epsilon_bmax_bot + epsilon_h;
    VectorXd pri_plus_h = epsilon_bmax_pri + epsilon_h;
    VectorXd epsilon_dr = pri_plus_h * ((1 - ni_foot) / 2) +
        (pri_plus_h.cwiseProduct(pri_plus_h) * (1 + ni_foot) / 2.0 * (1 + ni_foot) / 2.0 +
            epsilon_dmax.cwiseProduct(epsilon_dmax)).cwiseSqrt();

    VectorXd result(3);
    result(0) = epsilon_br_top.maxCoeff();
    result(1) = epsilon_br_bot.maxCoeff();
    result(2) = epsilon_dr.maxCoeff();
    return result;
}
#endif
