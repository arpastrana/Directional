#ifndef DIRECTIONAL_GET_SUBDIVISION_MATRIX_H
#define DIRECTIONAL_GET_SUBDIVISION_MATRIX_H
#include <igl/igl_inline.h>
#include <Eigen/Core>
#include <Eigen/Sparse>
#include "SubdivisionInternal/shm_edge_topology.h"
#include "SubdivisionInternal/Sv_triplet_provider.h"
#include "SubdivisionInternal/Sc_triplet_provider.h"
#include "SubdivisionInternal/Se_triplet_provider.h"
#include "SubdivisionInternal/Sc_directional_triplet_provider.h"
#include "SubdivisionInternal/Se_directional_triplet_provider.h"
#include "SubdivisionInternal/Sf_triplet_provider.h"
#include "SubdivisionInternal/build_directional_subdivision_operators.h"
#include "SubdivisionInternal/shm_halfcurl_coefficients.h"
#include "SubdivisionInternal/shm_oneform_coefficients.h"
#include "SubdivisionInternal/hbspline_coefficients.h"

namespace directional{
    IGL_INLINE void get_subdivision_matrix(
        const Eigen::MatrixXd& VCoarse,
        const Eigen::MatrixXi& FCoarse,
        const Eigen::MatrixXi& EVCoarse,
        int subdivisionLevel,
        int N,
        const Eigen::VectorXi& matchingCoarse,
        Eigen::SparseMatrix<double>& PCoarse,
        Eigen::MatrixXi& EVFine,
        Eigen::MatrixXi& FFine,
        Eigen::SparseMatrix<double>& S_Gamma,
        Eigen::SparseMatrix<double>& PInvFine,
        Eigen::VectorXi& matchingFine
    )
    {

    }
    IGL_INLINE void get_subdivision_matrix(
        const Eigen::MatrixXd& VCoarse,
        const Eigen::MatrixXi& FCoarse,
        const Eigen::MatrixXi& EVCoarse,
        int subdivisionLevel,
        Eigen::SparseMatrix<double>& PCoarse,
        Eigen::MatrixXi& EVFine,
        Eigen::MatrixXi& FFine,
        Eigen::SparseMatrix<double>& S_Gamma,
        Eigen::SparseMatrix<double>& PInvFine
    )
    {
    Eigen::MatrixXi EICoarse, SFECoarse, EFCoarse, EIFine, SFEFine, EFFine;
        shm_edge_topology(FCoarse, EVCoarse, std::ref(EFCoarse), EICoarse, SFECoarse);
        // Initial sizes of the subdivision operators
        std::vector<int> initialSizes = std::vector<int>({ (int)VCoarse.rows(), (int)(EVCoarse.rows()), (int)(EVCoarse.rows()), (int)FCoarse.rows() });
        std::vector<Eigen::SparseMatrix<double>> subdivisionOperators;

        using coeffProv = coefficient_provider_t;

        // Providers for all subdivision operators
        auto Sv_provider = triplet_provider_wrapper<coeffProv>(subdivision::loop_coefficients, subdivision::Sv_triplet_provider<coeffProv>);
        auto Sc_provider = triplet_provider_wrapper<coeffProv>(subdivision::shm_halfcurl_coefficients, subdivision::Sc_triplet_provider<coeffProv>);
        auto Se_provider = triplet_provider_wrapper<coeffProv>(subdivision::shm_oneform_coefficients, subdivision::Se_triplet_provider<coeffProv>);

        // Construct regular vertex subdivision
        build_subdivision_operators(VCoarse, FCoarse, EVCoarse, EFCoarse, EICoarse, SFECoarse, std::vector<int>({ (int)VCoarse.rows() }), subdivisionLevel,
            FFine, EVFine, EFFine, EIFine, SFEFine, subdivisionOperators, Sv_provider, Se_provider, Sc_provider, Sf_provider);
        Eigen::SparseMatrix<double> S_0 = subdivisionOperators[0];
        Eigen::SparseMatrix<double> S_1 = subdivisionOperators[1];
        Eigen::SparseMatrix<double> S_epsstar = subdivisionOperators[2];
        Eigen::SparseMatrix<double> S_2 = subdivisionOperators[3];

        Eigen::VectorXd VCoarseVec;
        VCoarseVec.resize(3 * VCoarse.rows());
        for (int v = 0; v < VCoarse.rows(); ++v) VCoarseVec.middleRows(3 * v, 3) = VCoarse.row(v).transpose();
        Eigen::VectorXd VFineVec = S_0 * VCoarseVec;
        Eigen::MatrixXd VFine;
        VFine.setConstant(VCoarseVec.rows() / 3, 3, -1);
        for (int v = 0; v < VFine.rows(); ++v)
        {
            VFine.row(v) = VFineVec.middleRows(3 * v, 3).transpose();
        }
        
        Eigen::SparseMatrix<double> subdivider,WCoarse,WInvFine;
        directional::block_diag({&S_0, &S_epsstar}, subdivider);
        

        get_W(FCoarse, EVCoarse, SFECoarse.leftCols(3), EFCoarse, WCoarse);
        get_P(VCoarse, FCoarse, EVCoarse, SFECoarse.leftCols(3),1, PCoarse);
        get_P_inverse(VFine, FFine, EVFine, SFEFine.leftCols(3),1, PInvFine);
        get_W_inverse(VFine, FFine, EVFine, SFEFine.leftCols(3), WInvFine);
        S_Gamma = WInvFine * subdivider * WCoarse;
    }
}
#endif