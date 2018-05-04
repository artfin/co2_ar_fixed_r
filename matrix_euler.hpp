#pragma once

#include <cmath>
#include <vector>
#include <Eigen/Dense>

#include "psp_pes.hpp"
#include "co2_ar_dipole.hpp"

void inertia_tensor( Eigen::Matrix3d & inertia_tensor, const double R, const double Theta );
void inertia_tensor_dR( Eigen::Matrix3d & i_dr, double R, double Theta);
void inertia_tensor_dTheta( Eigen::Matrix3d & i_dt, const double R, const double Theta);

void fill_a_matrix( Eigen::Matrix2d & a, const double R, const double Theta);
void fill_A_matrix( Eigen::Matrix<double, 3, 2> & A, const double R, const double Theta);

void W_matrix( Eigen::Matrix<double, 3, 3> & W, const double theta, const double psi );
void dW_dpsi( Eigen::Matrix<double, 3, 3> & dW, const double theta, const double psi );
void dW_dtheta( Eigen::Matrix<double, 3, 3> & dW, const double theta, const double psi );

void rhs( double * out, const double R, const double pR,
                        const double Theta, const double pT,
                        const double phi, const double pPhi,
                        const double theta, const double p_theta,
                        const double psi, const double p_psi );

void transform_dipole( std::vector<double> & output, const double R,
                                                     const double Theta,
                                                     const double phi,
                                                     const double theta,
                                                     const double psi );


