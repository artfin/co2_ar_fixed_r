#include "matrix_euler.hpp"

const double mu1 = 14579.0;
const double mu2 = 38183.0;
const double l = 4.398; 

void inertia_tensor_dR( Eigen::Matrix3d & i_dr, const double R, const double Theta)
{
    i_dr(0, 0) = i_dr(1, 1) = 2 * mu2 * R;
    
    i_dr(1, 0) = i_dr(2, 0) = 0;    
    i_dr(0, 2) = i_dr(1, 2) = i_dr(2, 2) = 0;
    i_dr(0, 1) = i_dr(2, 1) = 0;
}

void inertia_tensor_dTheta( Eigen::Matrix3d & i_dt, const double R, const double Theta)
{
    double l2 = l * l;
    double sin_t = std::sin(Theta);
    double cos_t = std::cos(Theta);

    i_dt(0, 0) = - 2 * mu1 * l2 * sin_t * cos_t;
    i_dt(1, 0) = 0;
    i_dt(2, 0) = i_dt(0, 2) = - mu1 * l2 * (cos_t * cos_t - sin_t * sin_t);

    i_dt(1, 2) = 0;
    i_dt(2, 2) = - i_dt(0, 0); 

    i_dt(0, 1) = i_dt(1, 1) = i_dt(2, 1) = 0;
}

void inertia_tensor( Eigen::Matrix3d & inertia_tensor, const double R, const double Theta)
{
    double sin_t = std::sin(Theta);
    double cos_t = std::cos(Theta);
    double l2 = l * l;
    double R2 = R * R;

    inertia_tensor(0, 0) = mu1 * l2 * cos_t * cos_t + mu2 * R2; 
    inertia_tensor(1, 0) = 0;
    inertia_tensor(2, 0) = -mu1 * l2 * sin_t * cos_t;

    inertia_tensor(0, 2) = inertia_tensor(2, 0);
    inertia_tensor(1, 2) = 0;
    inertia_tensor(2, 2) = mu1 * l2 * sin_t * sin_t;

    inertia_tensor(0, 1) = 0;
    inertia_tensor(1, 1) = inertia_tensor(0, 0) + inertia_tensor(2, 2);
    inertia_tensor(2, 1) = 0;
}

void fill_a_matrix( Eigen::Matrix2d & a, const double R, const double Theta)
{
    a(0, 0) = mu2;
    a(0, 1) = 0;
    a(1, 0) = 0;
    a(1, 1) = mu1 * l * l;
}

void fill_A_matrix( Eigen::Matrix<double, 3, 2> & A, const double R, const double Theta)
{
    A(0, 0) = A(0, 1) = 0;
    A(1, 0) = 0;
    A(1, 1) = mu1 * l * l;
    A(2, 0) = A(2, 1) = 0;
}

// J = (V^{-1})^\top pe = W pe
void W_matrix( Eigen::Matrix<double, 3, 3> & W, const double theta, const double psi )
{
    double sin_psi = std::sin( psi );
    double cos_psi = std::cos( psi );

    double sin_theta = std::sin( theta );
    double cos_theta = std::cos( theta );
	double ctg_theta = cos_theta / sin_theta;

	W(0, 0) = sin_psi / sin_theta;
	W(0, 1) = cos_psi;
	W(0, 2) = - sin_psi * ctg_theta;
	
	W(1, 0) = cos_psi / sin_theta;
	W(1, 1) = - sin_psi;
	W(1, 2) = - cos_psi * ctg_theta;

	W(2, 0) = 0;
	W(2, 1) = 0;
	W(2, 2) = 1;
}

// \frac{\partial W}{\partial \psi}
void dW_dpsi( Eigen::Matrix<double, 3, 3> & dW, const double theta, const double psi )
{
    double cos_psi = std::cos( psi );
    double sin_psi = std::sin( psi );

    double sin_theta = std::sin( theta );
    double cos_theta = std::cos( theta );
	double ctg_theta = cos_theta / sin_theta;

	dW(0, 0) = cos_psi / sin_theta;
   	dW(0, 1) = - sin_psi;
	dW(0, 2) = - cos_psi * ctg_theta;

	dW(1, 0) = - sin_psi / sin_theta;
	dW(1, 1) = - cos_psi;
	dW(1, 2) = sin_psi * ctg_theta;

	dW(2, 0) = 0;
	dW(2, 1) = 0;
	dW(2, 2) = 0;
}

// \frac{\partial W}{\partial \theta}
void dW_dtheta( Eigen::Matrix<double, 3, 3> & dW, const double theta, const double psi )
{
    double sin_psi = std::sin( psi );
    double cos_psi = std::cos( psi );

    double sin_theta = std::sin( theta );
    double cos_theta = std::cos( theta );

    dW(0, 0) = - sin_psi * cos_theta / sin_theta / sin_theta;
	dW(0, 1) = 0;
    dW(0, 2) = sin_psi / sin_theta / sin_theta;

    dW(1, 0) = - cos_psi * cos_theta / sin_theta / sin_theta;
	dW(1, 1) = 0;
    dW(1, 2) = cos_psi / sin_theta / sin_theta;

	dW(2, 0) = 0;
	dW(2, 1) = 0;
	dW(2, 2) = 0;
}

void transform_dipole( std::vector<double> &output, const double R,
                                                    const double Theta,
                                                    const double phi,
                                                    const double theta,
                                                    const double psi )
{
    Eigen::Matrix<double, 3, 3> S;
	
    double sin_phi = std::sin( phi );
    double cos_phi = std::cos( phi );

    double sin_theta = std::sin( theta );
    double cos_theta = std::cos( theta );

    double sin_psi = std::sin( psi );
    double cos_psi = std::cos( psi );

	S(0, 0) = cos_psi * cos_phi - cos_theta * sin_phi * sin_psi;
	S(0, 1) = - sin_psi * cos_phi - cos_theta * sin_phi * cos_psi;
	S(0, 2) = sin_theta * sin_phi;

	S(1, 0) = cos_psi * sin_phi + cos_theta * cos_phi * sin_psi;
	S(1, 1) = - sin_psi * sin_phi + cos_theta * cos_phi * cos_psi;
	S(1, 2) = - sin_theta * cos_phi;

	S(2, 0) = sin_theta * sin_psi;
	S(2, 1) = sin_theta * cos_psi;
	S(2, 2) = cos_theta;

	// vector of dipole in molecular frame
    Eigen::Vector3d mol_dipole;

	mol_dipole(0) = dipx( R, Theta );
	mol_dipole(1) = 0;
	mol_dipole(2) = dipz( R, Theta );

	// vector of dipole in laboratory frame
    Eigen::Vector3d lab_dipole = S * mol_dipole;

	output[0] = lab_dipole(0);
	output[1] = lab_dipole(1);
	output[2] = lab_dipole(2);
}

void rhs( double * out, const double R, const double pR,
                        const double Theta, const double pT,
                        const double phi, const double p_phi,
                        const double theta, const double p_theta,
                        const double psi, const double p_psi )
// input:
//     R, pR, Theta, pTheta, phi, pPhi, theta, p_theta, psi, p_psi
// output:
//    d(R)/dt, d(pR)/dt, d(Theta)/dt, d(pTheta)/dt, d(phi)/dt, d(pPhi)/dt, d(theta)/dt, d(pTheta)/dt, d(psi)/dt, d(pPsi)/dt
{
	// vector of euler impulses
    Eigen::Vector3d pe ( p_phi, p_theta, p_psi );

	// vector of impulses
    Eigen::Vector2d p_vector(pR, pT);
	
	// ##################################################################
    // filling matrices I, a, A and derivatives of I
    Eigen::Matrix<double, 3, 3> I;
    inertia_tensor(I, R, Theta); 
	//std::cout << "I: " << I << std::endl;

    Eigen::Matrix<double, 3, 3> I_dR;
    inertia_tensor_dR(I_dR, R, Theta);
    
    Eigen::Matrix<double, 3, 3> I_dTheta;
    inertia_tensor_dTheta(I_dTheta, R, Theta);
    
    Eigen::Matrix<double, 2, 2> a;
    fill_a_matrix(a, R, Theta);
    
    Eigen::Matrix<double, 3, 2> A;
    fill_A_matrix(A, R, Theta);

    Eigen::Matrix<double, 3, 3> W;
	W_matrix( W, theta, psi );
	// ##################################################################


	// ##################################################################
    // filling matrices G11, G22, G12
    Eigen::Matrix<double, 3, 3> I_inv = I.inverse();
    Eigen::Matrix<double, 2, 2> a_inv = a.inverse();

    Eigen::Matrix<double, 3, 3> G11;
    Eigen::Matrix<double, 2, 2> G22;
    Eigen::Matrix<double, 3, 2> G12;

    Eigen::Matrix<double, 3, 3> t1 = I;
    t1.noalias() -= A * a_inv * A.transpose();
    G11 = t1.inverse();
 
    Eigen::Matrix<double, 2, 2> t2 = a;
	t2.noalias() -= A.transpose() * I_inv * A;
    G22 = t2.inverse();

    G12.noalias() = - G11 * A * a.inverse();
	// ##################################################################
	

	// ##################################################################
    // derivatives \frac{\partial \mH}{\partial p} 
    Eigen::Vector2d dH_dp = G22 * p_vector + G12.transpose() * W * pe;
	//std::cout << "dH_dp: " << dH_dp << std::endl;
	// ##################################################################


	// auxiliary variables
   	double ang_term, kin_term, cor_term;
 
	// ##################################################################
   	// derivative \frac{\partial \mH}{\partial R}
    Eigen::Matrix<double, 3, 3> G11_dR;
    Eigen::Matrix<double, 2, 2> G22_dR;
    Eigen::Matrix<double, 3, 2> G12_dR;

   	G11_dR = - G11 * I_dR * G11;
   	G12_dR = - G11_dR * A * a_inv;
   	G22_dR = - G22 * A.transpose() * I_inv * I_dR * I_inv * A * G22; 
   
   	ang_term = 0.5 * pe.transpose() * W.transpose() * G11_dR * W * pe;
   	kin_term = 0.5 * p_vector.transpose() * G22_dR * p_vector;
	cor_term = pe.transpose() * W.transpose() * G12_dR * p_vector;

	double dH_dR = ang_term + kin_term + cor_term;
	// ##################################################################


	// ##################################################################
   	// derivative \frac{\partial \mH}{\partial Theta}
    Eigen::Matrix<double, 3, 3> G11_dTheta;
    Eigen::Matrix<double, 2, 2> G22_dTheta;
    Eigen::Matrix<double, 3, 2> G12_dTheta;

	G11_dTheta = - G11 * I_dTheta * G11;
	G22_dTheta = - G22 * A.transpose() * I_inv * I_dTheta * I_inv * A * G22;
	G12_dTheta = -G11_dTheta * A * a_inv;

    ang_term = 0.5 * pe.transpose() * W.transpose() * G11_dTheta * W * pe;
	kin_term = 0.5 * p_vector.transpose() * G22_dTheta * p_vector;
	cor_term = pe.transpose() * W.transpose() * G12_dTheta * p_vector;

	double dH_dTheta = ang_term + kin_term + cor_term;	
	// ##################################################################

	// auxiliary variables
	double ang_term1, ang_term2;
	// ##################################################################
	// derivative /frac{\partial \mH}{\partial \theta} (euler angle)
    Eigen::Matrix<double, 3, 3> W_dtheta;
   	dW_dtheta( W_dtheta, theta, psi );

	ang_term1 = 0.5 * pe.transpose() * W_dtheta.transpose() * G11 * W * pe;
	ang_term2 = 0.5 * pe.transpose() * W.transpose() * G11 * W_dtheta * pe;
	cor_term = pe.transpose() * W_dtheta.transpose() * G12 * p_vector;

	double dH_dtheta = ang_term1 + ang_term2 + cor_term;
	// ##################################################################

	// ##################################################################
	// derivative \frac{\partial \mH}{\partial \psi} (euler angle)
    Eigen::Matrix<double, 3, 3> W_dpsi;
	dW_dpsi( W_dpsi, theta, psi );

	ang_term1 = 0.5 * pe.transpose() * W_dpsi.transpose() * G11 * W * pe;
	ang_term2 = 0.5 * pe.transpose() * W.transpose() * G11 * W_dpsi * pe;
	cor_term = pe.transpose() * W_dpsi.transpose() * G12 * p_vector;

	double dH_dpsi = ang_term1 + ang_term2 + cor_term;
	// ##################################################################

	// ##################################################################
	// derivative \frac{\partial \mH}{\partial \phi} (euler angle)
	double dH_dphi = 0;
	// ##################################################################

	// ##################################################################
	// derivatives \frac{\partial \mH}{\partial pe} (euler impulses)
	
    Eigen::Vector3d dH_dpe = W.transpose() * G11 * W * pe + W.transpose() * G12 * p_vector;
	//std::cout << "dH_dpe: " << dH_dpe << std::endl;
	// ##################################################################

/*
	out[0] = dH_dp(0); 
	out[1] = dH_dp(1); 
	out[2] = - dH_dR - dpsp_pesdR( R, Theta ); 
	out[3] = - dH_dTheta - dpsp_pesdTheta( R, Theta );

	out[4] = dH_dpe(0);
	out[5] = dH_dpe(1);
   	out[6] = dH_dpe(2);

	out[7] = - dH_dphi; 
	out[8] = - dH_dtheta;
	out[9] = - dH_dpsi;
*/

    out[0] = dH_dp(0); // dR/dt
    out[1] = - dH_dR - dpsp_pesdR(R, Theta); // d(pR)/dt
    out[2] = dH_dp(1); // dTheta/dt
    out[3] = - dH_dTheta - dpsp_pesdTheta(R, Theta); // d(pTheta)/dt
    out[4] = dH_dpe(0); // d(phi)/dt
    out[5] = - dH_dphi; // d(pPhi)/dt
    out[6] = dH_dpe(1); // d(theta)/dt
    out[7] = - dH_dtheta; // d(pTheta)/dt
    out[8] = dH_dpe(2); // d(psi)/dt
    out[9] = - dH_dpsi; // d(psi)/dt
}
