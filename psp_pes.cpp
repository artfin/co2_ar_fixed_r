#include "psp_pes.hpp"
#include "leg_arr.hpp"

double psp_pes( const double R, const double Theta )
{
	double t5, t6, t7, t8, t11, t12, t17, t23, t24, t27, t28, t31, t32, t34, t35, t38, t39, t40,
	 t41, t42, t44, t46, t47, t49, t50, t51, t52, t53, t54, t56, t58, t59, t60, t62, t63, t64, t65, t66, t68, t71, t72, t74, t75, t78;
	
    double cosT = std::cos( Theta );

	double *legP = legendre_array(10, cosT );
	
	if (R < 6.329250)
	{
        t5 = std::exp(-0.16000000000000e-5 * R * (0.447680e6 + 0.54321e5 * R));
        t11 = std::exp(-0.10000000000000e-6 * R * (0.6504790e7 + 0.320299e6 * R));
        t17 = std::exp(-0.10000000000000e-5 * R * (0.811806e6 + 0.75313e5 * R));
        t24 = std::exp(-0.10000000000000e-6 * R * (0.6208770e7 + 0.310717e6 * R));
        t31 = std::exp(-0.20000000000000e-6 * R * (0.5878850e7 + 0.238569e6 * R));
        t38 = std::exp(-0.50000000000000e-5 * R * (0.147068e6 + 0.5935e4 * R));
        t42 = std::exp(-0.18813500000000e1 * R);
        t49 = std::exp(-0.10000000000000e-6 * R * (0.10554700e8 + 0.182219e6 * R));
        t53 = std::exp(-0.21459600000000e1 * R);
        t60 = std::exp(-0.10000000000000e-6 * R * (0.10394200e8 + 0.494781e6 * R));
        t64 = std::exp(-0.24461600000000e1 * R);
        t68 = std::exp(-0.20476500000000e1 * R);
		t71 = 0.224247e2 * t5 - 0.599056e0 * t11 + 0.635744e2 * legP[2] * t17 - 0.207391e0 * legP[2] * t24 + 0.991128e2 * legP[4] * t31 - 0.523497e-1 * legP[4] * t38 + 0.318652e3 * legP[6] * t42 - 0.374994e-1 * legP[6] * t49 + 0.332826e3 * legP[8] * t53 - 0.137376e-1 * legP[8] * t60 + 0.435837e3 * legP[10] * t64 - 0.108283e0 * legP[10] * t68;
		
		free( legP );
		return t71;
	}
	else if (R < 6.900260)
	{
        t5 = std::exp(-0.16000000000000e-5 * R * (0.447680e6 + 0.54321e5 * R));
		t6 = R * R;
		t7 = t6 * t6;
		t8 = t7 * t7;
		t12 = legP[2] * t8;
        t17 = std::exp(-0.10000000000000e-5 * R * (0.811806e6 + 0.75313e5 * R));
        t24 = std::exp(-0.10000000000000e-6 * R * (0.6208770e7 + 0.310717e6 * R));
		t27 = legP[4] * t8;
        t32 = std::exp(-0.20000000000000e-6 * R * (0.5878850e7 + 0.238569e6 * R));
        t39 = std::exp(-0.50000000000000e-5 * R * (0.147068e6 + 0.5935e4 * R));
		t42 = legP[6] * t8;
        t44 = std::exp(-0.18813500000000e1 * R);
        t51 = std::exp(-0.10000000000000e-6 * R * (0.10554700e8 + 0.182219e6 * R));
		t54 = legP[8] * t8;
        t56 = std::exp(-0.21459600000000e1 * R);
        t63 = std::exp(-0.10000000000000e-6 * R * (0.10394200e8 + 0.494781e6 * R));
		t66 = legP[10] * t8;
        t68 = std::exp(-0.24461600000000e1 * R);
        t72 = std::exp(-0.20476500000000e1 * R);
		t75 = 0.224247000e9 * t5 * t8 - 0.1145000000e10 * t6 - 0.23800000000e11 + 0.635744000e9 * t12 * t17 - 0.2073910e7 * t12 * t24 + 0.991128000e9 * t27 * t32 - 0.523497e6 * t27 * t39 + 0.3186520000e10 * t42 * t44 - 0.374994e6 * t42 * t51 + 0.3328260000e10 * t54 * t56 - 0.137376e6 * t54 * t63 + 0.4358370000e10 * t66 * t68 - 0.1082830e7 * t66 * t72;
		t78 = 0.10000000000000e-6 * t75 / t8;
		
		free( legP );
		return t78;
	}
	else if (R < 6.964500)
	{
        t5 = std::exp(-0.16000000000000e-5 * R * (0.447680e6 + 0.54321e5 * R));
		t6 = R * R;
		t7 = t6 * t6;
		t8 = t7 * t7;
        t17 = std::exp(-0.10000000000000e-5 * R * (0.811806e6 + 0.75313e5 * R));
		t23 = legP[4] * t8;
        t28 = std::exp(-0.20000000000000e-6 * R * (0.5878850e7 + 0.238569e6 * R));
        t35 = std::exp(-0.50000000000000e-5 * R * (0.147068e6 + 0.5935e4 * R));
		t38 = legP[6] * t8;
        t40 = std::exp(-0.18813500000000e1 * R);
        t47 = std::exp(-0.10000000000000e-6 * R * (0.10554700e8 + 0.182219e6 * R));
		t50 = legP[8] * t8;
        t52 = std::exp(-0.21459600000000e1 * R);
        t59 = std::exp(-0.10000000000000e-6 * R * (0.10394200e8 + 0.494781e6 * R));
		t62 = legP[10] * t8;
        t64 = std::exp(-0.24461600000000e1 * R);
        t68 = std::exp(-0.20476500000000e1 * R);
		t71 = 0.224247000e9 * t5 * t8 - 0.1145000000e10 * t6 - 0.23800000000e11 + 0.635744000e9 * legP[2] * t8 * t17 - 0.266000000e9 * legP[2] * t6 - 0.20800000000e11 * legP[2] + 0.991128000e9 * t23 * t28 - 0.523497e6 * t23 * t35 + 0.3186520000e10 * t38 * t40 - 0.374994e6 * t38 * t47 + 0.3328260000e10 * t50 * t52 - 0.137376e6 * t50 * t59 + 0.4358370000e10 * t62 * t64 - 0.1082830e7 * t62 * t68;
		t74 = 0.10000000000000e-6 * t71 / t8;
		
		free( legP );
		return t74;
	}
	else
	{
        t5 = std::exp(-0.16000000000000e-5 * R * (0.447680e6 + 0.54321e5 * R));
		t6 = R * R;
		t7 = t6 * t6;
		t8 = t7 * t7;
        t17 = std::exp(-0.10000000000000e-5 * R * (0.811806e6 + 0.75313e5 * R));
        t28 = std::exp(-0.20000000000000e-6 * R * (0.5878850e7 + 0.238569e6 * R));
		t32 = legP[6] * t8;
        t34 = std::exp(-0.18813500000000e1 * R);
        t41 = std::exp(-0.10000000000000e-6 * R * (0.10554700e8 + 0.182219e6 * R));
		t44 = legP[8] * t8;
        t46 = std::exp(-0.21459600000000e1 * R);
        t53 = std::exp(-0.10000000000000e-6 * R * (0.10394200e8 + 0.494781e6 * R));
		t56 = legP[10] * t8;
        t58 = std::exp(-0.24461600000000e1 * R);
        t62 = std::exp(-0.20476500000000e1 * R);
		t65 = 0.112123500e9 * t5 * t8 - 0.572500000e9 * t6 - 0.11900000000e11 + 0.317872000e9 * legP[2] * t8 * t17 - 0.133000000e9 * legP[2] * t6 - 0.10400000000e11 * legP[2] + 0.495564000e9 * legP[4] * t8 * t28 - 0.2050000000e10 * legP[4] + 0.1593260000e10 * t32 * t34 - 0.187497e6 * t32 * t41 + 0.1664130000e10 * t44 * t46 - 0.68688e5 * t44 * t53 + 0.2179185000e10 * t56 * t58 - 0.541415e6 * t56 * t62;
		t68 = 0.20000000000000e-6 * t65 / t8;

		free( legP );
		return t68;
	}
}

double dpsp_pesdR( const double R, const double Theta )
{
	double t5, t6, t7, t8, t9, t12, t13, t16, t21, t22, t24, t30, t31, t35, t36, t38, t39, t40, t41, t42, t44, t45, t48, t49, t50,
	 t51, t52, t54, t56, t57, t58, t59, t61, t62, t65, t66, t67, t71, t72, t73, t74, t75, t78,
	  t80, t81, t82, t84, t85, t86, t87, t88, t89, t92, t93, t95, t99;
	
    double cosT = std::cos( Theta );
	double *legP = legendre_array(10, cosT);
	
	if (R < 6.329250)
	{
        t5 = std::exp(-0.16000000000000e-5 * R * (0.447680e6 + 0.54321e5 * R));
        t13 = std::exp(-0.10000000000000e-6 * R * (0.6504790e7 + 0.320299e6 * R));
        t21 = std::exp(-0.10000000000000e-5 * R * (0.811806e6 + 0.75313e5 * R));
		t22 = legP[2] * t21;
        t30 = std::exp(-0.10000000000000e-6 * R * (0.6208770e7 + 0.310717e6 * R));
		t31 = legP[2] * t30;
        t39 = std::exp(-0.20000000000000e-6 * R * (0.5878850e7 + 0.238569e6 * R));
		t40 = legP[4] * t39;
        t48 = std::exp(-0.50000000000000e-5 * R * (0.147068e6 + 0.5935e4 * R));
		t49 = legP[4] * t48;
        t54 = std::exp(-0.18813500000000e1 * R);
        t61 = std::exp(-0.10000000000000e-6 * R * (0.10554700e8 + 0.182219e6 * R));
		t62 = legP[6] * t61;
        t67 = std::exp(-0.21459600000000e1 * R);
        t74 = std::exp(-0.10000000000000e-6 * R * (0.10394200e8 + 0.494781e6 * R));
		t75 = legP[8] * t74;
        t80 = std::exp(-0.24461600000000e1 * R);
        t84 = std::exp(-0.20476500000000e1 * R);
		t87 = -0.16062543513600e2 * t5 - 0.38980228118400e1 * t5 * R + 0.38967334782400e0 * t13 + 0.38375407548800e-1 * t13 * R - 0.51610079366400e2 * t22 - 0.95759575744000e1 * t22 * R + 0.12876430190700e0 * t31 + 0.12887981869400e-1 * t31 * R - 0.11653385685600e3 * t40 - 0.94580966332800e1 * t40 * R + 0.38494828398000e-1 * t49 + 0.31069546950000e-2 * t49 * R - 0.59949594020e3 * legP[6] * t54 + 0.39579491718000e-1 * t62 + 0.13666206337200e-2 * t62 * R - 0.71423128296e3 * legP[8] * t67 + 0.14279136192000e-1 * t75 + 0.13594206931200e-2 * t75 * R - 0.106612703592e4 * legP[10] * t80 + 0.22172568495e0 * legP[10] * t84;
		
		free( legP );
		return t87;
	}
	else if (R < 6.900260)
	{
        t5 = std::exp(-0.16000000000000e-5 * R * (0.447680e6 + 0.54321e5 * R));
		t6 = R * R;
		t7 = t6 * t6;
		t8 = t7 * t7;
		t9 = t8 * R;
		t12 = t8 * t6;
		t16 = legP[2] * t9;
        t21 = std::exp(-0.10000000000000e-5 * R * (0.811806e6 + 0.75313e5 * R));
		t24 = legP[2] * t12;
        t31 = std::exp(-0.10000000000000e-6 * R * (0.6208770e7 + 0.310717e6 * R));
		t36 = legP[4] * t9;
        t41 = std::exp(-0.20000000000000e-6 * R * (0.5878850e7 + 0.238569e6 * R));
		t44 = legP[4] * t12;
        t51 = std::exp(-0.50000000000000e-5 * R * (0.147068e6 + 0.5935e4 * R));
		t56 = legP[6] * t9;
        t58 = std::exp(-0.18813500000000e1 * R);
        t65 = std::exp(-0.10000000000000e-6 * R * (0.10554700e8 + 0.182219e6 * R));
		t71 = legP[8] * t9;
        t73 = std::exp(-0.21459600000000e1 * R);
        t80 = std::exp(-0.10000000000000e-6 * R * (0.10394200e8 + 0.494781e6 * R));
		t86 = legP[10] * t9;
        t88 = std::exp(-0.24461600000000e1 * R);
        t92 = std::exp(-0.20476500000000e1 * R);
		t95 = 0.40156358784000e15 * t5 * t9 + 0.97450570296000e14 * t5 * t12 - 0.17175000000000e17 * t6 - 0.47600000000000e18 + 0.12902519841600e16 * t16 * t21 + 0.23939893936000e15 * t24 * t21 - 0.3219107547675e13 * t16 * t31 - 0.322199546735e12 * t24 * t31 + 0.29133464214000e16 * t36 * t41 + 0.23645241583200e15 * t44 * t41 - 0.962370709950e12 * t36 * t51 - 0.77673867375e11 * t44 * t51 + 0.14987398505000e17 * t56 * t58 - 0.989487292950e12 * t56 * t65 - 0.34165515843e11 * legP[6] * t12 * t65 + 0.17855782074000e17 * t71 * t73 - 0.356978404800e12 * t71 * t80 - 0.33985517328e11 * legP[8] * t12 * t80 + 0.26653175898000e17 * t86 * t88 - 0.5543142123750e13 * t86 * t92;
		t99 = -0.40000000000000e-13 * t95 / t9;
		
		free( legP );
		return t99;
	}
	else if (R < 6.964500)
	{
        t5 = std::exp(-0.16000000000000e-5 * R * (0.447680e6 + 0.54321e5 * R));
		t6 = R * R;
		t7 = t6 * t6;
		t8 = t7 * t7;
		t9 = t8 * R;
		t12 = t8 * t6;
        t21 = std::exp(-0.10000000000000e-5 * R * (0.811806e6 + 0.75313e5 * R));
		t30 = legP[4] * t9;
        t35 = std::exp(-0.20000000000000e-6 * R * (0.5878850e7 + 0.238569e6 * R));
		t38 = legP[4] * t12;
        t45 = std::exp(-0.50000000000000e-5 * R * (0.147068e6 + 0.5935e4 * R));
		t50 = legP[6] * t9;
        t52 = std::exp(-0.18813500000000e1 * R);
        t59 = std::exp(-0.10000000000000e-6 * R * (0.10554700e8 + 0.182219e6 * R));
		t65 = legP[8] * t9;
        t67 = std::exp(-0.21459600000000e1 * R);
        t74 = std::exp(-0.10000000000000e-6 * R * (0.10394200e8 + 0.494781e6 * R));
		t80 = legP[10] * t9;
        t82 = std::exp(-0.24461600000000e1 * R);
        t86 = std::exp(-0.20476500000000e1 * R);
		t89 = 0.40156358784000e15 * t5 * t9 + 0.97450570296000e14 * t5 * t12 - 0.17175000000000e17 * t6 - 0.47600000000000e18 + 0.12902519841600e16 * legP[2] * t9 * t21 + 0.23939893936000e15 * legP[2] * t12 * t21 - 0.39900000000000e16 * legP[2] * t6 - 0.41600000000000e18 * legP[2] + 0.29133464214000e16 * t30 * t35 + 0.23645241583200e15 * t38 * t35 - 0.962370709950e12 * t30 * t45 - 0.77673867375e11 * t38 * t45 + 0.14987398505000e17 * t50 * t52 - 0.989487292950e12 * t50 * t59 - 0.34165515843e11 * legP[6] * t12 * t59 + 0.17855782074000e17 * t65 * t67 - 0.356978404800e12 * t65 * t74 - 0.33985517328e11 * legP[8] * t12 * t74 + 0.26653175898000e17 * t80 * t82 - 0.5543142123750e13 * t80 * t86;
		t93 = -0.40000000000000e-13 * t89 / t9;
	
		free( legP );
		return t93;
	}
	else
	{
        t5 = std::exp(-0.16000000000000e-5 * R * (0.447680e6 + 0.54321e5 * R));
		t6 = R * R;
		t7 = t6 * t6;
		t8 = t7 * t7;
		t9 = t8 * R;
		t12 = t8 * t6;
        t21 = std::exp(-0.10000000000000e-5 * R * (0.811806e6 + 0.75313e5 * R));
        t35 = std::exp(-0.20000000000000e-6 * R * (0.5878850e7 + 0.238569e6 * R));
		t42 = legP[6] * t9;
        t44 = std::exp(-0.18813500000000e1 * R);
        t51 = std::exp(-0.10000000000000e-6 * R * (0.10554700e8 + 0.182219e6 * R));
		t57 = legP[8] * t9;
        t59 = std::exp(-0.21459600000000e1 * R);
        t66 = std::exp(-0.10000000000000e-6 * R * (0.10394200e8 + 0.494781e6 * R));
		t72 = legP[10] * t9;
        t74 = std::exp(-0.24461600000000e1 * R);
        t78 = std::exp(-0.20476500000000e1 * R);
		t81 = 0.40156358784000e15 * t5 * t9 + 0.97450570296000e14 * t5 * t12 - 0.17175000000000e17 * t6 - 0.47600000000000e18 + 0.12902519841600e16 * legP[2] * t9 * t21 + 0.23939893936000e15 * legP[2] * t12 * t21 - 0.39900000000000e16 * legP[2] * t6 - 0.41600000000000e18 * legP[2] + 0.29133464214000e16 * legP[4] * t9 * t35 + 0.23645241583200e15 * legP[4] * t12 * t35 - 0.82000000000000e17 * legP[4] + 0.14987398505000e17 * t42 * t44 - 0.989487292950e12 * t42 * t51 - 0.34165515843e11 * legP[6] * t12 * t51 + 0.17855782074000e17 * t57 * t59 - 0.356978404800e12 * t57 * t66 - 0.33985517328e11 * legP[8] * t12 * t66 + 0.26653175898000e17 * t72 * t74 - 0.5543142123750e13 * t72 * t78;
		t85 = -0.40000000000000e-13 * t81 / t9;
	
		free( legP );
		return t85;
	}
}


double dpsp_pesdTheta( const double R, const double Theta )
{
	double t5, t6, t7, t8, t9, t14, t15, t23, t26, t27, t30, t33, t36, t37, t38, t40, t41,
	 t44, t45, t48, t49, t50, t55, t56, t57, t59, t63, t64, t66, t67, t73, t74, 
	 t75, t77, t81, t82, t85, t87, t89, t90, t91, t94, t95, t97, t99, t105; 
    double cosT = std::cos( Theta );
    double sinT = std::sin( Theta );
	double *legP = legendre_array( 10, cosT );

	if (R < 6.329250)
	{
        t5 = std::exp(-0.10000000000000e-5 * R * (0.811806e6 + 0.75313e5 * R));
		t6 = cosT * t5;
        t14 = std::exp(-0.10000000000000e-6 * R * (0.6208770e7 + 0.310717e6 * R));
		t15 = cosT * t14;
        t23 = std::exp(-0.20000000000000e-6 * R * (0.5878850e7 + 0.238569e6 * R));
        t33 = std::exp(-0.50000000000000e-5 * R * (0.147068e6 + 0.5935e4 * R));
        t40 = std::exp(-0.18813500000000e1 * R);
        t50 = std::exp(-0.10000000000000e-6 * R * (0.10554700e8 + 0.182219e6 * R));
        t57 = std::exp(-0.21459600000000e1 * R);
        t67 = std::exp(-0.10000000000000e-6 * R * (0.10394200e8 + 0.494781e6 * R));
        t74 = std::exp(-0.24461600000000e1 * R);
        t81 = std::exp(-0.20476500000000e1 * R);
		t87 = -0.317872000e9 * t6 + 0.317872000e9 * t6 * legP[2] + 0.1036955e7 * t15 - 0.1036955e7 * t15 * legP[2] - 0.991128000e9 * t23 * legP[3] + 0.991128000e9 * t23 * cosT * legP[4] + 0.523497e6 * t33 * legP[3] - 0.523497e6 * t33 * cosT * legP[4] - 0.4779780000e10 * t40 * legP[5] + 0.4779780000e10 * t40 * cosT * legP[6] + 0.562491e6 * t50 * legP[5] - 0.562491e6 * t50 * cosT * legP[6] - 0.6656520000e10 * t57 * legP[7] + 0.6656520000e10 * t57 * cosT * legP[8] + 0.274752e6 * t67 * legP[7] - 0.274752e6 * t67 * cosT * legP[8] - 0.10895925000e11 * t74 * legP[9] + 0.10895925000e11 * t74 * cosT * legP[10] + 0.2707075e7 * t81 * legP[9] - 0.2707075e7 * t81 * cosT * legP[10];
		t89 = cosT * cosT;
		t94 = -0.40000000000000e-6 * t87 * sinT / (t89 - 0.1e1);
	
		free( legP );
		return t94;
	}
	else if (R < 6.900260)
	{
        t5 = std::exp(-0.10000000000000e-5 * R * (0.811806e6 + 0.75313e5 * R));
		t6 = cosT * t5;
        t14 = std::exp(-0.10000000000000e-6 * R * (0.6208770e7 + 0.310717e6 * R));
		t15 = cosT * t14;
        t23 = std::exp(-0.20000000000000e-6 * R * (0.5878850e7 + 0.238569e6 * R));
        t33 = std::exp(-0.50000000000000e-5 * R * (0.147068e6 + 0.5935e4 * R));
        t40 = std::exp(-0.18813500000000e1 * R);
        t50 = std::exp(-0.10000000000000e-6 * R * (0.10554700e8 + 0.182219e6 * R));
        t57 = std::exp(-0.21459600000000e1 * R);
        t67 = std::exp(-0.10000000000000e-6 * R * (0.10394200e8 + 0.494781e6 * R));
        t74 = std::exp(-0.24461600000000e1 * R);
        t81 = std::exp(-0.20476500000000e1 * R);
		t87 = -0.317872000e9 * t6 + 0.317872000e9 * t6 * legP[2] + 0.1036955e7 * t15 - 0.1036955e7 * t15 * legP[2] - 0.991128000e9 * t23 * legP[3] + 0.991128000e9 * t23 * cosT * legP[4] + 0.523497e6 * t33 * legP[3] - 0.523497e6 * t33 * cosT * legP[4] - 0.4779780000e10 * t40 * legP[5] + 0.4779780000e10 * t40 * cosT * legP[6] + 0.562491e6 * t50 * legP[5] - 0.562491e6 * t50 * cosT * legP[6] - 0.6656520000e10 * t57 * legP[7] + 0.6656520000e10 * t57 * cosT * legP[8] + 0.274752e6 * t67 * legP[7] - 0.274752e6 * t67 * cosT * legP[8] - 0.10895925000e11 * t74 * legP[9] + 0.10895925000e11 * t74 * cosT * legP[10] + 0.2707075e7 * t81 * legP[9] - 0.2707075e7 * t81 * cosT * legP[10];
		t89 = cosT * cosT;
		t94 = -0.40000000000000e-6 * t87 * sinT / (t89 - 0.1e1);
		
		free( legP );
		return t94;
	}
	else if (R < 6.964500)
	{
        t5 = std::exp(-0.10000000000000e-5 * R * (0.811806e6 + 0.75313e5 * R));
		t6 = cosT * t5;
		t7 = R * R;
		t8 = t7 * t7;
		t9 = t8 * t8;
		t15 = cosT * t7;
        t26 = std::exp(-0.20000000000000e-6 * R * (0.5878850e7 + 0.238569e6 * R));
		t27 = t9 * t26;
		t30 = cosT * legP[4];
        t37 = std::exp(-0.50000000000000e-5 * R * (0.147068e6 + 0.5935e4 * R));
		t38 = t9 * t37;
        t44 = std::exp(-0.18813500000000e1 * R);
		t45 = t9 * t44;
		t48 = -0.317872000e9 * t6 * t9 + 0.317872000e9 * t6 * t9 * legP[2] + 0.133000000e9 * t15 - 0.133000000e9 * t15 * legP[2] + 0.10400000000e11 * cosT - 0.10400000000e11 * cosT * legP[2] - 0.991128000e9 * t27 * legP[3] + 0.991128000e9 * t27 * t30 + 0.523497e6 * t38 * legP[3] - 0.523497e6 * t38 * t30 - 0.4779780000e10 * t45 * legP[5];
		t49 = cosT * legP[6];
        t56 = std::exp(-0.10000000000000e-6 * R * (0.10554700e8 + 0.182219e6 * R));
		t57 = t9 * t56;
        t63 = std::exp(-0.21459600000000e1 * R);
		t64 = t9 * t63;
		t67 = cosT * legP[8];
        t74 = std::exp(-0.10000000000000e-6 * R * (0.10394200e8 + 0.494781e6 * R));
		t75 = t9 * t74;
        t81 = std::exp(-0.24461600000000e1 * R);
		t82 = t9 * t81;
		t85 = cosT * legP[10];
        t89 = std::exp(-0.20476500000000e1 * R);
		t90 = t9 * t89;
		t95 = 0.4779780000e10 * t45 * t49 + 0.562491e6 * t57 * legP[5] - 0.562491e6 * t57 * t49 - 0.6656520000e10 * t64 * legP[7] + 0.6656520000e10 * t64 * t67 + 0.274752e6 * t75 * legP[7] - 0.274752e6 * t75 * t67 - 0.10895925000e11 * t82 * legP[9] + 0.10895925000e11 * t82 * t85 + 0.2707075e7 * t90 * legP[9] - 0.2707075e7 * t90 * t85;
		t99 = cosT * cosT;
		t105 = -0.40000000000000e-6 * (t48 + t95) * sinT / t9 / (t99 - 0.1e1);
		
		free( legP );
		return t105;
	}
	else
	{
        t5 = std::exp(-0.10000000000000e-5 * R * (0.811806e6 + 0.75313e5 * R));
		t6 = cosT * t5;
		t7 = R * R;
		t8 = t7 * t7;
		t9 = t8 * t8;
		t15 = cosT * t7;
        t26 = std::exp(-0.20000000000000e-6 * R * (0.5878850e7 + 0.238569e6 * R));
		t27 = t9 * t26;
		t30 = cosT * legP[4];
        t36 = std::exp(-0.18813500000000e1 * R);
		t37 = t9 * t36;
		t40 = -0.317872000e9 * t6 * t9 + 0.317872000e9 * t6 * t9 * legP[2] + 0.133000000e9 * t15 - 0.133000000e9 * t15 * legP[2] + 0.10400000000e11 * cosT - 0.10400000000e11 * cosT * legP[2] - 0.991128000e9 * t27 * legP[3] + 0.991128000e9 * t27 * t30 + 0.4100000000e10 * legP[3] - 0.4100000000e10 * t30 - 0.4779780000e10 * t37 * legP[5];
		t41 = cosT * legP[6];
        t48 = std::exp(-0.10000000000000e-6 * R * (0.10554700e8 + 0.182219e6 * R));
		t49 = t9 * t48;
		t55 = exp(-0.21459600000000e1 * R);
		t56 = t9 * t55;
		t59 = cosT * legP[8];
        t66 = std::exp(-0.10000000000000e-6 * R * (0.10394200e8 + 0.494781e6 * R));
		t67 = t9 * t66;
        t73 = std::exp(-0.24461600000000e1 * R);
		t74 = t9 * t73;
		t77 = cosT * legP[10];
        t81 = std::exp(-0.20476500000000e1 * R);
		t82 = t9 * t81;
		t87 = 0.4779780000e10 * t37 * t41 + 0.562491e6 * t49 * legP[5] - 0.562491e6 * t49 * t41 - 0.6656520000e10 * t56 * legP[7] + 0.6656520000e10 * t56 * t59 + 0.274752e6 * t67 * legP[7] - 0.274752e6 * t67 * t59 - 0.10895925000e11 * t74 * legP[9] + 0.10895925000e11 * t74 * t77 + 0.2707075e7 * t82 * legP[9] - 0.2707075e7 * t82 * t77;
		t91 = cosT * cosT;
		t97 = -0.40000000000000e-6 * (t40 + t87) * sinT / t9 / (t91 - 0.1e1);
	
		free( legP );
		return t97;
	}
}
