#include "minperimeter.h"
#include <string>
#include <iostream>
#include <fstream>
#include <omp.h>
#include <LBFGS.h>

using namespace LBFGSpp;

//!Solve till viscosity becomes zero.
//!
//! Parameters:
//! -----------
//! a, b : x and y coordinates of the starting position
//! pointid : suffix for the output file
//! area: Area constraint value
//! ffft: surface equation to calculate first fundamental form with derivatives
//! numpts: number of points on the curve
//! S: system size
//! k: penalty coefficient
//! viscosity: the coefficient of the viscosity term
//! maxiter: maximum number of solver iterations
//! tol: when change in solution < tol, the solver exits.
void solve(scalar a, scalar b, const int pointid,
	const scalar area, generatedfunc ffft, const int numpts,
	const scalar S, const int maxiter=100, const scalar k=1e5,
	const scalar viscosity=2.0, const scalar viscoustol=1e-16,
	const scalar epsilon=1e-16){

    std::string filename = "path-" + std::to_string(pointid) + ".csv";
    std::ofstream outfile;
    outfile.open(filename);
    outfile.precision(digits);
    outfile << std::fixed;

    VectorXmp Rab(numpts + 2), Rold(numpts);
    Eigen::Map<VectorXmp> R(Rab.data(), numpts);

    const scalar initial_radius = 0.8*math::sqrt(area/PI);
    R.setConstant(initial_radius);
    Rold.setConstant(0.9*initial_radius);
    Rab(numpts) = a;
    Rab(numpts + 1) = b;

    LBFGSParam<scalar> param;
    param.past = 100;
    param.delta = 1e-5;
    //param.epsilon = epsilon;
    LBFGSSolver<scalar> solver(param);

    Assembler assemble(ffft, S, k, viscosity, area, numpts, Rold);

    int count = 0;
    while(count < maxiter){
	std::cout << "Count = " << count << std::endl;

	scalar L;
	int niter = solver.minimize(assemble, Rab, L);

	std::cout << "	niter = " << niter << std::endl;
	std::cout << "	L = " << L << std::endl;

	scalar viscousterm = 0.5*viscosity*((R - Rold).dot(R - Rold));
	std::cout << "	viscousterm = " << viscousterm << std::endl;
	if(viscousterm < viscoustol && count > 3){
	    for(int z = 0; z < numpts + 1; ++z){
		outfile << Rab(z) << ",";
	    }
	    outfile << Rab(numpts + 1) << std::endl;

	    outfile.close();
	    break;
	}

	Rold = R;
	a = Rab(numpts);
	b = Rab(numpts + 1);

	for(int z = 0; z < numpts + 1; ++z){
	    outfile << Rab(z) << ",";
	}

	outfile << Rab(numpts + 1) << std::endl;
	count += 1;
    }
    if(count == maxiter){
	std::cout << " Pointid = " << pointid
	    << "\n Maximum number of viscous iterations reached!\n"
	    << " The solver had difficulty converging." << std::endl;
    }

    outfile.close();
}

int main(){

#ifdef USEMPREAL
    mpreal::set_default_prec(mpfr::digits2bits(digits));
#endif

    const int pointid =  0;
    const int numpts  = 50;
    const int maxiter =  3;

    const scalar 	     k = 1e9;
    const scalar             S = 2.0;
    const scalar           b = -0.99;
    const scalar          a = -0.495;
    const scalar       area = 0.6139;
    const scalar     viscosity = 2.0;
    const scalar     epsilon = 1e-8;
    const scalar  viscoustol = 1e-8;

    generatedfunc ffft = firstfundaformterms();

    std::cout.precision(digits);
    std::cout << std::fixed;

    solve(a, b, pointid, area, ffft, numpts, S, maxiter, k, viscosity,
	    viscoustol, epsilon);
    return 0;
}
