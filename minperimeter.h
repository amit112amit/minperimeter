#ifndef __MINPERIMETER_H__
#define __MINPERIMETER_H__

#include <functional>
#include "settings.h"
#include "generatedcode.h"
#include "generatedEFGcode.h"

typedef std::function<tuple12(const scalar, const scalar, const scalar,
	const scalar)> generatedfunc;
typedef std::function<tuple3(const scalar, const scalar, const scalar,
	const scalar)> generatedEFGfunc;
typedef std::tuple<scalar, scalar> tuple2;
typedef std::tuple<scalar, scalar, scalar, scalar, scalar> tuple5;

tuple5 areawithjac(const scalar, const scalar, const scalar, const scalar,
	const scalar, const scalar, generatedfunc);

scalar totalarea(const VectorXmp&, const scalar, const scalar, generatedEFGfunc);
scalar elementarea(const scalar, const scalar, const scalar, const scalar rj,
	const scalar, const scalar, generatedEFGfunc);
tuple2 centroid(const VectorXmp&, const scalar, const scalar, generatedEFGfunc);

tuple5 lengthwithjac(const scalar, const scalar, const scalar, const scalar,
	const scalar, const scalar, generatedfunc);
scalar totallength(const VectorXmp &, const scalar, const scalar,
	generatedEFGfunc);

std::tuple<scalar, VectorXmp> assemble(VectorXmp&, const scalar, generatedfunc f,
		const scalar, const VectorXmp&, const scalar, const scalar);

//!Calculates both the objective function and the jacobian.
//!Notes:
//!    We will use Gauss-quadrature with 1 point over the angle t and
//!    Gauss-quadrature of degree with 5 points over the radial distance R.
class Assembler{
    public:
	Assembler(const generatedfunc func, const scalar systemsize,
		const scalar penaltycoeff, const scalar viscosity,
		const scalar areaconstraint, const int numtri, VectorXmp &R_old)
	    : Rold(R_old), ffft(func), S(systemsize), k(penaltycoeff),
	    v(viscosity), A0(areaconstraint), N(numtri), T(VectorXmp(numtri)){
		assert(N > 2);
		assert(Rold.size() == N);
		scalar dt = 2.0*PI/N;
		T.setLinSpaced(N, 0.0, 2.0*PI - dt);
	    }

	scalar operator()(VectorXmp& Rab, VectorXmp& grad){
	    Eigen::Map<VectorXmp> R(Rab.data(), N);

	    scalar a = Rab(N);
	    scalar b = Rab(N + 1);

	    scalar L = 0.0;                     //!< The perimeter
	    scalar A = 0.0;                     //!< The surface area

	    scalar dLda = 0.0;
	    scalar dLdb = 0.0;
	    scalar dAda = 0.0;
	    scalar dAdb = 0.0;

	    grad = VectorXmp::Zero(N + 2);
	    VectorXmp dAdr = VectorXmp::Zero(N);

	    for(int i=0; i < N; ++i){
		int j = (i < N-1)? i + 1 : 0;
		scalar ri = R(i);
		scalar rj = R(j);
		scalar ti = T(i);
		scalar tj = (i < N-1)? T(j) : 2.0*PI;

		scalar Li, dLdri, dLdrj, dLida, dLidb; 
		std::tie(Li, dLdri, dLdrj, dLida, dLidb) =
		    lengthwithjac(ti, tj, ri, rj, a, b, ffft);
		L += Li;
		grad(i) += dLdri;
		grad(j) += dLdrj;
		dLda += dLida;
		dLdb += dLidb;

		scalar Ai, dAdri, dAdrj, dAida, dAidb;
		std::tie(Ai, dAdri, dAdrj, dAida, dAidb) =
		    areawithjac(ti, tj, ri, rj, a, b, ffft);
		A += Ai;
		dAdr(i) += dAdri;
		dAdr(j) += dAdrj;
		dAda += dAida;
		dAdb += dAidb;
	    }

	    // The objective function to be minimized
	    scalar F = L/S + 0.5*(k/math::pow(S, 4))*(A - A0)*(A - A0) +
		    0.5*(v/N/S/S)*((R - Rold).dot(R - Rold));

	    // Add contributions from the area and viscous terms to the grad 
	    for(int i = 0; i < N; ++i){
		grad(i) /= S;
		grad(i) += (k/math::pow(S,4))*(A-A0)*dAdr(i) +
		    (v/N/S/S)*(R(i) - Rold(i));
	    }

	    grad(N) = dLda/S + (k/math::pow(S,4))*(A-A0)*dAda;
	    grad(N+1) = dLdb/S + (k/math::pow(S,4))*(A-A0)*dAdb;

	    return F;
	}

    private:
	int			   N;  //!< Number of triangles
	scalar 			   S;  //!< System size
	scalar 			   k;  //!< Penalty coefficient
	scalar			   v;  //!< Viscosity
	scalar 			  A0;  //!< Area constraint value
	generatedfunc 		ffft;  //!< Function to calculate E, F, G etc.
	VectorXmp	       &Rold;  //!< Previous radial positions
	VectorXmp	           T;  //!< Previous radial positions
};

#endif //__MINPERIMETER_H__
