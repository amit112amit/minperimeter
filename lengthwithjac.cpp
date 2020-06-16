#include "minperimeter.h"

//!
//!Calculate surface arc-length and its derivatives.
//!
tuple5 lengthwithjac(const scalar ti, const scalar tj, const scalar ri,
	const scalar rj, const scalar a, const scalar b,
	generatedfunc ffft){

    scalar r = 0.5*(ri + rj);
    scalar t = 0.5*(ti + tj);
    scalar E, F, G, Er, Fr, Gr, Ea, Fa, Ga, Eb, Fb, Gb;
    std::tie(E, F, G, Er, Fr, Gr, Ea, Fa, Ga, Eb, Fb, Gb) = ffft(r, t, a, b);

    scalar dt = tj - ti;
    scalar rp = (rj - ri)/dt;
    scalar dL = math::sqrt(E*rp*rp + 2.0*F*rp + G);
    scalar L = dt*dL;

    scalar fac = dt/dL;
    scalar term1 = (E*rp + F)/dt;
    scalar term2 = 0.25*(Er*rp*rp + 2.0*Fr*rp + Gr);

    scalar dLdri = fac*(-term1 + term2);
    scalar dLdrj = fac*(term1 + term2);
    scalar dLda = 0.5*fac*(Ea*rp*rp + 2.0*Fa*rp + Ga);
    scalar dLdb = 0.5*fac*(Eb*rp*rp + 2.0*Fb*rp + Gb);

    return {L, dLdri, dLdrj, dLda, dLdb};
}


//!Calculate surface arc-length and its derivatives.
//
//! Parameters:
//! -----------
//! R: radial distances
//! a, b: center of the domain
//! getEFG: function to calculate the first fundamental form terms
//! 
//! Output:
//! -------
//! totallength: the perimeter
//
scalar totallength(const VectorXmp &R, const scalar a, const scalar b,
	generatedEFGfunc efg){

    int N = R.size();

    scalar dt = 2.0*PI/N;
    VectorXmp T(N);
    T.setLinSpaced(N, 0.0, 2.0*PI - dt);

    scalar L = 0.0;

    for(int i=0; i < N; ++i){
	int j = (i < N-1)? i + 1 : 0;
	scalar ri = R(i);
	scalar rj = R(j);
	scalar ti = T(i);
	scalar tj = (i < N-1)? T(j) : 2.0*PI;

        scalar r = 0.5*(ri + rj);
        scalar t = 0.5*(ti + tj);
        scalar E, F, G;
	std::tie(E, F, G) = efg(r, t, a, b);

        scalar rp = (rj - ri)/dt;
        scalar dL = math::sqrt(E*rp*rp + 2.0*F*rp + G);
        L += dt*dL;
    }

    return L;
}
