#include "minperimeter.h"

//! Calculate the area of an element and its derivatives wrt ri and rj.
//
//! We will use 5-point Gauss quadrature.
tuple5 areawithjac(const scalar ti, const scalar tj, const scalar ri,
	const scalar rj, const scalar a, const scalar b,
	generatedfunc ffft){

    auto quad_obj = FivePointGaussQuad::getQuadObj();
    auto evalpoints = quad_obj.getEvalPoints();
    auto weights = quad_obj.getWeights();

    scalar dt = tj - ti;
    scalar t = 0.5*(ti + tj);
    scalar r = 0.25*(ri + rj);
    scalar A = 0.0;
    scalar jac_r = 0.0;
    scalar jac_a = 0.0;
    scalar jac_b = 0.0;

    for(int i=0; i < 5; ++i){
	scalar re = r*evalpoints[i];
	auto [E, F, G, Er, Fr, Gr, Ea, Fa, Ga, Eb, Fb, Gb] = ffft(re, t, a, b);
	scalar f = math::sqrt(E*G - F*F);
	A += weights[i]*f;
	jac_r += weights[i]*evalpoints[i]*(E*Gr - 2*F*Fr + G*Er)/f;
	jac_a += weights[i]*(E*Ga + G*Ea - 2*F*Fa)/f;
	jac_b += weights[i]*(E*Gb + G*Eb - 2*F*Fb)/f;
    }
    scalar dAdri = dt*(0.125*r*jac_r + 0.25*A);
    scalar dAdrj = dAdri;
    scalar dAda = 0.5*dt*r*jac_a;
    scalar dAdb = 0.5*dt*r*jac_b;
    A *= dt*r;
    return {A, dAdri, dAdrj, dAda, dAdb};
}


//! Get total area of all elements of the domain.
//
//  @param R: Nx1 array of the radial distances.
//  @param a, b: center of the domain.
//  @param getEFG: a function that calculates E, F, G.
scalar totalarea(const VectorXmp &R, const scalar a, const scalar b,
	generatedEFGfunc getEFG){

    auto quad_obj = FivePointGaussQuad::getQuadObj();
    auto evalpoints = quad_obj.getEvalPoints();
    auto weights = quad_obj.getWeights();

    int N = R.size();

    scalar dt = 2.0*PI/N;
    VectorXmp T(N);
    T.setLinSpaced(N, 0.0, 2.0*PI - dt);
    
    scalar A = 0.0;

    for(int i=0; i < N; ++i){
	int j = (i < N-1)? i + 1 : 0;
	scalar ri = R(i);
	scalar rj = R(j);
	scalar ti = T(i);
	scalar tj = (i < N-1)? T(j) : 2.0*PI;

	scalar t = 0.5*(ti + tj);
	scalar r = 0.25*(ri + rj);
	scalar Ai = 0.0;

        for(int z=0; z < 5; ++z){
            scalar re = r*evalpoints[z];
	    scalar E, F, G;
	    std::tie(E, F, G) = getEFG(re, t, a, b);
            scalar f = math::sqrt(E*G - F*F);
            Ai += weights[z]*f;
	}

        A += Ai*dt*r;
    }

    return A;
}


//!Calculates the area of a triangular element.

//! Parameters:
//! -----------
//! a, b: coordinates of one of the vertices of the triangle.
//! ri, rj: side-lengths from (a, b).
//! ti, tj: angles along which the two sides are located wrt (a, b).
//! getEFG: a function that calculates E, F, G.
//! 
//! Output:
//! -------
//! area: area of the triangle
scalar elementarea(const scalar ti, const scalar tj, const scalar ri,
	const scalar rj, const scalar a, const scalar b, generatedEFGfunc getEFG){

    auto quad_obj = FivePointGaussQuad::getQuadObj();
    auto evalpoints = quad_obj.getEvalPoints();
    auto weights = quad_obj.getWeights();

    scalar dt = tj - ti;
    scalar t = 0.5*(ti + tj);
    scalar r = 0.25*(ri + rj);
    scalar area = 0.0;

    for(int i=0; i < 5; ++i){
        scalar re = r*evalpoints[i];
	scalar E, F, G;
	std::tie(E, F, G) = getEFG(re, t, a, b);
        area += weights[i]*math::sqrt(E*G - F*F);
    }

    area *= dt*r;

    return area;
}


//! Find the centroid of the domain using weighted average of centroid
//! of the triangular elements.
//! 
//! Parameters:
//! -----------
//! a, b: x-y coordinates of origin of polar reference frame
//! R: N x 1 array of radial distances of N points along circumference
//!    of the domain
//! efg: a function that returns E, F and G i.e. the coefficients of
//!      the first fundamental form
//! 
//! Output:
//! -------
//! Xc, Yc: coordinates of the centroid.
std::tuple<scalar, scalar> centroid(const VectorXmp &R, const scalar a,
	const scalar b, generatedEFGfunc efg){

    int N = R.size();

    scalar dt = 2.0*PI/N;
    VectorXmp T(N);
    T.setLinSpaced(N, 0.0, 2.0*PI - dt);

    scalar area = 0.0;
    scalar Xc = 0.0;
    scalar Yc = 0.0;

    for(int i=0; i < N; ++i){
	int j = (i < N-1)? i + 1 : 0;
	scalar ri = R(i);
	scalar rj = R(j);
	scalar ti = T(i);
	scalar tj = (i < N-1)? T(j) : 2.0*PI;

        scalar currarea = elementarea(ti, tj, ri, rj, a, b, efg);
        area += currarea;

        //The coordinates of the vertices of the triangle
        scalar xi = a + ri*math::cos(ti);
        scalar yi = b + ri*math::sin(ti);
        scalar xj = a + rj*math::cos(tj);
        scalar yj = b + rj*math::sin(tj);

        Xc += currarea*(a + xi + xj);
        Yc += currarea*(b + yi + yj);
    }

    Xc /= 3*area;
    Yc /= 3*area;

    return {Xc, Yc};
}
