#ifndef __SETTINGS_H__
#define __SETTINGS_H__

#include <array>
#include <Eigen/Core>

#define USEMPREAL

#ifdef USEMPREAL
#include <mpreal.h>
#include <unsupported/Eigen/MPRealSupport>
#else
#define USE_MATH_DEFINES
#include <math.h>
#endif

#ifdef USEMPREAL
using mpfr::mpreal;
namespace math = mpfr;
typedef mpreal scalar;
inline int const digits = 32;
inline mpreal const PI = mpfr::const_pi(digits);
#else
namespace math = std;
typedef double scalar;
inline double const PI = 3.1415926535897932;
inline int const digits = 16;
#endif


typedef Eigen::Matrix<scalar, Eigen::Dynamic, Eigen::Dynamic>  MatrixXmp;
typedef Eigen::Matrix<scalar, Eigen::Dynamic, 1>	       VectorXmp;


//!
//! Five-point Gauss quadrature data for numerical integration.
//!
class FivePointGaussQuad{
    public:
	~FivePointGaussQuad() = default;
	static FivePointGaussQuad& getQuadObj(){
	    static FivePointGaussQuad fpgq;

	    // Set the points
	    fpgq.points[0] = -1.0/scalar(3.0) * math::sqrt(5 + 2 *
			    math::sqrt(10.0/scalar(7.0)));
	    fpgq.points[1] = -1.0/scalar(3.0) * math::sqrt(5 - 2 *
			    math::sqrt(10.0/scalar(7.0)));
	    fpgq.points[2] = 0.0;
	    fpgq.points[3] = -1.0 * fpgq.points[1];
	    fpgq.points[4] = -1.0 * fpgq.points[0];

	    // Set the weights
	    fpgq.weights[0] = (322.0 - 13.0*math::sqrt(70.0))/scalar(900.0);
	    fpgq.weights[1] = (322.0 + 13.0*math::sqrt(70.0))/scalar(900.0);
	    fpgq.weights[2] = 128.0/scalar(225.0);
	    fpgq.weights[3] = fpgq.weights[1];
	    fpgq.weights[4] = fpgq.weights[0];

	    // Set the evalpoints
	    for(int i=0; i < 5; ++i){
		fpgq.evalpoints[i] = 1 + fpgq.points[i];
	    }

	    return fpgq;
	}

	auto getEvalPoints(){ return evalpoints;}
	auto getWeights(){ return weights;}
    private:
	FivePointGaussQuad() = default;

	std::array<scalar, 5> weights;
	std::array<scalar, 5> points;
	std::array<scalar, 5> evalpoints;

}; // FivePointGaussQuad

#endif // __SETTINGS_H__
