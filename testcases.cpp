#include <list>
#include <iostream>
#include "minperimeter.h"

void test_areajac(){

    std::cout << "********** Running test_areajac() **********" << std::endl;

    scalar h = 1e-6;

    scalar ti = 0.5*PI;
    scalar tj = 0.5*PI + 0.01*PI;
    scalar ri = 1.1;
    scalar rj = 1.2;
    scalar a = 2.0;
    scalar b = 3.0;

    generatedfunc ffft = firstfundaformterms();

    scalar Ap, Am, fromfunc, fromdiff;

    std::cout.precision(digits);
    std::cout << std::fixed;

    // Test dAdri
    std::tie(std::ignore, fromfunc, std::ignore, std::ignore, std::ignore) = 
	areawithjac(ti, tj, ri, rj, a, b, ffft);
    std::tie(Ap, std::ignore, std::ignore, std::ignore, std::ignore) = 
	areawithjac(ti, tj, ri + h, rj, a, b, ffft);
    std::tie(Am, std::ignore, std::ignore, std::ignore, std::ignore) = 
	areawithjac(ti, tj, ri - h, rj, a, b, ffft);
    fromdiff = 0.5*(Ap - Am)/h;
    std::cout <<" From function, dAdri = " << fromfunc << std::endl;
    std::cout <<" From approx  , dAdri = " << fromdiff << std::endl;
    std::cout << std::endl;

    // Test dAdrj
    std::tie(std::ignore, std::ignore, fromfunc, std::ignore, std::ignore) = 
	areawithjac(ti, tj, ri, rj, a, b, ffft);
    std::tie(Ap, std::ignore, std::ignore, std::ignore, std::ignore) = 
	areawithjac(ti, tj, ri, rj + h, a, b, ffft);
    std::tie(Am, std::ignore, std::ignore, std::ignore, std::ignore) = 
	areawithjac(ti, tj, ri, rj - h, a, b, ffft);
    fromdiff = 0.5*(Ap - Am)/h;
    std::cout <<" From function, dAdrj = " << fromfunc << std::endl;
    std::cout <<" From approx  , dAdrj = " << fromdiff << std::endl;
    std::cout << std::endl;

    // Test dAda
    std::tie(std::ignore, std::ignore, std::ignore, fromfunc, std::ignore) = 
	areawithjac(ti, tj, ri, rj, a, b, ffft);
    std::tie(Ap, std::ignore, std::ignore, std::ignore, std::ignore) = 
	areawithjac(ti, tj, ri, rj, a + h, b, ffft);
    std::tie(Am, std::ignore, std::ignore, std::ignore, std::ignore) = 
	areawithjac(ti, tj, ri, rj, a - h, b, ffft);
    fromdiff = 0.5*(Ap - Am)/h;
    std::cout <<" From function, dAda = " << fromfunc << std::endl;
    std::cout <<" From approx  , dAda = " << fromdiff << std::endl;
    std::cout << std::endl;

    // Test dAdb
    std::tie(std::ignore, std::ignore, std::ignore, std::ignore, fromfunc) = 
	areawithjac(ti, tj, ri, rj, a, b, ffft);
    std::tie(Ap, std::ignore, std::ignore, std::ignore, std::ignore) = 
	areawithjac(ti, tj, ri, rj, a, b + h, ffft);
    std::tie(Am, std::ignore, std::ignore, std::ignore, std::ignore) = 
	areawithjac(ti, tj, ri, rj, a, b - h, ffft);
    fromdiff = 0.5*(Ap - Am)/h;
    std::cout <<" From function, dAdb = " << fromfunc << std::endl;
    std::cout <<" From approx  , dAdb = " << fromdiff << std::endl;
    std::cout << std::endl;
    std::cout << "*******************************************" << std::endl;

    return;
}


void test_totalarea(){
    std::cout << "******************** IMPORTANT ********************"
	<< std::endl;
    std::cout << "* To run `test_totalarea()` first generate code for "
	<< std::endl << "* `z = x^2 + y^2` using `makecppcode.py` "
	<< "Python script." << std::endl
	<< "***************************************************"
	<< std::endl;

    generatedEFGfunc efg = EFGFunc();

    VectorXmp R(100);
    R.setOnes();
    scalar calculated_area = totalarea(R, 0.0, 0.0, efg);
    scalar analytical_area = (5.0*math::sqrt(5.0) - 1.0)*PI/6.0;

    std::cout.precision(digits);
    std::cout << std::fixed;

    std::cout<< " Analytical area = " << analytical_area << std::endl;
    std::cout<< " Calculated area = " << calculated_area << std::endl;
    std::cout << std::endl;
}

void test_elementarea(){
    std::cout << "******************** IMPORTANT ********************"
	<< std::endl;
    std::cout << "* To run `test_elementarea()` first generate code for "
	<< std::endl << "* `z = x^2 + y^2` using `makecppcode.py` "
	<< "Python script." << std::endl
	<< "***************************************************"
	<< std::endl;

    generatedEFGfunc efg = EFGFunc();

    scalar ti = 0.0;
    scalar tj = PI/50.0;
    scalar ri = 1.0;
    scalar rj = 1.0;
    scalar a = 0.0;
    scalar b = 0.0;

    scalar analytical_area = (5*math::sqrt(5) - 1)*PI/600.0;
    scalar calculated_area = elementarea(ti, tj, ri, rj, a, b, efg);

    std::cout.precision(digits);
    std::cout << std::fixed;

    std::cout<< " Analytical element area = " << analytical_area << std::endl;
    std::cout<< " Calculated element area = " << calculated_area << std::endl;
    std::cout << std::endl;
}


void test_centroid(){
    std::cout << "* To run `test_centroid()` first generate code for "
	<< std::endl << "* `z = x^2 + y^2` using `makecppcode.py` "
	<< "Python script." << std::endl
	<< "***************************************************"
	<< std::endl;

    generatedEFGfunc efg = EFGFunc();
    VectorXmp R(100);
    R.setOnes();

    scalar Xc, Yc;
    std::tie(Xc, Yc) = centroid(R, 0.0, 0.0, efg);

    std::cout.precision(digits);
    std::cout << std::fixed;
    std::cout << "Analytical centroid = (0, 0)." << std::endl;
    std::cout << "Calculated centroid Xc = " << Xc << std::endl;
    std::cout << "Calculated centroid Yc = " << Xc << std::endl;
    std::cout << std::endl;
}

void test_lengthjac(){

    std::cout << "********** Running test_lengthjac() **********" << std::endl;

    scalar h = 1e-6;

    scalar ti = 0.5*PI;
    scalar tj = 0.5*PI + 0.01*PI;
    scalar ri = 1.1;
    scalar rj = 1.2;
    scalar a = 2.0;
    scalar b = 3.0;

    generatedfunc ffft = firstfundaformterms();

    scalar Lp, Lm, fromfunc, fromdiff;

    std::cout.precision(digits);
    std::cout << std::fixed;

    // Test dLdri
    std::tie(std::ignore, fromfunc, std::ignore, std::ignore, std::ignore) = 
	lengthwithjac(ti, tj, ri, rj, a, b, ffft);
    std::tie(Lp, std::ignore, std::ignore, std::ignore, std::ignore) = 
	lengthwithjac(ti, tj, ri + h, rj, a, b, ffft);
    std::tie(Lm, std::ignore, std::ignore, std::ignore, std::ignore) = 
	lengthwithjac(ti, tj, ri - h, rj, a, b, ffft);
    fromdiff = 0.5*(Lp - Lm)/h;
    std::cout <<" From function, dLdri = " << fromfunc << std::endl;
    std::cout <<" From approx  , dLdri = " << fromdiff << std::endl;
    std::cout << std::endl;

    // Test dLdrj
    std::tie(std::ignore, std::ignore, fromfunc, std::ignore, std::ignore) = 
	lengthwithjac(ti, tj, ri, rj, a, b, ffft);
    std::tie(Lp, std::ignore, std::ignore, std::ignore, std::ignore) = 
	lengthwithjac(ti, tj, ri, rj + h, a, b, ffft);
    std::tie(Lm, std::ignore, std::ignore, std::ignore, std::ignore) = 
	lengthwithjac(ti, tj, ri, rj - h, a, b, ffft);
    fromdiff = 0.5*(Lp - Lm)/h;
    std::cout <<" From function, dLdrj = " << fromfunc << std::endl;
    std::cout <<" From approx  , dLdrj = " << fromdiff << std::endl;
    std::cout << std::endl;

    // Test dLda
    std::tie(std::ignore, std::ignore, std::ignore, fromfunc, std::ignore) = 
	lengthwithjac(ti, tj, ri, rj, a, b, ffft);
    std::tie(Lp, std::ignore, std::ignore, std::ignore, std::ignore) = 
	lengthwithjac(ti, tj, ri, rj, a + h, b, ffft);
    std::tie(Lm, std::ignore, std::ignore, std::ignore, std::ignore) = 
	lengthwithjac(ti, tj, ri, rj, a - h, b, ffft);
    fromdiff = 0.5*(Lp - Lm)/h;
    std::cout <<" From function, dLda = " << fromfunc << std::endl;
    std::cout <<" From approx  , dLda = " << fromdiff << std::endl;
    std::cout << std::endl;

    // Test dLdb
    std::tie(std::ignore, std::ignore, std::ignore, std::ignore, fromfunc) = 
	lengthwithjac(ti, tj, ri, rj, a, b, ffft);
    std::tie(Lp, std::ignore, std::ignore, std::ignore, std::ignore) = 
	lengthwithjac(ti, tj, ri, rj, a, b + h, ffft);
    std::tie(Lm, std::ignore, std::ignore, std::ignore, std::ignore) = 
	lengthwithjac(ti, tj, ri, rj, a, b - h, ffft);
    fromdiff = 0.5*(Lp - Lm)/h;
    std::cout <<" From function, dLdb = " << fromfunc << std::endl;
    std::cout <<" From approx  , dLdb = " << fromdiff << std::endl;
    std::cout << std::endl;
    std::cout << "*********************************************" << std::endl;

    return;
}

void test_assembler(){

    int N = 4;
    VectorXmp Rab(N + 2);
    Rab.setConstant(3.0);
    Rab(N) = 1.0;
    Rab(N + 1) = 2.5;

    VectorXmp Rold(N);
    Rold.setConstant(2.8);
    
    const scalar A0 = 50.0;
    const scalar S = 20.0;
    const scalar k = 1e3;
    const scalar v = 2.0;

    generatedfunc ffft = firstfundaformterms();

    Assembler assemble(ffft, S, k, v, A0, N, Rold);

    std::cout << "Testing `assemble()...` " << std::endl;

    auto herr = std::list<scalar>{1e-16, 1e-14, 1e-12, 1e-10,
	1e-8, 1e-6, 1e-4, 1e-2, 1e+1};
    for(auto const &h : herr){ 
	scalar Fplus, Fminus;
	VectorXmp jacfromfunc(N + 2), jacfromdiff(N + 2), temp(N + 2);

	// Calculate the analytical derivative
	assemble(Rab, jacfromfunc);

	// Calculate the numerical derivative
	jacfromdiff.setZero();
	for(int i = 0; i < N + 2; ++i){
	    Rab(i) += h;
	    Fplus = assemble(Rab, temp);
	    Rab(i) -= 2.0*h;
	    Fminus = assemble(Rab, temp);
	    jacfromdiff(i) += 0.5*(Fplus - Fminus)/h;
	    Rab(i) += h;
	}
	
	std::cout.precision(digits);
	std::cout << std::fixed;
	std::cout << "	h = " << h
	    << " error norm = "
	    << ((jacfromfunc - jacfromdiff).norm())/jacfromdiff.norm()
	    << std::endl;
    }
    std::cout << std::endl;
}

int main(){

#ifdef USEMPREAL
    mpreal::set_default_prec(mpfr::digits2bits(digits));
#endif

    test_areajac();
    test_lengthjac();
    //test_totalarea();
    //test_elementarea();
    //test_centroid();
    test_assembler();
    return 0;

}
