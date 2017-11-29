#include <iostream>
#include <iomanip>
#include <cmath>
//#include <gsl/gsl_sf_hyperg.h>

using namespace std;

////////////////////////////////////////////////////////////////////////
////// Function Declarations
////////////////////////////////////////////////////////////////////////

double gammaSolve(double etaI, double delta, double m, double n, double tau);

int main()
{
	// Add fluid parameters
	double etaO = 10.0;
	double etaI = 0.01;
	double m = 1.0;
	double n = 0.666;
	double delta = etaO - etaI;
	
	// System inputs
	double H = 1.0;
	double L = 10.0;
	double delP = 50.05;
	double tauW = H * delP / L;
	
	double gamma = 0.0;
	gamma = gammaSolve(etaI, delta, m, n, tauW);

	cout << "Final gamma = " << gamma << endl;

	return 0;
}

////////////////////////////////////////////////////////////////////////
//////	gammaSolve Function
////////////////////////////////////////////////////////////////////////
//
//	Variabls are 
//	etaI - Infinite shear viscosity
//	delta - difference between infinite and zero shear viscosities
//	m - model parameter, from (m * gamma)^n
//	n - model parameter, from (m * gamma)^n
//	tau - shear rate, calculate this before function
//
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
double gammaSolve(double etaI, double delta, double m, double n, double tau)
{
	////////////////////////////////////////////////////////////////
	//////	Fixed Point Iteration
	////////////////////////////////////////////////////////////////
	
	// Numerical Input
	double tolerance = 1.0;
	int itLimit = 5;
	int iterate = 0;
	bool check = true;
	double error = 10000.0;

	// Initialize
	double gammaOld = 10.0;
	double gamma = 0.0;
	double mGammaOld = 0.0;
	double f = 0.0;
	double g = 0.0;

	// Carry out our fixed point iterations
	while (check) {
		mGammaOld = m * gammaOld;
		f = pow(mGammaOld, n );
		g = 1 + f;
		gamma = etaI + delta / g;
		gamma = tau * pow(gamma, -1 );

		// Checks
		iterate++;
		error = gamma - gammaOld;
		error = abs(error );
		if (error < tolerance ) check = false;
		if (iterate > itLimit ) check = false;

		gammaOld = gamma;
			
		cout << "gamma = " << gamma;
		cout << ";	iterate = " << iterate;
		cout << ";	error = " << error << endl;

		if (iterate == itLimit) { cout << "Warning: Iteration Limit Reached!\n\n"; }
	}

	////////////////////////////////////////////////////////////////
	////// Newton's Method Solve	
	////////////////////////////////////////////////////////////////

	// Numerical Input
	tolerance = 0.0000000001;
	itLimit = 10;
	iterate = 0;
	check = true;
	error = 10000.0;

	// Initialize
	double nFunc = 0.0;
	double ndFunc = 0.0;

	// Carry out Newton's method
	while (check) {
		mGammaOld = m * gammaOld;
		f = pow(mGammaOld, n );
		g = 1 + f;

		nFunc = ( etaI + delta / g ) * gammaOld - tau ;
		ndFunc = etaI + delta / g - n*delta*f/ pow(g, 2 ) ;

		gamma = gammaOld - nFunc/ ndFunc;
		
		// Checks
		iterate++;
		error = gamma - gammaOld;
		error = abs(error );
		if (error < tolerance ) check = false;
		if (iterate > itLimit ) check = false;

		gammaOld = gamma;	

		cout << "gamma = " << gamma;
		cout << ";	iterate = " << iterate;
		cout << ";	error = " << error << endl;

		if (iterate == itLimit) { cout << "Warning: Iteration Limit Reached!\n\n"; }
	}


	return gamma;
}














