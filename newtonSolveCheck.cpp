#include <iostream>
#include <iomanip>
#include <cmath>
//#include <gsl/gsl_sf_hyperg.h>

using namespace std;

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
	double tau = H * delP / L;
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
	}

	cout << "gamma = " << gamma << endl;
	cout << "iterate = " << iterate << endl;
	cout << "error = " << error << endl;

	return 0;
}

















