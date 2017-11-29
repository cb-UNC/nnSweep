#include <iostream>
#include <iomanip>
#include <cmath>
#include <gsl/gsl_sf_hyperg.h>

using namespace std;

////////////////////////////////////////////////////////////////////////
////// Function Declarations
////////////////////////////////////////////////////////////////////////

double gammaSolve(double etaI, double delta, double m, double n, double tau);
double velocitySolve(double etaI, double delta, double m, double n, double L, double delP, double gammaH, double gammaW);

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
	
	// Solve for gamma at the wall
	double gammaW = gammaSolve(etaI, delta, m, n, tauW);

	// Add distance gamma at new distance from centrum
	double h = H / 2.0;
	double tauH = h * delP / L;
	double gammaH = gammaSolve(etaI, delta, m, n, tauH);

	double vH = velocitySolve(etaI, delta, m, n, L, delP, gammaW, gammaH);




	return 0;
}

////////////////////////////////////////////////////////////////////////
//////	gammaSolve Function
////////////////////////////////////////////////////////////////////////
//
//	Variabls are:
//
//	etaI 	- Infinite shear viscosity
//	delta 	- difference between infinite and zero shear viscosities
//	m 	- model parameter, from (m * gamma)^n
//	n 	- model parameter, from (m * gamma)^n
//	tau 	- shear rate, calculate this before function
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

	cout << "\n\n\n";


	return gamma;
}


////////////////////////////////////////////////////////////////////////
//////	velocitySolve Function
////////////////////////////////////////////////////////////////////////
//
//	Variabls are:
//
//	etaI	-  Infinite shear viscosity
//	delta 	- difference between infinite and zero shear viscosities
//	m 	- model parameter, from (m * gamma)^n
//	n 	- model parameter, from (m * gamma)^n
//	L 	- length of the system
//	delP 	- change in pressure across the system
//	gammaH 	- shear rate at distance h from the centrum
//	gammaW 	- shear rate at distance H from the centrum (the wall)
//
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
double velocitySolve(double etaI, double delta, double m, double n, double L, double delP, double gammaH, double gammaW)
{
	// Generate inputs for I at distance h
	double gamma = gammaH;
	double mGamma = m * gamma;
	double f = pow(mGamma, n );
	double g = 1 + f; 

	// Calculate hypergeometric function for distance h
	double aF = 1.0;
	double bF = 2.0 / n;
	double cF = ( n + 1.0 ) / n;
	double dF = -f;

	double F = gsl_sf_hyperg_2F1 ( aF, bF, cF, dF );

	// Calculate I at distance h
	double iH = gamma*gamma * ( delta*g*F - 2*delta - etaI*g  ) / ( 2 * g );
	
	// Generate inputs for I at the wall
	gamma = gammaW;
	mGamma = m * gamma;
	f = pow(mGamma, n );
	g = 1 + f; 

	// Calculate hypergeometric function for distance h
	dF = -f;
	F = gsl_sf_hyperg_2F1 ( aF, bF, cF, dF );

	// Calculate I at distance h
	double iW = gamma*gamma * ( delta*g*F - 2*delta - etaI*g  ) / ( 2 * g );

	// Calculate I and velocity
	double I = iW - iH;
	double velocity = L * I / delP;

	cout << "I for velocity = " << I;
	cout << ";	velocity = " << velocity << "\n";

	return velocity;
}





























