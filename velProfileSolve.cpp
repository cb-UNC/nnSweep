#include <iostream>
#include <iomanip>
#include <cmath>
#include <gsl/gsl_sf_hyperg.h>

using namespace std;

/******************** Global Constants *********************************
 ***********************************************************************/ 



////////////////////////////////////////////////////////////////////////
////// Function Declarations
////////////////////////////////////////////////////////////////////////

double gammaSolve(double etaI, double delta, double m, double n, double tau);
double velocitySolve(double etaI, double delta, double m, double n, double L, double delP, double gammaH, double gammaW);
void velocityProfile(double etaI, double delta, double m, double n, double H,  double L, double delP, int DNP, double v[] );


////// Main ////////////////////////////////////////////////////////////
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
	
	// Profile inputs
	int DNP = 6;
	int vLength = 2*DNP - 1;
	double v[vLength];

	// Calculate the profile
	velocityProfile(etaI, delta, m, n, H,  L, delP, DNP, v );

	// Checking a couple of velocities
	cout << "Velocity at the wall is:	" << v[0] << endl;
	cout << "Velocity at the centrum is:	" << v[DNP-1] << endl;
	cout << "Velocity halfway between is:	" << v[(DNP-1)/2] << endl;

	return 0;
}
////// End Main ////////////////////////////////////////////////////////






////////////////////////////////////////////////////////////////////////
//////	gammaSolve Function
////////////////////////////////////////////////////////////////////////
//
//	Variables are:
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
			
//		cout << "gamma = " << gamma;
//		cout << ";	iterate = " << iterate;
//		cout << ";	error = " << error << endl;

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

//		cout << "gamma = " << gamma;
//		cout << ";	iterate = " << iterate;
//		cout << ";	error = " << error << endl;

		if (iterate == itLimit) { cout << "Warning: Iteration Limit Reached!\n\n"; }
	}

//	cout << "\n\n\n";


	return gamma;
}


////////////////////////////////////////////////////////////////////////
//////	velocitySolve Function
////////////////////////////////////////////////////////////////////////
//
//	Variables are:
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

//	cout << "I for velocity = " << I;
//	cout << ";	velocity = " << velocity << "\n";

	return velocity;
}


////////////////////////////////////////////////////////////////////////
//////	velocityProfile Function
////////////////////////////////////////////////////////////////////////
//
//	Variables are:
//
//	etaI	-  Infinite shear viscosity
//	delta 	- difference between infinite and zero shear viscosities
//	m 	- model parameter, from (m * gamma)^n
//	n 	- model parameter, from (m * gamma)^n
//	H 	- Distance to the wall from the centrum
//	L 	- length of the system
//	delP 	- change in pressure across the system
//	DNP 	- Discrete Number of Points across the cross section
//	v[]	- The array which will contain the profile
//
//	********* Notes! ***********************************************
//	The DNP will be the number of points from the centrum to the
//	wall. However, the array containing the profile will have the
//	full profile, from wall to wall; this is not as complicated as 
//	it sounds due to symmetry about the centrum. Thus, the length of 
//	the array that stores the profile will be 2*DNP-1. 
//	****************************************************************
//
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
void velocityProfile(double etaI, double delta, double m, double n, double H,  double L, double delP, int DNP, double v[] )
{
	// Calculate necessary lengths
	int N = DNP - 1;
	int vLength = 2*DNP - 1;

	// Calculate values at the wall
	double tauW = H * delP / L;
	double gammaW = gammaSolve(etaI, delta, m, n, tauW);

	// Initialize any re-used variables
	double h = 0.0;
	double tauH = 0.0;
	double gammaH = 0.0;
	double vH = 0.0;

	for (int i = 0; i <= N; i++) {
		h = (N - i) * H / N;
		tauH = h * delP / L;
		gammaH = gammaSolve(etaI, delta, m, n, tauH);
		vH = velocitySolve(etaI, delta, m, n, L, delP, gammaW, gammaH);

		v[i] = vH;
		if (i < N) { v[2*N-i] = vH;  }

		cout << "At distance " << h << " from the centrum, velocity equals:	";
		cout << vH << "\n";
	}


	// Output
	double x = 0.0;
	
	cout << "\nThe velocities across the cross section are:\n\n";
	cout << "x	v	\n";
	for (int i = 0; i < vLength; i++) {
		x = i * H / N;
		cout << x << "       " << v[i] << "\n";	
	}

	cout << "\n\n";
}























