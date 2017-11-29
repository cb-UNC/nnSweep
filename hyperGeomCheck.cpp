#include <iostream>
#include <iomanip>
#include <cmath>
#include <gsl/gsl_sf_hyperg.h>

using namespace std;

int main()
{
        // Add fluid parameters
        double etaO = 10.0;
        double etaI = 0.01;
        double m = 1.0;
        double n = 0.666;

        // Add known gamma
        double gamma = 1.0;
        double mgamma = m*gamma;
        double f = pow ( mgamma, n);

        // Inputs for hypergeometric function
        double aH = 1.0;
        double bH = 2.0/n;
        double cH = (n+1.0)/n;
        double dH = -f;

        double H = gsl_sf_hyperg_2F1 (aH, bH, cH, dH);

        cout << fixed;
        cout << "The hypergeometric function using GSL gives: " << H << endl;

        return 0;
}
