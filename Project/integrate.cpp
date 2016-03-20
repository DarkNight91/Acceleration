#ifndef INTEGRATE_CPP
#define INTEGRARTE_CPP

#include <math.h>
#include "integrate.h"
#include "iostream"

using namespace std;

double integrate(double x, double mean, double sigma){
    float step = 0.1;
    int iterations;
    double index;
    double lli;

    double sum=0.0;

    lli = -2.0+mean-3*sigma;

    if(x < lli){
        return 0.00;
    }

    iterations = (int)(ceil((x-lli)/step));

    for(int i=0; i < iterations; i++){
        index = pow( ((lli + i*step - step/2) - mean)/sigma, 2)/2;
        sum += step*pow(2.718281828, -index);
    }

    return sum/(sigma * sqrt(2*3.141592654));
}

#endif
