#include <iostream>
#include <cmath>
#include "Rate_param.h"
#include "Useful.h"
/*This function file contains the necessary functions to determine
the rates of the binding and unbinding for integrins from the cell 
membrane to the ECM */
//Function dictating the binding on the upper layers of the cell
double k_bind(double d)
{
    double k_b;
    double f_1 = 0.02;
    double h = 10.0;
    double y;
    if (d >= 0.0 && d <= h)
    {
        y = f_1 * (1.0 - (d/h)); 
    }
    else
    {
        y = 0.0;
    }
    k_b = y;
    return k_b;
}
//Function dictating the unbinding of the biomaterial and the cell
double k_unbind(double d)
{
    double s = 15.0;
    double g1 = 0.005;
    double g2 = 0.01;
    double g3 = 70.0; 
    
    double h = 10.0; 
    double k_u;
    double y;
    if (d >= 0.0 && d <= s)
    {
        y = g1 + (g2 * (d/h));
    }
    else
    {
        y = g3;
    }
    k_u = y;
    return k_u;
}
//Function dictating the unbinding between the cell and the AFM probe
double k_unbind_U(double d)
{
    double s = 15.0;
    double g1 = 0.02;
    double g2 = 0.01;
    double g3 = 1.0;

    double h = 10.0;
    double k_u;
    double y;
    if (d >= 0.0 && d <= s)
    {
        y = g1 + (g2 * (d / h)); 
    }
    else
    {
        y = g3;
    }
    k_u = y;
    return k_u;
}
//Function that simulates the biomaterial transitioning from the extended
//state to the contracted state
double k_e_c(int lambda)
{
    double y;
    if (lambda == 1)
    {
        y = 1.0e-04;
    }
    else
    {
        y = 1.0e-04;
    }
    return y;
}
//Function that simulates the biomaterial transitioning from the contracted
//state to the extended state
double k_c_e(int lambda)
{
    double y;
    if (lambda == 1)
    {
        y = 1.0e-04;
    }
    else
    {
        y = 1.0e-04;
    }
    return y;
}
