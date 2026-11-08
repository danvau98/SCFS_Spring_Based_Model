#pragma once
//Finds the maximum value in a 1D array
double abs_max_val_1d(double* vec, int N);
//Destroys the memory of a 2D array of the type integer
void destroyTwoDimenArray_int(int** ptr, int row, int col);
//Destroys the memory of a 2D array of the data type double
void destroyTwoDimenArray_double(double** ptr, int row, int col);
//sgn function 
double sgn(double v);
//Fills in the time domain for an array with uniform distribution
double* tDomain(double* t, double at, double bt, double dt, int Nt);
//Fills in the spatial domain for an array with unifrom distribution
double* xDomain(double* x, double ax, double bx, double dx, int Nx);
//Integrates a function using the trapezoidal method
double trapz(double* f, double dx, int N);
//Calculate a distance using the Pythagorus method
double two_dimensional_dist(double x, double L);
//Finds the length of an array
int getLengthFromFile(const std::string& filename);