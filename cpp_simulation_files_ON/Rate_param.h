#pragma once
//Function dictating the binding on the upper layers of the cell
double k_bind(double d);
//Function dictating the binding of the biomaterial in the extnded state
double k_bind_e(double d);
//Function dictating the binding of the biomaterial in the contracted state
double k_bind_c(double d);
//Function dictating the unbinding of the biomaterial and the cell
double k_unbind(double d);
//Function dictating the unbinding between the cell and the AFM probe
double k_unbind_U(double d);
//Function that simulates the biomaterial transitioning from the extended
//state to the contracted state
double k_e_c(int lambda1, double l);
//Function that simulates the biomaterial transitioning from the contracted
//state to the extended state
double k_c_e(int lambda2, double l);