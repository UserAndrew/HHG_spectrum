#ifndef GROUND_STATE_H
#define GROUND_STATE_H
#include <iostream>
#include <fstream>
#include <fftw3.h>
#include <cmath>
#include <vector>

double V(double x);
double Psi_solution(double x);
void ground_state(double func[], const int N,
                 const double dx, std::vector<double> coord);
#endif // GROUND_STATE_H
