//
//  integral.hpp
//  cp_4_13
//
//  Created by Nikolay Tikhonov on 09.05.16.
//  Copyright © 2016 Nikolay Tikhonov. All rights reserved.
//

#ifndef integral_hpp
#define integral_hpp

#include <iostream>
#include <cmath>
#include "matrix.hpp"
#include "SLE.h"

double Func(double x);
double Pol(double x, double A, double B, double C, double D);
Matrix Moments_NC(double a, double b);
//Matrix Moments_G(double a, double b);
double dg(double x, double A, double B, double C);
double NFR(double x, double a, double b, double A, double B, double C, double D);
void Newton(double* m, double a, double b, double A, double B, double C, double D);
double Newton_Cotes(double a, double b);
double Gauss(double a, double b);


#endif /* integral_hpp */
