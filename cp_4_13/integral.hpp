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

Matrix Moments(double a, double b);
double Newton_Cotes(double a, double b);
double Gauss(double a, double b);

#endif /* integral_hpp */
