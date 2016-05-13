//
//  integral.cpp
//  cp_4_13
//
//  Created by Nikolay Tikhonov on 09.05.16.
//  Copyright © 2016 Nikolay Tikhonov. All rights reserved.
//

#include "integral.hpp"
double Function(double x){
    return 2*cos(3.5*x)*exp(5*x/3)+3*sin(1.5*x)*exp(-4*x)+3;
}

Matrix Moments_NC(double a, double b){
    Matrix M(3,1);
    
    M.Add(0, 0, 5/4*(pow((b-1.5),(4/5)) - pow((a-1.5),(4/5))));
    M.Add(1, 0, 5/9*(pow((b-1.5),(9/5)) - pow((a-1.5),(9/5)))-3/2*M.Get(0, 0));
    M.Add(2, 0, pow((1.5),2*M.Get(0, 0)) - pow(5/3*(b-1.5),(9/5)) + pow(5/3*(a-1.5),(9/5)) +5/14*(pow((b-1.5),(14/5)) - pow((a-1.5),(14/5))));

    return M;
}

Matrix Moments_G(double a, double b){
    Matrix M(6,1);
    M.Add(0, 0, 5/4*(pow((b-1.5),(4/5)) - pow((a-1.5),(4/5))));
    M.Add(1, 0, 5/9*(pow((b-1.5),(9/5)) - pow((a-1.5),(9/5)))-3/2*M.Get(0, 0));
    M.Add(2, 0, pow((1.5),2*M.Get(0, 0)) - pow(5/3*(b-1.5),(9/5)) + pow(5/3*(a-1.5),(9/5)) +5/14*(pow((b-1.5),(14/5)) - pow((a-1.5),(14/5))));
    M.Add(3, 0, 1.1);
    M.Add(4, 0, 1.1);
    M.Add(5, 0, 1.1);
    return M;
}

double Newton_Cotes(double a, double b){
    Matrix m(Moments_NC(a,b));
    Matrix x(3, 1);
    
    x.Add(0,0,a);
    x.Add(1,0,(a + b)/2);
    x.Add(2,0,b);
    
    Matrix X(3), A(3, 1);
    
    for (int i = 0; i < 3; ++i){
        X.Add(i,0,pow(x.Get(0,0), i));
        X.Add(i,1,pow(x.Get(1,0), i));
        X.Add(i,2,pow(x.Get(2,0), i));
    }
    
    Matrix L(3), U(3), P1(3), P2(3);
    P1.insertDiag(1);
    P2.insertDiag(1);
    unsigned rank,swaps;
    P1P2LU(X, P1, P2, L, U, rank, swaps);
    
    if (rank==3){
        Matrix P1m(P1*m);
        SOLE(L, U, P1m, A, 3);
        
        A = P2*A;
        
        double S = 0;
        for (int i = 0; i < 3; ++i)
            S += A.Get(i,0) * Function(x.Get(i,0));
        return S;
    }
    else {
        throw("det == 0");
    }
}

double Gauss(double a, double b){
    return 1;
}