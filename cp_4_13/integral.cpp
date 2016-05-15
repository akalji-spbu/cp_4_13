//
//  integral.cpp
//  cp_4_13
//
//  Created by Nikolay Tikhonov on 09.05.16.
//  Copyright © 2016 Nikolay Tikhonov. All rights reserved.
//

#include "integral.hpp"

double Func(double x){
    return 2.5*cos(2 * x)*exp(2 * x / 3.) + 4 * sin(3.5*x)*exp(-3 * x) + 3 * x;
}

double g(double x, double A, double B, double C, double D){
    return A*pow(x, 3) + B*pow(x, 2) + C*x + D;
}

double dg(double x, double A, double B, double C){
    return 3 * A*pow(x, 2) + 2 * B*x + C;
}

double newton_r(double x, double a, double b, double A, double B, double C, double D){
    double sol = x;
    while (fabs(g(sol, A, B, C, D)) > 1e-13){
        sol = sol - g(sol, A, B, C, D) / dg(sol, A, B, C);
    }
    
    return sol;
}

void newton(double* m, double a, double b, double A, double B, double C, double D){
    double DIS = pow(2 * B, 2) - 4 * 3 * A*C;
    
    if (DIS <= 0)
        throw ("error1");
    
    //Точки максимума/минимума
    double x1 = (-2 * B + sqrt(DIS)) / (2 * 3 * A);
    double x2 = (-2 * B - sqrt(DIS)) / (2 * 3 * A);
    
    double c = -2 * B / 6 * A;
    
    
    if (x1>b || x1<a || x2>b || x2 < a)
        throw ("error2");
    
    if (x1 > x2){
        double tmp = x1;
        x1 = x2;
        x2 = tmp;
    }
    
    if (!((A > 0 && g(x1, A, B, C, D) >= 0 && g(x2, A, B, C, D) <= 0) || (A < 0 && g(x1, A, B, C, D) <= 0 && g(x2, A, B, C, D) >= 0)))
        throw ("error3");
    
    
    m[0] = newton_r(a, a, x1, A, B, C, D);
    m[1] = (g(c, A, B, C, D) > 0) ? newton_r(c, c, x2, A, B, C, D) : newton_r(c, x1, c, A, B, C, D);
    m[2] = newton_r(b, x2, b, A, B, C, D);
    
    
    if (m[0]<a || m[0]>b || m[1]<a || m[1]>b || m[2]<a || m[2]>b)
        throw ("error4");
    
}

double Newton_K(double a, double b){
    Matrix m(3, 1);
    Matrix x(3, 1);
    
    x.Add(0,0,a);
    x.Add(1,0, (a + b) / 2);
    x.Add(2,0, b);
    
    m.Add(0,0, 5. / 4 * pow(b - 0.1, 4. / 5) - 5. / 4 * pow(a - 0.1, 4. / 5));
    m.Add(1,0, 0.1*5. / 4 * pow(b - 0.1, 4. / 5) + 5. / 9 * pow(b - 0.1, 9. / 5) - (0.1*5. / 4 * pow(a - 0.1, 4. / 5) + 5. / 9 * pow(a - 0.1, 9. / 5)));
    m.Add(2,0, 5. / 14 * pow(b - 0.1, 14. / 5) + 2 * 0.1*5. / 9 * pow(b - 0.1, 9. / 5) + 0.1*0.1*5. / 4 * pow(b - 0.1, 4. / 5) - (5. / 14 * pow(a - 0.1, 14. / 5) + 2 * 0.1*5. / 9 * pow(a - 0.1, 9. / 5) + 0.1*0.1*5. / 4 * pow(a - 0.1, 4. / 5)));
    
    Matrix X(3), A(3, 1);
    
    for (int i = 0; i < 3; ++i){
        X.Add(i,0, pow(x.Get(0,0), i));
        X.Add(i,1, pow(x.Get(1,0), i));
        X.Add(i,2, pow(x.Get(2,0), i));
    }
    
    Matrix L(3), U(3), P(3, 3), G(3, 3);
    P.insertDiag(1);
    G.insertDiag(1);
    unsigned rank,swaps;
    P1P2LU(X, P, G, L, U, rank,swaps);
    if (rank==3){
        Matrix Pm(P*m);
        SOLE(L, U, Pm, A, 3);
        //SLAU_solution(L, U, P*m, 3, A);
        A = G*A;
        
        double S = 0;
        for (int i = 0; i < 3; ++i)
            S += A.Get(i,0) * Func(x.Get(i,0));
        return S;
    }
    else {
        throw("det == 0");
    }
}

double Gauss(double a, double b){
    Matrix m(6, 1);
    
    m.Add(0,0, 5. / 4 * pow(b - 0.1, 4. / 5) - 5. / 4 * pow(a - 0.1, 4. / 5));
    m.Add(1,0, 0.1*5. / 4 * pow(b - 0.1, 4. / 5) + 5. / 9 * pow(b - 0.1, 9. / 5) - (0.1*5. / 4 * pow(a - 0.1, 4. / 5) + 5. / 9 * pow(a - 0.1, 9. / 5)));
    m.Add(2,0 , 5. / 14 * pow(b - 0.1, 14. / 5) + 2 * 0.1*5. / 9 * pow(b - 0.1, 9. / 5) + 0.1*0.1*5. / 4 * pow(b - 0.1, 4. / 5) - (5. / 14 * pow(a - 0.1, 14. / 5) + 2 * 0.1*5. / 9 * pow(a - 0.1, 9. / 5) + 0.1*0.1*5. / 4 * pow(a - 0.1, 4. / 5)));
    m.Add(3,0, 5. / 19 * pow(b - 0.1, 19. / 5) + 3 * 0.1*5. / 14 * pow(b - 0.1, 14. / 5) + 3 * 0.1*0.1*5. / 9 * pow(b - 0.1, 9. / 5) + 0.1*0.1*0.1*5. / 4 * pow(b - 0.1, 4. / 5) - (5. / 19 * pow(a - 0.1, 19. / 5) + 3 * 0.1*5. / 14 * pow(a - 0.1, 14. / 5) + 3 * 0.1*0.1*5. / 9 * pow(a - 0.1, 9. / 5) + 0.1*0.1*0.1*5. / 4 * pow(a - 0.1, 4. / 5)));
    m.Add(4,0, 5. / 24 * pow(b - 0.1, 24. / 5) + 5. / 19 * 4 * 0.1* pow(b - 0.1, 19. / 5) + 6 * 0.1*0.1*5. / 14 * pow(b - 0.1, 14. / 5) + 4 * 0.1* 0.1*0.1*5. / 9 * pow(b - 0.1, 9. / 5) + 0.1*0.1*0.1*0.1*5. / 4 * pow(b - 0.1, 4. / 5) - (5. / 24 * pow(a - 0.1, 24. / 5) + 5. / 19 * 4 * 0.1* pow(a - 0.1, 19. / 5) + 6 * 0.1*0.1*5. / 14 * pow(a - 0.1, 14. / 5) + 4 * 0.1* 0.1*0.1*5. / 9 * pow(a - 0.1, 9. / 5) + 0.1*0.1*0.1*0.1*5. / 4 * pow(a - 0.1, 4. / 5)));
    m.Add(5,0, 5. / 29 * pow(b - 0.1, 29. / 5) + 5. / 24 * 5 * 0.1 * pow(b - 0.1, 24. / 5) + 5. / 19 * 10 * 0.1*0.1* pow(b - 0.1, 19. / 5) + 10 * 0.1*0.1*0.1*5. / 14 * pow(b - 0.1, 14. / 5) + 5 * 0.1 * 0.1* 0.1*0.1*5. / 9 * pow(b - 0.1, 9. / 5) + 0.1*0.1*0.1*0.1*0.1*5. / 4 * pow(b - 0.1, 4. / 5) - (5. / 29 * pow(a - 0.1, 29. / 5) + 5. / 24 * 5 * 0.1 * pow(a - 0.1, 24. / 5) + 5. / 19 * 10 * 0.1*0.1* pow(a - 0.1, 19. / 5) + 10 * 0.1*0.1*0.1*5. / 14 * pow(a - 0.1, 14. / 5) + 5 * 0.1 * 0.1* 0.1*0.1*5. / 9 * pow(a - 0.1, 9. / 5) + 0.1*0.1*0.1*0.1*0.1*5. / 4 * pow(a - 0.1, 4. / 5)));
    
    Matrix A(3), b1(3,1);
    
    for (int s = 0; s < 3; ++s){
        for (int i = 0; i < 3; ++i)
            A.Add(s,i,m.Get(i + s,0));
        b1.Add(s,0,-m.Get(3 + s,0));
    }
    
    Matrix a1(3, 1);
    Matrix L(3), U(3), P(3, 3), G(3, 3);
    P.insertDiag(1);
    G.insertDiag(1);
    unsigned rank, swaps;
    P1P2LU(A, P, G, L, U, rank, swaps);
    
    //SLAU_solution(L, U, P*b1, 3, a1);
    Matrix Pb1(P*b1);
    SOLE(L, U, Pb1, a1, 3);
    a1 = G*a1;
    
    double *x1 = new double[3];
    newton(x1, a, b, 1, a1.Get(2,0), a1.Get(1,0), a1.Get(0,0));
    
    /*double Q, R;
     double const pi = 3.14;
     double tmp = a1[2][0];
     Q = (pow(a1[2][0], 2.0) - 3.0 * a1[1][0]) / 9.0;
     R = (2.0 * pow(a1[2][0], 3.0) - 9.0 * a1[2][0] * a1[1][0] + 27.0 * a1[0][0]) / 54.0;
     //cout << "Q " << Q << endl << "R " << R << endl;
     double t = 0;
     if (pow(R, 2.0)<pow(Q, 3.0))						//Ф-ы Виета
     {
     t = acos(R / sqrt(pow(Q, 3.0))) / 3.0,
     x1[0] = -2.0 * sqrt(Q)*cos(t) - tmp / 3.0;
     x1[1] = -2.0 * sqrt(Q)*cos(t + (2.0 * pi / 3.0)) - tmp / 3.0;
     x1[2] = -2.0 * sqrt(Q)*cos(t - (2.0 * pi / 3.0)) - tmp / 3.0;
     //	cout << x1 << " " << x2 << " " << x3 << endl;
     }
     else
     cout << "net deistv. kornei v res3" << endl;*/
    //Проверка корней
    //cout << '[' << a << " ; " << b << ']' << endl;
    //cout << x1[0] << " " << x1[1] << " " << x1[2] << endl << endl;
    
    Matrix x(3, 1);
    
    x.Add(0,0,x1[0]);
    x.Add(1,0,x1[1]);
    x.Add(2,0,x1[2]);
    
    Matrix X(3), A1(3, 1);
    
    for (int i = 0; i < 3; ++i){
        X.Add(i,0, pow(x.Get(0,0), i));
        X.Add(i,1, pow(x.Get(1,0), i));
        X.Add(i,2, pow(x.Get(2,0), i));
    }
    
    Matrix m1(3, 1);
    
    m1.Add(0,0, m.Get(0,0));
    m1.Add(1,0, m.Get(1,0));
    m1.Add(2,0, m.Get(2,0));
    Matrix L1(3), U1(3), P1(3), G1(3);
    P.insertDiag(1);
    G.insertDiag(1);
    unsigned rank1, swaps1;
    P1P2LU(X, P1, G1, L1, U1, rank1, swaps1);
    Matrix P1m1 = P1*m1;
    SOLE(L1, U1, P1m1, A1, 3);
    A1 = G1*A1;
    
    
    double S = 0;
    for (int i = 0; i < 3; ++i)
        S += A1.Get(i,0) * Func(x.Get(i,0));
    //S += A1[i][0] * pow(x[i][0],5);
    //cout << S <<" "<<m[5][0]<< endl;
    delete[] x1;
    return S;
}
