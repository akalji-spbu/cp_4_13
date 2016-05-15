//
//  main.cpp
//  cp_4_13
//
//  Created by Nikolay Tikhonov on 03.05.16.
//  Copyright © 2016 Nikolay Tikhonov. All rights reserved.
//

#include <iostream>
#include "integral.hpp"
void a();
void b();
void c();
void d();
void e();

int main(int argc, const char * argv[]) {
    std::cout.setf(std::ios::fixed);
    std::cout.precision(13);
    std::cout << "--------------------" << std::endl << "N 3" << std::endl << "--------------------" << std::endl;
    a();
    std::cout << "--------------------" << std::endl << "N 4" << std::endl << "--------------------" << std::endl;
    b();
    std::cout << "--------------------" << std::endl << "N 5" << std::endl << "--------------------" << std::endl;
    c();
    std::cout << "--------------------" << std::endl << "N 6" << std::endl << "--------------------"<<std::endl;
    d();
    std::cout << "--------------------" << std::endl << "N 7" << std::endl << "--------------------" << std::endl;
    e();
    return 0;
}

void a(){
    //пункт 3
    double s1 = 0, s2 = 0, h = 0.0275, a = 0.1, b = 2.300001;
    int k = 1;
    
    while (a + k*h <= b) {
        s1 += Newton_K(a+(k-1)*h, a + k*h);
        s2 += Gauss(a+(k-1)*h, a + k*h);
        k++;
    }
    std::cout << "h:	" << h << std::endl;
    std::cout <<"NK:	"<< s1 << std::endl;
    std::cout <<"G:	"<<  s2 << std::endl;
}

void b(){
    //пункт 4
    double s1 = 0, s2 = 0, h1=2.2, a=0.1, b=2.3;
    
    s1 = Gauss(a, b);
    
    
    double h2 =(b-a)/ ceil(2 * (b - a) / h1);
    
    int k = 1;
    while (a + k*h2 <= b+1e-10) {
        s2 += Gauss(a + (k - 1)*h2, a + k*h2);
        k++;
    }
    
    
    double L = h1 / h2;
    double R = fabs((s2 - s1) / (pow(L, 6) - 1));
    
    while (R > 1e-11){
        s1 = s2;
        h1 = h2;
        h2 = (b - a) / ceil(2 * (b - a) / h1);
        
        s2 = 0;
        k = 1;
        while (a + k*h2 <= b + 1e-10) {
            s2 += Gauss(a + (k - 1)*h2, a + k*h2);
            k++;
            
        }
        
        
        L = h1 / h2;
        R = fabs((s2 - s1) / (pow(L, 6) - 1));
        double R2 = fabs((s2 - s1) / (-pow(L, -6) + 1));
        //Для проверки 6 пункта
        std::cout << " R2: " << R2 << std::endl << "h: " << h2 << " " << "R1: " << R;
    }
    
    std::cout <<std::endl<<"S: "<< s2 << std::endl;
    
}

void c(){
    //пункт 5
    double s1 = 0, s2 = 0, s3 = 0, h1 = 2.2, a = 0.1, b = 2.3;
    
    s1 = Gauss(a, b);
    
    
    double h2 = (b - a) / ceil(2 * (b - a) / h1);
    int k = 1;
    while (a + k*h2 <= b + 1e-10) {
        s2 += Gauss(a + (k - 1)*h2, a + k*h2);
        k++;
    }
    
    
    double h3 = (b - a) / ceil(2 * (b - a) / h2);
    k = 1;
    while (a + k*h3 <= b + 1e-10) {
        s3 += Gauss(a + (k - 1)*h3, a + k*h3);
        k++;
    }
    
    
    double L = h2 / h3;
    double m = -log(fabs((s3 - s2) / (s2 - s1))) / log(L);
    double R = fabs((s3 - s2) / (pow(L, m) - 1));
    std::cout <<"m: "<< m << std::endl;
    
    while (R > 1e-11){
        s1 = s2;
        s2 = s3;
        h2 = h3;
        h3 = (b - a) / ceil(2 * (b - a) / h2);
        
        s3 = 0;
        k = 1;
        while (a + k*h3 <= b + 1e-10) {
            s3 += Gauss(a + (k - 1)*h3, a + k*h3);
            k++;
        }
        
        
        L = h2 / h3;
        m = -log(fabs((s3 - s2) / (s2 - s1))) / log(L);
        R = fabs((s3 - s2) / (pow(L, m) - 1));
        std::cout << "m: " << m << std::endl;
    }
    
    std::cout << "S: " << s3 << std::endl;
}

void d(){
    double s1 = 0, s2 = 0, h1 = 1.1, a = 0.1, b = 2.3;
    
    
    int k = 1;
    while (a + k*h1 <= b + 1e-10) {
        s1 += Gauss(a + (k - 1)*h1, a + k*h1);
        k++;
    }
    
    
    double h2 = (b - a) / ceil(2 * (b - a) / h1);
    k = 1;
    while (a + k*h2 <= b + 1e-10) {
        s2 += Gauss(a + (k - 1)*h2, a + k*h2);
        k++;
    }
    
    double L = h1 / h2;
    double R = fabs((s2 - s1) / (pow(L, 6) - 1));
    
    double hopt = h2 * pow(1e-10 / R, 1. / 6);
    hopt = (b - a) / (ceil((b - a) / (0.95*hopt)));
    double s3 = 0;
    std::cout << "hopt:	" << hopt << std::endl;
    k = 1;
    while (a + k*hopt <= b + 1e-9) {
        s3 += Gauss(a + (k - 1)*hopt, a + k*hopt);
        k++;
    }
    double h4 = (b - a) / (ceil((b - a) / (2*hopt))), s4=0;
    k = 1;
    while (a + k*h4 <= b + 1e-10) {
        s4 += Gauss(a + (k - 1)*h4, a + k*h4);
        k++;
    }
    
    L = h4 / hopt;
    R = fabs((s3 - s4) / (pow(L, 6) - 1));
    double n = (b - a) / hopt;
    std::cout << "n:	" << n << std::endl;
    std::cout <<"R:	"<< R << std::endl;
    std::cout << "E:	" << 1e-10 << std::endl;
    std::cout <<"S:	"<< s3 << std::endl;
}

void e(){
    double s1 = 0, s2 = 0, h1 = 0.55, a = 0.1, b = 2.3;
    
    
    int k = 1;
    while (a + k*h1 <= b + 1e-10) {
        s1 += Gauss(a + (k - 1)*h1, a + k*h1);
        k++;
    }
    
    
    double h2 = (b - a) / ceil(2 * (b - a) / h1);
    k = 1;
    while (a + k*h2 <= b + 1e-10) {
        s2 += Gauss(a + (k - 1)*h2, a + k*h2);
        k++;
    }
    
    double h3 = (b - a) / ceil(2 * (b - a) / h2), s3 =0;
    k = 1;
    while (a + k*h3 <= b + 1e-10) {
        s3 += Gauss(a + (k - 1)*h3, a + k*h3);
        k++;
    }
    
    double L = h2 / h3;
    double m = -log(fabs((s3 - s2) / (s2 - s1))) / log(L);
    double R = fabs((s3 - s2) / (pow(L, m) - 1));
    
    double hopt = h2 * pow(1e-10 / R, 1. / 6);
    hopt = (b - a) / (ceil((b - a) / (0.95*hopt)));
    double s4 = 0;
    std::cout << "hopt:	" << hopt << std::endl;
    k = 1;
    while (a + k*hopt <= b + 1e-9) {
        s4 += Gauss(a + (k - 1)*hopt, a + k*hopt);
        k++;
    }
    double h4 = (b - a) / (ceil((b - a) / (2 * hopt)));
    double s5 = 0;
    k = 1;
    while (a + k*h4 <= b + 1e-10) {
        s5 += Gauss(a + (k - 1)*h4, a + k*h4);
        k++;
    }
    
    L = h4 / hopt;
    R = fabs((s4 - s5) / (pow(L, m) - 1));
    
    double n = (b - a) / hopt;
    std::cout << "n:	" << n << std::endl;
    std::cout << "m:	" << m << std::endl;
    std::cout << "R:	" << R << std::endl;
    std::cout << "E:	" << 1e-10 << std::endl;
    std::cout << "S:	" << s4 << std::endl;
}
