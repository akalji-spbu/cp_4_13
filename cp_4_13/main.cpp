//
//  main.cpp
//  cp_4_13
//
//  Created by Nikolay Tikhonov on 03.05.16.
//  Copyright © 2016 Nikolay Tikhonov. All rights reserved.
//

#include <iostream>
#include "integral.hpp"
void task3(double a, double b);
void task4(double a, double b);
void task5(double a, double b);
void task6(double a, double b);
void task7(double a, double b);

int main(int argc, const char * argv[]) {
    double a = 1.5;
    double b = 2.3;
    
    std::cout.setf(std::ios::fixed);
    std::cout.precision(13);
    
    std::cout << "///////////////////////////////////////////////////////////" << std::endl;
    std::cout<< "Task 3" << std::endl;
    std::cout<< "--------------------" << std::endl;
    task3(a,b);
    std::cout << std::endl;
    
    std::cout << "///////////////////////////////////////////////////////////" << std::endl;
    std::cout<< "Task 4" << std::endl;
    std::cout<< "--------------------" << std::endl;
    task4(a,b);
    std::cout << std::endl;
    
    std::cout << "///////////////////////////////////////////////////////////" << std::endl;
    std::cout<< "Task 5" << std::endl;
    std::cout<< "--------------------" << std::endl;
    task5(a,b);
    std::cout << std::endl;
    
    std::cout << "///////////////////////////////////////////////////////////" << std::endl;
    std::cout<< "Task 6" << std::endl;
    std::cout<< "--------------------" << std::endl;
    task6(a,b);
    std::cout << std::endl;
    
    std::cout << "///////////////////////////////////////////////////////////" << std::endl;
    std::cout<< "Task 7" << std::endl;
    std::cout<< "--------------------" << std::endl;
    task7(a,b);
    std::cout << std::endl;
    
    return 0;
}

void task3(double a, double b){
    //пункт 3
    double s1 = 0, s2 = 0, h = 0.01;
    int k = 1;
    b=b+0.00001;
    
    while (a + k*h <= b) {
        s1 += Newton_Cotes(a+(k-1)*h, a + k*h);
        s2 += Gauss(a+(k-1)*h, a + k*h);
        k++;
    }
    std::cout << "h:	" << h << std::endl;
    std::cout <<"NK:	"<< s1 << std::endl;
    std::cout <<"G:	"<<  s2 << std::endl;
}

void task4(double a, double b){
    //пункт 4
    double s1 = 0, s2 = 0, h1=0.8;
    
    s1 = Gauss(a, b);
    
    
    double h2 =(b-a)/ ceil(2 * (b - a) / h1);
    
    int k = 1;
    while (a + k*h2 <= b+1e-11) {
        s2 += Gauss(a + (k - 1)*h2, a + k*h2);
        k++;
    }
    
    
    double L = h1 / h2;
    double R = fabs((s2 - s1) / (pow(L, 6) - 1));
    
    //std::cout<<s2<<std::endl<<s1<<std::endl;
    
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

void task5(double a, double b){
    //пункт 5
    double s1 = 0, s2 = 0, s3 = 0, h1 = 0.8;
    
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
    
    while (R > 1e-12){
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
        //std::cout << s1<<"\t" << s2<<"\t" <<s3 << std::endl;
        //std::cout << "R:	"<< R << std::endl;
        std::cout << "m: " << m << std::endl;
    }
    
    std::cout << "S: " << s3 << std::endl;
}

void task6(double a, double b){
    double s1 = 0, s2 = 0, h1 = 0.8/2;
    
    
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

void task7(double a, double b){
    double s1 = 0, s2 = 0, h1 = 0.8/8;
    
    
    int k = 1;
    while (a + k*h1 <= b + 1e-5) {
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
    //std::cout << s1<<"\t" << s2<<"\t" <<s3 << std::endl;
    //std::cout << "R000:	\t" << R << std::endl;
    
    //std::cout<<L<<std::endl;
    //std::cout<<s1<<std::endl;
    //std::cout<<s2<<std::endl;
    //std::cout<<s3<<std::endl;
    
    double hopt = h3 * pow(1e-10 / R, 1. / 6);
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
    std::cout << "n: \t\t" << n << std::endl;
    std::cout << "m:	" << m << std::endl;
    std::cout << "R:	" << R << std::endl;
    std::cout << "E:	" << 1e-10 << std::endl;
    std::cout << "S:	" << s4 << std::endl;
}
