#define _USE_MATH_DEFINES
#include<iostream>
#include<math.h>
#include "math.h"
#include <fstream>
#include <vector>

using namespace std;

//44261 A = 1, B = 6, C=2

double m(double t) {
    double suma = 0;
    for (int i = 1; i <= 2; i++) {
        suma = suma + double(cos(12 * t * pow(i, 2)) + cos(16 * t * i)) / (pow(i, 2));
    }
    return suma;
}

double zA(double t, double kA, double fn) { return (kA * m(t) + 1) * cos(2 * M_PI * fn * t); }

double zp(double t, double kp, double fn) { return cos(2 * M_PI * fn * t + kp * m(t)); }

int main()
{       
    //zadanie 1
    double kA, kp, fn, t_P, t_K;

    //============TEST============//
    ofstream zad1_1("zad1_1.dat");
    while (t_P <= t_K) {
        zad1_1 << t_P << " " << zA(1,2,3) << endl;
        t_P = t_P + 1;
    }
    //============TEST============//


    return 0;
}
