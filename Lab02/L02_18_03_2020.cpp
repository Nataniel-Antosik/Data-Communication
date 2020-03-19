#define _USE_MATH_DEFINES
#include <iostream>
#include <math.h>
#include <fstream>

//44261 A = 1, B = 6, C=2s

using namespace std;

double ton_Prosty(double t, double A, double PHI, double f) { return double(A * sin(2 * M_PI * f * t + PHI)); }

int kwantyzacja_Sygnalu(double t, double A, double PHI, double f, int Q) { return (((pow(2, Q) / 2) - 0.5) * A * sin(2 * M_PI * f * t + PHI) + ceil((pow(2, Q) / 2) - 0.5)); }
//ceil - zaokrąglenie w górę
int main()
{
	double t_P = 0, t_K = 1, A = 1, f = 6, PHI = 2 * M_PI * f;
	//t_P - początek, t_K - koniec
	double FS = 0.001; //częstotliwość próbkowania
	//Zadanie 1
	ofstream zad1("zad1.dat");
	while (t_P <= t_K) {
		zad1 << t_P << " " << ton_Prosty(t_P, A, PHI, f) << endl;
		t_P = t_P + FS;
	}
	zad1.close();
	//Zadanie 2
	int Q = 16;
	t_P = 0;
	ofstream zad2("zad2.dat");
	while (t_P <= t_K) {
		zad2 << t_P << " " << kwantyzacja_Sygnalu(t_P, A, PHI, f, Q) << endl;
		t_P = t_P + FS;
	}
	zad2.close();
	//Zadanie 3
	Q = Q / 2;
	t_P = 0;
	FS = FS * 2;
	ofstream zad3("zad3.dat");
	while (t_P <= t_K) {
		zad3 << t_P << " " << kwantyzacja_Sygnalu(t_P, A, PHI, f, Q) << endl;
		t_P = t_P + FS;
	}
	zad3.close();
	return 0;
}
