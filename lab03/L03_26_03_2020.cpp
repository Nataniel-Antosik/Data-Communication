#define _USE_MATH_DEFINES
#include<iostream>
#include<math.h>
#include "math.h"
#include <fstream>

//44261 A = 1, B = 6, C=2

using namespace std;

struct zespolona {
	double RE, IM;
};

zespolona* DFT(double* tab, int N) { //N to rozmiar
	zespolona* tmp = new zespolona[N]();
	for (int i = 0; i < N; i++) {
		for (int k = 0; k < N; k++) {
			tmp[i].RE += tab[k] * cos((k * i * 2 * M_PI) / N);
			tmp[i].IM -= tab[k] * sin((k * i * 2 * M_PI) / N);
		}
	}
	return tmp;
}

double* widmo_Amplitudowe(zespolona* tab, int rozmiar) {
	double* tmp = new double[rozmiar]();
	for (int j = 0; j < rozmiar; j++) {
		tmp[j] = sqrt((tab[j].RE * tab[j].RE) + (tab[j].IM * tab[j].IM));
	}
	return tmp;
}

double* widmo_A_Dec(zespolona* tab, int rozmiar) {
	double* tmp = new double[rozmiar]();
	double* tmp2 = new double[rozmiar]();
	tmp = widmo_Amplitudowe(tab, rozmiar);
	for (int j = 0; j < rozmiar; j++) {
		tmp2[j] = (10 * log10(tmp[j]));
	}
	return tmp2;
}

double* skala_Czestotliwosci(double FS, int rozmiar) {
	double* tmp = new double[rozmiar]();
	for (int j = 0; j < rozmiar; j++) {
		tmp[j] = j * FS / rozmiar;
	}
	return tmp;
}

double ton_Prosty(double t, double A, double PHI, double f) { return double(A * sin(2 * M_PI * f * t + PHI)); }

int main()
{
	zespolona* dane;

	double t_P = 0, t_K = 162, A = 1.0, f = 6, PHI = (2 * M_PI);
	double FS = 81, *tmp, * decybel, * skala_C; //częstotliwość próbkowania, największy dzielnik końca przedziału
	tmp = new double[t_K](); //przechowuje wartości Y z tonu prostego
	int i = 0;
	//t_P - początek, t_K - koniec
	//Zadanie 1
	ofstream zad2_1("zad2_1.dat");
	while (t_P <= t_K) {
		zad2_1 << t_P << " " << ton_Prosty(1 / FS * t_P, A, PHI, f) << endl; //pomniejszamy wartość t_P próbki aby ładny wykres się wygenereował
		tmp[i] = ton_Prosty(1 / FS * t_P, A, PHI, f);
		i += 1;
		t_P = t_P + 1;
	}
	dane = new zespolona[t_K]();
	decybel = new double[t_K](); // wartość widmowa decybelowa
	skala_C = new double[t_K]();
	decybel = new double[t_K](); // wartość widmowa decybelowa

	dane = DFT(tmp, t_K); //dane z dft
	skala_C = skala_Czestotliwosci(FS, t_K);
	decybel = widmo_A_Dec(dane, t_K);
	i = 0;
	ofstream zad2_2("zad2_2.dat");
	while (t_P <= t_K) {
		zad2_2 << skala_C[i] << " " << decybel[i] << endl;
		i += 1;
		t_P = t_P + 1;
	}
	return 0;
}

