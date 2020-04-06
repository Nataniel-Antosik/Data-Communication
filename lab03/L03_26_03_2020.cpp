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

double a, b, c, suma = 0;

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

double* IDFT(zespolona* tab, int N) { //N to rozmiar
	double* tmp = new double[N]();
	double* wartosc_Funkcji = new double[N]();
	for (int i = 0; i < N; i++) {
		for (int k = 0; k < N; k++) {
			tmp[i] += (tab[k].RE * cos((k * i * 2 * M_PI) / N)) / N;
			tmp[i] -= (tab[k].IM * sin((k * i * 2 * M_PI) / N)) / N;
		}
		wartosc_Funkcji[i] = tmp[i];
	}
	return wartosc_Funkcji;
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
		//	if (tmp2[j] < 0) { tmp2[j] = 0; }
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

double x(double t) { return double(a * t * t + b * t + c); }
double y(double t) { return double(2 * pow(x(t), 2) + 12 * cos(t)); }
double z(double t) { return double(sin(2 * M_PI * 7 * t) * x(t) - 0.2 * log10(abs(y(t)) + M_PI)); }
double u(double t) { return double(sqrt(abs(y(t) * y(t) * z(t))) - 1.8 * sin(0.4 * t * z(t) * x(t))); }
double v(double t) {

	if (t < 0.22 && t >= 0) {
		return double((1 - 7 * t) * sin((2 * M_PI * t * 10) / (t + 0.04)));
	}
	if (t >= 0.22 && t < 0.7) {
		return double(0.63 * t * sin(125 * t));
	}
	if (t <= 1.0 && t >= 0.7) {
		return double(pow(t, (-0.662)) + 0.77 * sin(8 * t));
	}
}

double p(double t, int wybor) {
	//N należy do zbioru (2,4,16)
	switch (wybor) {
	case 1:
		for (int i = 1; i <= 2; i++) {
			suma = suma + double(cos(12 * t * pow(i, 2)) + cos(16 * t * i)) / (pow(i, 2));
		}
		break;
	case 2:
		for (int i = 1; i <= 4; i++) {
			suma = suma + double(cos(12 * t * pow(i, 2)) + cos(16 * t * i)) / (pow(i, 2));
		}
		break;
	case 3:
		for (int i = 1; i <= 16; i++) {
			suma = suma + double(cos(12 * t * pow(i, 2)) + cos(16 * t * i)) / (pow(i, 2));
		}
		break;
	}
	return suma;
}

int main()
{
	a = 1, b = 6, c = 2;
	zespolona* dane;

	double t_P = 0, t_K = 162, A = 1.0, f = 6, PHI = (2 * M_PI);
	double* x_T, * y_T, * z_T, * u_T, * v_T, * p_T_1, * p_T_2, * p_T_3;
	double FS = 81, *ton_P, * decybel, * skala_C, * idft; //częstotliwość próbkowania, największy dzielnik końca przedziału
	ton_P = new double[t_K](); //przechowuje wartości Y z tonu prostego
	int i = 0;
	//t_P - początek, t_K - koniec

	//Zadanie 1
	ofstream zad2_1("zad2_1.dat");
	while (t_P < t_K) {
		zad2_1 << 1 / FS * t_P << " " << ton_Prosty(1 / FS * t_P, A, PHI, f) << endl; //pomniejszamy wartość t_P próbki aby ładny wykres się wygenereował
		ton_P[i] = ton_Prosty(1 / FS * t_P, A, PHI, f);
		i += 1;
		t_P = t_P + 1;
	}
	//Zadanie 2
	dane = new zespolona[t_K]();
	skala_C = new double[t_K]();
	decybel = new double[t_K](); // wartość widmowa decybelowa

	dane = DFT(ton_P, t_K); //dane z dft
	skala_C = skala_Czestotliwosci(FS, t_K);
	decybel = widmo_A_Dec(dane, t_K);

	t_P = 0;
	i = 0;
	ofstream zad2_2("zad2_2.dat");
	while (t_P < t_K) {
		zad2_2 << skala_C[i] << " " << decybel[i] << endl;
		i = i + 1;
		t_P = t_P + 1;
	}
	//Zadanie 3
	int j_K = 0;
	x_T = new double[t_K]();
	y_T = new double[t_K]();
	z_T = new double[t_K]();
	u_T = new double[t_K]();
	v_T = new double[t_K]();
	p_T_1 = new double[t_K]();
	p_T_2 = new double[t_K]();
	p_T_3 = new double[t_K]();

	for (int j = 0; j < t_K; j++) { //pętla do zapisu danych z każdej funkcji 
		x_T[j] = x(j);
		y_T[j] = y(j);
		z_T[j] = z(j);
		u_T[j] = u(j);
		v_T[j] = x(j_K);
		p_T_1[j] = p(j, 1);
		p_T_2[j] = p(j, 2);
		p_T_3[j] = p(j, 3);
		j_K += 0.001;
	}

	i = 0;
	t_P = 0;
	ofstream x_t("x_t.dat");
	dane = DFT(x_T, t_K);
	skala_C = skala_Czestotliwosci(FS, t_K);
	decybel = widmo_A_Dec(dane, t_K);
	while (t_P < t_K) {
		x_t << skala_C[i] << " " << decybel[i] << endl;
		i += 1;
		t_P = t_P + 1;
	}

	i = 0;
	t_P = 0;
	ofstream y_t("y_t.dat");
	dane = DFT(y_T, t_K);
	skala_C = skala_Czestotliwosci(FS, t_K);
	decybel = widmo_A_Dec(dane, t_K);
	while (t_P < t_K) {
		y_t << skala_C[i] << " " << decybel[i] << endl;
		i += 1;
		t_P = t_P + 1;
	}

	i = 0;
	t_P = 0;
	ofstream z_t("z_t.dat");
	dane = DFT(z_T, t_K);
	skala_C = skala_Czestotliwosci(FS, t_K);
	decybel = widmo_A_Dec(dane, t_K);
	while (t_P < t_K) {
		z_t << skala_C[i] << " " << decybel[i] << endl;
		i += 1;
		t_P = t_P + 1;
	}

	i = 0;
	t_P = 0;
	ofstream u_t("u_t.dat");
	dane = DFT(u_T, t_K);
	skala_C = skala_Czestotliwosci(FS, t_K);
	decybel = widmo_A_Dec(dane, t_K);
	while (t_P < t_K) {
		u_t << skala_C[i] << " " << decybel[i] << endl;
		i += 1;
		t_P = t_P + 1;
	}

	i = 0;
	t_P = 0;
	ofstream v_t("v_t.dat");
	dane = DFT(v_T, t_K);
	skala_C = skala_Czestotliwosci(FS, t_K);
	decybel = widmo_A_Dec(dane, t_K);
	while (t_P < t_K) {
		v_t << skala_C[i] << " " << decybel[i] << endl;
		i += 1;
		t_P = t_P + 1;
	}

	i = 0;
	t_P = 0;
	ofstream p_t_1("p_t_1.dat");
	dane = DFT(p_T_1, t_K);
	skala_C = skala_Czestotliwosci(FS, t_K);
	decybel = widmo_A_Dec(dane, t_K);
	while (t_P < t_K) {
		p_t_1 << skala_C[i] << " " << decybel[i] << endl;
		i += 1;
		t_P = t_P + 1;
	}

	i = 0;
	t_P = 0;
	ofstream p_t_2("p_t_2.dat");
	dane = DFT(p_T_2, t_K);
	skala_C = skala_Czestotliwosci(FS, t_K);
	decybel = widmo_A_Dec(dane, t_K);
	while (t_P < t_K) {
		p_t_2 << skala_C[i] << " " << decybel[i] << endl;
		i += 1;
		t_P = t_P + 1;
	}

	i = 0;
	t_P = 0;
	ofstream p_t_3("p_t_3.dat");
	dane = DFT(p_T_3, t_K);
	skala_C = skala_Czestotliwosci(FS, t_K);
	decybel = widmo_A_Dec(dane, t_K);
	while (t_P < t_K) {
		p_t_3 << skala_C[i] << " " << decybel[i] << endl;
		i += 1;
		t_P = t_P + 1;
	}
	//zadanie 4

	idft = new double[t_K]();
	t_P = 0;
	i = 0;
	dane = DFT(ton_P, t_K);
	idft = IDFT(dane, t_K);
	ofstream IDFT_N("IDFT_N.dat");
	while (t_P < t_K) {
		IDFT_N << 1 / FS * t_P << " " << idft[i] << endl; //ustalamy częstotliwość probkowania, uzależniam rozmiar przedziału od wielkości próbkowania (FS) 
		i += 1;
		t_P = t_P + 1;
	}
	return 0;
}

