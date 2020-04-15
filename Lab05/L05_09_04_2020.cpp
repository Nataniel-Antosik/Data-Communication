#define _USE_MATH_DEFINES
#include<iostream>
#include<math.h>
#include "math.h"
#include <fstream>

using namespace std;

//44261 A = 1, B = 6, C=2

struct zespolona {
    double RE, IM;
};

//Sygnał informacyjny
double m(double t) {
    double suma = 0;
    for (int i = 1; i <= 2; i++) {
        suma = suma + double(cos(12 * t * pow(i, 2)) + cos(16 * t * i)) / (pow(i, 2));
    }
    return suma;
}
//Modulacja amplitudy
double zA(double t, double kA, double fn) { return (kA * m(t) + 1) * cos(2 * M_PI * fn * t); }
//Modulacja fazy
double zp(double t, double kp, double fn) { return cos(2 * M_PI * fn * t + kp * m(t)); }

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

double szerokosc_Pasma(double * wartosc_F, double * czestotliwosc, int rozmiar) {
    int max_X = 0;
    double n_Rozmiar = ceil(rozmiar / 2); //połowa tablicy
    for (int i = 0; i < n_Rozmiar; i++) { 
        if (wartosc_F[i] > wartosc_F[max_X]) {
            max_X = i;
        }
    }
    double pomniejszona_w_Funkcji = (wartosc_F[max_X] - 3); //wartosc pomniejszona o 3 db
    double roznica = pomniejszona_w_Funkcji - wartosc_F[0]; //różnica od lewej wartości
    double fmin = czestotliwosc[0];
    double fmax = czestotliwosc[max_X + 1];
    for (int i = 1; i < max_X; i++) {
        if (roznica > pomniejszona_w_Funkcji - wartosc_F[i]) { //jeżeli znajdę mniejszą różnicę to zamieniam różnicę i nową wartość minimalną
            roznica = pomniejszona_w_Funkcji - wartosc_F[i];
            fmin = czestotliwosc[i];
        }
    }
    roznica = pomniejszona_w_Funkcji - wartosc_F[max_X + 1]; //to samo co na górze tylko bierzemy wartość z prawej strony od naszego znalezionego maksa
    for (int i = max_X + 2; i < n_Rozmiar; i++) {
        if (roznica > pomniejszona_w_Funkcji - wartosc_F[i]) {
            roznica = pomniejszona_w_Funkcji - wartosc_F[i];
            fmax = czestotliwosc[i];
        }
    }
    return fmax - fmin;
}

int main()
{       
    //zadanie 1
    //a)
    double kA, kp, fn, t_P, t_K, fs, * skala_C;
    double* decybel1, * decybel2, * decybel3, * decybel4, * decybel5, * decybel6;
    double* zA_a, * zp_a, * zA_b, * zp_b, * zA_c, * zp_c;
    int j = 0;
    t_P = 0;
    t_K = 1;
    fn = 1000;
    fs = 5000;
    kA = 0.5;
    kp = 1;

    zA_a = new double[t_K * fs]();
    zp_a = new double[t_K * fs]();
    zA_b = new double[t_K * fs]();
    zp_b = new double[t_K * fs]();
    zA_c = new double[t_K * fs]();
    zp_c = new double[t_K * fs]();
    
    ofstream zad1_a_i("zad1_a_i.dat");
    while (t_P < t_K) {
        zad1_a_i << t_P << " " << m(t_P) << endl;
        t_P = t_P + 1/fs;
    }
    t_P = 0;
    ofstream zad1_a_zA("zad1_a_zA.dat");
    while (t_P < t_K) {
        zad1_a_zA << t_P << " " << zA(t_P, kA, fn) << endl;
        zA_a[j] = zA(t_P, kA, fn);
        j += 1;
        t_P = t_P + 1 / fs;
    }
    t_P = 0;
    j = 0;
    ofstream zad1_a_zp("zad1_a_zp.dat");
    while (t_P < t_K) {
        zad1_a_zp << t_P << " " << zp(t_P, kp, fn) << endl;
        zp_a[j] = zp(t_P, kp, fn);
        t_P = t_P + 1 / fs;
        j += 1;
    }

    //b)
    kA = 7;
    kp = 1.57;
    t_P = 0;
    j = 0;
    ofstream zad1_b_zA("zad1_b_zA.dat");
    while (t_P < t_K) {
        zad1_b_zA << t_P << " " << zA(t_P, kA, fn) << endl;
        zA_b[j] = zA(t_P, kA, fn);
        t_P = t_P + 1 / fs;
        j += 1;
    }
    t_P = 0;
    j = 0;
    ofstream zad1_b_zp("zad1_b_zp.dat");
    while (t_P < t_K) {
        zad1_b_zp << t_P << " " << zp(t_P, kp, fn) << endl;
        zp_b[j] = zp(t_P, kp, fn);
        t_P = t_P + 1 / fs;
        j += 1;
    }

    //c) BA = 61  , AB = 16 
    kA = 66;
    kp = 25;
    t_P = 0;
    j = 0;
    ofstream zad1_c_zA("zad1_c_zA.dat");
    while (t_P < t_K) {
        zad1_c_zA << t_P << " " << zA(t_P, kA, fn) << endl;
        zA_c[j] = zA(t_P, kA, fn);
        t_P = t_P + 1 / fs;
        j += 1;
    }
    t_P = 0;
    j = 0;
    ofstream zad1_c_zp("zad1_c_zp.dat");
    while (t_P < t_K) {
        zad1_c_zp << t_P << " " << zp(t_P, kp, fn) << endl;
        zp_c[j] = zp(t_P, kp, fn);
        t_P = t_P + 1 / fs;
        j += 1;
    }
    //zadanie 2
    zespolona* ZA_a, * Zp_a, * ZA_b, * Zp_b, * ZA_c, * Zp_c;
    ZA_a = new zespolona[t_K * fs]();
    Zp_a = new zespolona[t_K * fs]();
    ZA_b = new zespolona[t_K * fs]();
    Zp_b = new zespolona[t_K * fs]();
    ZA_c = new zespolona[t_K * fs]();
    Zp_c = new zespolona[t_K * fs]();

    skala_C = new double[t_K * fs]();

    decybel1 = new double[t_K * fs]();
    decybel2 = new double[t_K * fs]();
    decybel3 = new double[t_K * fs]();
    decybel4 = new double[t_K * fs]();
    decybel5 = new double[t_K * fs]();
    decybel6 = new double[t_K * fs]();

    skala_C = skala_Czestotliwosci(fs, t_K * fs);
    
    ZA_a = DFT(zA_a, t_K * fs);
    decybel1 = widmo_A_Dec(ZA_a, t_K * fs);
    ofstream zad2_a_zA("zad2_a_zA.dat");
    for (int i = 0; i < t_K * fs; i++) {
        zad2_a_zA << skala_C[i] << " " << decybel1[i] << endl;
    }
    cout << "Szerokosc pasma sygnalu zA(t) a) wynosi " << szerokosc_Pasma(decybel1, skala_C, t_K * fs) << endl;
    
    Zp_a = DFT(zp_a, t_K * fs);
    decybel2 = widmo_A_Dec(Zp_a, t_K * fs);
    ofstream zad2_a_zp("zad2_a_zp.dat");
    for (int i = 0; i < t_K * fs; i++) {
        zad2_a_zp << skala_C[i] << " " << decybel2[i] << endl;
    }
    cout << "Szerokosc pasma sygnalu zp(t) a) wynosi " << szerokosc_Pasma(decybel2, skala_C, t_K * fs) << endl;

    ZA_b = DFT(zA_b, t_K * fs);
    decybel3 = widmo_A_Dec(ZA_b, t_K * fs);
    ofstream zad2_b_zA("zad2_b_zA.dat");
    for (int i = 0; i < t_K * fs; i++) {
        zad2_b_zA << skala_C[i] << " " << decybel3[i] << endl;
    }
    cout << "Szerokosc pasma sygnalu zA(t) b) wynosi " << szerokosc_Pasma(decybel3, skala_C, t_K * fs) << endl;
 
    Zp_b = DFT(zp_b, t_K * fs);
    decybel4 = widmo_A_Dec(Zp_b, t_K * fs);
    ofstream zad2_b_zp("zad2_b_zp.dat");
    for (int i = 0; i < t_K * fs; i++) {
        zad2_b_zp << skala_C[i] << " " << decybel4[i] << endl;
    }
    cout << "Szerokosc pasma sygnalu zp(t) b) wynosi " << szerokosc_Pasma(decybel4, skala_C, t_K * fs) << endl;

    ZA_c = DFT(zA_c, t_K * fs);
    decybel5 = widmo_A_Dec(ZA_c, t_K * fs);
    ofstream zad2_c_zA("zad2_c_zA.dat");
    for (int i = 0; i < t_K * fs; i++) {
        zad2_c_zA << skala_C[i] << " " << decybel5[i] << endl;
    }
    cout << "Szerokosc pasma sygnalu zA(t) c) wynosi " << szerokosc_Pasma(decybel5, skala_C, t_K * fs) << endl;
  
    Zp_c = DFT(zp_c, t_K * fs);
    decybel6 = widmo_A_Dec(Zp_c, t_K * fs);
    ofstream zad2_c_zp("zad2_c_zp.dat");
    for (int i = 0; i < t_K * fs; i++) {
        zad2_c_zp << skala_C[i] << " " << decybel6[i] << endl;
    }
    cout << "Szerokosc pasma sygnalu zp(t) c) wynosi " << szerokosc_Pasma(decybel6, skala_C, t_K * fs) << endl;
   
    //zadanie 3

    //Szerokosc pasma sygnalu zA(t) a) wynosi 4
    //Szerokosc pasma sygnalu zp(t) a) wynosi 5
    //Szerokosc pasma sygnalu zA(t) b) wynosi 5
    //Szerokosc pasma sygnalu zp(t) b) wynosi 3
    //Szerokosc pasma sygnalu zA(t) c) wynosi 5
    //Szerokosc pasma sygnalu zp(t) c) wynosi 58

    return 0;
}
