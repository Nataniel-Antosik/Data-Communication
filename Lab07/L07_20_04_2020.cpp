#define _USE_MATH_DEFINES
#include <sstream>
#include <iostream>
#include <fstream>
#include <string>
#include <cstddef>
#include <bitset>

using namespace std;

//Kluczowanie amplitudy (ASK)
double zA(double A1, double A2, double t, double f, double phi, int dane) {
    if (dane == 0) { return (A1 * sin(2 * M_PI * f * t + phi)); }
    else { return  (A2 * sin(2 * M_PI * f * t + phi)); }
}

double zA1(double A2, double t, double f, double phi) {
    return  (A2 * sin(2 * M_PI * f * t + phi));
}
//Kluczowanie częstotliwości (FSK)
double zF(double A, double t, double f0, double f1, double phi, int dane) {
    if (dane == 0) { return (A * sin(2 * M_PI * f0 * t + phi)); }
    else { return (A * sin(2 * M_PI * f1 * t + phi)); }
}

double zF1(double A, double t, double f0, double phi) {
    return (A * sin(2 * M_PI * f0 * t + phi));
}

double zF2(double A, double t, double f1, double phi) {
    return (A * sin(2 * M_PI * f1 * t + phi));
}

//Kluczowanie fazy (PSK)
double zP(double A, double t, double f, double phi0, double phi1, int dane) {
    if (dane == 0) { return (A * sin(2 * M_PI * f * t + phi0)); }
    else { return (A * sin(2 * M_PI * f * t + phi1)); }
}

double zP1(double A, double t, double f, double phi1) {
    return (A * sin(2 * M_PI * f * t + phi1));
}

int* konwersja(string dane, int* wyjscie) {
    for (int i = 0; i <= dane.size(); i++) {
        if (dane[i] == '0') { wyjscie[i] = 0; }
        else if (dane[i] == '1') { wyjscie[i] = 1; }
    }
    return wyjscie;
}

int* rozszerzenie(int* wejscie, int* wyjscie, int czestotliwosc, int Tb, int bity) {
    int k = 0;
    for (int i = 0; i < bity; i++) {
        if (wejscie[i] == 0) {
            for (int j = 0; j < czestotliwosc * Tb; j++) {
                wyjscie[k] = 0;
                k += 1;
            }
        }
        else {
            for (int j = 0; j < czestotliwosc * Tb; j++) {
                wyjscie[k] = 1;
                k += 1;
            }
        }
    }
    return wyjscie;
}

string S2BS(char in[], bool Switch) {           //String to Binary Stream
    int dlugosc = strlen(in) - 1;               //długość łańcucha - 1
    ostringstream str;                          //Tworzymy obiekt który jest bufforem, przechowujemy po koleji bity
    string wynik;

    if (Switch == true) {                       //little endina
                                                //zapisu danych, w której najmniej znaczący bajt, umieszczony jest jako pierwszy
        for (int i = 0; i <= dlugosc; i++) {    //przechodzimy po całym łańcuchu
            bitset<8> x(in[i]);                 //zamiana liczby numerycznej na bity
            str << x;                           //dodanie do buffora bitów
        }

        wynik = str.str();                      //przesłanie buffora do wyniku
        return wynik;
    }
    else {                                      //big endian
                                                //zapisu danych, w której najbardziej znaczący bajt, umieszczony jest jako pierwszy
        for (int i = 0; i <= dlugosc; i++) {    //przechodzimy po całym łańcuchu
            bitset<8> x(in[i]);                 //zamiana liczby numerycznej na bity
            str << x;                           //dodanie do buffora bitów
        }

        wynik = str.str();
        string odwr;
        int n_Dlugosc = wynik.length();
        for (int i = n_Dlugosc - 1; i >= 0; i--) {
            odwr += wynik[i];
        }
        return odwr;
    }
}

double* x(double* wejscie1, double* wejscie2, double* wyjscie, double bity, double czestotliwosc) {
    for (int i = 0; i < czestotliwosc * bity; i++) {
        wyjscie[i] = wejscie1[i] * wejscie2[i];
    }
    return wyjscie;
}

double* p(double* wejscie1, double* wejscie2, double* wyjscie, double bity, double czestotliwosc) { //całka
    int j = 0, k = 0, n = czestotliwosc;
    double suma = 0;
    for (int i = 0; i < bity; i++) {
        for (j; j < n; j++) {
            suma += wyjscie[j] * 1 / czestotliwosc;
        }
        for (k; k < n; k++) {
            wyjscie[k] = suma;
        }
        n += czestotliwosc;
        suma = 0;
    }
    return wyjscie;
}

double* demodulator_ASK_PSK(double* wejscie1, double* wejscie2, double* wyjscie, double bity, double czestotliwosc, double h) {

    wyjscie = x(wejscie1, wejscie2, wyjscie, bity, czestotliwosc);
    wyjscie = p(wejscie1, wejscie2, wyjscie, bity, czestotliwosc);

    for (int i = 0; i < czestotliwosc * bity; i++)
        wyjscie[i] = (wyjscie[i] > h) ? 1 : 0;

    return wyjscie;
}

double* demodulator_FSK(double* wejscie1, double* wejscie2, double* wejscie3, double* wyjscie, double bity, double czestotliwosc, int h) {
    int rozmiar = czestotliwosc * bity;
    double* F1 = new double[rozmiar]();
    double* F2 = new double[rozmiar]();

    F1 = x(wejscie1, wejscie2, F1, bity, czestotliwosc);
    F2 = x(wejscie1, wejscie3, F2, bity, czestotliwosc);
    F1 = p(wejscie1, wejscie2, F1, bity, czestotliwosc);
    F2 = p(wejscie1, wejscie3, F2, bity, czestotliwosc);

    for (int i = 0; i < rozmiar; i++) {
        wyjscie[i] = F2[i] - F1[i];
    }

    for (int i = 0; i < czestotliwosc * bity; i++)
        wyjscie[i] = (wyjscie[i] > h) ? 1 : 0;

    return wyjscie;
}

int main()
{
    double A1, A2, A, F1, F0, F, Phi0, Phi1, Phi, Fs, Tb, N;
    double* zA_1, * zP_1, * zF_1, * zA_1_1, * zP_1_1, * zF_1_1, * zF_1_2, * zA_1_x, * zP_1_x, * zF_1_x_1, * zF_1_x_2, * zA_1_p, * zP_1_p, * zF_1_p_1, * zF_1_p_2, * zF_1_p;
    double* zA_1_m, * zP_1_m, * zF_1_m;

    char d[4] = "ABC";
    string d1 = S2BS(d, 1);

    cout << "Little endian:  " << d1 << endl;

    double bit10 = 0, bity = 0;
    //test//
    int rozmiar = d1.length() - 1;
    int* a1;
    a1 = new int[rozmiar]();
    konwersja(d1, a1);

    cout << endl;
    for (int i = 0; i <= rozmiar; i++) {
        cout << a1[i];
    }
    cout << endl << d1;
    cout << endl << "Rozmiar: " << rozmiar << endl;


    //k = rozszerzenie(a1, rozmiar, k);
    Tb = 1; //czas trwania jednego bitu
    N = 2;
    F0 = (N + 1) / Tb;
    F1 = (N + 2) / Tb;
    Phi0 = 0;
    Phi1 = M_PI;
    Phi = 4;
    Fs = 100;
    A1 = 1;
    A2 = 2;
    A = 4;
    F = N * Tb;

    zA_1 = new double[Fs * rozmiar];
    zP_1 = new double[Fs * rozmiar];
    zF_1 = new double[Fs * rozmiar];
    zA_1_1 = new double[Fs * rozmiar];
    zP_1_1 = new double[Fs * rozmiar];
    zF_1_1 = new double[Fs * rozmiar];
    zF_1_2 = new double[Fs * rozmiar];
    zA_1_x = new double[Fs * rozmiar];
    zP_1_x = new double[Fs * rozmiar];
    zF_1_x_1 = new double[Fs * rozmiar];
    zF_1_x_2 = new double[Fs * rozmiar];
    zA_1_p = new double[Fs * rozmiar];
    zP_1_p = new double[Fs * rozmiar];
    zF_1_p_1 = new double[Fs * rozmiar];
    zF_1_p_2 = new double[Fs * rozmiar];
    zA_1_m = new double[Fs * rozmiar];
    zP_1_m = new double[Fs * rozmiar];
    zF_1_m = new double[Fs * rozmiar];
    zF_1_p = new double[Fs * rozmiar];

    bity = rozmiar; //rozmiar od stringa bitowego
    int* k2 = new int[Fs * bity];
    k2 = rozszerzenie(a1, k2, Fs, Tb, bity);

    double* probka2 = new double[Fs * rozmiar](); //Próbkowanie
    double tmp2 = 0;

    for (int i = 0; i < Fs * rozmiar; i++) {
        probka2[i] = tmp2;
        tmp2 += 0.01;
        zA_1[i] = zA(A1, A2, probka2[i], F, Phi, k2[i]);
        zA_1_1[i] = zA1(A2, probka2[i], F, Phi);
        zP_1[i] = zP(A, probka2[i], F, Phi0, Phi1, k2[i]);
        zP_1_1[i] = zP1(A, probka2[i], F, Phi1);
        zF_1[i] = zF(A, probka2[i], F0, F1, Phi, k2[i]);
        zF_1_1[i] = zF1(A, probka2[i], F0, Phi);
        zF_1_2[i] = zF2(A, probka2[i], F1, Phi);
    }

    zA_1_m = demodulator_ASK_PSK(zA_1, zA_1_1, zA_1_m, bity, Fs, 1.5);
    zP_1_m = demodulator_ASK_PSK(zA_1, zP_1_1, zP_1_m, bity, Fs, 2);
    zF_1_m = demodulator_FSK(zF_1, zF_1_1, zF_1_2, zF_1_m, bity, Fs, 5);

    zA_1_x = x(zA_1, zA_1_1, zA_1_x, bity, Fs);
    cout << zA_1_x[1];
    for (int i = 0; i < Fs * rozmiar; i++) {
        zA_1_p[i] = zA_1_x[i];
    }  
    zA_1_p = p(zA_1, zA_1_1, zA_1_p, bity, Fs);

    zP_1_x = x(zP_1, zP_1_1, zP_1_x, bity, Fs);
    for (int i = 0; i < Fs * rozmiar; i++) {
        zP_1_p[i] = zP_1_x[i];
    }
    zP_1_p = p(zP_1, zP_1_1, zP_1_p, bity, Fs);

    zF_1_x_1 = x(zF_1, zF_1_1, zF_1_x_1, bity, Fs);
    for (int i = 0; i < Fs * rozmiar; i++) {
        zF_1_p_1[i] = zF_1_x_1[i];
    }
    zF_1_p_1 = p(zF_1, zF_1_1, zF_1_p_1, bity, Fs);

    zF_1_x_2 = x(zF_1, zF_1_2, zF_1_x_2, bity, Fs);
    for (int i = 0; i < Fs * rozmiar; i++) {
        zF_1_p_2[i] = zF_1_x_2[i];
    }
    zF_1_p_2 = p(zF_1, zF_1_2, zF_1_p_2, bity, Fs);

    for (int i = 0; i < Fs * rozmiar; i++) {
        zF_1_p[i] = zF_1_p_2[i] - zF_1_p_1[i];
    }
    //========================================================//
    ofstream zad_zA("zad_zA.dat");
    for (int i = 0; i < Fs * rozmiar; i++) {
        zad_zA << probka2[i] << " " << zA_1[i] << endl;
    }
    ofstream zad_zA_x("zad_zA_x.dat");
    for (int i = 0; i < Fs * rozmiar; i++) {
        zad_zA_x << probka2[i] << " " << zA_1_x[i] << endl;
    }
    ofstream zad_zA_p("zad_zA_p.dat");
    for (int i = 0; i < Fs * rozmiar; i++) {
        zad_zA_p << probka2[i] << " " << zA_1_p[i] << endl;
    }
    ofstream zad_zA_m("zad_zA_m.dat");
    for (int i = 0; i < Fs * bity; i++) {
        zad_zA_m << probka2[i] << " " << zA_1_m[i] << endl;
    }
    //========================================================//
    ofstream zad_zP("zad_zP.dat");
    for (int i = 0; i < Fs * rozmiar; i++) {
        zad_zP << probka2[i] << " " << zP_1[i] << endl;
    }
    ofstream zad_zP_x("zad_zP_x.dat");
    for (int i = 0; i < Fs * rozmiar; i++) {
        zad_zP_x << probka2[i] << " " << zP_1_x[i] << endl;
    }
    ofstream zad_zP_p("zad_zP_p.dat");
    for (int i = 0; i < Fs * rozmiar; i++) {
        zad_zP_p << probka2[i] << " " << zP_1_p[i] << endl;
    }
    ofstream zad_zP_m("zad_zP_m.dat");
    for (int i = 0; i < Fs * bity; i++) {
        zad_zP_m << probka2[i] << " " << zP_1_m[i] << endl;
    }
    //========================================================//
    ofstream zad_zF("zad_zF.dat");
    for (int i = 0; i < Fs * rozmiar; i++) {
        zad_zF << probka2[i] << " " << zF_1[i] << endl;
    }
    ofstream zad_zF_x_1("zad_zF_x_1.dat");
    for (int i = 0; i < Fs * rozmiar; i++) {
        zad_zF_x_1 << probka2[i] << " " << zF_1_x_1[i] << endl;
    }
    ofstream zad_zF_x_2("zad_zF_x_2.dat");
    for (int i = 0; i < Fs * rozmiar; i++) {
        zad_zF_x_2 << probka2[i] << " " << zF_1_x_2[i] << endl;
    }
    ofstream zad_zF_p("zad_zF_p.dat");
    for (int i = 0; i < Fs * rozmiar; i++) {
        zad_zF_p << probka2[i] << " " << zF_1_p[i] << endl;
    }
    ofstream zad_zF_m("zad_zF_m.dat");
    for (int i = 0; i < Fs * bity; i++) {
        zad_zF_m << probka2[i] << " " << zF_1_m[i] << endl;
    }

    return 0;
}
