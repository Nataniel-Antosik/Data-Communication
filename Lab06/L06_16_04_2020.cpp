#define _USE_MATH_DEFINES
#include <sstream>
#include <iostream>
#include <fstream>
#include <string>
#include <cstddef>
#include <bitset>

using namespace std;

//44261 A = 1, B = 6, C=2

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

double szerokosc_Pasma(double* wartosc_F, double* czestotliwosc, int rozmiar) {
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

//Kluczowanie amplitudy (ASK)
double zA(double A1, double A2, double t, double f,double phi ,int dane) {
        if (dane == 0) { return (A1 * sin(2 * M_PI * f * t + phi)); }
        else { return  (A2 * sin(2 * M_PI * f * t + phi)); }
}
//Kluczowanie częstotliwości (FSK)
double zF(double A, double t, double f0, double f1, double phi, int dane) {  
        if (dane == 0) { return (A * sin(2 * M_PI * f0 * t + phi)); }
        else { return (A * sin(2 * M_PI * f1 * t + phi)); }  
}
//Kluczowanie fazy (PSK)
double zP(double A, double t, double f, double phi0, double phi1, int dane) {
        if (dane == 0) { return (A * sin(2 * M_PI * f * t + phi0)); }
        else { return (A * sin(2 * M_PI * f * t + phi1)); }
}

int *konwersja(string dane,int * wyjscie) {
    for (int i = 0; i <= dane.size(); i++) {
        if (dane[i] == '0') { wyjscie[i] = 0; }
        else if (dane[i] == '1') { wyjscie[i] = 1; }
    }
    return wyjscie;
}

int* rozszerzenie(int * wejscie, int * wyjscie, int czestotliwosc, int Tb, int bity) {   
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

int main()
{
    char d[26] = "abcdefghchijklmnoprstuwyz";
    string d1 = S2BS(d, 1);
    
    cout << "Little endian:  " << d1 << endl;
    double A1, A2, A, F1, F0, F, Phi0, Phi1, Phi, Fs, Tb, N;
    int bit10 = 0, bit = 0;
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
    int a = 0;
    
    //k = rozszerzenie(a1, rozmiar, k);
    Tb = 1;
    N = 2;
    F0 = (N + 1) / Tb;
    F1 = (N + 2) / Tb;
    Phi0 = 0;
    Phi1 = M_PI;
    Phi = 4;
    Fs = 100;   
    A1 = 1; 
    A2 = 2;
    A = 3; 
    F = N * Tb;
    bit10 = 10;
    int* k = new int[Fs * Tb * bit10];
    k = rozszerzenie(a1, k, Fs, Tb, bit10);
    
      
    //zad 2
    int* k2 = new int[Fs * rozmiar];
    bit = rozmiar;
    k2 = rozszerzenie(a1, k, Fs, Tb, bit);

    double* probka2 = new double[Fs * rozmiar](); //Próbkowanie
    double tmp2 = 0;
    for (int i = 0; i < Fs * rozmiar; i++) {
        probka2[i] = tmp2;
        tmp2 += 0.01;
    }

    double tmp = 0;
    double* probka = new double[Fs * 10](); //Próbkowanie
    tmp = 0;
    for (int i = 0; i < Fs * 10; i++) {
        probka[i] = tmp;
        tmp += 0.01;
    }

    zespolona* zA_W, * zF_W, * zP_W;
    double* decybel, * skala_C, * zA_A, * zP_A, * zF_A;
    
    //=====Deklaracja_pamięci======//
    zA_W = new zespolona[Fs * bit]();
    zF_W = new zespolona[Fs * bit]();
    zP_W = new zespolona[Fs * bit]();

    zA_A = new double[Fs * bit]();
    zP_A = new double[Fs * bit]();
    zF_A = new double[Fs * bit]();

    decybel = new double[Fs * bit]();

    skala_C = skala_Czestotliwosci(Fs, Fs * bit);
    //=============================//
    
    ofstream zad2_zA("zad2_zA.dat");
    for (int i = 0; i < Fs * rozmiar; i++) {
        zad2_zA << probka2[i] << " " << zA(A1, A2, probka2[i], F, Phi, k2[i]) << endl;
        zA_A[i] = zA(A1, A2, probka2[i], F, Phi, k2[i]);
    }

    ofstream zad2_zF("zad2_zF.dat");
    for (int i = 0; i < Fs * rozmiar; i++) {
        zad2_zF << probka2[i] << " " << zF(A, probka2[i], F0, F1, Phi, k2[i]) << endl;
        zF_A[i] = zF(A, probka2[i], F0, F1, Phi, k2[i]);
    }

    ofstream zad2_zP("zad2_zP.dat");
    for (int i = 0; i < Fs * rozmiar; i++) {
        zad2_zP << probka2[i] << " " << zP(A, probka2[i], F, Phi0, Phi1, k2[i]) << endl;
        zP_A[i] = zP(A, probka2[i], F, Phi0, Phi1, k2[i]);
    }
    
    ofstream zad_inf("zad_inf.dat");
    for (int i = 0; i < Fs * 10; i++) {
        zad_inf << probka[i] << " " << k[i] << endl;
    }

    ofstream zad3_zA("zad3_zA.dat");
    for (int i = 0; i < Fs * 10; i++) {
        zad3_zA << probka[i] << " " << zA(A1, A2, probka[i], F, Phi, k[i]) << endl;      
    }

    ofstream zad3_zF("zad3_zF.dat");
    for (int i = 0; i < Fs * 10; i++) {
        zad3_zF << probka[i] << " " << zF(A, probka[i], F0, F1, Phi, k[i]) << endl;       
    }

    ofstream zad3_zP("zad3_zP.dat");
    for(int i = 0; i < Fs * 10; i++){
        zad3_zP << probka[i] << " " << zP(A, probka[i], F, Phi0, Phi1, k[i]) << endl;      
    }

    zA_W = DFT(zA_A, Fs * bit);
    decybel = widmo_A_Dec(zA_W, Fs * bit);
    ofstream zad4_zA("zad4_zA.dat");
    for (int i = 0; i < Fs * bit; i++) {
        zad4_zA << skala_C[i] << " " << decybel[i] << endl;
    }
    cout << "Szerokosc pasma sygnalu zA(t) a) wynosi " << szerokosc_Pasma(decybel, skala_C, Fs * bit) << endl;

    zF_W = DFT(zF_A, Fs * bit);
    decybel = widmo_A_Dec(zF_W, Fs * bit);
    ofstream zad4_zF("zad4_zF.dat");
    for (int i = 0; i < Fs * bit; i++) {
        zad4_zF << skala_C[i] << " " << decybel[i] << endl;
    }
    cout << "Szerokosc pasma sygnalu zF(t) a) wynosi " << szerokosc_Pasma(decybel, skala_C, Fs * bit) << endl;

    zP_W = DFT(zP_A, Fs * bit);
    decybel = widmo_A_Dec(zP_W, Fs * bit);
    ofstream zad4_zP("zad4_zP.dat");
    for (int i = 0; i < Fs * bit; i++) {
        zad4_zP << skala_C[i] << " " << decybel[i] << endl;
    }
    cout << "Szerokosc pasma sygnalu zP(t) a) wynosi " << szerokosc_Pasma(decybel, skala_C, Fs * bit) << endl;
    
    //Szerokość pasma sygnału zA(t) wynosi: 1
    //Szerokość pasma sygnału zF(t) wynosi: 0.6
    //Szerokość pasma sygnału zP(t) wynosi: 1.5

    return 0;
}

