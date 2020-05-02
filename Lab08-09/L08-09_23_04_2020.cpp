#define _USE_MATH_DEFINES
#include <sstream>
#include <iostream>
#include <fstream>
#include <string>
#include <cstddef>
#include <bitset>

using namespace std;

int* konwersja(string dane, int* wyjscie) {
    for (int i = 0; i <= dane.size(); i++) {
        if (dane[i] == '0') { wyjscie[i] = 0; }
        else if (dane[i] == '1') { wyjscie[i] = 1; }
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

int* CLK(int* wejscie, int czestotliwosc, int bity) {
    double n = czestotliwosc;
    for (int i = 0; i < bity; i++) {
        for (int j = 0; j < n; j++) {
            wejscie[j + (i * czestotliwosc)] = (j < n / 2.0) ? 1 : 0;
        }
    }
    return wejscie;
}

double TTL(double wejscie) { return wejscie; }

double Manchester() {
    return 0;
}

double NRZI() {
    return 0;
}

double BAMI() {
    return 0;
}

int main()
{
    char d[4] = "ABC";
    string d1 = S2BS(d, 1);
    cout << "Little endian:  " << d1 << endl;

    double FS = 100, TB = 1;
    int const bity = 16;
    int* zegar, * tab, * ABC;
    double tmp;
    double* probka;
    
    //test//
    int rozmiar = d1.length() - 1;
    int* a1;
    a1 = new int[rozmiar]();
    konwersja(d1, a1);

    tab = new int[bity]; //tworzenie sygnału dla zegara  
    zegar = new int[bity * FS]();
    ABC = new int[bity * FS];
    probka = new double[bity * FS];
    tmp = 0;
    
    for (int i = 0; i < bity * FS; i++) {
        probka[i] = tmp;
        tmp += 0.01;
    }
   
    zegar = CLK(zegar, FS, bity);

    ABC = rozszerzenie(a1, ABC, FS, TB, bity);
    for (int i = 0; i < FS * bity; i++) {
        cout << zegar[i] << " ";
    }


    ofstream CLK("CLK.dat");
    for (int i = 0; i < bity * FS; i++) {
        CLK << probka[i] << " " << zegar[i] << endl;
    }

    ofstream TTL("TTL.dat");
    for (int i = 0; i < bity * FS; i++) {
        TTL << probka[i] << " " << ABC[i] << endl;
    }
}
