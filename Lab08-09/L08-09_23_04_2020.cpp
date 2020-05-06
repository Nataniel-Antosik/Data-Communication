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

int * dekoder_TTL(int * wejscie, int czestotliwosc, int bity) {
    int* wyjscie = new int[bity];
    int n = 0;
    for (int i = 0; i < bity; i++) {
            wyjscie[i] = (wejscie[n] == 1) ? 1 : 0;
            n += czestotliwosc;
    }
    return wyjscie;
}

int * Manchester(int * TTL, int * CLK, int czestotliwosc, int bity) {
    int const rozmiar = czestotliwosc * bity;
    int* wyjscie = new int[rozmiar]();
    for (int i = 0; i < rozmiar; i++) {
        if (TTL[i] == 0 && CLK[i] == 0) {
            wyjscie[i] = -1;
        }
        else if (TTL[i] == 0 && CLK[i] == 1) {
            wyjscie[i] = 1;
        }
        else if (TTL[i] == 1 && CLK[i] == 0) {
            wyjscie[i] = 1;
        }
        else if (TTL[i] == 1 && CLK[i] == 1) {
            wyjscie[i] = -1;
        }
    }
    return wyjscie;
}

int* dekoder_Manchester(int* Manchester, int czestotliwosc, int bity) {
    int* wyjscie = new int[bity]();
    int k = 0;
    for (int i = 0; i < bity; i++) {
        wyjscie[i] = (Manchester[k] == -1) ? 1 : 0;
        k += czestotliwosc;
    }
    return wyjscie;
}

int* NRZI(int* sygnal, int czestotliwosc, int bity) { //sprawdzamy wartość sygnału
    int const rozmiar = czestotliwosc * bity;
    int* wyjscie = new int[rozmiar]();
    int NRZI_poprzedni = 1; //Ta zmienna sprawdza poprzedni sygnał
    double n = czestotliwosc;
    int k = 0; //do próbkowania
    for (int i = 0; i < bity; i++) {
        if (sygnal[i] == 1) { 
            if (NRZI_poprzedni == 1) NRZI_poprzedni = -1;
            else NRZI_poprzedni = 1;
        }
        for (int j = 0; j < n; j++) {
            wyjscie[k] = NRZI_poprzedni;
            k++;
        }    
    }
    return wyjscie;
}

int* dekoder_NRZI(int* NRZI, int czestotliwosc, int bity) {
    int n = czestotliwosc;
    int NRZI_poprzedni;
    int* wyjscie = new int[bity]();
    int k = 0;
    wyjscie[0] = 0;
    for (int i = 1; i < bity - 1; i++) {
        wyjscie[i] = (NRZI[k] ^ NRZI[k + n]) ? 1 : 0;
        k += n;
    }
    wyjscie[bity - 1] = 0;
    return wyjscie;
}

int* BAMI(int* sygnal, int czestotliwosc, int bity) {
    int const rozmiar = czestotliwosc * bity;
    int* wyjscie = new int[rozmiar]();
    int BAMI_poprzedni = 1;
    double n = czestotliwosc;
    int k = 0;
    int poprzedni_Sygnal = 1; //zmieniamy jak napotkamy 1 w sygnale
    for (int i = 0; i < bity; i++) {
        if (sygnal[i] == 0) {
            BAMI_poprzedni = 0;
        }
        else {
            if (poprzedni_Sygnal == 1) poprzedni_Sygnal = -1;
            else poprzedni_Sygnal = 1;
            BAMI_poprzedni = poprzedni_Sygnal;
        }  
        for (int j = 0; j < n; j++) {
            wyjscie[k] = BAMI_poprzedni;
            k++;
        }
    }
    return wyjscie;
}

int* dekoder_BAMI(int* BAMI, int czestotliwosc, int bity) {
    int n = czestotliwosc;
    int k = 0;
    int* wyjscie = new int[bity]();
    for (int i = 0; i < bity; i++) {
        wyjscie[i] = (BAMI[k] != 0) ? 1 : BAMI[k]; 
        k += czestotliwosc;
    }
    return wyjscie;
}

int main()
{
    char d[4] = "ABC";
    string d1 = S2BS(d, 1);
    cout << "Little endian:  " << d1 << endl;

    double FS = 100, TB = 1;
    int const bity = 16;
    int* zegar, * tab, * ABC, * Man, * nrzi ,* d_TTL, * bami, * d_BAMI, * d_NRZI, * d_Manchester;
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

    //===================Kodowanie==================//
    for (int i = 0; i < bity * FS; i++) { ABC[i] = TTL(ABC[i]); } // Sygnał TTL
    Man = Manchester(ABC, zegar, FS, bity);
    nrzi = NRZI(a1, FS, bity);
    bami = BAMI(a1, FS, bity);   
    //==============================================//

    //=============Sygnaly zdekodowane==============//
    d_TTL = dekoder_TTL(ABC, FS, bity);
    d_BAMI = dekoder_BAMI(bami, FS, bity);
    d_NRZI = dekoder_NRZI(nrzi, FS, bity);
    d_Manchester = dekoder_Manchester(Man, FS, bity);

    cout << endl << "Zdekodowany sygnal TTL:          ";
    for (int i = 0; i < bity; i++) {
        cout << d_TTL[i];
    }
    cout << endl << "Zdekodowany sygnal BAMI:         ";
    for (int i = 0; i < bity; i++) {
        cout << d_BAMI[i];
    }
    cout << endl << "Zdekodowany sygnal NRZI:         ";
    for (int i = 0; i < bity; i++) {
        cout << d_NRZI[i];
    }
    cout << endl << "Zdekodowany sygnal Manchester:   ";
    for (int i = 0; i < bity; i++) {
        cout << d_Manchester[i];
    }
    //==============================================//

    //================Zapis do pliku================//
    ofstream CLK("CLK.dat");
    for (int i = 0; i < bity * FS; i++) {
        CLK << probka[i] << " " << zegar[i] << endl;
    }

    ofstream TTL("TTL.dat");
    for (int i = 0; i < bity * FS; i++) {
        TTL << probka[i] << " " << ABC[i] << endl;
    }

    ofstream Manchester("Manchester.dat");
    for (int i = 0; i < bity * FS; i++) {
        Manchester << probka[i] << " " << Man[i] << endl;
    }

    ofstream NRZI("NRZI.dat");
    for (int i = 0; i < bity * FS; i++) {
        NRZI << probka[i] << " " << nrzi[i] << endl;
    }

    ofstream BAMI("BAMI.dat");
    for (int i = 0; i < bity * FS; i++) {
        BAMI << probka[i] << " " << bami[i] << endl;
    }
    //==============================================//
}
