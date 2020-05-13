#define _USE_MATH_DEFINES
#include <iostream>
#include <sstream>
#include <string>
#include <bitset>

using namespace std;

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

int* konwersja(string dane, int* wyjscie) {
    for (int i = 0; i <= dane.size(); i++) {
        if (dane[i] == '0') { wyjscie[i] = 0; }
        else if (dane[i] == '1') { wyjscie[i] = 1; }
    }
    return wyjscie;
}

int * kod_Hamminga(int * bity, bool wybor = true) { 
    int* h; //wynik
    int p4 = 0; //zmienna pomocnicza dla 4 bitu parzystości 
    h = (wybor != true) ? new int[7]() : new int[8](); //true - SECDED   
   
    h[0] = (bity[0] + bity[1] + bity[3]) % 2;
    h[1] = (bity[0] + bity[2] + bity[3]) % 2;
    h[2] = bity[0];
    h[3] = (bity[1] + bity[2] + bity[3]) % 2;
    h[4] = bity[1];
    h[5] = bity[2];
    h[6] = bity[3];

    if (wybor != true) { //SECDEC
        for (int i = 0; i < 7; i++) { p4 += h[i]; }
        p4 %= 2;
        h[7] = p4;
    }
    return h;
}

void negacja(int* h, int indeks) {
    if (indeks >= 0 && indeks < 8) { h[indeks] = !h[indeks]; }
    else { cout << endl << "Wprowadzono zly indeks" << endl; }
}

void indeks_Korekty(int* h) {
    int* p = new int[3];
    p[0] = (h[0] + h[2] + h[4] + h[6]) % 2;
    p[1] = (h[1] + h[2] + h[5] + h[6]) % 2;
    p[2] = (h[3] + h[4] + h[5] + h[6]) % 2;
    int n = p[0] * 1 + p[1] * 2 + p[2] * 4;
    if (n > 0) { negacja(h, n - 1); }
}

int* dekoder_Hamminga(int* h) {
    int* wynik = new int[4];
    indeks_Korekty(h);
    wynik[0] = h[2];
    wynik[1] = h[4];
    wynik[2] = h[5];
    wynik[3] = h[6];
    return wynik;
}

int* dekoder_SECDEO(int* h) {
    //Krok 1) 
    int p4 = 0; 
    for (int i = 0; i < 7; i++) p4 += h[i]; //Weryfikacja p4
    p4 %= 2;
    if (p4 != h[7]) { cout << endl << "Szansa na poprawny wynik: 50%" << endl; }
    //Krok 2), 3), 4)  
    indeks_Korekty(h); //sprawdzamy p1 p2 p3, wyznaczamy indeks korekty i negujemy wskazany bit
    //Krok 5)
    p4 = 0;
    for (int i = 0; i < 7; i++) p4 += h[i]; //Ponowna weryfikacja p4
    p4 %= 2;
    //Krok 6)
    if (p4 != h[7]) { 
        cout << endl << "Sa 2 bledy, wykonaj ponowna transmisje" << endl; 
        return NULL;
    }
    cout << endl << "Jest 1 blad, zostal juz naprawiony" << endl;
    return dekoder_Hamminga(h);
}

int main()
{
    char d[5] = "abc";
    string d1 = S2BS(d, 1);
    cout << "Little endian:  " << d1 << endl;
    //011000010110001001100011

    const int rozmiar = d1.length() - 1;
    int const bity = d1.length();

    int* Sygnal = new int[rozmiar]();
    int** Hamming = new int* [6]();      //wynik kodu Hamminga (7,4)
    int** Hamming8 = new int* [6]();    //wynik kodu Hamminga (8,4)
    int** dekodowanie = new int* [6];   //wynik dekodowania kodu Hamminga (7,4)
    int** dekodowanie8 = new int* [6];  //wynik dekodowania kodu Hamminga (8,4)

    konwersja(d1, Sygnal);
    cout << "Sygnal: " << endl;
    for (int i = 0; i < rozmiar + 1; i++) {
        cout << Sygnal[i];
    }
    
    cout << endl << endl << "Kod Hamminga (7,4): " << endl;
    cout << "(p1,p2,d1,p3,d2,d3,d4)" << endl;
    for (int i = 0; i < 6; i++) {
        Hamming[i] = kod_Hamminga(&Sygnal[i * 4]);
        for (int j = 0; j < 7; j++) {
            cout << Hamming[i][j];
        }
        cout << endl;
    }
    cout << "Dekodowanie kodu Hamminga (7,4): " << endl;
    for (int i = 0; i < 6; i++) {
        dekodowanie[i] = dekoder_Hamminga(Hamming[i]);
        for (int j = 0; j < 4; j++) {
            cout << dekodowanie[i][j];
        }
    }
    cout << endl;
    
    cout << endl << "Kod Hamminga (8,4): " << endl;
    cout << "(p1,p2,d1,p3,d2,d3,d4,p4)" << endl;
    for (int i = 0; i < 6; i++) {
        Hamming8[i] = kod_Hamminga(&Sygnal[i * 4], false);
        for (int j = 0; j < 8; j++) {
            cout << Hamming8[i][j];
        }
        cout << endl;
    }
    cout << "Dekodowanie kodu Hamminga (8,4): " << endl;
    for (int i = 0; i < 6; i++) {
        dekodowanie8[i] = dekoder_Hamminga(Hamming8[i]);
        for (int j = 0; j < 4; j++) {
            cout << dekodowanie8[i][j];
        }
    }
    cout << endl;
    //wynik: 
    //1100110
    //1101001
    //1100110
    //0101010
    //1100110
    //1000011
    return 0;
}
