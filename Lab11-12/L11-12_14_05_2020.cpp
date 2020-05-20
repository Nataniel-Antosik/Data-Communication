#define _USE_MATH_DEFINES
#include <iostream>
#include <sstream>
#include <string>
#include <bitset>
#include <fstream>
#include <string>

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

int* kod_Hamminga(int* bity, bool wybor = true) {
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

    if (wybor != true) { //SECDED
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

void indeks_Korekty(int* h, int n) {
    int* p = new int[3];
    p[0] = (h[0] + h[2] + h[4] + h[6]) % 2;
    p[1] = (h[1] + h[2] + h[5] + h[6]) % 2;
    p[2] = (h[3] + h[4] + h[5] + h[6]) % 2;
    n = p[0] * 1 + p[1] * 2 + p[2] * 4;
    if (n > 0) { negacja(h, n - 1); }
}

int* dekoder_Hamminga(int* h) {
    int* wynik = new int[4];
    int n = 0;
    indeks_Korekty(h, n);
    wynik[0] = h[2];
    wynik[1] = h[4];
    wynik[2] = h[5];
    wynik[3] = h[6];
    return wynik;
}

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

double* x(double* wejscie1, double* wejscie2, double* wyjscie, double bity, double czestotliwosc, double TB) {
    for (int i = 0; i < czestotliwosc * bity * TB; i++) {
        wyjscie[i] = wejscie1[i] *  wejscie2[i];
    }
    return wyjscie;
}

double* p(double* wyjscie, double bity, double czestotliwosc, double TB) { //całka
    int j = 0, k = 0, n = czestotliwosc * TB;
    double suma = 0;
    for (int i = 0; i < bity; i++) {
        for (j; j < n; j++) {
            suma += wyjscie[j] * 1 / czestotliwosc * TB;
            
        }
        for (k; k < n; k++) {
            wyjscie[k] = suma;
        }
        n += czestotliwosc * TB;
        suma = 0;
    }
    return wyjscie;
}

double* demodulator_ASK_PSK(double* wejscie1, double* wejscie2, double* wyjscie, double bity, double czestotliwosc, double TB, double h) {

    wyjscie = x(wejscie1, wejscie2, wyjscie, bity, czestotliwosc, TB);
    wyjscie = p(wyjscie, bity, czestotliwosc, TB);

    for (int i = 0; i < czestotliwosc * bity * TB; i++)
        wyjscie[i] = (wyjscie[i] >= h) ? 1 : 0;

    return wyjscie;
}

double* demodulator_FSK(double* wejscie1, double* wejscie2, double* wejscie3, double* wyjscie, double bity, double czestotliwosc, double TB, double h) {
    double rozmiar = czestotliwosc * bity * TB;
    double* F1 = new double[rozmiar]();
    double* F2 = new double[rozmiar]();

    F1 = x(wejscie1, wejscie2, F1, bity, czestotliwosc, TB);
    F2 = x(wejscie1, wejscie3, F2, bity, czestotliwosc, TB);
    F1 = p(F1, bity, czestotliwosc, TB);
    F2 = p(F2, bity, czestotliwosc, TB);

    for (int i = 0; i < rozmiar; i++) {
        wyjscie[i] = F2[i] - F1[i];
    }

    for (int i = 0; i < czestotliwosc * bity * TB; i++)
        wyjscie[i] = (wyjscie[i] >= h) ? 1 : 0;

    return wyjscie;
}

int* nowa_Tablica(int** Hamming, int rozmiar, int r_p_bitow, bool tryb = true) {
    int* nowa_Tablica = new int[rozmiar]();
    int r = 0;
    int wybor = (tryb == false) ? 8 : 7;
    for (int i = 0; i < r_p_bitow; i++) {
        for (int j = 0; j < wybor; j++) {
            nowa_Tablica[r] = Hamming[i][j];
            r++;
        }
    }
    return nowa_Tablica;
}

int** nowa_Tablica2(int* d_Sygnal, int rozmiar, int r_p_bitow, bool tryb = true) {
    int** Hamming = new int* [rozmiar]();
    int wybor = (tryb == false) ? 8 : 7;
    int r = 0;
    for (int i = 0; i < rozmiar; i++) {
        Hamming[i] = &d_Sygnal[r];
        r++;
    }
    return Hamming;
}

int* rozszerzenie(int* wejscie, int* wyjscie, double czestotliwosc, double Tb, int bity) {
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

int SB2S() {
    return 0;
}

//===============
int main()
{
    //Dane
    double A1, A2, A, F1, F0, F, Phi0, Phi1, Phi, Fs, Tb, N;
    double* zA_1, * zP_1, * zF_1, * zA_1_1, * zP_1_1, * zF_1_1, * zF_1_2, * zA_1_m, * zP_1_m, * zF_1_m, * probka2;
    int* d_zA, * d_zP, * d_zF;

    const int r_Slowa = 12;             //rozmiar słowa wchodzącego
    char d[r_Slowa] = "ALA MA KOTA";
    string d1 = S2BS(d, true);             //konwersja string na strumień bitów

    cout << "Little endian:  " << endl << d1 << endl; //wypis strumienia
    //test//
    const int rozmiar = d1.length() - 1;
    int* Sygnal = new int[rozmiar]();
    konwersja(d1, Sygnal);
    cout << "Nowa tablica typu int: " << endl;

    int parita_Bitow = 0;               //ilość parti do przekazania do Hamminga
    for (int i = 0, j = 1; i <= rozmiar; i++, j++) {
        cout << Sygnal[i];
        if (j % 4 == 0) parita_Bitow++;
    }
    cout << endl << "ilosc parti 4 bitow: " << parita_Bitow << endl;
    cout << endl << "Ilosc bitow sygnalu: " << rozmiar << endl;

    Tb = 0.1;                           //czas trwania jednego bitu
    N = 2;
    F0 = 1.0;
    F1 = 4.0;
    Phi0 = 0;
    Phi1 = M_PI;
    Phi = 0;
    Fs = 1000;
    A1 = 0.0;
    A2 = 1.0;
    A = 1.0;
    F = 2.0;
    //===========Kodowanie kanałowe===========//
    int** Hamming = new int* [r_Slowa * 7]();       //wynik kodu Hamminga (7,4)
    int rozmiar_H = 0;                              //rozmiar dla tablicy trzymającej sygnał Hamminga;
    int* n_Hamming;                                 //Nowa tablica (pojedynczy wskaźnik na sygnał po kodowaniu Hamminga)
    int wybor = 7;
    cout << endl << endl << "Kod Hamminga (7,4): " << endl;
    cout << "(p1,p2,d1,p3,d2,d3,d4)" << endl;
    for (int i = 0; i < parita_Bitow; i++) {        //pętla po ilości pakietow po 4 bity do hamminga
        Hamming[i] = kod_Hamminga(&Sygnal[i * 4]);
        for (int j = 0; j < wybor; j++) {
            cout << Hamming[i][j];
            rozmiar_H++;
        }
        cout << endl;
    }
    cout << endl << "Rozmiar kodu Hamminga: " << rozmiar_H << endl;
    n_Hamming = nowa_Tablica(Hamming, rozmiar_H, parita_Bitow);
    cout << "Nowa tablica: " << endl;
    for (int i = 0; i < rozmiar_H; i++) {
        cout << n_Hamming[i];
    }

    //test
    const double rozmiar_Rozszerzony = rozmiar_H * Fs * Tb;//rozszerzenie
    int* n_r_Hamming = new int[rozmiar_Rozszerzony]();
    n_r_Hamming = rozszerzenie(n_Hamming, n_r_Hamming, Fs, Tb, rozmiar_H);

    cout << endl << "Testowe rozszerzenie: " << endl;
    int licz = 0;
    bool test_rozmiaru_Tablicy;
    for (int i = 0; i < rozmiar_Rozszerzony; i++) {
        licz++;
    }
    test_rozmiaru_Tablicy = (licz == rozmiar_Rozszerzony) ? true : false; //test dla rozszerzonego rozmiaru
    cout << endl << "Rozmiar tablicy po rozszerzeniu: " << licz << " Udane (1 - true, 0 - false): " << test_rozmiaru_Tablicy << endl;
    //===========Modulacja===========//
    zA_1 = new double[rozmiar_Rozszerzony];
    zP_1 = new double[rozmiar_Rozszerzony];
    zF_1 = new double[rozmiar_Rozszerzony];
    probka2 = new double[rozmiar_Rozszerzony];
    double tmp2 = 0;
    for (int i = 0; i < rozmiar_Rozszerzony; i++) {
        probka2[i] = tmp2;
        tmp2 += 0.01;
        zA_1[i] = zA(A1, A2, probka2[i], F, Phi, n_r_Hamming[i]);
        zP_1[i] = zP(A, probka2[i], F, Phi0, Phi1, n_r_Hamming[i]);
        zF_1[i] = zF(A, probka2[i], F0, F1, Phi, n_r_Hamming[i]);
    }
    cout << endl;
    //===========Wykresy===========//
    ofstream zad_zA("zad_zA.dat");
    for (int i = 0; i < rozmiar_Rozszerzony; i++) {
        zad_zA << probka2[i] << " " << zA_1[i] << endl;
    }
    ofstream zad_zP("zad_zP.dat");
    for (int i = 0; i < rozmiar_Rozszerzony; i++) {
        zad_zP << probka2[i] << " " << zP_1[i] << endl;
    }
    ofstream zad_zF("zad_zF.dat");
    for (int i = 0; i < rozmiar_Rozszerzony; i++) {
        zad_zF << probka2[i] << " " << zF_1[i] << endl;
    }
    //===========Demodulacja===========//
    zA_1_1 = new double[rozmiar_Rozszerzony];
    zP_1_1 = new double[rozmiar_Rozszerzony];
    zF_1_1 = new double[rozmiar_Rozszerzony];
    zF_1_2 = new double[rozmiar_Rozszerzony];
    zA_1_m = new double[rozmiar_Rozszerzony];
    zP_1_m = new double[rozmiar_Rozszerzony];
    zF_1_m = new double[rozmiar_Rozszerzony];
    d_zA = new int[rozmiar_H];
    d_zP = new int[rozmiar_H];
    d_zF = new int[rozmiar_H];

    /*testowanie całek
    double* zA_1_p = new double[rozmiar_Rozszerzony];
    double* zP_1_p = new double[rozmiar_Rozszerzony];
    double* zF_1_p = new double[rozmiar_Rozszerzony];
    double* zF_2_p = new double[rozmiar_Rozszerzony];
    */

    for (int i = 0; i < rozmiar_Rozszerzony; i++) {
        zA_1_1[i] = zA1(A2, probka2[i], F, Phi);
        zP_1_1[i] = zP1(A, probka2[i], F, Phi1);
        zF_1_1[i] = zF1(A, probka2[i], F0, Phi);
        zF_1_2[i] = zF2(A, probka2[i], F1, Phi);
    }

    zA_1_m = demodulator_ASK_PSK(zA_1, zA_1_1, zA_1_m, rozmiar_H, Fs, Tb, 0.0000001);
    zP_1_m = demodulator_ASK_PSK(zP_1, zP_1_1, zP_1_m, rozmiar_H, Fs, Tb, 0.0001);
    zF_1_m = demodulator_FSK(zF_1, zF_1_1, zF_1_2, zF_1_m, rozmiar_H, Fs, Tb, 0.00001);
    cout << endl << "Demodulacja ASK: " << endl;
    int n = 0;
    for (int i = 0; i < rozmiar_H; i++) {
        cout << zA_1_m[n];
        d_zA[i] = zA_1_m[n];
        n += 100;
    }
    /*
    n = 0;
    cout << endl << "Demodulacja PSK: " << endl;
    for (int i = 0; i < rozmiar_H; i++) {
        cout << zP_1_m[n];
        d_zP[i] = zP_1_m[n];
        n += 100;
    }
    n = 0;
    cout << endl << "Demodulacja FSK: " << endl;
    for (int i = 0; i < rozmiar_H; i++) {
        cout << zF_1_m[n];
        d_zF[i] = zF_1_m[n];
        n += 100;
    }*/
    /*
    zA_1_p = x(zA_1, zA_1_1, zA_1_p, rozmiar_H, Fs, Tb);
    zA_1_p = p(zA_1_p, rozmiar_H, Fs, Tb);

    zP_1_p = x(zP_1, zP_1_1, zP_1_p, rozmiar_H, Fs, Tb);
    zP_1_p = p(zP_1_p, rozmiar_H, Fs, Tb);

    zF_2_p = x(zF_1, zF_1_2, zF_2_p, rozmiar_H, Fs, Tb);
    zF_1_p = x(zF_1, zF_1_1, zF_1_p, rozmiar_H, Fs, Tb);
    zF_1_p = p(zF_1_p, rozmiar_H, Fs, Tb);
    zF_2_p = p(zF_2_p, rozmiar_H, Fs, Tb);
    double* roz = new double[rozmiar_Rozszerzony];
    for (int i = 0; i < rozmiar_Rozszerzony; i++) {
        roz[i] = zF_2_p[i] - zF_1_p[i];
    }
    ofstream zad_zA_p("zad_zA_p.dat");
    for (int i = 0; i < rozmiar_Rozszerzony; i++) {
        zad_zA_p << probka2[i] << " " << zA_1_p[i] << endl;
    }

    ofstream zad_zP_p("zad_zP_p.dat");
    for (int i = 0; i < rozmiar_Rozszerzony; i++) {
        zad_zP_p << probka2[i] << " " << zP_1_p[i] << endl;
    }

    ofstream zad_zF_p("zad_zF_p.dat");
    for (int i = 0; i < rozmiar_Rozszerzony; i++) {
        zad_zF_p << probka2[i] << " " << roz[i] << endl;
    }
    */
    //===========Dekodowanie kanałowe===========//
    int** d_zA_Hamming;
    int** d_zP_Hamming;
    int** d_zF_Hamming;
    int** dekodowanie_zA_H = new int* [rozmiar_H];
    int** dekodowanie_zP_H = new int* [rozmiar_H];
    int** dekodowanie_zF_H = new int* [rozmiar_H];

    d_zA_Hamming = nowa_Tablica2(d_zA, rozmiar_H, parita_Bitow);
    
    cout << endl << "Test zapisu do nowej tablicy: " << endl;
    for (int i = 0; i < parita_Bitow; i++) {        //pętla po ilości pakietow po 4 bity do hamminga
        for (int j = 0; j < wybor; j++) {
            cout << d_zA_Hamming[i][j];
        }      
    }

    cout << endl << endl << "Oryginalny sygnal 'ALA MA KOTA': " << endl;
    for (int i = 0, j = 1; i <= rozmiar; i++, j++) { cout << Sygnal[i]; }

    cout << endl << "Dekodowanie kodu Hamminga (7,4) dla zA: " << endl;
    for (int i = 0; i < parita_Bitow; i++) {
        dekodowanie_zA_H[i] = dekoder_Hamminga(d_zA_Hamming[i]);
        for (int j = 0; j < 4; j++) {
            cout << dekodowanie_zA_H[i][j];
        }
    }
    cout << endl;
    //===========Dane===========//
        return 0;
}


