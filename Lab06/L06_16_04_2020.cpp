
#include<sstream>
#include <iostream>
#include <fstream>
#include <string>
#include <cstddef>
#include <bitset>
#include <vector>

typedef unsigned char byte;

using namespace std;

/*
double zA(double A1) {

}

double zF() {

}

double zP() {

}
*/

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
        str << endl; 
        wynik = str.str();                      //przesłanie buffora do wyniku
    }
    else {                                      //big endian
                                                //zapisu danych, w której najbardziej znaczący bajt, umieszczony jest jako pierwszy
        for (int i = dlugosc; i >= 0; i--) {    //przechodzimy po całym łańcuchu
            bitset<8> x(in[i]);                 //zamiana liczby numerycznej na bity
            str << x;                           //dodanie do buffora bitów
        }
        str << endl;
        wynik = str.str();
    }
    return wynik;
}


int main()
{
    char dupa[4] = "agp";

    cout << "Little endian:  " << S2BS(dupa, 0);
    cout << "Big endian:     " << S2BS(dupa, 1);
    return 0;
}

