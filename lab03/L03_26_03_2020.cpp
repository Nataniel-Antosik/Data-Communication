#define _USE_MATH_DEFINES
#include<iostream>
#include<math.h>
#include "math.h"

//44261 A = 1, B = 6, C=2

using namespace std;

struct zespolona {
	double RE, IM;
};

zespolona* DFT(double* tab, int N) { //N to rozmiar
	zespolona* tmp = new zespolona[N];
	for (int i = 0; i < N; i++) {
		for (int k = 0; k < N; k++) {
			tmp[i].RE += tab[k] * cos((k * i * 2 * M_PI) / N);
			tmp[i].IM -= tab[k] * sin((k * i * 2 * M_PI) / N);
		}
	}
	return tmp;
}

int main()
{
	zespolona* dane;
	double tab[] = { 1,2,3,4,5,6,7,8,9,10 };
	int rozmiar = 10;
	dane = DFT(tab, rozmiar);
	int k = 1;
	for (int i = 0; i < rozmiar; i++) {
		cout << "WARTOSC NR: " << k++ << ": " << "(" << dane[i].RE << ")" << " + " << "i(" << dane[i].IM << ")" << endl;
	}

	return 0;
}

