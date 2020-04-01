#define _USE_MATH_DEFINES
#include<iostream>
#include<math.h>

//44261 A = 1, B = 6, C=2

using namespace std;

struct zespolona {
	double RE, IM;
};

zespolona* DFT(double tab[], int N, zespolona dane[]) { //N to rozmiar
	for (int i = 0; i < N; i++) {
		for (int k = 0; k < N; k++) {
			dane[k].RE += tab[i] * cos(k * i * (M_PI / (N + 1)));
			dane[k].IM += tab[i] * sin(k * i * (M_PI / (N + 1)));
		}
	}
	return dane;
}

int main()
{
	zespolona dane[10];
	double tab[] = { 1,2,3,4,5,6,7,8,9,10 };
	int rozmiar = 10;
	DFT(tab, rozmiar, dane);
	for (int i = 0; i < rozmiar; i++) {
		cout << "WARTOSC NR: " << i + 1 << ": " << "(" << dane[i].RE << ")" << " + " << "i(" << dane[i].IM << ")" << endl;
	}
	
	return 0;
}

