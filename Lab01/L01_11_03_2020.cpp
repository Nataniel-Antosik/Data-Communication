#define _USE_MATH_DEFINES
#include <iostream>
#include <math.h>
#include <fstream>

//44261 A = 1, B = 6, C=2
using namespace std;
double a, b, c, x1, x2, delta, t_p = (-10), t_k = 10 , deltaT = 0.01, suma = 0;

void wyroznik_m_zerowe(){
	if (a == 0){
		cout << "To nie jest rownanie kwadratowe";
	}
	else{
		delta = (b * b) - (4 * a * c);
		if (delta == 0){
			x1 = -b / (2 * a);
			cout << "Miejsce zerowe wynosi: " << x1;
		}
		else{
			if (delta > 0){
				x1 = (-b - sqrt(delta)) / (2 * a);
				x2 = (-b + sqrt(delta)) / (2 * a);
				cout << "x1= " << x1 << endl << "x2= " << x2;
			}
			else
				cout << "Funkcja nie ma miejsc zerowych";
		}
	}
}

double x(double t) { return double(a * t * t + b * t + c); }
double y(double t) { return double(2 * pow(x(t), 2) + 12 * cos(t)); }
double z(double t) { return double(sin(2 * M_PI * 7 * t) * x(t) - 0.2 * log10(abs(y(t)) + M_PI)); }
double u(double t) { return double(sqrt(abs(y(t) * y(t) * z(t))) - 1.8 * sin(0.4 * t * z(t) * x(t))); }
double v(double t) {
	
	if (t < 0.22 && t >= 0) {
		return double((1 - 7 * t) * sin((2 * M_PI * t * 10) / (t + 0.04)));
	}
	if (t >= 0.22 && t < 0.7) {
		return double (0.63 * t * sin(125 * t));
	}
	if (t <= 1.0 && t >= 0.7) {
		return double (pow(t, (-0.662)) + 0.77 * sin(8 * t));
	}
}

double p(double t,int wybor) {
	//N należy do zbioru (2,4,16)
	switch (wybor) {
	case 1:
		for (int i = 1; i <= 2; i++) {
			suma = suma + double(cos(12 * t * pow(i, 2)) + cos(16 * t * i)) / (pow(i, 2));
		}
		break;
	case 2:
		for (int i = 1; i <= 4; i++) {
			suma = suma + double(cos(12 * t * pow(i, 2)) + cos(16 * t * i)) / (pow(i, 2));
		}
		break;
	case 3:
		for (int i = 1; i <= 16; i++) {
			suma = suma + double(cos(12 * t * pow(i, 2)) + cos(16 * t * i)) / (pow(i, 2));
		}
		break;
	}
	return suma;
}

int main(){
	a = 1, b = 6, c = 2;
	wyroznik_m_zerowe();
	//wykres x
	ofstream x_t("x_t.dat"); //tworzenie pliku
	double t1 = t_p; //zmienna pomocnicza dla wykresu i jest ona na początku
	while (t1 <= t_k){ //pętla do końca 
		x_t << t1 << " " << x(t1) << endl; //zapisujemy początek wykresu i naszego x o początku(wartość)
		t1 = t1 + deltaT; //dodajemy kolejny krok do początku
	}
	x_t.close();
	//wykres y
	ofstream y_t("y_t.dat");
	double t_p_2 = 0;
	double t_k_2 = 1;
	double deltaT_2 = (double)1 / (double)22050;
	while (t_p_2 <= t_k_2){
		y_t << t_p_2 << " " << y(t_p_2) << endl;
		t_p_2 = t_p_2 + deltaT_2;
	}
	y_t.close();

	ofstream z_t("z_t.dat");
	t_p_2 = 0;
	while (t_p_2 <= t_k_2) {
		z_t << t_p_2 << " " << z(t_p_2) << endl;
		t_p_2 = t_p_2 + deltaT_2;
	}
	z_t.close();

	ofstream u_t("u_t.dat");
	t_p_2 = 0;
	while (t_p_2 <= t_k_2) {
		u_t << t_p_2 << " " << u(t_p_2) << endl;
		t_p_2 = t_p_2 + deltaT_2;
	}
	z_t.close();
	
	ofstream v_t("v_t.dat");
	t_p_2 = 0;
	while (t_p_2 <= t_k_2) {
		v_t << t_p_2 << " " << v(t_p_2) << endl;
		t_p_2 = t_p_2 + deltaT_2; 
	}
	z_t.close();

	ofstream p_t_1("p_t_1.dat");
	t_p_2 = 0;
	while (t_p_2 <= t_k_2) {
		p_t_1 << t_p_2 << " " << p(t_p_2,1) << endl;
		t_p_2 = t_p_2 + deltaT_2; 
	}
	z_t.close();

	ofstream p_t_2("p_t_2.dat");
	t_p_2 = 0;
	while (t_p_2 <= t_k_2) {
		p_t_2 << t_p_2 << " " << p(t_p_2, 2) << endl;
		t_p_2 = t_p_2 + deltaT_2;
	}
	z_t.close();

	ofstream p_t_3("p_t_3.dat");
	t_p_2 = 0;
	while (t_p_2 <= t_k_2) {
		p_t_3 << t_p_2 << " " << p(t_p_2, 3) << endl;
		t_p_2 = t_p_2 + deltaT_2;
	}
	z_t.close();
	return 0;
}
