#define _USE_MATH_DEFINES
#include <stdio.h>
#include <math.h>
#include "parameter.h"
#include <time.h>
#include<complex>
#pragma warning(disable : 4996)

int main() {
	clock_t start, end;
	FILE* fp;
	//FILE* fp2;
	FILE* fp3;
	char filename[20];
	double* psi = new double[NX]();
	double* psi_n = new double[NX]();
	double* psi_t = new double[NX]();
	double* psi_o = new double[NX]();
	start = clock();
	//èâä˙îg
#pragma omp parallel for
	for (int i = 0; i < NX; i++) {
		double r = i2r(i);
		psi[i] = psi_n[i] = psi_o[i] = gauss(r);
	}
	//hRämîF
	fp3 = fopen("hR.txt", "w");
	for (int i = 0; i < NX; i++) {
		double r = i2r(i);
		fprintf(fp3, "%e\t%e\n", r * 197 / alpha, hR(r, psi));
	}
	fclose(fp3);
	//éûä‘î≠ìW
	for (int k = 0; k < looplim; k++) {
		//ç≈èâÇ…ëOÇÃéûä‘ÇÃä÷êîÇï€ë∂
		save(psi_n, psi);
		save(psi_n, psi_t);
		//ÉKÉEÉXÉUÉCÉfÉãîΩïú
		for (int j = 0; j < looplim; j++) {
			gsmeth1(psi, psi_n, psi_t);
			if (fabs(error(psi, psi_n)) < epsgz) { if (k % 20 == 0) { printf("i = %d\tGSloop = %d\tmu = %e\terror %e\n", k, j, mu, fabs((1 - nrm(psi_n)) / 2 / DT)); } break; }
			else { save(psi_n, psi); }
		}
		//muÇÃä…òa
		mu += (1 - nrm(psi_n)) / 2 / DT;
		//é˚ë©ÇµÇΩÇÁó£íE&ãKäiâª
		if (fabs((1 - nrm(psi_n)) / 2 / DT) < eps) { normaliz(psi_n, nrm(psi_n) / N); break; }
		normaliz(psi_n, nrm(psi_n));
	}
	//gpèëÇ´çûÇ›
	fp = fopen("gp.txt", "w");
	for (int i = 0; i < NX; i++) {
		double r = i2r(i);
		fprintf(fp, "%e\t%e\n", r * 197 / alpha, psi_n[i] * pow(alpha, 1.5)*pow(197, -1.5));
	}
	fclose(fp);
	end = clock();
	printf("\nTotal time: %.0fmin %fsec\n", floor((double)(end - start) / CLOCKS_PER_SEC / 60), fmod((double)(end - start) / CLOCKS_PER_SEC, 60));
	delete[] psi;
	delete[] psi_n;
	delete[] psi_t;
	delete[] psi_o;
	return 0;
}