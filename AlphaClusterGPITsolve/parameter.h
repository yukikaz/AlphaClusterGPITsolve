#pragma once

const double N = 3;
const int NX = 1000;
const double sigma = 1.0;
const double m = 3733;
const double omega = 2.58;
const double alpha = sqrt(omega*m);
double mu = 1;
const double mur = 0.7 / sqrt(m*omega) * 197;
const double mua = 0.475 / sqrt(m*omega) * 197;
const double Vr = 400 / omega;
const double Va = -130 / omega;
const double eps = 2.3e-9;
const double epsgz = 1e-9;
const int looplim = 50000;

const double R = 8;
const double DR = R / (NX + 1);
const double DT = DR*DR;

inline double i2r(int i) { return (i + 1)*DR; }

inline double gauss(double x) {
	return 2 / sqrt(sqrt(M_PI)*sigma*sigma*sigma)*exp(-x*x / (2 * sigma*sigma));
}

inline double hR(double x, double* p) {
	double temp = 0;
#pragma omp parallel for reduction(+:temp)
	for (int i = 0; i < NX; i++) {
		double r = i2r(i);
		temp += r / x*p[i] * p[i] * DR* (Vr / (2 * mur*mur)*(exp(-mur*mur*(x - r)*(x - r)) - exp(-mur*mur*(x + r)*(x + r))) + Va / (2 * mua*mua)*(exp(-mua*mua*(x - r)*(x - r)) - exp(-mua*mua*(x + r)*(x + r))));
	}
	return mu - 0.5*x*x - N*0.5 * temp;
}

inline double nrm(double* p) {
	double temp = 0;
#pragma omp parallel for reduction(+:temp)
	for (int i = 0; i < NX; i++) {
		double r = i2r(i);
		temp += pow(p[i], 2)*DR * r*r;
	}
	return temp;
}

inline double error(double* p, double* p2) {
	double temp = 0;
#pragma omp parallel for reduction(+:temp)
	for (int i = 0; i < NX; i++) {
		double r = i2r(i);
		temp += pow(p[i] - p2[i], 2)*DR * r*r;
	}
	return temp;
}

double* normaliz(double* p, double nrm) {
#pragma omp parallel for
	for (int i = 0; i < NX; i++) {
		p[i] /= sqrt(nrm);
	}
	return p;
}

void gsmeth1(double* p, double* p_n, double* p_t) {
	for (int i = 0; i < NX; i++) {
		double r = i2r(i);
		if (i == 0) { p_n[0] = ((1 + 1 / (1 + i))*(p[1 + i] + p_t[i + 1]) - 2 * (1 - DR*DR*hR(r, p) - 2 * DR*DR / DT)*p_t[i]) / (2 * (1 - DR*DR*hR(r, p) + 2 * DR*DR / DT)); }
		else if (i == NX - 1) { p_n[i] = ((1 - 1 / (i + 1))*(p_n[i - 1] + p_t[i - 1]) - 2 * (1 - DR*DR*hR(r, p) - 2 * DR*DR / DT)*p_t[i]) / (2 * (1 - DR*DR*hR(r, p) + 2 * DR*DR / DT)); }
		else p_n[i] = ((1 + 1 / (1 + i))*(p[1 + i] + p_t[i + 1]) + (1 - 1 / (i + 1))*(p_n[i - 1] + p_t[i - 1]) - 2 * (1 - DR*DR*hR(r, p) - 2 * DR*DR / DT)*p_t[i]) / (2 * (1 - DR*DR*hR(r, p) + 2 * DR*DR / DT));
	}
	return;
}

void gsmeth2(double* p, double* p_n, double* p_t) {
	for (int i = 0; i < NX; i++) {
		double r = i2r(i);
		if (i == 0) { p_n[0] = ((1 + 1 / (1 + i))*(p[1 + i] + p_t[i + 1]) - 2 * (1 - DR*DR*hR(r, p_t) - 2 * DR*DR / DT)*p_t[i]) / (2 * (1 - DR*DR*hR(r, p_t) + 2 * DR*DR / DT)); }
		else if (i == NX - 1) { p_n[i] = ((1 - 1 / (i + 1))*(p_n[i - 1] + p_t[i - 1]) - 2 * (1 - DR*DR*hR(r, p_t) - 2 * DR*DR / DT)*p_t[i]) / (2 * (1 - DR*DR*hR(r, p_t) + 2 * DR*DR / DT)); }
		else p_n[i] = ((1 + 1 / (1 + i))*(p[1 + i] + p_t[i + 1]) + (1 - 1 / (i + 1))*(p_n[i - 1] + p_t[i - 1]) - 2 * (1 - DR*DR*hR(r, p_t) - 2 * DR*DR / DT)*p_t[i]) / (2 * (1 - DR*DR*hR(r, p_t) + 2 * DR*DR / DT));
	}
	return;
}

void save(double* psave, double* pout) {
#pragma omp parallel for
	for (int i = 0; i < NX; i++) {
		pout[i] = psave[i];
	}
	return;
}