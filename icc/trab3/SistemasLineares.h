// GRR20197152 Guilherme Costa Pateiro
// Universidade Federal do Parana , Introducao a computacao cientifica
// ultima edicao: 11/06/2021 02:20
#ifndef __SISLINEAR_H__
#define __SISLINEAR_H__

// Parâmetros para teste de convergência

#define MAXIT   50  // Número máximo de iterações em métodos iterativos

#define M_PI 3.14159265358979311600

typedef double real_t;

typedef struct {
	real_t *D, *Di2, *Di1, *Ds1, *Ds2, *B;
	int n;
}SistLinear_t;

typedef struct{
	int n;			//tamanho
	double a, b;	//intervalo
	double ya, yb;	//condicoes de contorno
	double (* p) (double), (* q) (double), (* r) (double);	
} Edo;

typedef struct{
	int nx;			//tamanho x
	int ny;			//tamanho y
	double ax, bx;	//intervalo x
	double ay, by;  //intervalo y
	double (* ya) (double,double), (* yb) (double,double);	//condicoes de contorno y 
	double (* xa) (double,double), (* xb) (double,double);	//condicoes de contorno x 
	double (* q) (double,double), (* r) (double,double);	
} Mdf;

// Alocaçao e desalocação de memória
SistLinear_t* alocaSistLinear (int n,int m);
void liberaSistLinear (SistLinear_t *SL);

Edo* alocaEdo(unsigned int n, double a, double b, double ya, double yb);
Mdf* alocaMdf(int nx, int ny, double ax, double bx, double ay, double by);
real_t** alocaMatriz(int m, int n);

// Impressão de sistemas lineares
void prnSL(SistLinear_t *SL);

// Retorna a normaL2 do resíduo. 
real_t normaL2Residuo(SistLinear_t *SL, real_t *x,int m);

// Método de Gauss-Seidel. Gera um sistema linear tri ou pentadiagonal em SL, e resultados do sistema 'Y' ou 'y' 
int gaussSeidelEDO (Edo *edoeq,SistLinear_t *SL, double *Y,double *tTotal);
int gaussSeidelMDF (Mdf *mdfeq,SistLinear_t* SL, double **Y,double *tTotal);

#endif // __SISLINEAR_H__

