// GRR20197152 Guilherme Costa Pateiro
// Universidade Federal do Parana , Introducao a computacao cientifica
// ultima edicao: 11/06/2021 02:20
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "utils.h"
#include "SistemasLineares.h"
#include "Myfunc.h"

#define L 6 
#define W 8   

//-----------------------------------------------//
//              elementos da EDO A 
// ----------------------------------------------//
double Ap(double x){
	return 0;
}

double Aq(double x){
	return 0;
}

double Ar(double x){
	double w;

	w = 6*x - 0.5*(x*x);
	return w;
}

//-----------------------------------------------//
//              elementos da MFP B
// ----------------------------------------------//
double Bq(double x,double y){
	return -1;
}

double Br(double x,double y){
	double w;
	w = sin(x)*sin(x); 
	return w;
}
double xaB(double x,double y){
	return 20;
}

double xbB(double x,double y){
	return 45;
}

double yaB(double x,double y){
	return 0;
}
double ybB(double x,double y){
	return 100;
}


//-----------------------------------------------//
//              elementos da EDO C
// ----------------------------------------------//
double Cp(double x){
	return 0;
}

double Cq(double x){
	return 1;
}

double Cr(double x){
	return 0;
}
//-----------------------------------------------//
//              elementos da MFP D
// ----------------------------------------------//
double Dq(double x,double y){
	return 0;
}

double Dr(double x,double y){
	double w;
	w = -cos(x+y) -cos(x-y); 
	return w;
}
double xaD(double x,double y){
	return cos(y);
}

double xbD(double x,double y){
	return -cos(y);
}

double yaD(double x,double y){
	return cos(x);
}
double ybD(double x,double y){
	return 0;
}


int main (){
	double	tempo;
	int cc;
	real_t norma;
	SistLinear_t* SL;
	real_t *y;
	real_t **Y;
	//--------------------------------------------//
	//
	//					Funcao A 
	//				   Tamanho 5 
	//--------------------------------------------//
	Edo* edoA = alocaEdo(5,0,12,0,0);          	//aloca memoria para a Edo
	if (edoA == NULL)  return -1;
	edoA->p = &Ap;								//insere as funcoes 
	edoA->q = &Aq;
	edoA->r = &Ar;
	y = alocaVetor(edoA->n);					//aloca o vetor resposta
	if (y == NULL)  return -1;
	SL = alocaSistLinear(edoA->n,1);			//aloca o sistema linar
	if (SL == NULL)  return -1;
	cc = gaussSeidelEDO (edoA, SL, y, &tempo);	//faz o gauss-seidel
	printf("**** letra A  tamanho 5    -->%i interacoes\n",cc);  //imrpime resultados
	printf("    -->tempo = %0.6f ms\n\n",tempo);
	prnSL(SL);
	prnVetor (y,edoA->n);
	norma = normaL2Residuo(SL,y,edoA->n);		//faz a norma L2 do residuo
  	printf("    -->Norma L2 do residuo: %0.6e\n\n",norma);
  	liberaSistLinear(SL);						//libera memoria
	free(y);
	//-----------------------------------------------//
	//                  Tamanho 10 
	// ----------------------------------------------//
	edoA->n = 10;								//modifica o tamanho 
	y = alocaVetor(edoA->n);
	if (y == NULL)  return -1;
	SL = alocaSistLinear(edoA->n,1);
	if (SL == NULL)  return -1;
	cc = gaussSeidelEDO (edoA, SL, y, &tempo);
	printf("**** letra A  tamanho 10    -->%i interacoes\n",cc);
	printf("    -->tempo = %0.6f ms\n\n",tempo);
	prnSL(SL);
	prnVetor (y,edoA->n);
	norma = normaL2Residuo(SL,y,edoA->n);
  	printf("    -->Norma L2 do residuo: %0.6e\n\n",norma);
  	liberaSistLinear(SL);
	free(y);

	//--------------------------------------------//
	//
	//					Funcao B 
	//				 Tamanho 5 x 3 
	//--------------------------------------------//
	Mdf* MdfB = alocaMdf(5,3,0,L,0,W);			//aloca memoria para a MDF
	if (MdfB == NULL)  return -1;
	MdfB->xa = &xaB;							//completa as funcoes 
	MdfB->xb = &xbB;	
	MdfB->ya = &yaB; 
	MdfB->yb = &ybB;
	MdfB->q = &Bq;
	MdfB->r = &Br;
	SL = alocaSistLinear(MdfB->nx,3);			//aloca o sitema linear
	if (SL == NULL)  return -1;
	Y = alocaMatriz(MdfB->ny,MdfB->nx);			//aloca a matriz de resposta
	if (Y == NULL)  return -1;
	cc = gaussSeidelMDF (MdfB, SL, Y, &tempo);	// faz o gauss-seidel 
	printf("**** letra B  tamanho 5    -->%i interacoes\n",cc); //imprime resultados
	printf("    -->tempo = %0.6f ms\n\n",tempo);
	prnSL(SL);
	prnMatriz(Y,MdfB->nx,MdfB->ny);
	y = alocaVetor(MdfB->ny*MdfB->nx);      	//aloca um vetor 
	if (y == NULL)  return -1;				
	converteMatrizVetor(Y,y,MdfB->nx,MdfB->ny); //para a matriz resposta para o vetor 
	norma = normaL2Residuo(SL, y,MdfB->nx);		//faz a norma L2 do resuduo com o vetor 
	printf("    -->Norma L2 do residuo: %0.6e \n\n",norma);
	liberaMatriz(Y,MdfB->nx);					//libera memoria 
	liberaSistLinear(SL);
	free(y);
	//-----------------------------------------------//
	//               Tamanho 10 x 3 
	// ----------------------------------------------//
	MdfB->nx = 10;								//modifica o tamanho 
	SL = alocaSistLinear(MdfB->nx,3);
	if (SL == NULL)  return -1;
	Y = alocaMatriz(MdfB->ny,MdfB->nx);
	if (Y == NULL)  return -1;
	cc = gaussSeidelMDF (MdfB, SL, Y, &tempo);
	printf("**** letra B  tamanho 10    -->%i interacoes\n",cc);
	printf("    -->tempo = %0.6f ms\n\n",tempo);
	prnSL(SL);
	prnMatriz(Y,MdfB->nx,MdfB->ny);
	y = alocaVetor(MdfB->ny*MdfB->nx);
	if (y == NULL)  return -1;
	converteMatrizVetor(Y,y,MdfB->nx,MdfB->ny);
	norma = normaL2Residuo(SL, y,MdfB->nx);
	liberaMatriz(Y,MdfB->nx);
	printf("    -->Norma L2 do residuo: %0.6e \n\n",norma);
	liberaSistLinear(SL);
	free(y);

	//--------------------------------------------//
	//
	//					Funcao C 
	//				   Tamanho 5 
	//--------------------------------------------//
	Edo* edoC = alocaEdo(5,0,1,0,1);
	if (edoC == NULL)  return -1;
	edoC->p = &Cp;
	edoC->q = &Cq;
	edoC->r = &Cr;
	y = alocaVetor(edoC->n);
	if (y == NULL)  return -1;
	SL = alocaSistLinear(edoC->n,1);
	if (SL == NULL)  return -1;
	cc = gaussSeidelEDO (edoC, SL, y, &tempo);
	printf("**** letra C  tamanho 5    -->%i interacoes\n",cc);
	printf("    -->tempo = %0.6f ms\n\n",tempo);
	prnSL(SL);
	prnVetor (y,edoC->n);
	norma = normaL2Residuo(SL,y,edoC->n);
  	printf("    -->Norma L2 do residuo: %0.6e\n\n",norma);
  	liberaSistLinear(SL);
	free(y);
	//-----------------------------------------------//
	//                  Tamanho 10 
	// ----------------------------------------------//
	edoC->n = 10;
	y = alocaVetor(edoC->n);
	if (y == NULL)  return -1;
	SL = alocaSistLinear(edoC->n,1);
	if (SL == NULL)  return -1;
	cc = gaussSeidelEDO (edoC, SL, y, &tempo);
	printf("**** letra C  tamanho 10    -->%i interacoes\n",cc);
	printf("    -->tempo = %0.6f ms\n\n",tempo);
	prnSL(SL);
	prnVetor (y,edoC->n);
	norma = normaL2Residuo(SL,y,edoC->n);
  	printf("    -->Norma L2 do residuo: %0.6e\n\n",norma);
  	liberaSistLinear(SL);
	free(y);

	//--------------------------------------------//
	//
	//					Funcao D 
	//				 Tamanho 5 x 3 
	//--------------------------------------------//
	Mdf* MdfD = alocaMdf(5,3,0,M_PI,0,M_PI/2);
	if (MdfD == NULL)  return -1;
	MdfD->xa = &xaD;
	MdfD->xb = &xbD;	
	MdfD->ya = &yaD; 
	MdfD->yb = &ybD;
	MdfD->q = &Dq;
	MdfD->r = &Dr;
	SL = alocaSistLinear(MdfD->nx,3);
	if (SL == NULL)  return -1;
	Y = alocaMatriz(MdfD->ny,MdfD->nx);
	if (Y == NULL)  return -1;
	cc = gaussSeidelMDF (MdfD, SL, Y, &tempo);
	printf("**** letra D  tamanho 5    -->%i interacoes\n",cc);
	printf("    -->tempo = %0.6f ms\n\n",tempo);
	prnSL(SL);
	prnMatriz(Y,MdfD->nx,MdfD->ny);
	y = alocaVetor(MdfD->ny*MdfD->nx);
	if (y == NULL)  return -1;
	converteMatrizVetor(Y,y,MdfD->nx,MdfD->ny);
	norma = normaL2Residuo(SL, y,MdfD->nx);
	liberaMatriz(Y,MdfD->nx);
	printf("    -->Norma L2 do residuo: %0.6e \n\n",norma);
	liberaSistLinear(SL);
	free(y);
	//-----------------------------------------------//
	//               Tamanho 10 x 3 
	// ----------------------------------------------//
	MdfD->nx = 10;
	SL = alocaSistLinear(MdfD->nx,3);
	if (SL == NULL)  return -1;
	Y = alocaMatriz(MdfD->ny,MdfD->nx);
	if (Y == NULL)  return -1;
	cc = gaussSeidelMDF (MdfD, SL, Y, &tempo);
	printf("**** letra D  tamanho 10    -->%i interacoes\n",cc);
	printf("    -->tempo = %0.6f ms\n\n",tempo);
	prnSL(SL);
	prnMatriz(Y,MdfD->nx,MdfD->ny);
	y = alocaVetor(MdfD->ny*MdfD->nx);
	if (y == NULL)  return -1;
	converteMatrizVetor(Y,y,MdfD->nx,MdfD->ny);
	norma = normaL2Residuo(SL, y,MdfD->nx);
	liberaMatriz(Y,MdfD->ny);
	printf("    -->Norma L2 do residuo: %0.6e \n\n",norma);
	liberaSistLinear(SL);
	free(y);

	//-----------------------------------------------//
	//               libera memoria 
	// ----------------------------------------------//
	free(edoA);
	free(edoC);
	free(MdfB);
	free(MdfD);
	return 1;
}	
