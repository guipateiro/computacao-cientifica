// GRR20197152 Guilherme Costa Pateiro
// Universidade Federal do Parana , Introducao a computacao cientifica
// ultima edicao: 11/06/2021 02:20
#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#include "SistemasLineares.h"
#include "utils.h"
#include "Myfunc.h"


real_t Vdiff(real_t *x , real_t *xnovo,unsigned int n){
	real_t result[n];
	for (int i = 0; i < n; ++i) {
		result[i] = fabs(x[i] - xnovo[i]);
	}
	int maior = 0;       //o maior número contido no vetor
	for (int i = 0; i < n; ++i){
		if(result[i] > result[maior])
			maior = i;
	}
	return result[maior];
}


void copiaVetor(real_t *x , real_t *xnovo,unsigned int n){
	for (int i = 0; i < n; ++i){
		x[i] = xnovo[i];
	}
}

void copiaMatriz(real_t **x,real_t **xnovo,int nx,int ny){
	for (int i = 0; i < nx; ++i)
		for (int j = 0; j < ny; ++j)
			xnovo[i][j] = x[i][j];
}


real_t* alocaVetor(unsigned int n){
	real_t *novovetor;
	novovetor = (real_t*) calloc (n , sizeof (real_t)) ;
	return novovetor;
}

real_t Mdiff(real_t **xnovo, real_t **x,int nx,int ny){
	real_t result[nx][ny];
	int i,j;
	for (i = 0; i < nx; ++i)
		for (j = 0; j < ny; ++j)
			result[i][j] = fabs(x[i][j] - xnovo[i][j]);	

	int n = 0;
	int m = 0;       //o maior número contido no vetor
	for (i = 0; i < nx; ++i)
		for (j = 0; j < ny; ++j)
			if(result[i][j] > result[n][m]){
				n = i;
				m = j;
			}	
	return result[n][m];		
}

void liberaMatriz(real_t **y ,int n){
 	for (int i = 0 ; i < n ; ++i)
    	free(y[i]);
  	free(y);
} 

void prnVetor (real_t *v, unsigned int n){
  printf("    -->X: ");
  for (int i = 0; i < n; ++i){
    printf("X%u: %0.6lf ",i+1,v[i]);
  }
 printf("\n"); 
}

void prnMatriz (real_t **v,int n ,int m){
  for (int i = 0; i < n; ++i){
    printf("    -->X%i: ",i+1);
    for (int j = 0; j < m; ++j){
      printf("Y%u: %0.6lf ",j+1,v[i][j]);
    }
  printf("\n");  
  }
}

void converteMatrizVetor(real_t **matriz, real_t *vetor,int n, int m){
    for (int j = 0; j < m; ++j){
      for (int i = 0; i < n; ++i)
        vetor[i+n*j] = matriz[i][j];
    }  
}

real_t** alocaMatriz(int m, int n){
  real_t **novo;
    novo = malloc (n * sizeof (real_t*)) ;
  for (int i=0; i < n; i++)
    novo[i] = (real_t*) malloc (m * sizeof (real_t)) ;
  return novo;
}
