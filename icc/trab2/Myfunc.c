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

void limpaVetor(real_t *x,unsigned int n){
	for (int i = 0; i < n; ++i){
		x[i] = 0;
	}
}

real_t* alocaVetor(unsigned int n){
	real_t *novovetor;
	novovetor = (real_t*) malloc (n * sizeof (real_t)) ;
	return novovetor;
}

SistLinear_t* copiaSL(SistLinear_t *SL){
	SistLinear_t *novo;
	int i,j;
	novo = alocaSistLinear(SL->n);
	if (novo == NULL){
    	fprintf(stderr,"erro ao alocar a matriz\n");
    	return NULL;
  	}
	novo->n = SL->n;
	novo->erro = SL->erro;
	for ( i = 0; i < novo->n; i++)
    	for (j = 0; j < novo->n; j++)
      		novo->A[i][j]= SL->A[i][j];
  	for(i = 0; i < novo->n ; i++)
    	novo->b[i] = SL->b[i];  
  	return(novo);
}

/* Seja um S.L. de ordem ‘n’
*/
void retrossubs(real_t **SL, real_t *b, real_t *x, unsigned int n) {
	for (int i = n-1; i >= 0; --i) {
		x[i] = b[i];
		for (int j = i+1; j < n; ++j)
			x[i] -= SL[i][j] * x[j];
		x[i] /= SL[i][i];
 	}
}

void pivoteamento(real_t **SL, real_t *b, unsigned int n){
	real_t aux;
	unsigned int maior;
	for (int i = 0; i < n; ++i){
    	maior = i; 
    	for (int k = i+1; k < n; ++k)
      		if (fabs(SL[k][i]) > fabs(SL[maior][i]))
        		maior = k;
    	if (maior != i ){
      		for (int k = i; k < n; ++k){
        		aux = SL[i][k];
        		SL[i][k] = SL[maior][k];
        		SL[maior][k] = aux;
      		}
      		aux = b[i];
      		b[i] = b[maior];
      		b[maior] = aux;
    	}  
	}
}
