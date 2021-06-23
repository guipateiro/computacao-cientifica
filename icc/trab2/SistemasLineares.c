#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#include "utils.h"
#include "SistemasLineares.h"
#include "Myfunc.h"

/*!
  \brief Essa função calcula a norma L2 do resíduo de um sistema linear 

  \param SL Ponteiro para o sistema linear
  \param x Solução do sistema linear
  \param res Valor do resíduo

  \return Norma L2 do resíduo.
*/
real_t normaL2Residuo(SistLinear_t *SL, real_t *x, real_t *res){
  int i,j; 
  for (i = 0; i < SL->n; ++i){
    res[i] = 0;
    for (j = 0; j < SL->n; ++j)
      res[i] += (SL->A[i][j] * x[j]);    
    res[i] -= SL->b[i]; 
  }  
  real_t norma;
  for (i = 0; i < SL->n; ++i){
    norma += (res[i]*res[i]);
  }
  norma = sqrt(norma);
  return norma;
}


/*!
  \brief Método da Eliminação de Gauss

  \param SL Ponteiro para o sistema linear
  \param x ponteiro para o vetor solução
  \param tTriangulariza tempo gasto na triangularização
  \param tRetroSubst tempo gasto na retrosubstituição

  \return código de erro. 0 em caso de sucesso.
*/

int eliminacaoGauss (SistLinear_t *SL, real_t *x, double *tTotal){
  tTotal[0] = timestamp();
  real_t aux;
  int i,j,k;
  double m;
  unsigned int maior;

  //pivoteamento(SL->A,SL->b,SL->n);
  for (i = 0; i < SL->n ; ++i){
    maior = i; 
    for (k = i+1; k < SL->n; ++k)
      if (fabs(SL->A[k][i]) > fabs(SL->A[maior][i]))
        maior = k;
    if (maior != i ){
      for (k = i; k < SL->n; ++k){
        aux = SL->A[i][k];
        SL->A[i][k] = SL->A[maior][k];
        SL->A[maior][k] = aux;
      }
      aux = SL->b[i];
      SL->b[i] = SL->b[maior];
      SL->b[maior] = aux;
    }
    for (k = i+1; k < SL->n; ++k){
      m = SL->A[k][i] / SL->A[i][i];
      SL->A[k][i] = 0.0;
      for (j = i+1; j < SL->n; ++j)
        SL->A[k][j] -=  SL->A[i][j] * m;
      SL->b[k] -= SL->b[i] * m;
    }
  }  
  retrossubs(SL->A, SL->b,x,SL->n);

  for (i = 0; i < SL->n; ++i){
    if(x[i] ==  INFINITY)
      return -1;
    if (__isnanf(x[i])){
      return -2;
    }
  }
  tTotal[0] = timestamp() - tTotal[0];
  return 0; 
}

/*!
  \brief Método de Gauss-Jacobi

  \param SL Ponteiro para o sistema linear
  \param x ponteiro para o vetor solução. Ao iniciar função contém
            valor inicial
  \param erro menor erro aproximado para encerrar as iterações
  \param tTotal tempo gasto pelo método
  \param tIteração tempo gasto em cada iteração

  \return código de erro. Um nr positivo indica sucesso e o nr
          de iterações realizadas. Um nr. negativo indica um erro:
          -1 (não converge) -2 (sem solução)
*/
int gaussJacobi (SistLinear_t *SL, real_t *x, double *tTotal){
  tTotal[0] = timestamp();
  real_t xnovo[SL->n];
  real_t soma;
  int k = 0;
  int i,j;
  
  limpaVetor(x,SL->n);
  for (i = 0; i < SL->n; ++i) {
    soma = 0;
    for (j = 0; j < SL->n; ++j) {
      if (j != i) {
        soma += abs(SL->A[i][j]);
      }
    }
    if(soma >= abs(SL->A[i][i])){
      tTotal[0] = timestamp() - tTotal[0]; 
      return(-1);
    }
  }
  real_t erro =INFINITY; 
  while (erro >= SL->erro && k < MAXIT) {
    for (i = 0; i < SL->n; ++i) {
      soma = 0;
      for (j = 0; j < SL->n; ++j) {
        if (j != i) {
          soma += SL->A[i][j] * x[j] / SL->A[i][i];
        }
      xnovo[i] = (SL->b[i] / SL->A[i][i]) - soma;
      }
    }
    erro = Vdiff(x,xnovo,SL->n);
    copiaVetor(x,xnovo,SL->n);
    k++;
  }  
  for(i = 0; i < SL->n; ++i){
     if (__isnanf(x[i])){
      tTotal[0] = timestamp() - tTotal[0];
      return -2;
    }
    if(x[i] == INFINITY || x[i] == -INFINITY){
       tTotal[0] = timestamp() - tTotal[0];
      return -2;
    }
  }
  tTotal[0] = timestamp() - tTotal[0];  
  return k;
}

/*!
  \brief Método de Gauss-Seidel

  \param SL Ponteiro para o sistema linear
  \param x ponteiro para o vetor solução. Ao iniciar função contém
            valor inicial
  \param erro menor erro aproximado para encerrar as iterações
  \param tTotal tempo gasto pelo método
  \param tIteração tempo gasto em cada iteração

  \return código de erro. Um nr positivo indica sucesso e o nr
          de iterações realizadas. Um nr. negativo indica um erro:
          -1 (não converge) -2 (sem solução)
  */
int gaussSeidel (SistLinear_t *SL, real_t *x, double *tTotal){
  tTotal[0] = timestamp();
  real_t xnovo[SL->n];
  real_t soma;
  int i,j,l;
  int k = 0;
  limpaVetor(xnovo,SL->n);
  limpaVetor(x,SL->n);
  pivoteamento(SL->A,SL->b,SL->n);
  for (i = 0; i < SL->n; ++i) {
    soma = 0;
    for (j = 0; j < SL->n; ++j) {
      if (j != i) {
        soma += abs(SL->A[i][j]);
      }
    }
    if(soma >= abs(SL->A[i][i])){
      tTotal[0] = timestamp() - tTotal[0]; 
      return(-1);
    }
  }
  real_t erro = INFINITY; 
  while (erro > SL->erro && k < MAXIT){
    for (i = 0; i < SL->n; ++i) {
      soma = 0;
      for (j = 0; j < i; ++j) {
        soma += SL->A[i][j] * xnovo[j];
      }
      for (l = i + 1; l < SL->n; l++) {
        soma += SL->A[i][l] * x[l];
      }
      xnovo[i] = (SL->b[i] - soma) / SL->A[i][i];
    }  
    erro = Vdiff(x,xnovo,SL->n);
    copiaVetor(x,xnovo,SL->n);
    k++;
  }
  for(i = 0; i < SL->n; ++i){
     if (__isnanf(x[i])){
      tTotal[0] = timestamp() - tTotal[0];
      return -2;
    }
    if(x[i] == INFINITY || x[i] == -INFINITY){
       tTotal[0] = timestamp() - tTotal[0];
      return -3;
    }
  }
  tTotal[0] = timestamp() - tTotal[0];  
  return k;
}


/*!
  \brief Método de Refinamento

  \param SL Ponteiro para o sistema linear
  \param x ponteiro para o vetor solução. Ao iniciar função contém
            valor inicial para início do refinamento
  \param erro menor erro aproximado para encerrar as iterações

  \return código de erro. Um nr positivo indica sucesso e o nr
          de iterações realizadas. Um nr. negativo indica um erro:
          -1 (não converge) -2 (sem solução)
  */
int refinamento (SistLinear_t *SL, real_t *x, double *tTotal){
  tTotal[0] = timestamp();
  real_t res[SL->n];
  real_t w[SL->n];
  int cont = 0;
  real_t aux;
  int i,j,k;
  double m;
  unsigned int maior;

  limpaVetor(w,SL->n);
  limpaVetor(res,SL->n);
  real_t norma  = 10;//= normaL2Residuo(SL, x, res);

  if (norma < 5.0){
    tTotal[0] = timestamp() - tTotal[0];
    return 0;
  }

  while(cont < MAXIT){
    real_t norma = normaL2Residuo(SL, x, res);

    if (norma < 5.0){
      tTotal[0] = timestamp() - tTotal[0];
      return cont;
    }
    if (norma == INFINITY){
      
      tTotal[0] = timestamp() - tTotal[0];
      return -1;
    }

    for(i = 0; i < SL->n ; ++i){
      for (k = i+1; k < SL->n; ++k){
        m = SL->A[k][i] / SL->A[i][i];
        SL->A[k][i] = 0.0;
        for (j = i+1; j < SL->n; ++j)
          SL->A[k][j] -=  SL->A[i][j] * m;
        res[k] -= res[i] * m;
      }
    }
    w[SL->n-1] = res[SL->n-1]/SL->A[SL->n-1][SL->n-1];
    for (i = SL->n-2; i >= 0; --i) {
      w[i] = 0;
      for (j = i+1; j < SL->n; ++j)
        w[i] += SL->A[i][j] * w[j];
      w[i] = (res[i] - w[i]) /SL->A[i][i];
    }

    for (i = 0 ; i < SL->n ; ++i){
      x[i] = x[i] + w[i];
    }
    cont++;
  }
  tTotal[0] = timestamp() - tTotal[0];
  return cont;
}

/*!
  \brief Alocaçao de memória 

  \param n tamanho do SL

  \return ponteiro para SL. NULL se houve erro de alocação
  */
SistLinear_t* alocaSistLinear (unsigned int n){
  SistLinear_t *novo;
  novo = (SistLinear_t*) malloc (sizeof(SistLinear_t));
  novo->A = malloc (n * sizeof (real_t*)) ;
  for (int i=0; i < n; i++)
    novo->A[i] = (real_t*) malloc (n * sizeof (real_t)) ;
  novo->b = (real_t*) malloc (n * sizeof (real_t)) ;
  return(novo);  
}

/*!
  \brief Liberaçao de memória 

  \param sistema linear SL
  */
void liberaSistLinear (SistLinear_t *SL){
  free(SL->b);
  for (int i = 0 ; i < SL->n ; i++)
      free(SL->A[i]);
  free(SL->A);
  free(SL);
  return;
}

/*!
  \brief Leitura de SL a partir de Entrada padrão (stdin).

  \return sistema linear SL. NULL se houve erro (leitura ou alocação)
  */
SistLinear_t *lerSistLinear (){
  SistLinear_t *novo;
  int i,j;
  unsigned int n;
  scanf("%u",&n);
  if(n == 0 || n == EOF || n > 100)
    return NULL;
  novo = alocaSistLinear(n);
  if (novo == NULL){
    fprintf(stderr,"erro ao alocar a matriz\n");
    return NULL;
  }
  novo->n = n;
  scanf("%f",&novo->erro);
  for ( i = 0; i < novo->n; i++)
    for (j = 0; j < novo->n; j++)
      scanf("%f ",&novo->A[i][j]);
  for(i = 0; i < novo->n ; i++)
    scanf("%f ",&novo->b[i]);  
  return(novo);
}


// Exibe SL na saída padrão
void prnSistLinear (SistLinear_t *SL){
  for (int i = 0; i < SL->n; i++){
    for(int j = 0; j < SL->n; j++)
      printf("%0.2e ",SL->A[i][j]);
    printf(" || %0.2e \n",SL->b[i]);
  }
  printf("\n"); 
}

// Exibe um vetor na saída padrão
void prnVetor (real_t *v, unsigned int n){
  printf("    -->X: ");
  for (int i = 0; i < n; ++i){
    printf("%0.6f ",v[i]);
  }
 printf("\n"); 
}

