// GRR20197152 Guilherme Costa Pateiro
// Universidade Federal do Parana , Introducao a computacao cientifica
// ultima edicao: 11/06/2021 02:20
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
  \param m Tamaho do m de uma matriz m*n

  \return Norma L2 do resíduo.
*/
real_t normaL2Residuo(SistLinear_t *SL, real_t *x,int m){
  int i,j; 
  int n = SL->n/m;
  real_t res[SL->n];
  real_t norma = 0;
  for (j = 0; j < n; ++j)
    for (i = 0; i < m; ++i){
      res[i+m*j] = 0;
    
      if (i == 0 && j == 0)        
        res[i+m*j] = SL->Ds1[i+m*j]*x[i+1+m*j] + SL->Ds2[i+m*j]*x[i+m*(j+1)] + SL->D[i+m*j]*x[i+m*j]; 
      else if (i == 0 && j == n-1)
        res[i+m*j] = SL->Di2[i+m*j]*x[i+m*(j-1)] + SL->Ds1[i+m*j]*x[i+1+m*j] + SL->D[i+m*j]*x[i+m*j] ;
      else if (i == m-1 && j == 0)
        res[i+m*j] = SL->Di1[i+m*j]*x[i-1+m*j] + SL->Ds2[i+m*j]*x[i+m*(j+1)] + SL->D[i+m*j]*x[i+m*j] ;
      else if (i == m-1 && j == n-1) 
        res[i+m*j] = SL->Di2[i+m*j]*x[i+m*(j-1)] + SL->Di1[i+m*j]*x[i-1+m*j] + SL->D[i+m*j]*x[i+m*j] ;
      else if (i == 0 && j != 0)
        res[i+m*j] = SL->Ds1[i+m*j]*x[i+1+m*j] + SL->Ds2[i+m*j]*x[i+m*(j+1)] + SL->Di2[i+m*j]*x[i+m*(j-1)] + SL->D[i+m*j]*x[i+m*j];
      else if (i != 0 && j == 0)    
        res[i+m*j] = SL->Ds1[i+m*j]*x[i+1+m*j] + SL->Ds2[i+m*j]*x[i+m*(j+1)] + SL->Di1[i+m*j]*x[i-1+m*j] + SL->D[i+m*j]*x[i+m*j];
      else if (i == m-1 && j != n-1) 
        res[i+m*j] = SL->Di2[i+m*j]*x[i+m*(j-1)] + SL->Di1[i+m*j]*x[i-1+m*j] + SL->Ds2[i+m*j]*x[i+m*(j+1)] + SL->D[i+m*j]*x[i+m*j];
      else if (i != m-1 && j == n-1) 
        res[i+m*j] = SL->Di2[i+m*j]*x[i+m*(j-1)] + SL->Di1[i+m*j]*x[i-1+m*j] + SL->Ds1[i+m*j]*x[i+1+m*j] + SL->D[i+m*j]*x[i+m*j];
      else               
        res[i+m*j] = SL->Ds2[i+m*j]*x[i+m*(j+1)] + SL->Ds1[i+m*j]*x[i+1+m*j] + SL->Di1[i+m*j]*x[i-1+m*j] + SL->Di2[i+m*j]*x[i+m*(j-1)] + SL->D[i+m*j]*x[i+m*j];

      res[i+m*j] -= SL->B[i+m*j];

      norma += (res[i+m*j]*res[i+m*j]);      
    }  
  norma = sqrt(norma);
  return norma;
}

/*!
  \brief Método de Gauss-Seidel EDO

  \param edoeq Ponteiro para o sistema EDO
  \param SL Ponteiro para o sistema linear
  \param y ponteiro para o vetor solução. 
  \param tTotal tempo gasto pelo método

  \return  K numero de iteracoes realizadas 
  */
int gaussSeidelEDO (Edo *edoeq,SistLinear_t *SL, double *Y,double *tTotal){
  *tTotal = timestamp();
  int k = 0, i;
  int n = edoeq->n;
  real_t *Yantigo = alocaVetor(n);
  double h, xi, bi, d, di, ds, diff = 1;
  h = (edoeq->b - edoeq->a)/(n+1);   
  while(k < 50 && diff > 0.0001){         //comente essa linha para mudar para o metodo de seomente interacoes
  //while(k < 50 ){                       //descomente essa linha para midar apara o metodo de somente interacoes
    for (i=0; i < n; ++i) {           
      xi = edoeq->a + (i+1)*h;        
      bi = h*h * edoeq->r(xi);       
      di = 1 - h * edoeq->p(xi)/2.0;  
      d = -2 + h*h * edoeq->q(xi);   
      ds = 1 + h * edoeq->p(xi)/2.0;  

      SL->D[i] = d;

      if (i == 0){
        bi -= edoeq->ya * (1 - h*edoeq->p(edoeq->a+h)/2.0);
        SL->Di1[i] = 0;
        SL->Di2[i] = 0;
        SL->Ds1[i] = ds;
        SL->Ds2[i] = 0;
        SL->B[i] = bi;
      }
      else if (i == n-1){
        bi -= edoeq->yb * (1 + h*edoeq->p(edoeq->b-h)/2.0);
        SL->Di1[i] = di;
        SL->Di2[i] = 0;
        SL->Ds1[i] = 0;
        SL->Ds2[i] = 0;
        SL->B[i] = bi;
      }
      else{
        SL->Di1[i] = di;
        SL->Di2[i] = 0;
        SL->Ds1[i] = ds;
        SL->Ds2[i] = 0;
        SL->B[i] = bi;
      }


      if (i == 0)
        bi -= ds*Y[i+1] ;
      else if (i == n-1)
        bi -= di*Y[i-1] ;
      else
        bi -= ds*Y[i+1] + di*Y[i-1];
      Y[i] = bi / d;    
    }    
    diff = Vdiff(Y,Yantigo,n);     //comente essas duas linha se desejar utilizar 
    copiaVetor(Yantigo,Y,n);       //somente o criterio de interacoe
    k++;
  }
  *tTotal = timestamp() - *tTotal;
  return k;
}


/*!
  \brief Método de Gauss-Seidel MDF

  \param mdfeq Ponteiro para o sistema MDF
  \param SL Ponteiro para o sistema linear
  \param Y ponteiro para a matriz solução.
  \param tTotal tempo gasto pelo método

  \return  K numero de iteracoes realizadas 
  */
int gaussSeidelMDF (Mdf *mdfeq,SistLinear_t *SL, double **Y,double *tTotal){
  *tTotal = timestamp();
  int k = 0, i, j;
  real_t **Yantigo = alocaMatriz(mdfeq->ny,mdfeq->nx);
  double hx, hy, xi, bi, yi, d, di1, di2, ds1, ds2, diff = 1;
  hx = (mdfeq->bx - mdfeq->ax)/(mdfeq->nx+1);   
  hy = (mdfeq->by - mdfeq->ay)/(mdfeq->ny+1);
  while(k < 50 && diff > 0.0001){             //comente essa linha para mudar para o metodo de seomente interacoes
  //while(k < 50 ){                           //descomente essa linha para midar apara o metodo de somente interacoes
    for (j = 0; j < mdfeq->ny; ++j){
      yi = mdfeq->ay + (j+1)*hy;
      for (i = 0; i < mdfeq->nx; ++i) {           
        xi = mdfeq->ax + (i+1)*hx;        
        bi = hx*hx * hy*hy * mdfeq->r(xi,yi);  
        di2 = hx *hx;
        di1 = hy *hy;
        d = (-2 * ((hx*hx)+(hy*hy)) + mdfeq->q(xi,yi));   
        ds1 = hy *hy;
        ds2 = hx *hx;  

        SL->D[i+mdfeq->nx*j] = d;

        if (i == 0 && j == 0){        
          bi -= mdfeq->xa(mdfeq->ax+hx,yi) * (hy*hy) + mdfeq->ya(xi,mdfeq->ay+hy) * (hx*hx) ;
          SL->Di1[i+mdfeq->nx*j] = 0;
          SL->Di2[i+mdfeq->nx*j] = 0;
          SL->Ds1[i+mdfeq->nx*j] = ds1;
          SL->Ds2[i+mdfeq->nx*j] = ds2;
          SL->B[i+mdfeq->nx*j] = bi;
        }  
         else if (i == 0 && j == mdfeq->ny-1){
          bi -= mdfeq->yb(xi,mdfeq->by-hy) * (hx*hx) + mdfeq->xa(mdfeq->ax+hx,yi) * (hy*hy);
          SL->Di1[i+mdfeq->nx*j] = 0;
          SL->Di2[i+mdfeq->nx*j] = di2;
          SL->Ds1[i+mdfeq->nx*j] = ds1;
          SL->Ds2[i+mdfeq->nx*j] = 0;
          SL->B[i+mdfeq->nx*j] = bi;
        }  
        else if (i == mdfeq->nx-1 && j == 0){
          bi -= mdfeq->xb(mdfeq->bx-hx,yi) * (hy*hy) + mdfeq->ya(xi,mdfeq->ay+hy) * (hx*hx);
          SL->Di1[i+mdfeq->nx*j] = di1;
          SL->Di2[i+mdfeq->nx*j] = 0;
          SL->Ds1[i+mdfeq->nx*j] = 0;
          SL->Ds2[i+mdfeq->nx*j] = ds2;
          SL->B[i+mdfeq->nx*j] = bi;
        }

        else if (i == mdfeq->nx-1 && j == mdfeq->ny-1){
          bi -= mdfeq->yb(xi,mdfeq->by-hy) * (hx*hx) + mdfeq->xb(mdfeq->bx-hx,yi) * (hy*hy);
          SL->Di1[i+mdfeq->nx*j] = di1;
          SL->Di2[i+mdfeq->nx*j] = di2;
          SL->Ds1[i+mdfeq->nx*j] = 0;
          SL->Ds2[i+mdfeq->nx*j] = 0;
          SL->B[i+mdfeq->nx*j] = bi;
        }

        else if (i == 0 && j != 0){
          bi -= mdfeq->xa(mdfeq->ax+hx,yi) * (hy*hy);
          SL->Di1[i+mdfeq->nx*j] = 0;
          SL->Di2[i+mdfeq->nx*j] = di2;
          SL->Ds1[i+mdfeq->nx*j] = ds1;
          SL->Ds2[i+mdfeq->nx*j] = ds2;
          SL->B[i+mdfeq->nx*j] = bi;
        }
        else if (i != 0 && j == 0){    
          bi -= mdfeq->ya(xi,mdfeq->ay+hy) * (hx*hx);
          SL->Di1[i+mdfeq->nx*j] = di1;
          SL->Di2[i+mdfeq->nx*j] = 0;
          SL->Ds1[i+mdfeq->nx*j] = ds1;
          SL->Ds2[i+mdfeq->nx*j] = ds2;
          SL->B[i+mdfeq->nx*j] = bi;
        }  
        else if (i == mdfeq->nx-1 && j != mdfeq->ny-1){ 
          bi -= mdfeq->xb(mdfeq->bx-hx,yi) * (hy*hy);
          SL->Di1[i+mdfeq->nx*j] = di1;
          SL->Di2[i+mdfeq->nx*j] = di2;
          SL->Ds1[i+mdfeq->nx*j] = 0;
          SL->Ds2[i+mdfeq->nx*j] = ds2;
          SL->B[i+mdfeq->nx*j] = bi;
        }  
        else if (i != mdfeq->nx-1 && j == mdfeq->ny-1){
          bi -= mdfeq->yb(xi,mdfeq->by-hy) * (hx*hx);
          SL->Di1[i+mdfeq->nx*j] = di1;
          SL->Di2[i+mdfeq->nx*j] = di2;
          SL->Ds1[i+mdfeq->nx*j] = ds1;
          SL->Ds2[i+mdfeq->nx*j] = 0;
          SL->B[i+mdfeq->nx*j] = bi;
        }
        else{
          SL->Di1[i+mdfeq->nx*j] = di1;
          SL->Di2[i+mdfeq->nx*j] = di2;
          SL->Ds1[i+mdfeq->nx*j] = ds1;
          SL->Ds2[i+mdfeq->nx*j] = ds2;
          SL->B[i+mdfeq->nx*j] = bi;
        }

        if (i == 0 && j == 0)        
          bi -= ds1*Y[i+1][j] + ds2*Y[i][j+1]; 
        else if (i == 0 && j == mdfeq->ny-1)
          bi -= di2*Y[i][j-1] + ds1*Y[i+1][j] ;
        else if (i == mdfeq->nx-1 && j == 0)
          bi -= di1*Y[i-1][j] + ds2*Y[i][j+1] ;
        else if (i == mdfeq->nx-1 && j == mdfeq->ny-1) 
          bi -= di2*Y[i][j-1] + di1*Y[i-1][j] ;
        else if (i == 0 && j != 0)
          bi -= ds1*Y[i+1][j] + ds2*Y[i][j+1] + di2*Y[i][j-1];
        else if (i != 0 && j == 0)    
          bi -= ds1*Y[i+1][j] + ds2*Y[i][j+1] + di1*Y[i-1][j];
        else if (i == mdfeq->nx-1 && j != mdfeq->ny-1) 
          bi -= di2*Y[i][j-1] + di1*Y[i-1][j] + ds2*Y[i][j+1];
        else if (i != mdfeq->nx-1 && j == mdfeq->ny-1) 
          bi -= di2*Y[i][j-1] + di1*Y[i-1][j] + ds1*Y[i+1][j];
        else               
          bi -= ds2*Y[i][j+1] + ds1*Y[i+1][j] + di1*Y[i-1][j] + di2*Y[i][j-1]; 

        Y[i][j] = bi / d; 
        
      }
    }      
    diff = Mdiff(Y,Yantigo,mdfeq->nx,mdfeq->ny);        //comente essas duas linha se desejar utilizar 
    copiaMatriz(Y,Yantigo,mdfeq->nx,mdfeq->ny);         //somente o criterio de interacoes 
    k++;
  }
  liberaMatriz(Yantigo,mdfeq->ny);
  *tTotal = timestamp() - *tTotal;
  return k;
}

/*!
  \brief Alocaçao de memória 

  \param n e m tamanho da estrutura     

  \return ponteiro para SL. NULL se houve erro de alocação
  */
SistLinear_t* alocaSistLinear (int n, int m){
  SistLinear_t *novo;
  novo = (SistLinear_t*) malloc (sizeof(SistLinear_t));
  if (novo == NULL)
    return NULL;
  novo->D = (real_t*) malloc (n * m * sizeof (real_t));
  if (novo->D == NULL)
    return NULL;
  novo->Di2 = (real_t*) malloc (n * m * sizeof (real_t));
  if (novo->Di2 == NULL)
    return NULL;
  novo->Di1 = (real_t*) malloc (n * m * sizeof (real_t));
  if (novo->Di1 == NULL)
    return NULL;
  novo->Ds1 = (real_t*) malloc (n * m * sizeof (real_t));
  if (novo->Ds1 == NULL)
    return NULL;
  novo->Ds2 = (real_t*) malloc (n * m * sizeof (real_t));
  if (novo->Ds2 == NULL)
    return NULL;
  novo->B = (real_t*) malloc (n * m * sizeof (real_t)) ;
  if (novo->B == NULL)
    return NULL;
  novo->n = n*m;
  return(novo);  
}

Edo* alocaEdo(unsigned int n, double a, double b, double ya, double yb){
  Edo *novo;
  novo = (Edo*) malloc (sizeof (Edo));
  if (novo == NULL)
    return NULL;
  novo->n = n;
  novo->a = a;
  novo->b = b;
  novo->ya = ya;
  novo->yb = yb;
  return(novo);
}


Mdf* alocaMdf(int nx, int ny, double ax, double bx, double ay, double by){
  Mdf *novo;
  novo = (Mdf*) malloc (sizeof (Mdf));
  if (novo == NULL)
    return NULL;
  novo->nx = nx;
  novo->ny = ny;
  novo->ax = ax;
  novo->bx = bx;
  novo->ay = ay;
  novo->by = by;
  return(novo);
}
/*!
  \brief Liberaçao de memória 

  \param sistema linear SL
  */
void liberaSistLinear (SistLinear_t *SL){
  free(SL->D);
  free(SL->Di1);
  free(SL->Di2);
  free(SL->Ds1);
  free(SL->Ds2);
  free(SL->B);
  free(SL);
  return;
}

/*!
  \brief Impressão de sistema linerar 

  \imprime um SL na saida padrão , caso seja um sistem tridiagonal Di2 e Ds2 ficarão em branco 
  */

void prnSL(SistLinear_t *SL){
  printf("   DS2    ||    DS1    ||      D     ||    DI1    ||    DI2    ||      B        ||\n");
  for (int i = 0; i < SL->n; ++i){
    if(SL->Ds2[i] == 0.0)
      printf("          || ");
    else 
      printf(" %.2e || ",SL->Ds2[i]);
    if(SL->Ds1[i] == 0)
      printf("          || ");
    else 
      printf(" %.2e || ",SL->Ds1[i]);

    if(SL->D[i] == 0)
      printf("          || ");
    else 
      printf(" %.2e || ",SL->D[i]);
    if(SL->Di1[i] == 0.0)
      printf("          || ");
    else 
      printf(" %.2e || ",SL->Di1[i]);
    if(SL->Di2[i] == 0)
      printf("          || ");
    else 
      printf(" %.2e || ",SL->Di2[i]);
    if (SL->B[i] < 0)
      printf(" %.6f \t||\n",SL->B[i]);
    else 
      printf("  %.6f \t||\n",SL->B[i]);
  }  
  printf("\n");

}

