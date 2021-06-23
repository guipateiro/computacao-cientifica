// GRR20197152 Guilherme Costa Pateiro
// Universidade Federal do Parana , Introducao a computacao cientifica
// ultima edicao: 11/06/2021 02:20
#ifndef __MYFUNC_H__
#define __MYFUNC_H__

// Compara 2 estruturas devolvendo a maior diferenca entre elementos
real_t Vdiff(real_t *x , real_t *xnovo,unsigned int n);
real_t Mdiff(real_t **xnovo, real_t **x,int nx,int ny);

//copia todos os elementos de uma estrutura em outra estrutura igual 
void copiaVetor(real_t *x , real_t *xnovo,unsigned int n);
void copiaMatriz(real_t **x,real_t **xnovo,int nx,int ny);
// copia todos os elementos de uma matriz em um vetor
void converteMatrizVetor(real_t **matriz, real_t *vetor,int n, int m);

//aloca memoria para uma estrutura
real_t* alocaVetor(unsigned int n);
real_t** alocaMatriz(int m, int n);

//libera memoria de uma estrutura
void liberaMatriz(real_t **y ,int n);

//imprime uma estrutura
void prnVetor (real_t *v, unsigned int n);
void prnMatriz (real_t **v,int n ,int m);


#endif //__MYFUNC_H__