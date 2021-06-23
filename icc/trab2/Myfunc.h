#ifndef __MYFUNC_H__
#define __MYFUNC_H__


real_t Vdiff(real_t *x , real_t *xnovo,unsigned int n);

void copiaVetor(real_t *x , real_t *xnovo,unsigned int n);

void limpaVetor(real_t *x,unsigned int n);

real_t* alocaVetor(unsigned int n);

SistLinear_t* copiaSL(SistLinear_t *SL);

void retrossubs(real_t **SL, real_t *b, real_t *x, unsigned int n);

void pivoteamento(real_t **SL, real_t *b, unsigned int n);

#endif