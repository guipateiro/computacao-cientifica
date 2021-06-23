// GRR20197152 Guilherme Costa Pateiro
// Universidade Federal do Parana , Introducao a computacao cientifica
// ultima edicao: 19/05/2021 00:01
#include <stdio.h>
#include <math.h>
#include <float.h>
#include <stdlib.h>

//declaracao de struct de apoio para leitura da entrada
struct funcaot{
	char carac;
	int Xresp;
	int Xin1;
	int Xin2;
	char op;
} typedef funcao;

int compararfloatN(float A, float B){			//funcao para comparar igualdade entre 2 floats 
	float diff = fabs(A - B);
	if ((A <= 0 && B > 0) || (A > 0 && B <= 0))
		return 0;
	if ((A == INFINITY && B != INFINITY) || (B == INFINITY && A != INFINITY) || (A == -INFINITY && B != -INFINITY) || (B == -INFINITY && A != -INFINITY))
		return 0;
	if (A == NAN || B == NAN || B == -INFINITY || A == INFINITY)
		return 1;
	A = fabs(A);
	B = fabs(B);
	float maior = (B > A) ? B : A ;
	if (diff <= maior*FLT_EPSILON)
		return 1;
	return 0;
}

int compararfloatM(float A, float B){			//funcao para comparar igualdade entre 2 floats 
	float diff = fabs(A - B);
	if ((A <= 0 && B > 0) || (A > 0 && B <= 0))
		return 0;
	if ((A == INFINITY && B != INFINITY) || (B == INFINITY && A != INFINITY) || (A == -INFINITY && B != -INFINITY) || (B == -INFINITY && A != -INFINITY))
		return 0;
	if (A == NAN || B == NAN || B == -INFINITY || A == INFINITY)
		return 1;
	A = fabs(A);
	B = fabs(B);
	float maior = (B > A) ? B : A ;
	if (diff < maior*FLT_EPSILON)
		return 1;
	return 0;
}

int verificadiv0(float A, float B){ 					// verifica se existe um 0 no espaco entre A e B 
	if (compararfloatN(A+B,0.0) || compararfloatN(A+B,-0.0))
		return 1;
	return 0;
}

//main
int main(){
	int m,n,i,j;					//variaveis para contagem 

	scanf("%i %i",&m , &n);			//leitura de M e N 
	float entrada;					//numero de entrada 
	float min[m+n+1];				//vetor de minimos
	float max[m+n+1];				// vetor de maximos
	float aux[4];					//vetor auxiliar para calculos da multiplicacao
	funcao func = {0,0,0,0,0};		//struct para leitura do teclado
	
	//printf("%0.8e\n",FLT_EPSILON );
	for (i = 1; i <= m; i++){		//insercao dos M primeiros numeros 
		func.carac = 0;
		while (func.carac != 'X' && func.carac != 'x')
			scanf("%c",&func.carac);
		scanf("%i",&func.Xresp);
		scanf("%f", &entrada);
		min[i] = nextafterf(entrada,-INFINITY);			//criacao dos intervalos usando como base o arredondamento do numeros 
		max[i] = nextafterf(entrada,INFINITY);	

		if(compararfloatM(min[i],max[i])){ 				// verificacao se nenhum espaco invalido foi criado durante a execusao ou se o [a,b] é unitario
			printf("erro nos valores de X%i\n", i);
			exit(-1);
		}	
	}
	for(i = 0; i < n; i++){								//leitura das funcoes de N 
		func.carac = 0;
		while (func.carac != 'X' && func.carac != 'x')
			scanf("%c",&func.carac);
		scanf("%i",&func.Xresp);
		func.carac = 0;
		while (func.carac != 'X' && func.carac != 'x')
			scanf("%c",&func.carac);
		scanf("%i",&func.Xin1);
		func.carac = 0;
		while (func.carac != '*' && func.carac != '/' && func.carac != '+' && func.carac != '-')
			scanf("%c",&func.carac);
		func.op = func.carac;
		func.carac = 0;
		while (func.carac != 'X' && func.carac != 'x')
			scanf("%c",&func.carac);
		scanf("%i",&func.Xin2);

		if(func.Xresp > m+i+1){
			printf("Erro, X%i nao pode ser definido antes de definir todos os X antecessores\n", m+i+1);	//tratamento de erro caso as varaivies nao existam ou o numero delas nao esteja correto
			exit(1);
		}
		if(func.Xin1 <= 0 || func.Xin1 >= func.Xresp ){
			printf("Erro [Xin1] variavel nao existente\n");
			exit(1);
		}
		if(func.Xin2 <= 0 || func.Xin2 >= func.Xresp ){
			printf("Erro [Xin2] variavel nao existente\n");
			exit(1);
		}
		switch(func.op){												//area de execusao algebrica 
			case '+':													//soma + correcao de intervalos
				min[m+i+1]= min[func.Xin1] + min[func.Xin2];
				max[m+i+1]= max[func.Xin1] + max[func.Xin2];
				min[m+i+1]= nextafterf(min[m+i+1],-INFINITY);
				max[m+i+1]= nextafterf(max[m+i+1],INFINITY);

			break;

			case '-':													//subtracao + correcao de intervalos
				min[m+i+1]= min[func.Xin1] - max[func.Xin2];
				max[m+i+1]= max[func.Xin1] - min[func.Xin2];
				min[m+i+1]= nextafterf(min[m+i+1],-INFINITY);
				max[m+i+1]= nextafterf(max[m+i+1],INFINITY);
			break;

			case '*':													//multiplicacao + correcao de intervalos
	    		aux[0] = min[func.Xin1] * max[func.Xin2];
	  			aux[1] = max[func.Xin1] * min[func.Xin2];
   				aux[2] = max[func.Xin1] * max[func.Xin2];
  			    aux[3] = min[func.Xin1] * min[func.Xin2];
  			    max[m+i+1] = aux[0];
			    min[m+i+1] = aux[0];
			    for(j=0; j<=3; j++){
			        if(aux[j]>max[m+i+1]){
			            max[m+i+1] = aux[j];
        			}
			        if(aux[j]<min[m+i+1]){
			            min[m+i+1] = aux[j];
        			}
    			}
    			min[m+i+1]= nextafterf(min[m+i+1],-INFINITY);
				max[m+i+1]= nextafterf(max[m+i+1],INFINITY);
			break;

			case '/':													//divisao + correcao de intervalos
				min[m+i+1]= min[func.Xin1] * (1/ max[func.Xin2]);
				max[m+i+1]= max[func.Xin1] * (1/ min[func.Xin2]);
				if (verificadiv0(min[func.Xin2],max[func.Xin2])){		//caso especial caso tenha um 0 no intervalo
					min[m+i+1]= -INFINITY;
					max[m+i+1]= INFINITY;
				}
				else {
					min[m+i+1]= nextafterf(min[m+i+1],-INFINITY);
					max[m+i+1]= nextafterf(max[m+i+1],INFINITY);
				}	
			break;
		}
	}	
	for (i = m+1; i <= m+n; i++){					// verificacao se nenhum espaco invalido foi criado durante a execusao ou se o [a,b] é unitario
		if(compararfloatN(min[i],max[i])){
			printf("erro nos valores de X%i\n", i);
			exit(-1);
		}	
	}

	for (i = 1; i <= m+n; i++){
		printf("X%i = [ %1.8e\t, %0.8e]\n",i, min[i], max[i]);		// impressao de todos os elementos 
	}
	printf("\nNão unitários:\n");
	for (i = m+1; i <= m+n; i++){
		printf("X%i = [ %1.8e\t, %0.8e]\n",i, min[i], max[i]);		//impressao dos elementos nao unitarios
	}

	return 0;
}


