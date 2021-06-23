#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "utils.h"
#include "SistemasLineares.h"
#include "Myfunc.h"

int main (){
	SistLinear_t *sistema = lerSistLinear();  	//sistema linear
	SistLinear_t *sgauss;					  	// sistema modificado pela eliminacao de gauss
	real_t *x;								  	//vetor resposta X
	real_t norma;								//variavel para guardar a norma L2 de um vetor
	real_t *res;								//vetor com a diferenca entre X e o ressultado
	int interacoes = 0;							//contador de interacoes para os metodos interativos 
	int i;										//guardar os codigos de saida das funcoes
	double tempo[1];							//varive para guardar o tempo
	int count = 1;								//contador


	/*enquanto a variavel sistema nao for nula  o programa continuara, 
	encontrar o EOF ou algum erro de alocacao de memoria gera um NULL*/ 
	while(sistema != NULL ){					
		x = alocaVetor(sistema->n);				//aloca memoria para um vetor
		res = alocaVetor(sistema->n);			
		limpaVetor(x,sistema->n);				//preenche o vetor com 0
		sgauss = copiaSL(sistema);

		//codigo da eliminacao de gauss 
		interacoes = eliminacaoGauss(sgauss,x,tempo);	

		//imprime  o cabecalho do sistema  	
		printf("***** Sistema %i --> n = %i, erro: %0.3f\n",count, sistema->n, sistema->erro);
		
		//imprime o cabecalho do metodo da eliminacao de gauss 
		printf("===> Eliminação Gauss: %0.10f ms\n",tempo[0]);

		//caso o vetor tenha encontrado um elemento inf, significando um overflow ou uma diviao por zero 
		if(interacoes == -1){
			printf("    -->erro ao encontrar resultado\n\n");
        	fprintf(stderr,"elemento infinito encontrado, overflow ou erro no calculo\n");
		}
		//caso o vetor resposta tenha encontrado um elemento NAN, significando que houve um erro de calculo
		else if(interacoes == -2){
			printf("    -->erro ao encontrar resultado\n\n");
			fprintf(stderr,"elemento NAN encontrado em vetor de respostas,sistema sem solucao ou erro de calculo\n");
		}

		//caso nao tenha acontecido nenhum erro	
		else{

			//imprime o vetor na saida padrao
			prnVetor(x,sistema->n);

			//calculo da norma L2 + impressao dos resultados 
			norma = normaL2Residuo(sistema, x, res);
  			printf("    -->Norma L2 do residuo: %0.6e\n\n",norma);

  			//caso i maior que 0 signifia que houve pelo menos 1 interacao, 
			//i negativos significam erros, 
			//0 significa que nao houve interacoes e o programa nao fez alteracoes
			i = refinamento(sistema,x,tempo);
			if (i > 0){
				//imprime o cabecalho do metodo de refinamento com os resultados
				printf("===>refinamento: %0.6f ms --> %i iteracoes \n", tempo[0], i);
  				prnVetor(x,sistema->n);
  				norma = normaL2Residuo(sistema, x, res);
  				printf("    -->Norma L2 do residuo: %0.6e\n\n",norma);
  			}

  			//imprime em caso de erro 
  			else if (i == -1){
  				printf("===>refinamento: %0.6f ms --> erro \n", tempo[0]);
  				fprintf(stderr,"ERRO: overflow no calculo de da normal\n\n");
  			}
  		}		

  		//codigo interativo de jacobi
		interacoes = gaussJacobi(sistema,x,tempo);

		//imprime o cabecalho do metodo de jacobi e o numero de interacoes , caso haja um erro ele devolve o numero do erro em "interacoes"
		printf("===> Jacobi: %0.10f ms --> %i interacoes\n",tempo[0],interacoes);
		
		//caso a matriz nao tenha atendido o criterio de convergencia
		if( interacoes == -1)
			printf("    --> matriz nao atende requisitos de convergencia\n\n");

		//caso o vetor resposta tenha encontrado um elemento NAN, significando que houve um erro de calculo
		else if (interacoes == -2){
			printf("    -->erro ao encontrar resultado\n\n");
			fprintf(stderr,"elemento NAN encontrado em vetor de respostas,sistema sem solucao ou erro de calculo\n");
		}
		//caso o vetor tenha encontrado um elemento inf, significando um overflow ou uma diviao por zero 
		else if (interacoes == -3){
			printf("    -->erro ao encontrar resultado\n\n");
        	fprintf(stderr,"elemento infinito encontrado, overflow ou erro no calculo\n");
		}

		//caso nao tenha acontecido nenhum erro	
		else {
			//imprime o vetor na saida padrao
			prnVetor(x,sistema->n); 		

			//calculo da norma L2 + impressao dos resultados 
			norma = normaL2Residuo(sistema, x, res);
  			printf("    -->Norma L2 do residuo: %0.6e\n\n",norma);

  			//refinamento interativo do resultado
			i = refinamento(sgauss,x,tempo);

			//caso i maior que 0 signifia que houve pelo menos 1 interacao, 
			//i negativos significam erros, 
			//0 significa que nao houve interacoes e o programa nao fez alteracoes
			if (i > 0){

				//imprime o cabecalho do metodo de refinamento com os resultados
				printf("===>refinamento: %0.6f ms --> %i iteracoes \n", tempo[0], i);
  				prnVetor(x,sistema->n);
  				norma = normaL2Residuo(sistema, x, res);
  				printf("    -->Norma L2 do residuo: %0.6e\n\n",norma);
  			}

  			//imprime em caso de erro 
  			else if (i == -1){
  				printf("===>refinamento: %0.6f ms --> erro \n", tempo[0]);
  				fprintf(stderr,"ERRO: overflow no calculo de da normal\n\n");
  			}
  		}	


  		//codigo interativo de seidel
		interacoes = gaussSeidel(sistema,x,tempo);

		//imprime o cabecalho do metodo de jacobi e o numero de interacoes , caso haja um erro ele devolve o numero do erro em "interacoes"
		printf("===> Gauss-Seidel: %0.10f ms --> %i interacoes\n",tempo[0],interacoes);
			
		//caso a matriz nao tenha atendido o criterio de convergencia
		if( interacoes == -1)
			printf("    -->Matriz nao atende requisitos de convergencia\n\n");

		//caso o vetor resposta tenha encontrado um elemento NAN, significando que houve um erro de calculo
		else if (interacoes == -2){
			printf("    -->Erro ao encontrar resultado\n\n");
			fprintf(stderr,"elemento NAN encontrado em vetor de respostas,sistema sem solucao ou erro de calculo\n");
		}

		//caso o vetor tenha encontrado um elemento inf, significando um overflow ou uma diviao por zero 
		else if (interacoes == -2){
			printf("    -->Erro ao encontrar resultado\n\n");
        	fprintf(stderr,"elemento infinito encontrado, overflow ou erro no calculo\n");
		}

		//caso nao tenha acontecido nenhum erro	
		else {
			//imprime o vetor resposta do metodo na saida padrao
			prnVetor(x,sistema->n);

			//calculo da norma L2 + impressao dos resultados 
			norma = normaL2Residuo(sistema, x, res);
  			printf("    -->Norma L2 do residuo: %0.6e\n\n",norma);

  			//refinamento interativo do resultado
			i = refinamento(sgauss,x,tempo);

			//caso i maior que 0 signifia que houve pelo menos 1 interacao, 
			//i negativos significam erros, 
			//0 significa que nao houve interacoes e o programa nao fez alteracoes
			if (i > 0){

				//imprime o cabecalho do metodo de refinamento com os resultados
				printf("===>refinamento: %0.6f ms --> %i iteracoes \n", tempo[0], i);
  				prnVetor(x,sistema->n);
  				norma = normaL2Residuo(sistema, x, res);
  				printf("    -->Norma L2 do residuo: %0.6e\n\n",norma);
  			}

  			//imprime em caso de erro 
  			else if (i == -1){
  				printf("===>refinamento: %0.6f ms --> erro \n", tempo[0]);
  				fprintf(stderr,"ERRO: overflow no calculo de da normal\n\n");
  			}
  		}	
  		//incrementa o contador 
		count++;

		//desaloca os vetore e o sistema alocados durante a execusao  
		free(x);
		free(res);
		liberaSistLinear(sgauss);
		liberaSistLinear(sistema);

		//aloca o novo sistema, caso nao haja nenhum sistema, retorna NULL o que parara o programa
		sistema = lerSistLinear();
	}
}

