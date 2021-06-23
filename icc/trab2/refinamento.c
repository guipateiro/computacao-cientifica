//Objetivo do programa : Encontrar a aproximacao xx = (xx1,...,xxn) ou a
// mensagem de que o numero de iteracoes foi excedido
// e uma aproximacao COND para K(A)
// O programa foi desenvolvido em Linguagem C
//Dificuldades : 1)Implementar a resolucao do sistema Ay = r, usando o metodo
// de eliminacao de Gauss, devido ao numero de iteracoes;
// 2)Implementar o calculo do maximo valor absoluto dos vetores
// y, xx e do vetor denominado "sb", que e a subtracao do vetor
// xx de x.
#include <stdio.h>
//#include <conio.h>
#include <math.h>
#include <stdlib.h>

double max(double value1,double value2){
return ( (value1 > value2) ? value1 : value2);
}

void gauss(); //Metodo de eliminacao de gauss
	static char temp[10];
	double A[11][11];
	double X,Y,XX,SB,COND,TOL,B[11],xx[11],x[11],y[11],r[11],sb[11];
	int aux,i,j,n,t,s,k,K,N;
//programa principal
void main(){
	int i,j;
	 system("cls");
	printf ("Digite a ordem n da matriz quadrada A (max. 10): ");
	gets(temp);
	printf ("\n\n");
	n = atoi(temp);
	for (i=1; i<=n;i++){
		for (j=1; j<=n; j++){
			printf ("Entre com valor de A[%d][%d]==> ",i,j);
			gets (temp);
			A[i][j]=atof(temp);
		}
		printf ("Entre com o valor de b[%d]===> ",i);
		gets (temp);
		B[i]=atof(temp);
	}
	printf ("Entre com o numero de digitos de precisao t:");
	scanf ("%d",&t);
	printf ("Entre com o numero maximo de iteracoes:");
	scanf ("%d",&N);
	printf ("Entre com a tolerancia:");
	scanf ("%d",&TOL);
	gauss(); //para resolver o sistema Ax = b
	K=1;
	//faz todo esse procedimento enquanto numero de iteracoes for maior ou igual a K
	do{
		for (i=1;i<=n;i++){ //calculo de r{
			s=0;
			for (j=1;j<=n;j++) 
				s = s + (A[i][j]*x[j]);
			r[i] = B[i] - s;
		}
		aux = 1; //variavel de controle
		gauss(); //para resolver Ay=r
		for (i=1;i<=n;i++) 
			xx[i]=x[i]+y[i]; //calculo da aproximacao xx
		if (K==1){
			Y = y[1];
			XX = xx[1];
			for (j=1;j<=n;j++){
				sb[j] = x[j] - xx[j]; //subtrai xx de x
				Y = max(abs(Y),abs(y[j+1])); //pega o maximo valor, em modulo, do vetor y
				XX = max(abs(XX),abs(xx[j+1])); //pega o maximo valor, em modulo, de xx
			}
			COND = (Y/XX)*(pow(10,t)); //calcula o valor da aproximacao COND para K(A)
		}
		SB = sb[1];
		for (j=1;j<=n;j++)
			SB = max(abs(SB),abs(sb[j+1]));
			//calcula o maximo valor, em modulo, da subtracao x - xx
		if (SB < TOL){ //se tal valor for menor que a tolerancia,
			//retorna o valor da aproximacao xx e da aproximacao COND para K(A)
			 system("cls");
			printf ("Os valores de xx sao: \n");
			for (j=1; j<=n; j++){
				printf ("\nxx[%d]= %.4f",j,xx[j]);
			}
			printf("\nO valor de COND e %d. ",&COND);
			getchar();
		}
		K++; //incrementa K
		for (j=1;j<=n;j++) 
			x[j]=xx[j];
	}
	while(N>=K);

	if (SB >= TOL){ //se tal valor for maior que a tolerancia,
 //retorna mensagem de que numero de iteracoes foi excedido e
		system("cls"); //o valor da aproximacao COND para K(A)
		printf("\nMaximo numero de iteracoes foi excedido!");
		printf("\nO valor de COND e %d. ",&COND);
		getchar();
	}
}
//resolve o sistema Ax=b e Ay=r, pelo metodo de eliminacao de Gauss
void gauss(){
	int l;
	double m;
	for (k=1;k<=n-1;k++){
		for (i=k+1;i<=n;i++){
			m=A[i][k]/A[k][k];
			A[i][k]=0;
			for (j=k+1;j<=n;j++){
				A[i][j]=A[i][j]-(m*A[k][j]);
			}
			if (aux==1) 
				r[i]=r[i]-m*r[k]; //se for o sistema Ay=r
			else 
				B[i]=B[i]-m*B[k]; //se for o sistema Ax=b
		}
	}
	if (aux==1){ //se for o sistema Ay=r
		y[n]=r[n]/A[n][n];
		for (l=n-1;l>=1;l--){
			y[l]=0;
			for (j=l+1;j<=n;j++){
				y[l]=y[l]+A[l][j]*y[j];
			}
			y[l]=(r[l]-y[l])/A[l][l];
		}
	}
	else{ //se for o sistema Ax=b
		x[n]=B[n]/A[n][n];
		for (l=n-1;l>=1;l--){
			x[l]=0;
			for (j=l+1;j<=n;j++){
				x[l]=x[l]+A[l][j]*x[j];
			}
			x[l]=(B[l]-x[l])/A[l][l];
		}
	}
}