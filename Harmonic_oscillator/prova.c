#include<math.h>
#include<stdio.h>
#include<stdlib.h>
#include<time.h>
#include<string.h>



int main(){
	int N=3, L=4, k=2;
	double **matrix, *vector;
	FILE *file1;

	file1 = fopen("prova.txt", "r");
    if(file1==NULL){
        perror("Errore in apertura del file");
        exit(1);
    }

	matrix = calloc(L, sizeof(double*));
	for(int i=0; i<L; i++){
		matrix[i] = calloc(N, sizeof(double));
	}

	for(int i=0; i<L; i++){
		for(int j=0; j<N; j++){
			fscanf(file1, "%lf", &matrix[i][j]);
		}
	}

	vector = calloc(L, sizeof(double));

	for(int i=0; i<L; i++){
		vector[i] = matrix[i][k];
	}

	for(int i=0; i<L; i++){
		printf("%lf\n", vector[i]);
	}
	printf("\n");

	for(int i=0; i<L; i++){
		for(int j=0; j<N; j++){
			printf("%lf  ", matrix[i][j]);
		}
		printf("\n");
	}

	return 0;
}