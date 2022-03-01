// ------------------------------------------------------------------
// This program compute the mean and std of a physical quantity
// given in input from file using the bootstrap resampling technique
// ------------------------------------------------------------------


#include<math.h>
#include<stdio.h>
#include<stdlib.h>
#include<time.h>
#include<string.h>

void make_binned_resampling(int L, int len, double *sample, double *resampled_chain);


int main(){

	int len, M=50, N, L;
	double *sample, *resampled_chain, *MEAN, sum=0, final_mean=0, std=0, bh, omega, eta;
	char sample_name[70];
	FILE *file_sample_name, *mean_sigma_len_file, *mean_res_chains_file;

	srand(time(NULL));

	// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!! //
	/////////// CAMBIARE QUI ///////////
	N = 1300;
	omega = 1;
	eta = 0.01;
	L = 1000000;
	sprintf(sample_name, "../results/output/kinetic_energy/eta_%.2lf_omega_%.0lf/energy_N_%d.txt", eta, omega, N);
    ////////////////////////////////////
    // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!! //

    printf("sample_name=%s\n", sample_name);

    // open file with data
    file_sample_name = fopen(sample_name, "r");
    if(file_sample_name==NULL){
        perror("Error opening data file");
        exit(1);
    }

    // file in which writing mean,std,correlation_len
	mean_sigma_len_file = fopen("bootstrap_mean_sigma_len_file.txt","w");
	if(mean_sigma_len_file==NULL){
        perror("Error opening file mean_sigma_len_file.txt");
        exit(1);
    }

    // file in which writing the means of the resampled chains
	mean_res_chains_file = fopen("bootstrap_mean_res_chains_file.txt","w");
	if(mean_res_chains_file==NULL){
        perror("Error opening file mean_res_chains_file.txt");
        exit(1);
    }


    // read from file the data array
	sample = calloc(L, sizeof(double));

	for(int i=0; i<L; i++){
		fscanf(file_sample_name, "%lf", &sample[i]);
	}

	//// CYCLE OVER DIFFERENT CORRELATION LENGTH ////

	MEAN = calloc(M, sizeof(double)); // array with means of each singol chain
	resampled_chain = calloc(L, sizeof(double)); // array of reseampled chain
	
	for (int j=1; j<12; j++){
		len = pow(2,j);

		// let's make the resampled chains
		for(int m=0; m<M; m++){
			//m-resampled chain
			make_binned_resampling(L, len, sample, resampled_chain);

			// compute the mean of the resampled chain
			for(int i=0; i<L; i++){
				sum += resampled_chain[i];
			}
			sum = sum/L;

			// write it over the MEAN[] array
			MEAN[m] = sum;
			fprintf(mean_res_chains_file, "%lf\n", sum);

			sum = 0;
		}

		// compute the final mean...
		for(int i=0; i<M; i++){
				final_mean += MEAN[i];
			}
		final_mean = final_mean/M;

		// ..and its std
		for(int i=0; i<M; i++){
				std += pow((MEAN[i]),2);
			}
		std = sqrt(std/M - pow(final_mean,2));

		// print in stdout...
		printf("your quantity = %.10lf +/- %.10lf\n", final_mean, std);

		// ..and write it over file
		fprintf(mean_sigma_len_file, "%lf   %lf   %d\n", final_mean, std, len);

		final_mean = 0;
		std = 0;
	
	}

	free(sample);
	free(resampled_chain);
	free(MEAN);
	
	fclose(file_sample_name);
	fclose(mean_sigma_len_file);
	fclose(mean_res_chains_file);

	return 0;
}



//-------------------------------------------------//
// FUNCTIONS DEFINITION
//-------------------------------------------------//

void make_binned_resampling(int L, int len, double *sample, double *resampled_chain){
	int *index, start;
	int n=L/len;

	index = calloc(n, sizeof(int));

	// create an array with random index in [0,n-1]
	for(int i=0; i<n; i++){
		index[i] = rand()%n;
	}


	for(int i=0; i<n; i++){
		start = index[i]*len;
		for(int j=0; j<len; j++){
			resampled_chain[len*i+j] = sample[start+j];
		}
	}

	free(index);

	return ;
}