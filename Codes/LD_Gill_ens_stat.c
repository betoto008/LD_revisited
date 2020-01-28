#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <sys/stat.h>

double evolution_population_i(int n_pop, double Total_T, double L, double k_mean_0, double mu, double s, double C, double N0, double *K, double *F, double *K_mean_ensemble, double *K_div_ensemble, double *N_mean_ensemble);
double initial_population(double mu, double s, double L, double k_mean_0, double N0, double *N);
double step_gillespie(double mu, double L, double C, double *N, double *F, double *events, double *events_cumulative, double *T_i, int step, double sum_pop);
void set_fitness(double y, double s, double k0, double L, double *F, FILE *out_fitness);
double update_events(double mu, double L, double C, double *N, double *F, double *events, double *events_cumulative, double sum_pop);
void mean_div_trait(double *N, double *K, double L, double sum_pop, double *K_mean_i, double *K_div_i);
void linspace_R(double *array, double n_i, double n_f, int how_many);
int cfileexists(const char* filename);

int main(int argc, char **argv){ //Parameters from console: mu, s, C, N0, G_0

	clock_t begin = clock();
	//Initialize seed
	srand48(time(NULL));
	//Initialize parameters
  	double k_mean_0;//Initial pop_mean
  	sscanf(argv[5], "%lf", &k_mean_0);
  	double mu; //Mutation rate
	sscanf(argv[1], "%lf", &mu);
  	//printf("mu=%2.2e\n", mu);
  	double s; // Selection coeff.
  	sscanf(argv[2], "%lf", &s);
	//printf("s=%2.2e\n", s);
 	double k0 = 8.0; //lanscape threshold 
	double y = 1; //temperature
  	double C; //Carrying capacity
	sscanf(argv[3], "%lf", &C);
	//printf("C=%.1e\n", C);
  	double N0;
  	sscanf(argv[4], "%lf", &N0);
  	//printf("N0=%.1e\n", N0);
  	double L=1; //Seq. lenght
  	//sscanf(argv[5], "%lf", &L);
  	//printf("L=%f\n", L);
	double Total_T = 400.0; //Total time of simulations 
	int n_pop = 500; //Populations

	//For arrays to be printed
	double *K_mean_ensemble = malloc(sizeof(double)*(n_pop));
	double *K_div_ensemble = malloc(sizeof(double) * (n_pop));
	double *N_mean_ensemble = malloc(sizeof(double) * (n_pop));

	// Initializing arrays K_mean_ensemble and K_div_ensemble
	for(int i=1 ; i<=n_pop; i++){
		K_mean_ensemble[i] = 0.0;
		K_div_ensemble[i] = 0.0;
		N_mean_ensemble[i] = 0.0;
	}

  double *K = malloc(sizeof(double)*(L+2)); //Array with labels of error classes
  double *F = malloc(sizeof(double)*(L+2)); //Array with fitness landscape
  linspace_R(K,0,L,(int)L+1); //Fill the array K of labels
  FILE *out_fitness=fopen("Text_files/Stats/fitness.txt","w+");
  set_fitness(y, s, k0, L, F, out_fitness); // fill the array F with the values of fitness
  fclose(out_fitness);

  // maximum time founded in the evolution of a population.
  double T_max = Total_T;

  //Ensemble simulation
  for(int i = 1; i<=n_pop ; i++){
    // Run the evolution of a population. 
    double T_i_max = evolution_population_i(i, Total_T , L, k_mean_0, mu, s, C, N0, K, F, K_mean_ensemble, K_div_ensemble, N_mean_ensemble);

    //fprintf(out_times, "%f\n", T_i_max);
  	if(T_i_max < T_max){
  		T_max = T_i_max;
  	}
  }
  //fclose(out_times);

  int adapted = 0;
  for(int i=1; i<= n_pop; i++){
    if(K_mean_ensemble[i]<(k0+2)){
      adapted ++;
    }
  }

  //print array of data
  int rep = 0;
  char buf[0x100];
  snprintf(buf, sizeof(buf), "Text_files/Stats/K_ensemble_stat_mu%.1e_s%.1e_C%.1e_N0%.1e_L%.0f_G0%.0f_rep%d.txt", mu, s, C, N0, L, k_mean_0, rep); //create name of new output file
  int exist = cfileexists(buf);

  while(exist==1){
  	rep++;
  	snprintf(buf, sizeof(buf), "Text_files/Stats/K_ensemble_stat_mu%.1e_s%.1e_C%.1e_N0%.1e_L%.0f_G0%.0f_rep%d.txt", mu, s, C, N0, L, k_mean_0, rep); //create name of new output file
  	exist = cfileexists(buf);
  }

  FILE *out_ensemble = fopen(buf,"w+");
  for(int i=1; i<= n_pop; i++){
  	fprintf(out_ensemble, "%f %f %f %f %d\n", Total_T, K_mean_ensemble[i], K_div_ensemble[i], N_mean_ensemble[i], adapted);
  }
	fclose(out_ensemble);

 	//print in console T_max
  printf("Maximum time: %f\n", T_max);

	clock_t end = clock();
	double time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
	// print in console time_spent in minutes
	printf("Ensemble stats Running time: %f minutes.\n\n", time_spent/60.0);
} 
double evolution_population_i(int n_pop, double Total_T, double L, double k_mean_0, double mu, double s, double C, double N0, double *K, double *F, double *K_mean_ensemble, double *K_div_ensemble, double *N_mean_ensemble){

	int iter = 10000000;
	double *N = malloc(sizeof(double)*(L+2)); //Arrays with current population
	double *T_i = malloc(sizeof(double)*(iter+1)); //Arrays with current time
  	double K_mean_i;
  	double K_div_i;

  	T_i[1] = 0.0;

  	//current pop size
  	double sum_pop = 0.0;

  	sum_pop = initial_population(mu, s, L, k_mean_0, N0, N); // fill the array N with the initial popul. Return current pop size

  	// Arrays for saving the events and the cumulative events
  	double *events = malloc(sizeof(double)*(4*L+2+1));
  	double *events_cumulative = malloc(sizeof(double)*(4*L+2));

  	//Counter for current simulation
  	int i = 1;
  	while(T_i[i]<Total_T){
  		if (i>iter){ //check if it is necessary to allocate more memory
  			iter = iter*2;
  			T_i = realloc(T_i, sizeof(double)*(iter+1));
			}
  		double sum_pop_i = step_gillespie(mu, L, C, N, F, events, events_cumulative, T_i, i, sum_pop); //run the Gillespie step. Return the current pop size
		sum_pop = sum_pop_i;		
  		i++;
  	}
  	mean_div_trait(N, K, L, sum_pop, &K_mean_i, &K_div_i); // return the last mean and div trait.
  	K_mean_ensemble[n_pop] = K_mean_i;
 	K_div_ensemble[n_pop] = K_div_i;
 	N_mean_ensemble[n_pop] = sum_pop;
  	// return the last time of the evolution of this population
  return T_i[i];
}
double initial_population(double mu, double s, double L, double k_mean_0, double N0, double *N){

	double sum_pop = 0.0;

    N[1] = N0;
    N[2] = 0;


    return N0;
}
double step_gillespie(double mu, double L, double C, double *N, double *F, double *events, double *events_cumulative, double *T_i, int step, double sum_pop){

	double sum_events = update_events(mu, L, C, N, F, events, events_cumulative, sum_pop);

	//Here we create two random numbers to use in the iteration
	double r1 = drand48();
	double r2 = drand48();

	//Here we choose when will the next event occur.
	double t = -log(r1)/sum_events;

	double delta_N=0.0;

	T_i[step+1] = T_i[step] + t;

	//Now we decide what will occur.

	int j = 1;

	while(r2 > events_cumulative[j] && j < (4*L)+2){
		j++;
	}
    // a division occurs
	if(j <= L+1){
		N[j]++;
		delta_N++;
	}
    // a mutation occurs
	else{
		j = j - (int)(L+2);

		int mut_type = (j/((int)L));
		int mut_id = j % (int)L;
        // a beneficial mutation occurs
		if(mut_type == 0){
			N[mut_id+1]++;
			N[mut_id+2]--;
			//delta_N++;
		}
        // a deleterious mutation occurs
        else if(mut_type == 1){
			N[mut_id+2]++;
			N[mut_id+1]--;
			//delta_N++;
		}
		// a death occurs
		else if(mut_type == 2){
			N[(mut_id+1)]--;
			delta_N--;
		}
		else{
			N[(int)L+1]--;
			delta_N--;
		}
	}
    return sum_pop + delta_N;
}
void set_fitness(double y, double s, double k0, double L, double *F, FILE *out_fitness){
  	F[1] = 0.1;
  fprintf(out_fitness, "%f\n", F[1]);
  F[2] = F[1] + s;
  fprintf(out_fitness, "%f\n", F[2]);
}
double update_events(double mu, double L, double C, double *N, double *F, double *events, double *events_cumulative, double sum_pop){

	
	double sum = 0.0;

	double f_bar = 0.0;

	for(int i = 1; i<=L+1; i++){
		f_bar += F[i]*(N[i]/sum_pop);
	}

	// Birth
	for(int i=1;i<=L+1;i++){
		events[i] = F[i] * N[i];//*(1-(sum_pop/C));
		sum += events[i];
	}
	// Beneficial mutatios
	for(int i=L+2;i<=2*L+1;i++){
		events[i] = mu * N[(int)(i-(L))] * (i-(L+1)); // * F[(int)(i-(L))] * (1-(sum_pop/C));
		sum += events[i];
	}
	// Deletirious mutatios
	for(int i=2*L+2;i<=3*L+1;i++){
		events[i] = mu * N[(int)(i-(2*L+1))] * (L-(i-(2*L+2)));// * F[(int)(i-(2*L+1))] * (1-(sum_pop/C));
		sum += events[i];
	}
	//Death
	for(int i=3*L+2;i<=4*L+2;i++){
		//events[i] = (N[i-((3*(int)L)+1)]*sum_pop) / C;
		//events[i] = (N[i-((3*(int)L)+1)]*sum_pop*F[i-((3*(int)L)+1)]) / C;
    	events[i] = (N[i-((3*(int)L)+1)]*sum_pop*f_bar) / C;
		sum += events[i];
	}
	
    events_cumulative[1] = events[1]/sum;
    //printf("%f ",events_cumulative[1]);

	for(int i=2;i<=4*L+2;i++){
		events_cumulative[i] = events_cumulative[i-1] + events[i]/sum;
	}

	return sum;
}
void mean_div_trait(double *N, double *K, double L, double sum_pop, double *K_mean_i, double *K_div_i){
	// fill the array X
	double mean = 0.0;
	double div = 0.0;
	double x_i = 0.0;
	for(int i=1; i<= L+1 ; i++){
		x_i = N[i]/sum_pop;
		mean += x_i*K[i];
	}
	for(int i=1; i<= L+1 ; i++){
		x_i = N[i]/sum_pop;
		div += x_i*pow((K[i]-mean),2);
	}
	*K_mean_i = mean;
	*K_div_i = div;
}
void linspace_R(double *array, double n_i, double n_f, int how_many){
  for(int i=0;i<how_many;i++){
    array[i+1] = n_i + i*((n_f-n_i)/(how_many-1));
  }
}
int cfileexists(const char* filename){
    struct stat buffer;
    int exist = stat(filename,&buffer);
    if(exist == 0)
        return 1;
    else // -1
        return 0;
}
