#include <random>
#include <iostream>
#include <vector>
#include <fstream>
#include <map>
#include <algorithm>
#include <math.h>
#include "SIS_AR1.h"

using namespace std;

void print_matrix( vector < vector < size_t > > m );
void print_matrix( vector < vector < double > > M );
void print_vector( vector < double > v );




//This is the code for the method:
//Firstly I calculate the lower triangular 3-dimentional matrices called for the sampled events
//and for the samples with substitutions for every particle.
//I then calculate the lower triangular 3-dimensional matrices for the sampled observations
//and for the sampled observations with substitutions
//Finally I calculate the weights.

int main(){

  random_device rd;
  mt19937 generator( rd() );

  ar1();

  //DEFINITIONS
  
  //define the container for the sampled events and the sampled observations (0s and 1s)
  vector < vector < vector < double > > > sample;
  vector < vector < vector < size_t > > > sam_obs;
  //define the container for the new sampled events and the new sampled observations (0s and 1s)
  vector < vector < vector < size_t > > > new_sam_obs;
  vector < vector < vector < double > > > new_sample;
  //define the containter for the unnormalised weights
  vector < vector < vector < double > > > un_weights;
  //define the container for the normalised weights
  vector < vector < vector < double > > > weights;
  //define the number of particles
  double n = 1000;

  //Sampling from a normal distribution with mean 0
  //and variance sigma^2/(1-phi^2) for every particle
  //and store this in a vector clalled "vector_y0"
  //This vector will be used as the starting point to sample all other vectors of events
  //that will populate the matrix "sample"
  vector < double > vector_y0;
  for (unsigned j = 0; j < n; j++){
    normal_distribution < double > normalDist(0, sigmasq / (1 - phi * phi));
    vector_y0.push_back(normalDist (generator));
  }
  
  //Sampling for every particle from a normal distribution centred in the previous event times phi
  //and with variance sigma^2, filling the container "sample".
  //Making the substitiution every time I have an observation in real life,
  //filling the container for the new updated events "new_sample".
  for (size_t j = 0; j < n; j++){
    vector < vector < double > > matrix_sample;
    vector < vector < double > > matrix_new_sample;
    vector < double > row_matrix_sample;
    vector < double > row_matrix_new_sample;
    double y;
    y = vector_y0[j];
    row_matrix_sample.push_back(y);
    row_matrix_new_sample.push_back(y);
    matrix_sample.push_back(row_matrix_sample);
    matrix_new_sample.push_back(row_matrix_new_sample);
    for (size_t i = 1; i < N; i++){
      vector < size_t > row_obs;
      row_obs = obs[i];
      normal_distribution < double > normalDist( phi * row_matrix_new_sample[i - 1], sigmasq );
      double gen = normalDist (generator);
      row_matrix_new_sample.push_back(gen);
      row_matrix_sample.push_back(gen);
      for (size_t k = 0; k < i; k++){
	if (row_obs[k] == 1){
	  row_matrix_sample[k] = row_matrix_new_sample[k];
	}
      }
      for (size_t k = 0; k < i + 1; k++){
	if (row_obs[k] == 1){
	  row_matrix_new_sample[k] = X[k];
	}
      }
      matrix_sample.push_back(row_matrix_sample);
      matrix_new_sample.push_back(row_matrix_new_sample);
    }
    row_matrix_sample.clear();
    row_matrix_new_sample.clear();
    sample.push_back(matrix_sample);
    new_sample.push_back(matrix_new_sample);
    matrix_sample.clear();
    matrix_new_sample.clear();
  }

  //Sampling for every particle from a bernoulli distribution with probability p
  //filling the container of the sampled observations "sam_obs".
  //Substituting a0 with a 1 every time I have an observation in real life,
  //filling the matrix of updated sampled observations "new_sam_obs"
  for (size_t j = 0; j < n; j++){
    vector < vector < size_t > > matrix_obs;
    vector < vector < size_t > > matrix_new_obs;
    vector < size_t > row_matrix_obs;
    vector < size_t > row_matrix_new_obs;
    for (size_t i = 0; i < N; i++){
      vector < size_t > row_obs;
      row_obs = obs[i];
      bernoulli_distribution BerDist(p);
      double gen = BerDist (generator);
      row_matrix_obs.push_back(0);
      row_matrix_new_obs.push_back(0);
      for (size_t k = 0; k < i+1; k++){
	if (row_obs[k] == 1){
	  row_matrix_new_obs[k] = 1;
	}
	else if (row_obs[k] == 0){
	  row_matrix_obs[k] = gen;
	}
      }
      for (size_t k = 0; k < i+1; k++){
	if (row_obs[k] == 0 && row_matrix_obs[k] == 1){
	  row_matrix_new_obs[k] = 0;
	}
	else if (row_obs[k] == 0 && row_matrix_obs[k] == 0){
	  row_matrix_new_obs[k] = 0;
	}
      }
      matrix_obs.push_back(row_matrix_obs);
      matrix_new_obs.push_back(row_matrix_new_obs);
    }
    for (size_t i = 0; i < N-1; i++){
      for (size_t k = 0; k < i+1; k++){
	if (matrix_new_obs[i][k] == 1){
	matrix_obs[i+1][k] = matrix_new_obs[i][k];
	}
      }
    }
    row_matrix_obs.clear();
    row_matrix_new_obs.clear();
    sam_obs.push_back(matrix_obs);
    new_sam_obs.push_back(matrix_new_obs);
    matrix_obs.clear();
    matrix_new_obs.clear();
  }

  /* from here few sanity checks to print the relevant matrices for the last (present) time */
  //Sanity check populating and printing "sam_obs_last_time" with the sampled observations for the last time
  vector < vector < size_t > > sam_obs_last_time;
  for (size_t i = 0; i < n; i++){
    vector < vector < size_t > > temp_matrix;
    temp_matrix = sam_obs[i];
    vector < size_t > temp_vector;
    for (size_t j = 0; j < N; j++){
      temp_vector = temp_matrix[N-1];
    }
    sam_obs_last_time.push_back(temp_vector);
  }
  cout << "printing matrix sam_obs for the last time" << endl;
  print_matrix(sam_obs_last_time);

  //Sanity check populating and printing "new_sam_obs_last_time" with new sampled observations for the last time
  vector < vector < size_t > > new_sam_obs_last_time;
  for (size_t i = 0; i < n; i++){
    vector < vector < size_t > > temp_matrix;
    temp_matrix = new_sam_obs[i];
    vector < size_t > temp_vector;
    for (size_t j = 0; j < N; j++){
      temp_vector = temp_matrix[N-1];
    }
    new_sam_obs_last_time.push_back(temp_vector);
  }
  cout << "printing matrix new_sam_obs for the last time" << endl;
  print_matrix(new_sam_obs_last_time);
  
  //Sanity check populating and printing "sam_last_time" with new samples for the last time
  vector < vector < double > > sam_last_time;
  for (size_t i = 0; i < n; i++){
    vector < vector < double > > temp_matrix;
    temp_matrix = sample[i];
    vector < double > temp_vector;
    for (size_t j = 0; j < N; j++){
      temp_vector = temp_matrix[N-1];
    }
    sam_last_time.push_back(temp_vector);
  }
  cout << "printing matrix sam_last_time" << endl;
  print_matrix(sam_last_time);

  //Sanity check populating and printing "new_sam_last_time" with new samples for the last time
  vector < vector < double > > new_sam_last_time;
  for (size_t i = 0; i < n; i++){
    vector < vector < double > > temp_matrix;
    temp_matrix = new_sample[i];
    vector < double > temp_vector;
    for (size_t j = 0; j < N; j++){
      temp_vector = temp_matrix[N-1];
    }
    new_sam_last_time.push_back(temp_vector);
  }
  cout << "printing matrix new_sam_last_time" << endl;
  print_matrix(new_sam_last_time);

  //Transposing 
  vector < vector < vector < double > > > tsample;
  vector < vector < double > > matrix_tsample;
  vector < double > vector_tsample;
  vector < vector < vector < size_t > > > tsam_obs;
  vector < vector < size_t > > matrix_tsam_obs;
  vector < size_t > vector_tsam_obs;
  vector < vector < vector < size_t > > > tnew_sam_obs;
  vector < vector < size_t > > matrix_tnew_sam_obs;
  vector < size_t > vector_tnew_sam_obs;
  vector < vector < vector < double > > > tnew_sam;
  vector < vector < double > > matrix_tnew_sam;
  vector < double > vector_tnew_sam;
  for ( size_t j = 0; j < N; j++ ){
    double sam{0};
    double sam_ob{0};
    double new_sam{0};
    double new_sam_ob{0};
    for ( size_t i = 0; i < n; i++ ){
      for (size_t k = 0; k < j+1; k++ ){
	vector_tsample.push_back( sam );
	vector_tsam_obs.push_back( sam_ob );
	vector_tnew_sam_obs.push_back( new_sam );
	vector_tnew_sam.push_back( new_sam_ob );
      }   
      matrix_tsample.push_back( vector_tsample );
      vector_tsample.clear();
      matrix_tsam_obs.push_back( vector_tsam_obs );
      vector_tsam_obs.clear();
      matrix_tnew_sam_obs.push_back( vector_tnew_sam_obs );
      vector_tnew_sam_obs.clear();
      matrix_tnew_sam.push_back( vector_tnew_sam );
      vector_tnew_sam.clear();
    }
    tsample.push_back( matrix_tsample );
    matrix_tsample.clear();
    tsam_obs.push_back( matrix_tsam_obs );
    matrix_tsam_obs.clear();
    tnew_sam_obs.push_back( matrix_tnew_sam_obs );
    matrix_tnew_sam_obs.clear();
    tnew_sam.push_back( matrix_tnew_sam );
    matrix_tnew_sam.clear();
  }
  for ( size_t j = 0; j < N; j++ ){
    for (size_t i = 0; i < n; i++ ){
      for (size_t k = 0; k < j+1; k++ ){
	tsample[j][i][k]=sample[i][j][k];
	tsam_obs[j][i][k]=sam_obs[i][j][k];
	tnew_sam_obs[j][i][k]=new_sam_obs[i][j][k];
	tnew_sam[j][i][k]=new_sample[i][j][k];
      }
    }
  }
  cout << "tsample" <<endl;
  print_matrix(tsample[9]);
  cout << "tnewsample" <<endl;
  print_matrix(tnew_sam[9]);

  //Finding the unnormalised weights (using log then exponentiating)
  //filling the container "un_weights"
  //This is an important part of the code, should be always sure it is correct.

  //vector < double > vector_un_weights;
  vector < vector < double > > matrix_un_weights;
  for ( size_t j = 0; j < N; j++ ){
    vector < vector < double > > temp_matrix_sample;
    temp_matrix_sample = tsample[j];
    vector < vector < size_t > > temp_matrix_obs;
    temp_matrix_obs = tsam_obs[j];
    vector < vector < double > > temp_matrix_new_sample;
    temp_matrix_new_sample = tnew_sam[j];
    vector < vector < size_t > > temp_matrix_new_obs;
    temp_matrix_new_obs = tnew_sam_obs[j];
    for (size_t i = 0; i < n; i++){
      vector < double > vector_w;
      vector < double > vector_log_num;
      vector < double > vector_log_den;
      vector < double > row_sample;
      row_sample = temp_matrix_sample[i];
      vector < size_t > row_obs;
      row_obs = temp_matrix_obs[i];
      vector < double > row_new_sample;
      row_new_sample = temp_matrix_new_sample[i];
      vector < size_t > row_new_obs;
      row_new_obs = temp_matrix_new_obs[i];
      double w{1};
      for (size_t k = 0; k < j+1; k++){
	double num{0};
	double den{0};
	double exp_sum_num{1};
	double exp_sum_den{1};
	if (k==0){w=1;}
	else if (((row_new_sample[k] == row_sample[k]) && (row_new_sample[k-1] == row_sample[k-1]))){}
	else {
	  num = ((row_new_sample[k] - phi * row_new_sample[k-1]) * (row_new_sample[k] - phi * row_new_sample[k-1]));
	  den = ((row_sample[k] - phi * row_sample[k-1]) * (row_sample[k] - phi * row_sample[k-1]));
	}
	if (k==0){}
	else if (row_new_obs[k] == row_obs[k]){}
	else if (row_new_obs[k] == 1 && row_obs[k] == 0){ num = num + log(p) ; den = den + log(1-p); }
	else { num = num + log(1-p) ; den = den + log((p)); }
	vector_log_num.push_back(num);
	vector_log_den.push_back(den);
	double sum_num = accumulate(vector_log_num.begin(), vector_log_num.end(), 0.0);
	double sum_den = accumulate(vector_log_den.begin(), vector_log_den.end(), 0.0);
	exp_sum_num = exp(sum_num);
	exp_sum_den = exp(sum_den);
	if (exp_sum_den != 0 ){w = exp_sum_num/exp_sum_den;}
	vector_w.push_back(w);
      }
      matrix_un_weights.push_back(vector_w);
      vector_w.clear();
    }
    un_weights.push_back(matrix_un_weights);
    matrix_un_weights.clear();
  }

  
  //Transposing againthe matrix of the weights
  vector < vector < vector < double > > > tweights;
  vector < vector < double > > matrix_tweights;
  vector < double > vector_tweights;
  for ( size_t j = 0; j < n; j++ ){
    double twe{0};
    for ( size_t i = 0; i < N; i++ ){
      for (size_t k = 0; k < i+1; k++ ){ 
	vector_tweights.push_back( twe );
      }   
      matrix_tweights.push_back( vector_tweights );
      vector_tweights.clear();
    }
    tweights.push_back( matrix_tweights );
    matrix_tweights.clear();
    }
  for ( size_t j = 0; j < n; j++ ){
    for (size_t i = 0; i < N; i++ ){
      for (size_t k = 0; k < i+1; k++ ){
	tweights[j][i][k]=un_weights[i][j][k];
      }
    }
  }
  cout << "tweights" <<endl;
  print_matrix(tweights[0]);
  
  
  //some sanity check printing for the first particle
  //(can be done on any particle changing the value of "part_num")
  int part_num = 0;
  //Printing one matrix of the 3 dimensional vector un_weights of the unnormalised weights
  cout << "printing one of the un_weights matrix" << endl;
  print_matrix( un_weights[part_num]);

  //Creating a container of 0s of the correct size (lower triangular NxNxn) called "weights"
  //for the normalised weights
  vector < vector < double > > matrix_w;
  vector < double > vector_w;
  for ( size_t j = 0; j < n; j++ ){
    double elem;
    for ( size_t i = 0; i < N; i++ ){
      for (size_t k = 0; k < i + 1; k++ ){
	elem = 0;  
	vector_w.push_back( elem );
      }   
      matrix_w.push_back( vector_w );
      vector_w.clear();
    }
    weights.push_back( matrix_w );
    matrix_w.clear();
  }


  //Normalising the importance weights and puting them in matrix Weights
  for ( size_t i = 0; i < N; i++ ){
    for (size_t k = 0; k < i + 1; k++){
      double sum{0};
      for (size_t l = 0; l < n; l++ ){
	sum += tweights[l][i][k];
      }
      for (size_t j = 0; j < n; j++ ){
	weights[j][i][k] = tweights[j][i][k] / sum;
      }
    }
  }

  //Resampling (every time)
  for (size_t l = 0; l < N; l++ ){
    vector < vector < double > > weights_each_time;
    for ( size_t i = 0; i < n; i++ ){
      vector < vector < double > > temp_matrix;
      temp_matrix = weights[i];
      vector < double > temp_vector;
      for ( size_t j = 0; j < N; j++ ){
	temp_vector = temp_matrix[l];
      }
      weights_each_time.push_back(temp_vector);
    }
    for ( size_t i = 0; i < l + 1; i++ ){
      vector < double > column_vec;
      for ( size_t k = 0; k < n; k++ ){
	column_vec.push_back( weights_each_time[k][i] );
      }
      vector < double > new_temp;
      discrete_distribution< int > discrete( column_vec.begin(), column_vec.end() );
      for ( size_t k = 0; k < n; k++ ){
	new_temp.push_back( new_sample[discrete(generator)][l][i] );
      }
      for ( size_t j = 0; j < n; j++ ){
	new_sample[j][l][i] = new_temp[j];
      }
      new_temp.clear();
    }
  }
  vector < vector < double > > new_sample_last_time2;
  for ( size_t i = 0; i < n; i++ ){
    vector < vector < double > > temp_matrix;
    temp_matrix = new_sample[i];
    vector < double > temp_vector;
    for ( size_t j = 0; j < N; j++ ){
      temp_vector = temp_matrix[N-1];
    }
    new_sample_last_time2.push_back(temp_vector);
  }
  cout << "printing matrix new_sample_last_time2" << endl;
  print_matrix(new_sample_last_time2);
  cout << "printing matrix new_sample_last_time" << endl;
  print_matrix(new_sam_last_time);



  //sanity checks for the weights printing the last (current) time N for all particle
  
  //Finding the matrix of the unnormalised weights for the last time "un_weights_last_time"
  vector < vector < double > > un_weights_last_time;
  for ( size_t i = 0; i < n; i++ ){
    vector < vector < double > > temp_matrix;
    temp_matrix = tweights[i];
    vector < double > temp_vector;
    for ( size_t j = 0; j < N; j++ ){
      temp_vector = temp_matrix[N-1];
    }
    un_weights_last_time.push_back(temp_vector);
  }
  cout << "printing matrix un_weights_last_time" << endl;
  print_matrix(un_weights_last_time);
  //Finding the matrix of the weights for the last time "weights_last_time"
  vector < vector < double > > weights_last_time;
  for ( size_t i = 0; i < n; i++ ){
    vector < vector < double > > temp_matrix;
    temp_matrix = weights[i];
    vector < double > temp_vector;
    for ( size_t j = 0; j < N; j++ ){
      temp_vector = temp_matrix[N-1];
    }
    weights_last_time.push_back(temp_vector);
  }
  cout << "printing matrix weights_last_time" << endl;
  print_matrix(weights_last_time);


  //Create all the dat files for the plots


  //Create a dat file with the values of the sam_last_time
  //to craete boxplots
  ofstream outFile4( "./sam_last_time.dat" );
  outFile4 << endl;
  for ( size_t lin = 0; lin < n; lin++ ){
    for ( size_t col = 0; col < N; col++ ){
      outFile4 << sam_last_time[lin][col] << " ";
    }
    outFile4 << endl;
  }
  outFile4.close();

  //Create a dat file with the values of the W_N
  //to craete boxplots
  ofstream outFile5( "./weights_last_time.dat" );
  outFile5 << endl;
  for ( size_t lin = 0; lin < n; lin++ ){
    for ( size_t col = 0; col < N; col++ ){
      outFile5 << weights_last_time[lin][col] << " ";
    }
    outFile5 << endl;
  }
  outFile5.close();

  //Create a dat file with the resampled particles
  //to craete boxplots
  ofstream outFile6( "./resample.dat" );
  outFile6 << endl;
  for ( size_t lin = 0; lin < n; lin++ ){
    for ( size_t col = 0; col < N; col++ ){
      outFile6 << new_sample_last_time2[lin][col] << " ";
    }
    outFile6 << endl;
  }
  outFile6.close();
  


  /* This is the code to find the analytic estimates of expectation and variance
     of the missing values in the AR(1) model. These estimates will then be used as 
     gold standard to evaluate our method.*/
  
  
  for ( size_t j = 0; j < N; j++) {cout << vect_obs_N[j] << endl;}
  cout << "end" << endl;
  //case 1: I have one observation at the beginning and one at the end 
  if ( (vect_obs_N[0] == 1) && (vect_obs_N[N-1] == 1) ) {
    vector < size_t > missing{};
    vector < double > expectations{X[0]};
    vector < double > variances{0};
    for ( size_t i = 1; i < N; i++ ){
      if ( (vect_obs_N[i] == vect_obs_N[i-1]) && (vect_obs_N[i] == 1)) {
	expectations.push_back(X[i]);
	variances.push_back(0);
      }
      else if ((vect_obs_N[i] == vect_obs_N[i-1]) && (vect_obs_N[i] == 0)){}
      else {
        missing.push_back(i);
	//these are all the calculations for variances and expectations
	if (missing.size() == 2) {
	  size_t tau = missing[0];
	  size_t m = missing[1];
	  //beginning calculating the denomiator for variances and expectations
	  double den{};
	  vector < double > vecsum{};
	  double power{};
	  for ( size_t j = 0; j < (m-tau+1); j++) {
	    power = pow(phi, 2*j);
	    vecsum.push_back(power); 
	  }
	  for (auto& n : vecsum)
	    den += n;
	  //ending calculating the denomiator for variances and expectations
    	  double power1{};
	  double power2{};
	  vector < double > vecsum1{};
	  vector < double > vecsum2{};
	  double numexp1{};
	  double numexp2{};
	  double expect{};
	  double variance{};
	  //begining calculating the numerator of expectations and variances
	  for (size_t t = tau; t < m; t++){
	    for (size_t k = 0; k < (t-tau+1); k++) {
	      power1 = pow (phi, 2*k);
	      vecsum1.push_back(power1);
	    }
	    double sum1{};
	    for (auto& n : vecsum1)
	      sum1 += n;
	    numexp1 =  sum1 * pow(phi, (m-t)) * X[m];
	    vecsum1.clear();
	    for (size_t j = 0; j < (m-t); j++) {
	      power2 = pow (phi, 2*j);
	      vecsum2.push_back(power2);
	    }
	    double sum2{};
	    for (auto& n : vecsum2)
	      sum2 += n;
	    numexp2 = sum2 * pow(phi, (t-tau+1)) * X[tau-1];
	    vecsum2.clear();
	    variance = sigmasq*sum1*sum2 / den*den;
	    expect = (numexp1 + numexp2) / den;
	    expectations.push_back(expect);
	    variances.push_back(variance);
	  }
	  expectations.push_back(X[m]);
	  variances.push_back(0);
	  //endinging calculating the numerator of expectations and variances
	  missing.clear();
	}
	else {continue;}
      }
    }
    
    //Create a dat file with the values of the Expectations
    ofstream outFile7( "./AR1_interp_exp.dat" );
    outFile7 << endl;
    for ( double n : expectations ){
      outFile7 << n << endl;
    }
    outFile7.close();
    
    //Create a dat file with the values of the Variances
    ofstream outFile8( "./AR1_interp_var.dat" );
    outFile8 << endl;
    for ( double n : variances ){
      outFile8 << n << endl;
    }
    outFile8.close();
  }
  
  //case 2: I have one observation at the beginning and none at the end 
  else if ( (vect_obs_N[0] == 1) && (vect_obs_N[N-1] == 0) ) {
    vector < size_t > missing{};
    vector < double > expectations{X[0]};
    vector < double > variances{0};
    size_t lastvalue{};
    for ( int i = N-1; i >= 0; --i ){
      if ( vect_obs_N[i] == 1 ) {
	lastvalue = i+1;
	break;
      }
    }
    for ( size_t i = 1; i < lastvalue; i++ ){
      if ( (vect_obs_N[i] == vect_obs_N[i-1]) && (vect_obs_N[i] == 1)) {
	expectations.push_back(X[i]);
	variances.push_back(0);
      }
      else if ((vect_obs_N[i] == vect_obs_N[i-1]) && (vect_obs_N[i] == 0)){}
      else {
        missing.push_back(i);
	//these are all the calculations for variances and expectations
	if (missing.size() == 2) {
	  size_t tau = missing[0];
	  size_t m = missing[1];
	  //beginning calculating the denomiator for variances and expectations
	  double den{};
	  vector < double > vecsum{};
	  double power{};
	  for ( size_t j = 0; j < (m-tau+1); j++) {
	    power = pow(phi, 2*j);
	    vecsum.push_back(power); 
	  }
	  for (auto& n : vecsum)
	    den += n;
	  //ending calculating the denomiator for variances and expectations
	  double power1{};
	  double power2{};
	  vector < double > vecsum1{};
	  vector < double > vecsum2{};
	  double numexp1{};
	  double numexp2{};
	  double expect{};
	  double variance{};
	  //begining calculating the numerator of expectations and variances
	  for (size_t t = tau; t < m; t++){
	    for (size_t k = 0; k < (t-tau+1); k++) {
	      power1 = pow (phi, 2*k);
	      vecsum1.push_back(power1);
	    }
	    double sum1{};
	    for (auto& n : vecsum1)
	      sum1 += n;
	    numexp1 =  sum1 * pow(phi, (m-t)) * X[m];
	    vecsum1.clear();
	    for (size_t j = 0; j < (m-t); j++) {
	      power2 = pow (phi, 2*j);
	      vecsum2.push_back(power2);
	    }
	    double sum2{};
	    for (auto& n : vecsum2)
	      sum2 += n;
	    numexp2 = sum2 * pow(phi, (t-tau+1)) * X[tau-1];
	    vecsum2.clear();
	    variance = sigmasq*sum1*sum2 / den*den;
	    expect = (numexp1 + numexp2) / den;
	    expectations.push_back(expect);
	    variances.push_back(variance);
	  }
	  expectations.push_back(X[m]);
	  variances.push_back(0);
	  //endinging calculating the numerator of expectations and variances
	  missing.clear();
	}
	else {continue;}
      }
    }
    // these are the estimations of the last missing points. Step needed because
    // I don't know the value of the last time
    for (size_t i = 0; i < N-lastvalue; i++) {
    double expect{};
    double variance{};
    double power{};
    vector < double > vecsum{};
    expect = pow(phi,(i+1)) * X[lastvalue-1];
    power = pow (phi,2*i);
    expectations.push_back(expect);
    vecsum.push_back(power);
    double sum{};
    for (auto& n : vecsum)
      sum += n;
    variance = sum*sigmasq;
    variances.push_back(variance);
    }
    // end of the estimations of the last missing points.
    print_vector(expectations);
    print_vector(variances);
    
    //Create a dat file with the values of the Expectations
    ofstream outFile7( "./AR1_interp_exp.dat" );
    outFile7 << endl;
    for ( double n : expectations ){
      outFile7 << n << endl;
    }
    outFile7.close();
    
    //Create a dat file with the values of the Variances
    ofstream outFile8( "./AR1_interp_var.dat" );
    outFile8 << endl;
    for ( double n : variances ){
      outFile8 << n << endl;
    }
    outFile8.close();
  }
  
  //case 3: I have no observation at the beginning but one at the end

  else if ( (vect_obs_N[0] == 0) && (vect_obs_N[N-1] == 1) ) {
    // these are estimations of the first missing point. Step needed because
    // I don't know the value of the first time
    vector < size_t > missing{};
    vector < double > expectations{};
    vector < double > variances{};
    size_t firstvalue{};
    for ( size_t i = 0; i < N; i++ ){
      if ( vect_obs_N[i] == 1 ) {
	firstvalue = i;
	break;
      }
    }
    for ( size_t i = 0; i < firstvalue; i++ ){
      double expect{};
      double variance{};
      double power{};
      vector < double > vecsum{};
      expect = pow(phi,(firstvalue-i)) * X[firstvalue];
      power = pow(phi,2*(firstvalue-i-1));
      expectations.push_back(expect);
      vecsum.push_back(power);
      double sum{};
      for (auto& n : vecsum)
	sum += n;
      variance = sum*sigmasq;
      variances.push_back(variance);
    }
    // end of the estimations of the first missing points.
    expectations.push_back(X[firstvalue]);
    for ( size_t i = firstvalue+1; i < N; i++ ){
      if ( (vect_obs_N[i] == vect_obs_N[i-1]) && (vect_obs_N[i] == 1)) {
	expectations.push_back(X[i]);
	variances.push_back(0);
      }
      else if ((vect_obs_N[i] == vect_obs_N[i-1]) && (vect_obs_N[i] == 0)){}
      else {
        missing.push_back(i);
	//these are all the calculations for variances and expectations
	if (missing.size() == 2) {
	  size_t tau = missing[0];
	  size_t m = missing[1];
	  //beginning calculating the denomiator for variances and expectations
	  double den{};
	  vector < double > vecsum{};
	  double power{};
	  for ( size_t j = 0; j < (m-tau+1); j++) {
	    power = pow(phi, 2*j);
	    vecsum.push_back(power); 
	  }
	  for (auto& n : vecsum)
	    den += n;
	  //ending calculating the denomiator for variances and expectations
	  double power1{};
	  double power2{};
	  vector < double > vecsum1{};
	  vector < double > vecsum2{};
	  double numexp1{};
	  double numexp2{};
	  double expect{};
	  double variance{};
	  //begining calculating the numerator of expectations and variances
	  for (size_t t = tau; t < m; t++){
	    for (size_t k = 0; k < (t-tau+1); k++) {
	      power1 = pow (phi, 2*k);
	      vecsum1.push_back(power1);
	    }
	    double sum1{};
	    for (auto& n : vecsum1)
	      sum1 += n;
	    numexp1 =  sum1 * pow(phi, (m-t)) * X[m];
	    vecsum1.clear();
	    for (size_t j = 0; j < (m-t); j++) {
	      power2 = pow (phi, 2*j);
	      vecsum2.push_back(power2);
	    }
	    double sum2{};
	    for (auto& n : vecsum2)
	      sum2 += n;
	    numexp2 = sum2 * pow(phi, (t-tau+1)) * X[tau-1];
	    vecsum2.clear();
	    variance = sigmasq*sum1*sum2 / den*den;
	    expect = (numexp1 + numexp2) / den;
	    expectations.push_back(expect);
	    variances.push_back(variance);
	  }
	  expectations.push_back(X[m]);
	  variances.push_back(0);
	  //endinging calculating the numerator of expectations and variances
	  missing.clear();
	}
	else {continue;}
      }
    }

    //Create a dat file with the values of the Expectations
    ofstream outFile7( "./AR1_interp_exp.dat" );
    outFile7 << endl;
    for ( double n : expectations ){
      outFile7 << n << endl;
    }
    outFile7.close();
    
    //Create a dat file with the values of the Variances
    ofstream outFile8( "./AR1_interp_var.dat" );
    outFile8 << endl;
    for ( double n : variances ){
      outFile8 << n << endl;
    }
    outFile8.close();
  }

  //case 4: I have no observation at the beginning or at the end 
  else {
    // these are estimations of the first missing points. Step needed because
    // I don't know the value of the first time
    vector < size_t > missing{};
    vector < double > expectations{};
    vector < double > variances{};
    size_t firstvalue{};
    for ( size_t i = 0; i < N; i++ ){
      if ( vect_obs_N[i] == 1 ) {
	firstvalue = i;
	break;
      }
    }
    for ( size_t i = 0; i < firstvalue; i++ ){
      double expect{};
      double variance{};
      double power{};
      vector < double > vecsum{};
      expect = pow(phi,(firstvalue-i)) * X[firstvalue];
      power = pow(phi,2*(firstvalue-i-1));
      expectations.push_back(expect);
      vecsum.push_back(power);
      double sum{};
      for (auto& n : vecsum)
	sum += n;
      variance = sum*sigmasq;
      variances.push_back(variance);
    }
    cout << "first value " << firstvalue << endl;
    expectations.push_back(X[firstvalue]);
    // end of the estimations of the first missing points.
    size_t lastvalue{};
    for ( int i = N-1; i >= 0; --i ){
      if ( vect_obs_N[i] == 1 ) {
	lastvalue = i+1;
	break;
      }
    }
    cout << "last value " << lastvalue << endl;
    for ( size_t i = firstvalue+1; i < lastvalue; i++ ){
      if ( (vect_obs_N[i] == vect_obs_N[i-1]) && (vect_obs_N[i] == 1)) {
	expectations.push_back(X[i]);
	variances.push_back(0);
      }
      else if ((vect_obs_N[i] == vect_obs_N[i-1]) && (vect_obs_N[i] == 0)){}
      else {
        missing.push_back(i);
	//these are all the calculations for variances and expectations
	if (missing.size() == 2) {
	  size_t tau = missing[0];
	  size_t m = missing[1];
	  cout << "tau" << tau << "m" << m << endl;
	  //beginning calculating the denomiator for variances and expectations
	  double den{};
	  vector < double > vecsum{};
	  double power{};
	  for ( size_t j = 0; j < (m-tau+1); j++) {
	    power = pow(phi, 2*j);
	    vecsum.push_back(power); 
	  }
	  for (auto& n : vecsum)
	    den += n;
	  cout << " denominator " << den << endl;
	  //ending calculating the denomiator for variances and expectations
	  double power1{};
	  double power2{};
	  vector < double > vecsum1{};
	  vector < double > vecsum2{};
	  double numexp1{};
	  double numexp2{};
	  double expect{};
	  double variance{};
	  //begining calculating the numerator of expectations and variances
	  for (size_t t = tau; t < m; t++){
	    for (size_t k = 0; k < (t-tau+1); k++) {
	      power1 = pow (phi, 2*k);
	      vecsum1.push_back(power1);
	    }
	    print_vector(vecsum1);
	    double sum1{};
	    for (auto& n : vecsum1)
	      sum1 += n;
	    numexp1 =  sum1 * pow(phi, (m-t)) * X[m];
	    vecsum1.clear();
	    cout << "numexp1 " << numexp1 << endl;
	    cout << "sum1 " << sum1 << endl;
	    for (size_t j = 0; j < (m-t); j++) {
	      power2 = pow (phi, 2*j);
	      vecsum2.push_back(power2);
	    }
	    print_vector(vecsum2);
	    double sum2{};
	    for (auto& n : vecsum2)
	      sum2 += n;
	    numexp2 = sum2 * pow(phi, (t-tau+1)) * X[tau-1];
	    vecsum2.clear();
	    variance = sigmasq*sum1*sum2 / den*den;
	    cout << "numexp2 " << numexp2 << endl;
	    cout << "sum2 " << sum2 << endl;
	    expect = (numexp1 + numexp2) / den;
	    expectations.push_back(expect);
	    variances.push_back(variance);
	  }
	  expectations.push_back(X[m]);
	  variances.push_back(0);
	  //endinging calculating the numerator of expectations and variances
	  missing.clear();
	}
	else {continue;}
      }
    }
    print_vector(expectations);
    print_vector(variances);
    // this is the estimation for the last missing points. It is needed because
    // I don't know the value of the last time
    for (size_t i = 0; i < N-lastvalue; i++) {
    double expect{};
    double variance{};
    double power{};
    vector < double > vecsum{};
    expect = pow(phi,(i+1)) * X[lastvalue-1];
    power = pow (phi,2*i);
    expectations.push_back(expect);
    vecsum.push_back(power);
    double sum{};
    for (auto& n : vecsum)
      sum += n;
    variance = sum*sigmasq;
    variances.push_back(variance);
    }  
    print_vector(expectations);
    print_vector(variances);

    //Create a dat file with the values of the Expectations
    ofstream outFile7( "./AR1_interp_exp.dat" );
    outFile7 << endl;
    for ( double n : expectations ){
      outFile7 << n << endl;
    }
    outFile7.close();
    
    //Create a dat file with the values of the Variances
    ofstream outFile8( "./AR1_interp_var.dat" );
    outFile8 << endl;
    for ( double n : variances ){
      outFile8 << n << endl;
    }
    outFile8.close();
    
  }
  
  return 0;
  
}



//functions definitions

//this function prints a matrix of unsigned size_t
void print_matrix ( vector < vector < size_t > > m ){
  for ( const vector < size_t > v : m ){
    for  ( size_t x : v ) cout << x << ' ';
    cout << endl;
  }
}

//this function prints a matrix of signed doubles
void print_matrix ( vector < vector < double > > M ){
  for ( const vector < double > v : M ){
    for  ( double x : v ) cout << x << ' ';
    cout << endl;
  }
}

//this function prints a vector of doubles
void print_vector ( vector < double > v ){
  for ( const double x : v ) cout << x << ' ';
  cout << endl;
}



