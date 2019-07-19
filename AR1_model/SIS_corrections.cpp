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
  
  //define the container for the sampled events and the new sampled events
  vector < vector < vector < double > > > sample;
  vector < vector < vector < size_t > > > sam_obs;
  //define the container for the sampled observations (of 0s and 1s)
  //and the new sampled observations
  vector < vector < vector < size_t > > > new_sam_obs;
  vector < vector < vector < double > > > new_sample;
  //define the containter for the unnormalised weights
  vector < vector < vector < double > > > un_weights;
  //define the container for the normalised weights
  vector < vector < vector < double > > > weights;
  //define the number of particles
  double n = 10;

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

  //Finding the unnormalised weights (using log then exponentiating)
  //filling the container "un_weights"
  //This is an important part of the code, should be always sure it is correct.
  vector < vector < double > > matrix_un_weights;
  vector < double > vector_un_weights;
  const double constant = ( 1 / ( 2 * sigmasq ) );
  for ( size_t j = 0; j < n; j++ ){
    vector < vector < double > > temp_matrix_sample;
    temp_matrix_sample = sample[j];
    vector < vector < size_t > > temp_matrix_obs;
    temp_matrix_obs = sam_obs[j];
    vector < vector < double > > temp_matrix_new_sample;
    temp_matrix_new_sample = new_sample[j];
    vector < vector < size_t > > temp_matrix_new_obs;
    temp_matrix_new_obs = new_sam_obs[j];
    vector < double > vector_w{1};
    double log_weight;
    double w;
    for (size_t i = 1; i < N; i++){
      vector < double > vector_log_weights;
      vector_log_weights.push_back(1);
      vector < double > row_sample;
      row_sample = temp_matrix_sample[i];
      vector < size_t > row_obs;
      row_obs = temp_matrix_obs[i];
      vector < double > row_new_sample;
      row_new_sample = temp_matrix_new_sample[i];
      vector < size_t > row_new_obs;
      row_new_obs = temp_matrix_new_obs[i];
      double ys;
      double xs;
      for (size_t k = 1; k < row_sample.size(); k++){
	if (row_new_sample[k-1] == row_sample[k-1] && row_new_sample[k] == row_sample[k]){}
	else {
	  ys = - (row_new_sample[k] - phi * row_new_sample[k-1]) * (row_new_sample[k] - phi * row_new_sample[k-1]);
	  xs = (row_sample[k] - phi * row_sample[k-1]) * (row_sample[k] - phi * row_sample[k-1]);
	}
	if (row_new_obs[k] == row_obs[k]){}
	else if (row_new_obs[k] == 1 && row_obs[k] == 0){ys = ys * (1 - p) ; xs = xs * p;}
	else {ys = ys * p ; xs = xs * (1 -p); }
	log_weight = constant * ( ys + xs );
	vector_log_weights.push_back(log_weight);
      }
      double sum = accumulate(vector_log_weights.begin(), vector_log_weights.end(), 0.0);
      w = exp(sum);
      vector_w.push_back(w);
      }
    for (size_t i = 0; i < N; i++){
      for (size_t k = 0; k < i+1; k++){
	vector_un_weights.push_back(vector_w[k]);
      }
      matrix_un_weights.push_back(vector_un_weights);
      vector_un_weights.clear();
    }
    un_weights.push_back(matrix_un_weights);
    matrix_un_weights.clear();
  }
  
  //some sanity check printing for the first particle
  //(can be done on any particle changing the value of "part_num")
  int part_num = 0;
  //Printing one matrix of the 3 dimensional vector "sample"
  cout << "printing one of the sample matrix" << endl;
  print_matrix( sample[part_num] );
  //Printing one matrix the 3 dimensional vector "new_sample"
  cout << "printing one of the new_sample matrix" << endl;
  print_matrix( new_sample[part_num] );
  cout << "printing the matrix of observations" << endl;
  print_matrix(obs);
  //Printing one matrix of the 3 dimensional vector "sam_obs"
  cout << "printing one of the sam_obs matrix " << endl;
  print_matrix( sam_obs[part_num] );
  //Printing one matrix the 3 dimensional vector "new_sam_obs"
  cout << "printing one of the new_sam_obs matrix" << endl;
  print_matrix( new_sam_obs[part_num] );
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
	sum += un_weights[l][i][k];
      }
      for (size_t j = 0; j < n; j++ ){
	weights[j][i][k] = un_weights[j][i][k] / sum;
      }
    }
  }

  //Resampling (no effective sample size)
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
    temp_matrix = un_weights[i];
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

  //Creating a container of 0s of the correct size (lower triangular NxNxn) 
  //called "sam_times_weights" multiplying respectively simulated events and weights
  vector < vector < vector < double > > > sam_times_weights;
  vector < vector < double > > matrix_sam_wei;
  vector < double > vector_sam_wei;
  for ( size_t j = 0; j < n; j++ ){
    double elem;
    for ( size_t i = 0; i < N; i++ ){
      for (size_t k = 0; k < i + 1; k++ ){
	elem = 0;  
	vector_sam_wei.push_back( elem );
      }   
      matrix_sam_wei.push_back( vector_sam_wei );
      vector_sam_wei.clear();
    }
    sam_times_weights.push_back( matrix_sam_wei );
    matrix_sam_wei.clear();
  }
  
  //calculating the 3 dimensional matrix sim_time_weights for simulation and weights
  for ( size_t j = 0; j < n; j++ ){
    for ( size_t i = 0; i < N; i++ ){
      for ( size_t k = 0; k < i + 1; k++ ){
 	sam_times_weights[j][i][k] = new_sample[j][i][k] * weights[j][i][k];
      }
    }
  }

  //calculating the expectation E[y_i] at the last (current) time N
  vector < double > E;
  vector < vector < double > > multiplication_matrix;
  for ( size_t j = 0; j < n; j++ ){
    multiplication_matrix.push_back( sam_times_weights[j][N - 1] );
  }
  for ( size_t i = 0; i < N; i++ ){
    vector < double > multiplication_vector;
    for ( size_t j = 0; j < n; j++ ){
      multiplication_vector.push_back( multiplication_matrix[j][i] );
    }
    double sum;
    for ( auto & n : multiplication_vector)
      sum += n;
    E.push_back( sum );
  }
  cout << "printing vector of expectations " << endl;
  print_vector(E);
  print_vector(X);



  

  //from here I create all the dat files for the plots
  
  //Create a dat file with the values of E
  ofstream outFile3( "./vector_E.dat" );
  outFile3 << endl;
  for ( double n : E ){
    outFile3 << n << endl;
  }
  outFile3.close();

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



