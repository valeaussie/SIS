#include <random>
#include <iostream>
#include <vector>
#include <fstream>
#include <map>
#include <algorithm>
#include <math.h>
#include <tuple>

using namespace std;

/* this is the code to simulate the AR(1) model and find the vector of the observations.
Sampling from a normal distribution I populate a vector X of events.
Sampling form a geometric distribution I populate a vector "vector_ti" for the times of observations
for the event that happened at time i (event x_i will be observed at time t_i from when it happened).
Then creating a matrix of observations "Obs" that will have on each row 
the events that have been observed up to the time corresponding to the row number. 
So, at time 0 I will have 1 element in the row that might or might not have been observed,
at time 1 I will have two elements on the row, some observed, some not, and so on
I will have 0s whenever the element have have not been yet obeserved.
The values of the parameters at this stage are fixed and are sigma^2 = 1, 
phi = 0.5, p = 0.4 */

//DEFINITIONS

const double sigmasq = 1;
const float phi = 0.5;
const float p = 0.3;
const double N = 5;
vector < double > X;
vector < vector < size_t > > obs;
vector < size_t > vect_obs_N;

void print_matrix(vector < vector < size_t > > m);
void print_matrix(vector < vector < double > > M);
void print_vector(vector < double > v);
void print_vector(vector < int > V);
void print_vector(vector < size_t > V);



random_device rd;
mt19937 generator(rd());


int main(){

  //Sampling from a normal distribution. Put the values in a vector "X"
  normal_distribution < double > normalDist(0, sigmasq / (1 - phi * phi));
  X.push_back(normalDist(generator));
  
  for (size_t i = 1; i < N; i++){  
    normal_distribution < double > normalDist(phi * X[i - 1], sigmasq);
    X.push_back( normalDist(generator));
  }  

  //Sampling from a geometric distribution the values ti
  //and create a vector called "vector_ti".
  vector < double > vector_ti;
  for (size_t i = 0; i < N; i++){
    geometric_distribution <> geoDist(p);
    size_t ti = geoDist (generator);
    vector_ti.push_back(ti);
  }
  
  
  //Populating the matrix of observations callled "obs"
  //then creating a vector "vect_obs_N" for the final time
  for (size_t j = 0; j < N; j++){
    vector < size_t > tempvec{};
    for (size_t i = 0; i < j + 1; i++){
      if (vector_ti[i] <= j - i){
	tempvec.push_back(1);
      }
      else  {
	tempvec.push_back(0);
      }
    }
    obs.push_back(tempvec);
  }
  vect_obs_N = obs [N - 1];
  
  /*  cout << "vector x \n";
  print_vector(X);
  cout << "vector ti \n";
  print_vector(vector_ti);
  cout << "matrix Obs \n";
  print_matrix_double(Obs); */

  //Creating a dat file with the values of the vector of the observed events "vect_obs_N"
  //at the current time N calling it "real_data.dat"
  ofstream outFile1("./real_data.dat");
  outFile1 << endl;
  for (double n : vect_obs_N){
    outFile1 << n << endl;
  }
  outFile1.close();

  //Printing on the screen the vector "X" of all samples 
  cout << "printing the vector X of data \n";
  print_vector(X);

  //Creating a dat file with the values of X
  //calling it "vector_X.dat""
  std::ofstream outFile2("./vector_X.dat");
  outFile2 << endl;
  for (double n : X){
    outFile2 << n << endl;
  }
  outFile2.close();

  
  //This is the code for the method:
  //Firstly I calculate the lower triangular 3-dimentional matrices called for the sampled events
  //and for the samples with substitutions for every particle.
  //I then calculate the lower triangular 3-dimensional matrices for the sampled observations
  //and for the sampled observations with substitutions
  //Finally I calculate the weights.

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

//this function prints a vector of integers
void print_vector ( vector < int > v ){
  for ( const int x : v ) cout << x << ' ';
  cout << endl;
}

//this function prints a vector of integers
void print_vector ( vector < size_t > v ){
  for ( const int x : v ) cout << x << ' ';
  cout << endl;
}
