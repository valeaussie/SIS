#include <random>
#include <iostream>
#include <vector>
#include <fstream>
#include <map>
#include <algorithm>
#include <math.h>

using namespace std;

/* this is the code to sample from my model and find the vector of the observations.
Sampling from the AR1 model I populate a vector X. This will be the vector of the events.
Sampling form a geometric distribution I populate a vector "vector_ti" for the times of observations 
for the event that happened at time i (event x_i will be observed at time t_i from when it happened).
The matrix "Obs" will have on each row the events that have been observed up to time i
wich is the number of the row. So at time 0 I will have 1 element in the row, 
that might or might not have been observed,
at time 1 I will have two elements on the row, some observed, some not, and so on
I will have 0s whenever the element have have not been yet obeserved at the time
corresponding to the row
the values of the parameters at this stage are fixed and are sigma^2 = 1, 
phi = 0.5, p = 0.4
*/

const double sigmasq = 1;
const float phi = 0.5;
const float p = 0.2;
const double N = 10;
vector < double > X{};
vector < vector < size_t > > obs{};
vector < size_t > vect_obs_N{};
vector < size_t > vect_obs_half{};
vector < size_t > vect_obs_3quarter{};


void print_matrix( vector < vector < size_t > > m );
void print_matrix( vector < vector < double > > M );
void print_vector( vector < double > v);


random_device rd;
mt19937 generator( 0 );


int main(){

  // Sample from a normal distribution. Put the values in a vecotr X
  normal_distribution < double > normalDist( 0, sigmasq / ( 1 - phi * phi ) );
  X.push_back( normalDist ( generator ) );
  
  for ( size_t i = 1; i < N; i++ ){  
    normal_distribution < double > normalDist( phi * X[i - 1], sigmasq );
    X.push_back( normalDist ( generator ) );
  }  

  // Sample from a geometric distribution the values ti
  // and create a vector called "vector_ti".
  vector < double > vector_ti;
  for ( size_t i = 0; i < N; i++ ){
    geometric_distribution <> geoDist(p);
    size_t ti = geoDist ( generator );
    vector_ti.push_back(ti);
  }
  
  
  // Populate the matrix of observations callled obs
  // then create a vector vect_obs_N for the final time
  // vect_obs_quarter, 
  for ( size_t j = 0; j < N; j++ ){
    vector < size_t > tempvec{};
    for ( size_t i = 0; i < j + 1; i++ ){
      if ( vector_ti[i] <= j - i ){
	tempvec.push_back(1);
      }
      else  {
	tempvec.push_back(0);
      }
    }
    obs.push_back( tempvec );
  }
  vect_obs_N = obs [N - 1];

  cout << "printing the matrix of observations" << endl;
  print_matrix(obs);

  //create a vector for half time vect_obs_half
  for ( size_t j = 0; j < N; j++ ){
    vector < size_t > tempvec{};
    for ( size_t i = 0; i < j + 1; i++ ){
      if ( vector_ti[i] <= j - i ){
	tempvec.push_back(1);
      }
      else  {
	tempvec.push_back(0);
      }
    }
    obs.push_back( tempvec );
  }
  vect_obs_half = obs [round( (N - 1)/2 )];

  //create a vector for 3/4 of time vect_obs_3quarter
  for ( size_t j = 0; j < N; j++ ){
    vector < size_t > tempvec{};
    for ( size_t i = 0; i < j + 1; i++ ){
      if ( vector_ti[i] <= j - i ){
	tempvec.push_back(1);
      }
      else  {
	tempvec.push_back(0);
      }
    }
    obs.push_back( tempvec );
  }
  vect_obs_3quarter = obs [ round( 3*(N - 1)/4 )];

  
  
  /*  cout << "vector x \n";
  print_vector(X);
  cout << "vector ti \n";
  print_vector(vector_ti);
  cout << "matrix Obs \n";
  print_matrix_double(Obs); */

  // Create a dat file with the values of the vector of the observed events at time N
  // call it real_data.dat
  ofstream outFile1( "./real_data.dat" );
  outFile1 << endl;
  for ( double n : vect_obs_N ){
    outFile1 << n << endl;
  }
  outFile1.close();

  // Print that vector X on the screen
  cout << "printing the vector X of data \n";
  for ( size_t i = 0; i < N; i++ ){
    cout << X[i] << endl;
    }

  // Create a dat file with the values of X
  std::ofstream outFile2( "./vector_X.dat" );
  outFile2 << endl;
  for ( double n : X ){
    outFile2 << n << endl;
  }
  outFile2.close();
  
  /* this is the code for the method
     i here is the index for the current time that goes from 1 to N, 
     j is the index for the particles that goes from 1 to n. */

  // Here I calculate the 3 dimensional vector called "sample" that will store the new samples
  // for every particle.
  // I also calculate the 3 dimensional vector of the times of observations
  // Tese vectors are stacked from time zero to time N forming a matrix
  // and on the third ax I list the particles  
  
  vector < vector < vector < double > > > sample{};
  vector < vector < vector < size_t > > > sam_obs{};
  vector < vector < vector < size_t > > > new_sam_obs{};
  vector < vector < vector < double > > > new_sample{};
  //define the unnormalised weights 3 dimensional matrix
  vector < vector < vector < double > > > V{};
  //define the normalised weights 3 dimensional matrix
  vector < vector < vector < double > > > W{};
  // define the 3 dimensional matrix Y_W for simulation and weights
  vector < vector < vector < double > > > Y_W{};
  //number of particles
  double n = 5;

  // First I sample from a normal distribution with mean 0
  // and variance sigma^2/(1-phi^2) for every particle
  // and store this in a vector clalled "vector_y0"
  // This vector will be used as the starting point to simulate all other vectors of events
  // that will populate the matrix "sample"
  vector < double > vector_y0{};
  for ( unsigned j = 0; j < n; j++){
    normal_distribution < double > normalDist( 0, sigmasq / ( 1 - phi * phi ) );
    vector_y0.push_back( normalDist ( generator ) );
  }

  /* cout << "vector_y0 \n";
     print_vector(vector_y0); */
  
  // Sampling for every particle from a normal distribution centred in the previous event time phi
  // and with variance sigma^2, find the matrix of events "matrix_sample".
  // Making the substitiution every time I have an observation in real life,
  // find the matrix of updated events "matrix_new_sample"
  for ( size_t j = 0; j < n; j++ ){
    vector < vector < double > > matrix_sample{};
    vector < vector < double > > matrix_new_sample{};
    vector < double > row_matrix_sample{};
    vector < double > row_matrix_new_sample{};
    double y{};
    y = vector_y0[j];
    row_matrix_sample.push_back(y);
    row_matrix_new_sample.push_back(y);
    matrix_sample.push_back( row_matrix_sample );
    matrix_new_sample.push_back( row_matrix_new_sample );
    for ( size_t i = 1; i < N; i++ ){
      vector < size_t > row_obs{};
      row_obs = obs[i];
      normal_distribution < double > normalDist( phi * row_matrix_new_sample[i - 1], sigmasq );
      double gen = normalDist ( generator );
      row_matrix_new_sample.push_back( gen );
      row_matrix_sample.push_back( gen );
      for (size_t k = 0; k < i; k++){
	if (row_obs[k] != 0 ){
	  row_matrix_sample[k] = row_matrix_new_sample[k];
	}
      }
      for (size_t k = 0; k < i + 1; k++){
	if (row_obs[k] != 0 ){
	  row_matrix_new_sample[k] = X[k];
	}
      }
      matrix_sample.push_back( row_matrix_sample );
      matrix_new_sample.push_back( row_matrix_new_sample );
    }
    row_matrix_sample.clear();
    row_matrix_new_sample.clear();
    sample.push_back( matrix_sample );
    new_sample.push_back( matrix_new_sample );
    matrix_sample.clear();
    matrix_new_sample.clear();
  }

  // Sampling for every particle from a bernoulli distribution with probability p
  // find the matrix of observations "matrix_obs"
  // Making the substitiution every time I have an observation in real life,
  // find the matrix of updated observations "matrix_new_obs"
  for ( size_t j = 0; j < n; j++ ){
    vector < vector < size_t > > matrix_obs{};
    vector < vector < size_t > > matrix_new_obs{};
    vector < size_t > row_matrix_obs{};
    vector < size_t > row_matrix_new_obs{};
    for ( size_t i = 0; i < N; i++ ){
      vector < size_t > row_obs{};
      row_obs = obs[i];
      bernoulli_distribution BerDist(p);
      double gen = BerDist ( generator );
      row_matrix_new_obs.push_back( gen );
      row_matrix_obs.push_back( gen );
      for (size_t k = 0; k < i; k++){
	if (row_obs[k] != 0 ){
	  row_matrix_obs[k] = row_matrix_new_obs[k];
	}
      }
      for (size_t k = 0; k < i + 1; k++){
	if (row_obs[k] != 0 ){
	  row_matrix_new_obs[k] = 1;
	}
      }
      matrix_obs.push_back( row_matrix_obs );
      matrix_new_obs.push_back( row_matrix_new_obs );
    }
    row_matrix_obs.clear();
    row_matrix_new_obs.clear();
    sam_obs.push_back( matrix_obs );
    new_sam_obs.push_back( matrix_new_obs );
    matrix_obs.clear();
    matrix_new_obs.clear();
  }
  



  // Print one matrix of the 3 dimensional vector "sample"
  vector < vector < double > > printer1{};
  printer1 = sample[0];
  cout << "printing one of the sample matrix \n";
  print_matrix( printer1 );
  // Print one matrix the 3 dimensional vector "new_sample"
  vector < vector < double > > printer2{};
  printer2 = new_sample[0];
  cout << "printing one of the new_sample matrix \n";
  print_matrix( printer2 );

  // Print one matrix of the 3 dimensional vector "sam_obs"
  vector < vector < size_t > > printer3{};
  printer3 = sam_obs[0];
  cout << "printing one of the sam_obs matrix \n";
  print_matrix( printer3 );
  // Print one matrix the 3 dimensional vector "new_sam_obs"
  vector < vector < size_t > > printer4{};
  printer4 = new_sam_obs[0];
  cout << "printing one of the new_sam_obs matrix \n";
  print_matrix( printer4 );
  
  // find the matrix Y_N of new sample for the last time
  vector < vector < double > > Y_N{};
  for ( size_t i = 0; i < n; i++ ){
    vector < vector < double > > temp_matrix{};
    temp_matrix = new_sample[i];
    vector < double > temp_vector;
    for ( size_t j = 0; j < N; j++ ){
      temp_vector = temp_matrix[N-1];
    }
    Y_N.push_back(temp_vector);
  }
  
  cout << " printing matrix Y_N " << endl;
  print_matrix(Y_N);
  
  // Creating a 3 dimentional matrix "W" of zeroes for the normalised weights
  vector < vector < double > > matrix_w{};
  vector < double > vector_w{};
  for ( size_t j = 0; j < n; j++ ){
    double elem{};
    for ( size_t i = 0; i < N; i++ ){
      for (size_t k = 0; k < i + 1; k++ ){
	elem = 0;  
	vector_w.push_back( elem );
      }   
      matrix_w.push_back( vector_w );
      vector_w.clear();
    }
    W.push_back( matrix_w );
    matrix_w.clear();
  }

  // Creating a 3 dimentional matrix "Y_W" of zeroes
  vector < vector < double > > matrix_y_w{};
  vector < double > vector_y_w{};
  for ( size_t j = 0; j < n; j++ ){
    double elem{};
    for ( size_t i = 0; i < N; i++ ){
      for (size_t k = 0; k < i + 1; k++ ){
	elem = 0;  
	vector_y_w.push_back( elem );
      }   
      matrix_y_w.push_back( vector_y_w );
      vector_y_w.clear();
    }
    Y_W.push_back( matrix_y_w );
    matrix_y_w.clear();
  }

  
  //Finding the unnormalised weights (using log then exponentiating)
  const double constant = ( 1 / ( 2 * sigmasq ) );
  for ( size_t j = 0; j < n; j++ ){
    vector < vector < double > > matrix_V{0};
    vector < double > temp_vec{1};
    matrix_V.push_back( temp_vec );
    vector < vector < double > > temp_matrix_sample{};
    temp_matrix_sample = sample[j];
    vector < vector < double > > temp_matrix_new_sample{};
    temp_matrix_new_sample = new_sample[j];
    vector < double > vector_V{};
    double log_weight{};
    for ( size_t i = 1; i < N; i++ ){
      vector_V.push_back(1);
      vector < double > vector_log_weights{};
      vector_log_weights.push_back(1);
      vector < double > row_sample{};
      row_sample = temp_matrix_sample[i];
      vector < double > row_new_sample{};
      row_new_sample = temp_matrix_new_sample[i];
      for ( size_t k = 1; k < row_sample.size(); k++){
	if ( row_new_sample[k-1] == row_sample[k-1] && row_new_sample[k] == row_sample[k] ){}
	else {
	  double ys = - (row_new_sample[k] - phi * row_new_sample[k-1]) * (row_new_sample[k] - phi * row_new_sample[k-1]);
	  double xs = (row_sample[k] - phi * row_sample[k-1]) * (row_sample[k] - phi * row_sample[k-1]);
	  log_weight = constant * ( ys + xs );
	  vector_log_weights.push_back(log_weight);
	}
      double v{};
      double sum = accumulate(vector_log_weights.begin(), vector_log_weights.end(), 0.0);
      v = exp(sum);
      vector_V.push_back(v);
      }
    matrix_V.push_back(vector_V);
    vector_V.clear();
    }
    V.push_back( matrix_V );
  }

  
  /* // Print one matrix of the 3 dimensional vector V
  vector < vector < double > > printer3{};
  printer3 = V[0];
  cout << "printing one of the V matrix \n";
  print_matrix_double( printer3 ); */

  // normalise the importance weights and put it in matrix W
  for ( size_t i = 0; i < N; i++ ){
    for (size_t k = 0; k < i + 1; k++){
      double sum{0};
      for (size_t l = 0; l < n; l++ ){
	sum += V[l][i][k];
      }
      for (size_t j = 0; j < n; j++ ){
	W[j][i][k] = V[j][i][k] / sum;
      }
    }
  }
  /*
  // Resample
  double tresh = n/2;
  for ( size_t i = 0; i < N; i++ ){
    for (size_t k = 0; k < i + 1; k++){
      double ess{0};
      double sumsq{0};
      for (size_t l = 0; l < n; l++ ){
	// calculate the effective sample size, called ess
	sum += W[l][i][k] *  W[l][i][k];
	ess = 1 / sumsq;
      }
      // draw particles from the current particle set with probabilities
      // proportional to their weights and replace the current particles
      for (size_t j = 0; j < n; j++ ){
	if (ess < tresh){
	  W[j][i][k];
	}
      }
    }
  }
      
  */

  // Find the matrix of the weights for the last time
  // and call it W_N
  vector < vector < double > > W_N{};
  for ( size_t i = 0; i < n; i++ ){
    vector < vector < double > > temp_matrix{};
    temp_matrix = W[i];
    vector < double > temp_vector;
    for ( size_t j = 0; j < N; j++ ){
      temp_vector = temp_matrix[N-1];
    }
    W_N.push_back(temp_vector);
  }

  cout << " printing matrix W_N " << endl;
  print_matrix(W_N);
  

  // calculate the 3 dimensional matrix Y_W for simulation and weights
  for ( size_t j = 0; j < n; j++ ){
    for ( size_t i = 0; i < N; i++ ){
      for ( size_t k = 0; k < i + 1; k++ ){
 	Y_W[j][i][k] = new_sample[j][i][k] * W[j][i][k];
      }
    }
  }  
  
  //calculate the expectation E[y_i] for the last observation
  vector < double > E{};
  vector < vector < double > > Y_W_matrix{};
  for ( size_t j = 0; j < n; j++ ){
    Y_W_matrix.push_back( Y_W[j][N - 1] );
  }
  for ( size_t i = 0; i < N; i++ ){
    vector < double > multiplication_vector{};
    for ( size_t j = 0; j < n; j++ ){
      multiplication_vector.push_back( Y_W_matrix[j][i] );
    }
    double sum{};
    for ( auto & n : multiplication_vector)
      sum += n;
    E.push_back( sum );
  }
  
  // Create a dat file with the values of E
  ofstream outFile3( "./vector_E.dat" );
  outFile3 << endl;
  for ( double n : E ){
    outFile3 << n << endl;
  }
  outFile3.close();

  // Create a dat file with the values of the Y_N
  // to craete boxplots
  ofstream outFile4( "./Y_N.dat" );
  outFile4 << endl;
  for ( size_t lin = 0; lin < n; lin++ ){
    for ( size_t col = 0; col < N; col++ ){
      outFile4 << Y_N[lin][col] << " ";
    }
    outFile4 << endl;
  }
  outFile4.close();

  // Create a dat file with the values of the W_N
  // to craete boxplots
  ofstream outFile5( "./W_N.dat" );
  outFile5 << endl;
  for ( size_t lin = 0; lin < n; lin++ ){
    for ( size_t col = 0; col < N; col++ ){
      outFile5 << W_N[lin][col] << " ";
    }
    outFile5 << endl;
  }
  outFile5.close();
  
  return 0;
}




// functions definitions

// this function prints a matrix of unsigned size_t
void print_matrix( vector < vector < size_t > > m ){
  for ( const vector < size_t > & v : m ){
    for  ( size_t x : v ) cout << x << ' ';
    cout << endl;
  }
}

// this function prints a matrix of signed
void print_matrix( vector < vector < double > > M ){
  for ( const vector < double > & v : M ){
    for  ( double x : v ) cout << x << ' ';
    cout << endl;
  }
}

// this function prints a vector
void print_vector( vector < double > v){
  for ( size_t i = 0; i < v.size(); i++){
    cout << v[i] << ' ';
  }
  cout << endl;
}
