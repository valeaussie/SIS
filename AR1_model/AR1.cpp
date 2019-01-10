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

const int sigmasq = 1;
const float phi = 0.5;
const float p = 0.4;
const double N = 5;
vector < double > X{};
vector < vector < double > > Obs{};

void print_matrix_sizet( vector < vector < size_t > > m );
void print_matrix_double( vector < vector < double > > M );
void transpose( vector < vector < size_t > > &b );
void print_vector( vector < double > v);


random_device rd;
mt19937 generator( rd () );


int main(){

  // Sample from a normal distribution. Put the values in a vecotr X.
  
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
  
  
  // Populate the matrix of observations callled Obs 
  for ( size_t j = 0; j < N; j++ ){
    vector < double > tempvec{};
    for ( size_t i = 0; i < j + 1; i++ ){
      if ( vector_ti[i] <= j - i ){
	tempvec.push_back( X[i]);
      }
      else  {
	tempvec.push_back(0);
      }
    }
    Obs.push_back( tempvec );
  }
  
  cout << "vector x \n";
  print_vector(X);
  cout << "vector ti \n";
  print_vector(vector_ti);
  cout << "matrix Obs \n";
  print_matrix_double(Obs);  
  
  /* this is the code for the method
     i here is the index for the current time that goes from 1 to N, 
     j is the index for the particles that goes from 1 to n. */

  // Here I calculate the 3 dimensional vector called "sample" that will store the new samples
  // for every particle.
  // I also calculate the 3 dimensional vector of the times of observations
  // Tese vectors are stacked from time zero to time N forming a matrix
  // and on the third ax I list the particles  
  
  vector < vector < vector < double > > > sample{};
  vector < vector < vector < size_t > > > times{};
  vector < vector < vector < double > > > new_sample{};
  // vector < vector < double > > V{};
  // vector < vector < double > > W{};
  double n =3;

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
  
  // Sampling from a geometric distribution with p = 0.4
  // find the matrix_ti of the times of observations for each particle
  // Sampling for every particle from a normal distribution centred in the previous event time phi
  // and with variance sigma^2, I find the matrix of events "matrix_y"
  vector < vector < size_t > > matrix_ti{};
  vector < vector < double > > matrix_y{};
  for ( size_t j = 0; j < n; j++ ){
    vector < size_t > temp_ti{};
    vector < double > temp_y{};
    double y{};
    y = vector_y0[j];
    temp_y.push_back(y);
    for (size_t i = 1; i < N + 1; i++ ){
      normal_distribution < double > normalDist( phi * temp_y[i - 1], sigmasq );
      temp_y.push_back( normalDist ( generator ) );
      geometric_distribution < size_t > geoDist(p);
      size_t ti = geoDist ( generator );
      temp_ti.push_back( ti );
    }
    temp_y.pop_back();
    matrix_ti.push_back( temp_ti );
    matrix_y.push_back( temp_y );
    }

  // Print matrix_y and matrix_ti
  cout << "matrix_y \n";
  print_matrix_double( matrix_y );
  cout << "matrix_t \n";
  print_matrix_sizet( matrix_ti );

  // create the 3 dimensional vectors for observations and events  
  for ( size_t k = 0; k < n; k++){
    vector < double > row_matrix_y{};
    vector < size_t > row_matrix_ti{}; 
    row_matrix_y = matrix_y[k];
    row_matrix_ti = matrix_ti[k];
    vector < vector < double > > sim_matrix{};
    for ( size_t j = 0; j < N; j++ ){
      vector < double > tempvec{};
      for ( size_t i = 0; i < j + 1; i++ ){
      	if ( row_matrix_ti[i] <= j - i ){
	  tempvec.push_back( row_matrix_y[i]);
	}
	else {
	  tempvec.push_back(0);
	}
      }
      sim_matrix.push_back( tempvec );
    }
    sample.push_back( sim_matrix );
  }

  // Print one matrix of the 3 dimensional vector "sample"
  vector < vector < double > > printer1{};
  printer1 = sample[1];
  cout << "printing one of the sample matrix \n";
  print_matrix_double( printer1 );

  // make the substitiution and find the 3 dimensional vector "new_sample"
  
  for ( size_t k = 0; k < n; k++ ){
    vector < vector < double > > matrix_sample{};
    matrix_sample = sample[k];
    vector < vector < double > > new_matrix{};
    for ( size_t i = 0; i < N; i++ ){
      vector < double > row_matrix_sample{};
      vector < double > row_matrix_Obs{};
      row_matrix_sample = matrix_sample[i];
      row_matrix_Obs = Obs[i];
      for ( size_t j = 0; j < i + 1; j++ ){
	if ( row_matrix_Obs[j] != 0 ){
	  row_matrix_sample[j] = row_matrix_Obs[j];
	}
      }
      new_matrix.push_back( row_matrix_sample );
    }
    new_sample.push_back( new_matrix );
  }

  // Print one matrix the 3 dimensional vector "new_sample"
  vector < vector < double > > printer2{};
  printer2 = new_sample[1];
  cout << "printing one of the new_sample matrix \n";
  print_matrix_double( printer2 );
  
  // Find the weights


  /* If i for vector L is 0 and the element S[j][i] is 1 
     (no observation in either real life or simulation) 
     the importance weights will be w_(i_j) = w_[(i - 1)_j] */
  /* If i for vector L is 1 and the element S[j][i] is 0
     (observation in real life, no observation simulated), 
     make the substitutions putting the elements of the vector X_(k_i) in column i of Y, 
     then calculate the weights w_i^(j) = w_[(i - 1)_j] * (1 - p)^(t_i + i - N - 1) */
  /* If i for vector L is 1 and the element S[j][i] is 0
     (observation in real life, no observation simulated)
     make the substitutions putting the elements of the vector X_(k_i) in column i of Y,
     then calculate the weights w_(i_j) = w_[(i - 1)_j] * (1 - p)^(N - t_i - i + 1) */
  /* If i for vector L is 1 and the element S[j][i] is 1
     (observation in real life and in the simulation) 
     the importance weights will be w_(i_j) = w_[(i - 1)_j]*/
  /*
  for ( size_t i = 1; i < N; i++ ){
    vector < double > w{};
    
    bool contains_value {false};
    for ( size_t k = 0; k < L.size(); k++){
      if ( L[k] == i ){
	contains_value = true;
	break;
      }
    }
    for ( size_t j = 0; j < n; j++ ){
      if ( ( S[j][i] == 0 ) && ( contains_value == false ) ){
	V[j][i] = V[j][i - 1];
      }
      
      else if ( ( S[j][i] == 1 ) && ( contains_value == true ) ){
	
	// Calculate the power
	double pow1{};
	if ( ( vector_ti[i] + i - N - 1 ) == 0){
	  pow1 = 1;
	}
	else if ( (vector_ti[i] + i - N - 1) == 1){
	  pow1 = (1 - p);
	}
	else {
	  pow1 = pow ( (1 - p), ( vector_ti[i] + i - N - 1 ) );
	}
	
	V[j][i] = V[j][i - 1] * pow1 ;
      }
      
      else if ( ( S[j][i] == 1 ) && ( contains_value == false ) ){
	// Calculate the power
	double pow2{};
	if ( (N - vector_ti[i] - i + 1) == 0){
	  pow2 = 1;
	}
	else if ( (N - vector_ti[i] - i + 1) == 1){
	  pow2 = (1 - p);
	}
	else {
	  pow2 = pow ( (1 - p), ( N - vector_ti[i] - i + 1 ) );
	}
	
	V[j][i] = V[j][i - 1] * pow2 ;
      }
      else if ( ( S[j][i] == 0 ) && ( contains_value == true ) ){
	V[j][i] = V[j][i - 1];
      }
    }
  }
  // normalise the importance weights and put it in matrix W
  for ( size_t i = 1; i < N; i++ ){
    vector < double > column_vector{};
    double sum{};
    for ( size_t k = 0; k < n; k++ ){
      column_vector.push_back(V[k][i]);
    }
    for ( auto & n : column_vector)
      sum += n;
    for ( size_t j = 0; j < n; j++ ){
    W[j][i] = V[j][i] / sum;
    }
  }
  
  // calculate the expectation E[y_i] = sum for j from 2 to N y_(i_j) W_(i_j)
  vector < double > E{};
  for ( size_t i = 0; i < N; i++ ){
    vector < double > moltiplication_vector{};
    for ( size_t j = 0; j < n; j++ ){
      moltiplication_vector.push_back( Y[j][i] * W[j][i] );
    }
    double sum{};
    for ( auto & n : moltiplication_vector)
      sum += n;
    E.push_back( sum );
  }
  // Create a dat file with the values of E, this is useful to graph with gnuplot
    std::ofstream outFile2( "./vector_E.dat" );
  outFile2 << endl;
  for ( double n : E ){
    outFile2 << n << endl;
  }
  outFile2.close();
  for ( size_t i = 0; i < N; i++ ){
    cout << E[i] << endl;
  }
  */
  
  return 0;
}




// functions definitions

// this function prints a matrix of unsigned size_t
void print_matrix_sizet( vector < vector < size_t > > m ){
  for ( const vector < size_t > & v : m ){
    for  ( size_t x : v ) cout << x << ' ';
    cout << endl;
  }
}

// this function prints a matrix of signed
void print_matrix_double( vector < vector < double > > M ){
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
  cout << "\n";
}
