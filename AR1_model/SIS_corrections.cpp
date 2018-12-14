#include <random>
#include <iostream>
#include <vector>
#include <fstream>
#include <map>
#include <algorithm>
#include <math.h>
#include "SIS_AR1.h"

using namespace std;

void print_matrix_sizet( vector < vector < size_t > > m );
void print_matrix_double( vector < vector < double > > M );
void transpose( vector < vector < size_t > > &b );
int ar1();


random_device rd;
mt19937 generator( rd () );


int main(){
  
  vector < vector < double > > sample{};
  vector < vector < double > > new_sample{};
  vector < vector < double > > V{};
  vector < vector < double > > W{};
  double n =3;

  // This is the matrix S of the observations
  // The columns are indexd as j and are for  the particles (n)
  // The rows are indexed as i and are for the number of events (N)
  // It is a matrix of 0s (observation) and 1s (no observation)

  // This is the calculations matrix Y of the simulated events.
  // In the first raw of the matrix we have the values for f(x_1) for all the j particles

  // Creating a vector y simulating from a normal distribution with mean 0
  // and variance sigma^2/(1-phi^2)
  // This vector will be used as a starting point to simulate all other n vectors y of lenght n
  // that will populate the matrix Y
  
  vector < double > y{};
  for ( unsigned j = 0; j < n; j++){
    normal_distribution < double > normalDist( 0, sigmasq / ( 1 - phi * phi ) );
    y.push_back( normalDist ( generator ) );
  }
  sample.push_back(y);
  new_sample.push_back(y);

  // This is the vector v of the unnormalised weights
  vector < double > v{};
  for ( size_t j = 0; j < n; j++) {
    v.push_back(0);
  }
  V.push_back(v);
  v.clear();

  // This is the vector w of the normalised weights
  vector < double > w{};
  for ( size_t j = 0; j < n; j++) {
    w.push_back(1 / n);
  }
  W.push_back(w);
  w.clear();
  
  // Create a vector L with the first column of Z
  // These are the observations in real life
  vector < size_t > L;
  for ( size_t i = 0; i < Z.size(); i++ ){
    L.push_back(Z[i][0]);
  }
  
  for ( size_t i = 0; i < n; i++ ){
  }
    
  
  // Sampling from a geometric distribution with p = 0.4
  // find t_(i_j) for particle j
  // (time of the next observation of the event that happened at time i).
  // If t_(i_j) > N - i multiply y_i by (1 - p) and find the new y_i in the vectors Y_j.
  // Else multiply y_i by p(1 - p) and find the new y_i.
  // and in the vector S_(i_j) in position k_i = i + t_i, substitute the existing value with a 1.

 
  vector < vector < size_t > > matrix_ti{};
  //size_t ki(n);
  for ( size_t i = 1; i < N; i++ ){
    vector < size_t > vector_ti{};
    vector < size_t > vector_i{};
    vector < double > temp_y{};
 
    for ( size_t j = 0; j < n; j++ ){
      normal_distribution < double > normalDist( phi * y[i - 1], sigmasq );
      temp_y.push_back( normalDist ( generator ) );
      geometric_distribution < size_t > geoDist(p);
      size_t ti = geoDist ( generator );

      if ( ti >= N - i ){
	temp_y[j] = temp_y[j] * (1 - p);
      }

      else {
	temp_y[j] = temp_y[j] * p * (1 - p);
      }
    vector_ti.push_back(ti);

    }
    
    matrix_ti.push_back(vector_ti);
    sample.push_back(temp_y);
    new_sample.push_back(temp_y);
    y.clear();
    y = temp_y;

    // Make substitiution
    for ( size_t j = 0; j < n; j++){
      bool contains_value {false};
      for ( size_t k = 0; k < L.size(); k++){
	if ( L[k] == i ){
	  contains_value = true;
    	  break;
   	}
      }
      if (contains_value == true ){
    	new_sample[i][j] = X[i];
      }
    }
  }
  
  vector < vector < size_t > > S (n, vector < size_t > (N - 1) );
  transpose(matrix_ti);
  
  for ( size_t i = 0; i < matrix_ti.size(); i++ ){
    for ( size_t j = 0; j < matrix_ti[i].size(); j++ ){
      if ( matrix_ti[i][j] + j >= matrix_ti[i].size() ){}
      else if ( matrix_ti[i][j] == 0 ){}
      else {
	S[i][ j + matrix_ti[i][j] ] = 1;
      }
    }
  }

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



