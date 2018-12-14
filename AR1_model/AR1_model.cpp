#include <random>
#include <iostream>
#include <vector>
#include <fstream>
#include <map>
#include <algorithm>
#include <math.h>
#include "SIS_AR1.h"

using namespace std;


/* this is the code to sample from my model and find the vector of the observations
I am going to generate a vector X with the x_i as elements which are found sampling 
from the AR1 model I have specified. This will be the vector of the events.
Then I am going to generate N, R_(k_i) vectors of the times of the events 
that happened at time i and have been observed at time k_i.
Finally, I will generate N, Z_(k_i) vectors of observations of events 
that happened at time i and have been observed at time k_i.
the values of the parameters at this stage are fixed and are sigma^2 = 1, 
phi = 0.5, p = 0.4, I also choose the value of N = 1000
*/

const int sigmasq = 1;
const float phi = 0.5;
const float p = 0.4;
const double N = 5;


void print_matrix_sizet( vector < vector < size_t > > m );
void print_matrix_double( vector < vector < double > > M );
void transpose( vector < vector < size_t > > &b );
vector < double > X;
vector < vector < double > > Z;


int ar1(){
  

  random_device rd;
  mt19937 generator( rd () );


  // Sample from a normal distribution. Put the values in a vecotr X.
  
  normal_distribution < double > normalDist( 0, sigmasq / ( 1 - phi * phi ) );
  X.push_back( normalDist ( generator ) );
  
  for ( size_t i = 1; i < N; i++ ){  
    normal_distribution < double > normalDist( phi * X[i - 1], sigmasq );
    X.push_back( normalDist ( generator ) );
  }
  
  // Create a dat file with the values of X, this is useful to graph with gnuplot
  std::ofstream outFile( "./vector_X.dat" );
  outFile  << endl;
  for ( double n : X ){
    outFile << n << endl;
  }
  outFile.close();
  

  vector < vector < size_t > > r{};
  
  // Sample from a geometric distribution the values ti
  for ( size_t i = 0; i < N; i++ ){
    vector < size_t > vector_ti{};
    vector < size_t > vector_Ki{};
    vector < size_t > vector_i{};
    geometric_distribution <> geoDist(p);
    size_t ti = geoDist ( generator );
    vector_ti.push_back(ti);
    // Find Ki as ti + i.
    //Create the RKi vectors and put Ki as first value then i.
    //Store all this vectors in a vector of vectors (matrix) r
    if ( ti >= N - i ){}
    else {
      size_t ki = i + ti;
      vector_i.push_back(i);
      vector_Ki.push_back(ki);
      vector < size_t > Rki;
      Rki.push_back(ki);
      Rki.push_back(i);
      r.push_back(Rki);
    }
  }
  
  // If two lines of the matrix r have the same ki (first element)
  // add the second element of the second vector to the first vector
  // and set the second vector to be all made of 0
  for ( size_t i = 0; i < r.size(); i++ ){
    for ( size_t j = i + 1 ; j < r.size(); j++ ){
      if ( r[j][0] == r[i][0] ){
	r[i].push_back( r[j][1] );
	r[j][0] = 0;
	r[j][1] = 0;
      }
    }
  }
  
  // Delete all zeroes from r and call it R
  vector < vector < size_t > > R{};
  for ( size_t i = 0; i < r.size(); i++ ){
    if ( r[i][0] != 0 ){
      R.push_back( r[i] );
    }
  }
  
  // Sort the matrix R by the first column
  sort( R.begin(), R.end() );
  
  // Create a matrix z stacking all the times of observations up to time ki
  // (specified in the fisrt column). The times are also sorted.
  vector < vector < size_t > > z{};
  vector < size_t > temp_vector{};
  for ( size_t i = 0; i < R.size(); i++ ){
    for ( size_t j = 1; j < R[i].size(); j++ ){
      temp_vector.push_back( R[i][j] );
      sort( temp_vector.begin(), temp_vector.end() );
    }
    z.push_back( temp_vector);
  }
  for ( size_t i = 0; i < R.size(); i++ ){
    z[i].insert(z[i].begin(), R[i][0]);
  }
  
  // Create a matrix Z stacking all the values of the observations up to time ki
  vector < vector < double > > Z{};
  for ( size_t i = 0; i < z.size(); i++ ){
    vector < double > temp_vect_double{};
    for ( size_t j = 1; j < z[i].size(); j++ ){
      temp_vect_double.push_back( X[ ( z[i][j] ) ] );
    }
    Z.push_back( temp_vect_double );
    temp_vect_double.clear();
  }
  for ( size_t i = 0; i < z.size(); i++ ){
    Z[i].insert( Z[i].begin(), z[i][0] );
  }

  
  return 0;
}


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

// this function transposes a matrix
void transpose( vector < vector < size_t > > &b ){
  if (b.size() == 0 )
    return;
  vector < vector < size_t > > trans_vec(b[0].size(), vector < size_t > () );
  for ( size_t i = 0; i < b.size(); i++ ){ 
    for ( size_t j = 0; j < b[i].size(); j++ ){
      trans_vec[j].push_back(b[i][j]);
    }
  }
  b = trans_vec;
}


