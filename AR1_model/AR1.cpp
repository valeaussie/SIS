#include <random>
#include <iostream>
#include <vector>
#include <fstream>
#include <map>
#include <algorithm>
#include <math.h>

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
const double N = 10;
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
  // and create a vector of 0s and 1s called vector_Ki.
  // We will have a 0 when ti >= N - i and 1 otherwise
  vector < double > vector_ti;
  vector < size_t > vector_Ki{};
  for ( size_t i = 0; i < N; i++ ){
    geometric_distribution <> geoDist(p);
    size_t ti = geoDist ( generator );
    vector_ti.push_back(ti);
    if ( ti >= N - i ){
      vector_Ki.push_back(0);
    }
    else {
      vector_Ki.push_back(1);
    }
  }
  vector_Ki[0] = 0;

  
  
  // Populate the matrix of observations callled Obs
 
  for ( size_t j = 0; j < N + 1; j++ ){
    vector < double > tempvec{};
    for ( size_t i = 0; i < j; i++ ){
      if ( vector_ti[i] <= j - i - 1 ){
	tempvec.push_back( X[i]);
      }
      else  {
	tempvec.push_back(0);
      }
    }
    Obs.push_back( tempvec );
  }
  
  cout << " vector x \n ";
  print_vector(X);
  cout << " vector ti \n ";
  print_vector(vector_ti);
  cout << " matrix Obs \n ";
  print_matrix_double(Obs);
  

  // Create a dat file with the values of X, this is useful to graph with gnuplot
  //std::ofstream outFile( "./vector_Obs.dat" );
  //outFile  << endl;
  //for ( double n : X ){
  //  outFile << n << endl;
  //}
  //outFile.close();
  
  /* this is the code for the method
     i here is the index for the current time that goes from 1 to N, 
     j is the index for the particles that goes from 1 to n. I choose n to be 10 */
  
  
  vector < vector < double > > sample{};
  vector < vector < double > > new_sample{};
  vector < vector < double > > V{};
  vector < vector < double > > W{};
  double n =3;
  
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

    /* Make substitiution
    for ( size_t j = 0; j < n; j++){
      bool contains_value {false};
      for ( size_t k = 0; k < Obs.size(); k++){
	if ( L[k] == i ){
	  contains_value = true;
    	  break;
   	}
      }
      if (contains_value == true ){
    	new_sample[i][j] = X[i];
	}
	} */
    }
  
  vector < vector < size_t > > S (n, vector < size_t > (N - 1) );

  transpose(matrix_ti);
  
  //for ( size_t i = 0; i < matrix_ti.size(); i++ ){
  //  for ( size_t j = 0; j < matrix_ti[i].size(); j++ ){
  //    if ( matrix_ti[i][j] + j >= matrix_ti[i].size() ){
  //	cout << matrix_ti[i][j] << endl;
  //    }
  //    else if ( matrix_ti[i][j] == 0 ){}
  //    else {
  //	cout << matrix_ti[i][j] << endl;
  //	S[i][ j + matrix_ti[i][j] ] = 1;
  //    }
  //  }
  //}
  
  
  //print_matrix_sizet(S);
  transpose(S);
  

  vector < vector < size_t > > col_matrix_S{};
  //for ( size_t j = 0; j < n; j++ ){
  for ( size_t i = 0; i < N - 1; i++ ){
    //col_S[ i + col_ki[ i ] - 1 ] = 1;
    //cout << col_ti[ i ] << endl;
  }
    // col_matrix_S.push_back(col_S);
    // }
  
  //cout << " matrix transposed S " << endl;
  //print_matrix_sizet(S);


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

// this function transpose a matrix
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

// this function prints a vector
void print_vector( vector < double > v){
  for ( size_t i = 0; i < v.size(); i++){
    cout << v[i] << ' ';
  }
  cout << " \n ";
}
