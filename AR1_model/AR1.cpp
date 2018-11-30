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
const double N = 5;
vector < double > X{};

void print_matrix_sizet( vector < vector < size_t > > m );
void print_matrix_double( vector < vector < double > > M );

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
  
  // Create a dat file with the values of X, this is useful to graph with gnuplot
  std::ofstream outFile( "./vector_X.dat" );
  outFile  << endl;
  for ( double n : X ){
    outFile << n << endl;
  }
  outFile.close();
   
  vector < size_t > vector_ti{};
  vector < size_t > vector_Ki{};
  vector < size_t > vector_i{};
  vector < vector < size_t > > r{};

  // Sample from a geometric distribution the values ti
    for ( size_t i = 0; i < N; i++ ){
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

  /* this is the code for the method
     i here is the index for the current time that goes from 1 to N, 
     j is the index for the particles that goes from 1 to n. I choose n to be 10 */

  vector < vector < size_t > > S{};
  vector < vector < double > > sample{};
  vector < vector < double > > new_sample{};
  vector < vector < double > > V{};
  vector < vector < double > > W{};
  double n =3;

  // This is the matrix S of the observations
  // The columns are indexd as j and are for  the particles (n)
  // The rows are indexed as i and are for the number of events (N)
  // It is a matrix of 0s (observation) and 1s (no observation)
  /*   for ( size_t i = 0; i < N; i++ ){
    vector < size_t > s{};
    for ( size_t j = 0; j < n; j++ ){
      s.push_back(0);
    }
    S.push_back(s);
  }

    print_matrix_sizet(S); */

  // This is the calculations matrix Y of the simulated events.
  // In the first raw of the matrix we have the values for f(x_1) for all the j particles

  // Creating a vector y simulating from a normal distribution with mean 0 and variance sigma^2/(1-phi^2)
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
    cout << y[i] << endl;
  }
    
  
  // Sampling from a geometric distribution with p = 0.4
  // find t_(i_j) for particle j (time of the next observation of the event that happened at time i).
  // If t_(i_j) > N - i multiply y_i by (1 - p) and find the new y_i in the vectors Y_j.
  // Else multiply y_i by p(1 - p) and find the new y_i.
  // and in the vector S_(i_j) in position k_i = i + t_i, substitute the existing value with a 1.

  vector < vector < size_t > > matrix_ki{}; 
  for ( size_t i = 1; i < N; i++ ){
    vector < double > temp_y{};
    vector < size_t > vec_ki{};
 
    for ( size_t j = 0; j < n; j++ ){
      normal_distribution < double > normalDist( phi * y[i - 1], sigmasq );
      temp_y.push_back( normalDist ( generator ) );
      geometric_distribution < size_t > geoDist(p);
      size_t ti = geoDist ( generator );
      cout << " ti " << ti << endl;

      if ( ti >= N - i ){
	temp_y[j] = temp_y[j] * (1 - p);
	vec_ki.clear();
      }

      else {
	temp_y[j] = temp_y[j] * p * (1 - p);
      size_t ki{};
      ki = ti + i;;
      cout << " ki " << ki << endl;
      vec_ki.push_back(ki);
      }
    }

    matrix_ki.push_back(vec_ki);
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
  
  print_matrix_double(sample);
  cout << " new sample " << endl;
  print_matrix_double(new_sample);
  cout << " matrix Ki " << endl;
  print_matrix_sizet(matrix_ki);


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






