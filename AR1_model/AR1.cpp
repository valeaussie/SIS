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
vector < double > X;

void print_matrix_unsigned( vector < vector < unsigned > > m );
void print_matrix_double( vector < vector < double > > M );

random_device rd;

int main(){

  // Sample from a normal distribution. Put the values in a vecotr X.
  mt19937 generator( rd () );
  normal_distribution < double > normalDist( 0, sigmasq / ( 1 - phi * phi ) );
  X.push_back( normalDist ( generator ) );
  
  for ( unsigned i = 1; i < N; i++ ){
    normal_distribution < double > normalDist( phi * X[i - 1], sigmasq );
    X.push_back( normalDist ( generator ) );
  }

  cout << X.size() << endl;
  
  // Create a dat file with the values of X, this is useful to graph with gnuplot
    std::ofstream outFile( "./vector_X.dat" );
  outFile << "values of X" << endl;
  for ( double n : X ){
    outFile << n << endl;
  }
  outFile.close();
   
  vector < unsigned > vector_ti{};
  vector < unsigned > vector_Ki{};
  vector < unsigned > vector_i{};
  vector < vector < unsigned > > r{};

  // Sample from a geometric distribution the values ti
    for ( unsigned i = 0; i < N; i++ ){
    geometric_distribution <> geoDist(p);
    unsigned ti = geoDist ( generator );
    vector_ti.push_back(ti);
    cout << " vector ti " << ti << endl;
    // Find Ki as ti + i.
    //Create the RKi vectors and put Ki as first value then i.
    //Store all this vectors in a vector of vectors (matrix) r
    if ( ti >= N - i ){}
    else {
      unsigned ki = i + ti;
      vector_i.push_back(i);
      vector_Ki.push_back(ki);
      vector < unsigned > Rki;
      Rki.push_back(ki);
      Rki.push_back(i);
      r.push_back(Rki);
    }
  }

      
  // Print out the matrix r
  cout << "Matrix r" << endl;
  print_matrix_unsigned(r);
  
  // If two lines of the matrix r have the same ki (first element)
  // add the second element of the second vector to the first vector
  // and set the second vector to be all made of 0
    for ( unsigned i = 0; i < r.size(); i++ ){
    for ( unsigned j = i + 1 ; j < r.size(); j++ ){
      if ( r[j][0] == r[i][0] ){
	r[i].push_back( r[j][1] );
	r[j][0] = 0;
	r[j][1] = 0;
      }
    }
  }

    // Print out the matrix r
    cout << "New Matrix r" << endl;
    print_matrix_unsigned(r);

    // Delete all zeroes from r and call it R
    vector < vector < unsigned > > R{};
    for ( unsigned i = 0; i < r.size(); i++ ){
      if ( r[i][0] != 0 ){
	R.push_back( r[i] );
      }
    }

     // Print out the matrix R
    cout << "Matrix R" << endl;
    print_matrix_unsigned(R);


  // Sort the matrix R by the first column
    sort( R.begin(), R.end() );

    
     // Print out the matrix R
    cout << "New Matrix R" << endl;
    print_matrix_unsigned(R);

  // Create a matrix z stacking all the times of observations up to time ki
  // (specified in the fisrt column). The times are also sorted.
    vector < vector < unsigned > > z{};
    vector < unsigned > temp_vector{};
    for ( unsigned i = 0; i < R.size(); i++ ){
      for ( unsigned j = 1; j < R[i].size(); j++ ){
	temp_vector.push_back( R[i][j] );
	sort( temp_vector.begin(), temp_vector.end() );
      }
      z.push_back( temp_vector );
    }
    for ( unsigned i = 0; i < R.size(); i++ ){
      z[i].insert(z[i].begin(), R[i][0]);
    }

    
     // Print out the matrix z
    cout << "Matrix z" << endl;
    print_matrix_unsigned(z);
   
  // Create a matrix Z stacking all the values of the observations up to time ki
  vector < vector < double > > Z{};
  for ( unsigned i = 0; i < z.size(); i++ ){
    vector < double > temp_vect_double{};
    for ( unsigned j = 1; j < z[i].size(); j++ ){
      temp_vect_double.push_back( X[ ( z[i][j] ) ] );
    }
    Z.push_back( temp_vect_double );
    temp_vect_double.clear();
  }
  for ( unsigned i = 0; i < z.size(); i++ ){
    Z[i].insert( Z[i].begin(), z[i][0] );
  }

  // Print out the matrix Z
  cout << "Matrix Z" << endl;
  print_matrix_double(Z);

  /* this is the code for the method
   i here is the index for the current time that goes from 1 to N, 
   j is the index for the particles that goes from 1 to n. I choose n to be 10 */

  vector < vector < unsigned > > S{};
  vector < vector < double > > Y{};
  vector < vector < double > > V{};
  vector < vector < double > > w{};
  vector < vector < double > > W{};
  unsigned n = 3;

  // This is the matrix S of the observations
  // The rows are indexd as j and are for  the particles (n)
  // The columns are indexed as i and are for the number of events (N)
  // It is a matrix of 0s (observation) and 1s (no observation)
    for ( unsigned j = 0; j < n; j++ ){
    vector < unsigned > s{};
    for ( unsigned i = 0; i < N; i++ ){
      s.push_back(0);
    }
    S.push_back(s);
  }

  // This is the calculations matrix Y of the simulated events.
  // In the first column of the matrix we have the values for f(x_1) for all the j particles

  // Creating a vector temp_y simulating from a normal distribution
  // This vector will be used as a starting point to simulate all other n vectors y of lenght N
  // that will populate the matrix Y
  vector < double > y{};
  vector < double > temp_y{};
  for ( unsigned j = 0; j < n; j++){
    mt19937 generator( rd () );
    normal_distribution < double > normalDist( 0, sigmasq / ( 1 - phi * phi ) );
    temp_y.push_back( normalDist ( generator ) );
  }

  // This is the matrix V of the unnormalised weights
  for ( unsigned j = 0; j < n; j++){
    vector < double > v{};
    for ( unsigned i = 0; i < N; i++) {
      v.push_back(1);
    }
    V.push_back(v);
  }

  // This is the matrix W of the normalised weights
  for ( unsigned j = 0; j < n; j++ ){
    vector < double > w{};
    for ( unsigned i = 0; i < N; i++) {
      w.push_back( 1 / N );
    }
    W.push_back(w);
  }

  // This is the calculation of Y with n rows and N columns
  for ( unsigned j = 0; j < n; j++ ){
    for ( unsigned i = 1; i < N; i++ ){
      normal_distribution < double > normalDist( phi * temp_y[i - 1], sigmasq );
      y.push_back( normalDist ( generator ) );
    }
    Y.push_back(y);
    y.clear();
  }
  
   // Putting back the vector temp_y as first column for Y
  for ( unsigned i = 0; i < n; i++ ){
    Y[i].insert( Y[i].begin(), temp_y[i]);
  }
  
  // Print out the vector Y
  //cout << "Matrix Y" << endl;
  //print_matrix_double(Y);

  
  // Sampling from a geometric distribution with p = 0.4
  // find t_(i_j) for particle j (time of the next observation of the event that happened at time i).
  // If t_(i_j) > N - i multiply y_i by (1 - p)^(N - i) and find the new y_i in the vectors Y_j.
  // Else multiply y_i by p(1 - p)^(t_i - 1) and find the new y_i
  // and in the vector S_(i_j) in position k_i = i + t_i, substitute the existing value with a 1.
  for ( unsigned i = 0; i < N; i++ ){
    for ( unsigned j = 0; j < n; j++ ){
      geometric_distribution <> geoDist(p);
      unsigned ti = geoDist ( generator );
      double exp1;
      double exp2;
      if ( ti >= N - i ){
	if ( (N - i) == 1){
	  exp1 = (1 - p);
	}
	else {
	  for ( unsigned k = 1; k < (N - i); k++ ){
	    exp1 = (1 - p) * (1 - p);
	  }
	}
	Y[j][i] = Y[j][i] * exp1;
      }
      else {
	if ( ti == 0 ){
	  exp2 = 1;
	    }
	else if ( ti == 1 ){
	  exp2 =(1 - p);
	}
	else {
	  for ( unsigned k = 1; k < ti - 1; k++ ){
	    exp2 = (1 - p) * (1 - p);
	  }
	}
	Y[j][i] = Y[j][i] * p * exp2;
	unsigned ki = ti + i;
	S[j][ki] = 1;
      }
    }
  }

  // Print out the vector S
  cout << "Matrix S" << endl;
  print_matrix_unsigned(S);

  // Create a vector L with the first column of Z
  // These are the observations in real life
  vector < unsigned > L;
  for ( unsigned i = 0; i < z.size(); i++ ){
    L.push_back(Z[i][0]);
  }
  cout << "vector L" <<  endl; 
  for (unsigned i = 0; i < z.size(); i++){
    cout << L[i] << endl;
  }

  cout << " " << endl;

  /* if i for vector L is 0 and the element S[j][i] is 1 
     (no observation in either real life or simulation) 
     the importance weights will be w_(i_j) = w_[(i - 1)_j] */
  /* if i for vector L is 1 and the element S[j][i] is 0
     (observation in real life, no observation simulated), 
     make the substitutions putting the elements of the vector Z_(k_i) in Y_i^(j) in positions i, 
     then calculate the weights w_i^(j) = w_[(i - 1)_j] * (1 - p)^(t_i + i - N - 1) */
  /* if i for vector L is 1 and the element S[j][i] is 0
     (observation in real life, no observation simulated)
     make the substitutions putting the elements of the vector Z_(k_i) in Y_(i_j) in positions i,
     then calculate the weights w_(i_j) = w_[(i - 1)_j] * (1 - p)^(N - t_i - i + 1) */
  /* if i for vector L is 1 and the element S[j][i] is 1
     (observation in real life and in the simulation) 
     the importance weights will be w_(i_j) = w_[(i - 1)_j]*/
  
  for ( unsigned i = 1; i < N; i++ ){
    vector < double > w{};
    
    bool contains_value {false};
    for ( unsigned k = 0; k < L.size(); k++){
      if ( L[k] == i ){
	contains_value = true;
	break;
      }
    }
    cout << contains_value << endl;
    for ( unsigned j = 0; j < n; j++ ){
      if ( ( S[j][i] == 0 ) && ( contains_value == false ) ){
	cout << " a " << endl;
	V[j][i] = V[j][i - 1];
      }
      
      else if ( ( S[j][i] == 0 ) && ( contains_value == true ) ){
	cout << " b " << endl;
	
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
	
	cout << "pow1 " << pow1 << endl;
	V[j][i] = V[j][i - 1] * pow1 ;
      }
      
      else if ( ( S[j][i] == 1 ) && ( contains_value == false ) ){
	cout << " c " << endl;
	
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
	
	cout << " pow2 " << pow2 << endl;
	V[j][i] = V[j][i - 1] * pow2 ;
      }
      else if ( ( S[j][i] == 1 ) && ( contains_value == true ) ){
	cout << " d " << endl;
	V[j][i] = V[j][i - 1];
      }
    }
  }

  // Print out the vector V
  cout << "Matrix V" << endl;
  print_matrix_double(V);
  

  

  /* for j = 1, ..., n */
  /* normalise the importance weights W_(i_j) = w_(i_j) / sum for j from 2 to N of w_(t_j)*/

  /* calculate the expectation E[y_i] = sum for j from 2 to N y_(i_j) W_(i_j) */
  return 0;
}

// functions definitions

// this function prints a matrix of unsigned
void print_matrix_unsigned( vector < vector < unsigned > > m ){
  for ( const vector < unsigned > & v : m ){
    for  ( unsigned x : v ) cout << x << ' ';
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






