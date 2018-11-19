#include <random>
#include <iostream>
#include <vector>
#include <fstream>
#include <map>
#include <algorithm>

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
  
  // Create a dat file with the values of X, this is useful to graph with gnuplot
    std::ofstream outFile( "./vector_X.dat" );
  outFile << "values of X" << endl;
  for ( double n : X ){
    outFile << n << endl;
  }
  outFile.close();
   

  vector< unsigned > vector_Ki{};
  vector< unsigned > vector_i{};
  vector< vector < unsigned > > r{};

  // Sample from a geometric distribution the values ti
    for ( unsigned i = 1; i < N + 1; i++ ){
    geometric_distribution <> geoDist(p);
    unsigned ti = geoDist ( generator );

    // Find Ki as ti + i.
    //Create the RKi vectors and put Ki as first value then i.
    //Store all this vectors in a vector of vectors (matrix) r
    if ( ti > N - i ){}
    else {
      unsigned ki = i + ti;
      vector_i.push_back(i);
      /* cout << "i " << vector_i.back() << endl; */
      vector_Ki.push_back(ki);
      /*  cout << "Ki " << vector_Ki.back() << endl; */
      vector < unsigned > Rki;
      Rki.push_back(ki);
      Rki.push_back(i);
      r.push_back(Rki);
    }
  }
  
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

  // Delete all zeroes from r and call it R
    vector < vector < unsigned > > R{};
    for ( unsigned i = 0; i < r.size(); i++ ){
      if ( r[i][0] != 0 ){
	R.push_back( r[i] );
      }
    }

  // Sort the matrix R by the first column
    sort( R.begin(), R.end() );

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
    
  // Create a matrix Z stacking all the values of the observations up to time ki
  vector < vector < double > > Z{};
  vector < double > temp_vect_double{};
  for ( unsigned i = 0; i < z.size(); i++ ){
    for ( unsigned j = 1; j < z[i].size(); j++ ){
      temp_vect_double.push_back( X[ ( z[i][j] - 1 ) ] );
    }
    Z.push_back( temp_vect_double );
    temp_vect_double.clear();
  }
  for ( unsigned i = 0; i < z.size(); i++ ){
    Z[i].insert( Z[i].begin(), z[i][0] );
  }
  
  // Print out the matrix R
  cout << "Matrix R" << endl;
  print_matrix_unsigned(R);
  
  // Print out the matrix z
  cout << "Matrix z" << endl;
  print_matrix_unsigned(z);

  // Print out the vector Z
  cout << "Matrix Z" << endl;
  print_matrix_double(Z);

 
  /* this is the code for the method
   i here is the index for the current time that goes from 1 to N, 
   j is the index for the particles that goes from 1 to n. I choose n to be 10 */

  vector < vector < unsigned > > S{};
  vector < vector < double > > Y{};
  vector < vector < double > > V{};
  vector < vector < double > > W{};
  unsigned n = 10;

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
    v.push_back(1);
    V.push_back(v);
  }

  // This is the matrix W of the normalised weights
  for ( unsigned j = 0; j < n; j++ ){
    vector < double > w{};
    w.push_back( 1 / N );
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

  // Print out the vector temp_y
  cout << " vector temp_y " << endl;
  for ( unsigned i = 0; i < n; i++ ){
    cout << temp_y[i] << endl; 
  }

  // Sampling from a geometric distribution with p = 0.4
  // find t_(i_j) for particle j (time of the next observation of the event that happened at time i).
  // If t_(i_j) > N - i multiply y_i by (1 - p)^(N - i) and find the new y_i in the vectors Y_j.
  // Else multiply y_i by p(1 - p)^(t_i - 1) and find the new y_i
  // and in the vector S_(i_j) in position k_i = i + t_i, substitute the existing value with a 1.
  for ( unsigned i = 0; i < N; i++ ){
    for ( unsigned j = 0; j < n; j++ ){
      geometric_distribution <> geoDist(p);
      unsigned ti = geoDist ( generator );
      // cout << ti << endl;
      double exp1;
      double exp2;
      if ( ti >= N - i ){
	// cout <<  " ti >= N - i " << endl;
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
	// cout <<  " ti infinity " << endl;
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


      /* create j vectors V_(i_j) taking the first i elements of S_(i_j) for each particle j. These are the vectors that will list all the observations (1) or missed observation (0) for the particle j up to time i*/
      /* if k_i for vector V_(i_j) is 0 and R_(k_i) is empty (no observation in either real life or simulation) the importance weights will be w_(i_j) = w_[(i - 1)_j] else continue*/
      /* if k_i for vector V_(i_j) is 0 and R_(k_i) is not empty (observation in real life, no observation simulated), 
         make the substitutions putting the elements of the vector Z_(k_i) in Y_i^(j) in positions i, 
         then calculate the weights w_i^(j) = w_[(i - 1)_j] * (1 - p)^(t_i + i - N - 1) else continue */
      /* if k_i for vector V_(i_j) is 1 and R_(k_i) is not empty (observation in real life and in the simulation) the importance weights will be w_(i_j) = w_[(i - 1)_j] else continue*/
      /* if k_i for vector V_(i_j) is 0 and R_(k_i) is empty (observation in real life, no observation simulated
         make the substitutions putting the elements of the vector Z_(k_i) in Y_(i_j) in positions i,
         then calculate the weights w_(i_j) = w_[(i - 1)_j] * (1 - p)^(N - t_i - i + 1) else continue */

   /* for j = 1, ..., n */
   /* normalise the importance weights W_(i_j) = w_(i_j) / sum for j from 2 to N of w_(t_j)*/

/* calculate the expectation E[y_i] = sum for j from 2 to N y_(i_j) W_(i_j) */
  return 0;
}

// functions definitions

void print_matrix_unsigned( vector < vector < unsigned > > m ){
  for ( const vector < unsigned > & v : m ){
    for  ( unsigned x : v ) cout << x << ' ';
    cout << endl;
  }
}

void print_matrix_double( vector < vector < double > > M ){
  for ( const vector < double > & v : M ){
    for  ( double x : v ) cout << x << ' ';
    cout << endl;
  }
}





