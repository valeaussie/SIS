#include <random>
#include <iostream>
#include <vector>
#include <fstream>
#include <map>
#include <algorithm>

using namespace std;

/* this is the code to sample from my model and find the vector of the observations
I am going to generate a vector X with the x_i as elements which are found sampling from the AR1 model I have specified. This will be the vector of the events.
Then I am going to generate N, R_(k_i) vectors of the times of the events that happened at time i and have been observed at time k_i.
Finally, I will generate N, Z_(k_i) vectors of observations of events that happened at time i and have been observed at time k_i.
the values of the parameters at this stage are fixed and are sigma^2 = 1, phi = 0.5, p = 0.4, I also choose the value of N = 1000
*/


const int sigmasq = 1;
const float phi = 0.5;
const float p = 0.4;
const int N = 10;
vector< double > X;


random_device rd;

int main(){

  // Sample from a normal distribution. Put the values in a vecotr X.
    
  mt19937 generator( rd () );
  normal_distribution < double > normalDist( 0,sigmasq / ( 1 - phi * phi ) );

  X.push_back( normalDist ( generator ) );
  
  for ( int i = 2; i < N + 1; i++ ){
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

    // Find Ki as ti + i. Create the RKi vectors and put Ki as first value then i. Store all this vectors in a vector of vectors (matrix) r
    
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
  
  // If two lines of the matrix r have the same ki (first element) add the second element of the second vector to the first vector and set the second vector to be all made of 0
  
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
      //int temp_var = z[i][j];
      //cout << temp_var << endl;
      temp_vect_double.push_back(X[(z[i][j] - 1)]);
    }
    Z.push_back( temp_vect_double );
    temp_vect_double.clear();
  }
  
  
  // Print out the matrix R

  cout << "Matrix R" << endl;
  
  for ( const vector < unsigned > & v : R ){
    for  ( int x : v ) cout << x << ' ';
    cout << endl;
  }
  
  // Print out the matrix z

  cout << "Matrix z" << endl;


  for ( const vector < unsigned > & v : z ){
    for  ( int x : v ) cout << x << ' ';
    cout << endl;
  }

  // Print out the vector Z

  cout << "Matrix Z" << endl;

  for ( const vector < double > & v : Z ){
    for  ( double x : v ) cout << x << ' ';
    cout << endl;
  }
 
  return 0;
}

/* calculate the expectation of the x_i */


/* this is the code for the method
   i here is the index for the current time that goes from 1 to N, j is the index for the particles that goes from 1 to n. I choose n to be 1000 */



/* create j vectors S_j of N 0s. These are the vectors that will list if it has been (1) or not (0) an observations at time i for the particle j */
/* create j empty vectors Y_j for the sampled events */
/* create j empty vectors w_j for the unnormalised weights */
/* create j empty vectors W_j for the normalised weights */


/* for i = 1 */
/* for j = 1, ..., n */
/* sample j values for y_(1_j) and put each value as first element of the corresponding vector Y_j. We define the fist observations to be 0 for every j (we don't need to change anything in the vecotrs S_j) */
      /* set j unnormalised importance weights w_(1_j) = 1, and save this as first elements of all w_j */
      /* set j normalised importance weight W_(1_j) = 1/N, and save this as first elements of all W_j */
/* for i = 2, ..., N */
   /* for j = 1, ..., n */
      /* sample from a normal distribution with mean = phi * y_(i - 1) and var = sigma^2 and find the values for y_(i_j), put these as ith elements of the vectors Y_j */
      /* sample from a geometric distribution with p = 0.4 and find t_(i_j) for particle j, which is the time of the next observation of the event that happened at time i */
         /* if t_(i_j) > N - i multiply y_i by (1 - p)^(N - i) and find the new y_i in the vectors Y_j then do nothing else*/
         /* else multiply y_i by p(1 - p)^(t_i - 1) and find the new y_i then in the vector S_(i_j) in position k_i = i + t_i, substitute the existing value with a 1 */

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
