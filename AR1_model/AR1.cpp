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
const float p = 0.4;
const double N = 5;
vector < double > X{};
vector < vector < double > > Obs{};
vector < double > vector_Obs{};

void print_matrix_sizet( vector < vector < size_t > > m );
void print_matrix_double( vector < vector < double > > M );
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
  
  vector_Obs = Obs [N - 1];
  
  cout << "vector x \n";
  print_vector(X);
  cout << "vector ti \n";
  print_vector(vector_ti);
  cout << "matrix Obs \n";
  print_matrix_double(Obs);

  // Create a dat file with the values of X, this is useful to graph with gnuplot
  std::ofstream outFile1( "./vector_X.dat" );
  outFile1 << endl;
  for ( double n : X ){
    outFile1 << n << endl;
  }
  outFile1.close();
  for ( size_t i = 0; i < N; i++ ){
    cout << X[i] << endl;
  }
  
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
  //define the unnormalised weights vector
  vector < vector < double > > V{};
  //define the normalised weights vector
  vector < vector < double > > W{};
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
  // and with variance sigma^2, find the matrix of events "matrix_sample"

  vector < vector < double > > matrix_sample{};
  vector < vector < double > > matrix_new_sample{};
  for ( size_t j = 0; j < n; j++ ){
    vector < double > temp_y{};
    vector < double > row_matrix_new_sample{};
    double y{};
    y = vector_y0[j];
    temp_y.push_back(y);
    for ( size_t i = 1; i < N; i++ ){
      vector < double > row_matrix_sample{};
      normal_distribution < double > normalDist( phi * temp_y[i - 1], sigmasq );
      temp_y.push_back( normalDist ( generator ) );
    }
    for ( size_t l = 0; l < N; l++ ){
      row_matrix_new_sample.push_back( temp_y[l] );
    }
    for ( size_t k = 0; k < N; k++ ){
      if ( vector_Obs[k] != 0 ){
	row_matrix_new_sample[k] = vector_Obs[k];
      }
    }
    matrix_sample.push_back( temp_y );
    matrix_new_sample.push_back( row_matrix_new_sample );
    row_matrix_new_sample.clear();
  }

  vector < vector < size_t > > matrix_ti{};
  for ( size_t j = 0; j < n; j++ ){
    vector < size_t > temp_ti{}; 
    for (size_t i = 0; i < N; i++ ){ 
      geometric_distribution < size_t > geoDist(p);
      size_t ti = geoDist ( generator );
      temp_ti.push_back( ti );
    }
    matrix_ti.push_back( temp_ti );
  }

  // Print matrix_sample, matrix_new_sample and matrix_ti
  cout << "matrix_sample \n";
  print_matrix_double( matrix_sample );
  cout << "matrix_new_sample \n";
  print_matrix_double( matrix_new_sample );
  cout << "matrix_t \n";
  print_matrix_sizet( matrix_ti );

  // create the 3 dimensional vectors for observations and events  
  for ( size_t k = 0; k < n; k++){
    vector < double > row_matrix_sample{};
    vector < double > row_matrix_new_sample{};
    row_matrix_sample = matrix_sample[k];
    row_matrix_new_sample = matrix_new_sample[k];
    vector < vector < double > > sim_matrix{};
    vector < vector < double > > new_sim_matrix{};
    for ( size_t j = 0; j < N; j++ ){
      vector < double > tempvec_sim{};
      vector < double > tempvec_new_sim{};
      for ( size_t i = 0; i < j + 1; i++ ){
	tempvec_sim.push_back( row_matrix_sample[i]);
	tempvec_new_sim.push_back( row_matrix_new_sample[i]);
      }
      sim_matrix.push_back( tempvec_sim );
      new_sim_matrix.push_back( tempvec_new_sim );
    }
    sample.push_back( sim_matrix );
    new_sample.push_back( new_sim_matrix );
  }

  // Print one matrix of the 3 dimensional vector "sample"
  vector < vector < double > > printer1{};
  printer1 = sample[0];
  cout << "printing one of the sample matrix \n";
  print_matrix_double( printer1 );

  // Print one matrix the 3 dimensional vector "new_sample"
  vector < vector < double > > printer2{};
  printer2 = new_sample[0];
  cout << "printing one of the new_sample matrix \n";
  print_matrix_double( printer2 );
  
  // find the vector Y of new sample for all the times
  vector < vector < double > > Y{};
  vector <vector < double > > temp_matrix{};
  for ( size_t i = 0; i < n; i++ ){
      temp_matrix = new_sample[i];
      vector < double > temp_vector;
      for ( size_t j = 0; j < N; j++ ){
         temp_vector = temp_matrix[N-1];
      }
      Y.push_back(temp_vector);
  }
  
  // Creating a matrix "W" of zeroes for the normalised weights
  // and add the first element for each particle = 1/n 
  
  vector < double > w{};
  for ( size_t j = 0; j < n; j++ ){
      double elem{};
   for ( size_t i = 0; i < N; i++ ){
    elem = 0;  
   w.push_back(elem);
   }   
   W.push_back(w);
   w.clear();
  }

  for ( size_t j = 0; j < n; j++ ){
  W[j][0] =  1 / n;
  }
  
  //Finding the unnormalised weights
  
  for ( size_t k = 0; k < n; k++ ){
    vector < vector < double > > matrix_sample{};
    matrix_sample = sample[k];
    vector < vector < double > > matrix_new_sample{};
    matrix_new_sample = new_sample[k];
    double weight(1);
    vector < double > vector_weights{};
    vector_weights.push_back(1);
    for ( size_t i = 0; i < N-1; i++ ){
      vector < double > row_sample{};
      row_sample = matrix_sample[N - 1];
      vector < double > row_new_sample{};
      row_new_sample = matrix_new_sample[N - 1];
      double num = (row_new_sample[i+1] - phi * row_new_sample[i]) * (row_new_sample[i+1] - phi * row_new_sample[i]);
      double den = (row_sample[i+1] - phi * row_sample[i]) * (row_sample[i+1] - phi * row_sample[i]);
      double con = - ( 1 / ( 2 * sigmasq ) );
      if ( row_sample[i] == row_new_sample[i] ){
	weight = weight * 1;
      }
      else if ( row_new_sample[i] != 0 && row_sample[i] == 0 ){
	weight = weight * ( (exp (con * num)) * p ) / ( (exp (con * den)) * (1 - p));
      }
      else if ( row_new_sample[i] == 0 && row_sample[i] != 0 ){
	weight = weight * ( (exp (con * num)) * (1 - p)) / ( (exp (con * den)) * p);
	}
	  else {
	weight = weight * ( (exp (con * num)) ) / ( (exp (con * den)) );
	}
      vector_weights.push_back(weight);
      }
    V.push_back(vector_weights);
  }

  cout << "printing matrix V" << endl;
  print_matrix_double(V);

  
  // normalise the importance weights and put it in matrix W
  for ( size_t i = 0; i < N; i++ ){
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
    
    cout << "printing matrix W" << endl;
    print_matrix_double(W);
    
    cout << "printing matrix Y" << endl;
    print_matrix_double(Y);
  
  // calculate the expectation E[y_i] = sum for j from 2 to N y_(i_j) W_(i_j)
  vector < double > E{};
  for ( size_t i = 0; i < N; i++ ){
    vector < double > multiplication_vector{};
    for ( size_t j = 0; j < n; j++ ){
      multiplication_vector.push_back( Y[j][i] * W[j][i] );
    }
    double sum{};
    for ( auto & n : multiplication_vector)
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
