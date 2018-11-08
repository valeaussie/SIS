#include <stdio.h>
#include <stdlib.h>
#include <random>
#include <iostream>
#include <string>

/* this is the code to sample from my model and find the vector of the observations
I am going to generate a vector X with the x_i as elements which are found sampling from the AR1 model I have specified. This will be the vector of the events.
Then I am going to generate N, R_(k_i) vectors of the times of the events that happened at time i and have been observed at time k_i.
Finally, I will generate N, Z_(k_i) vectors of observations of events that happened at time i and have been observed at time k_i.
the values of the parameters at this stage are fixed and are sigma^2 = 1, phi = 0.5, p = 0.4, I also choose the value of N = 10
*/


int sigmasq = 1;
float phi = 0.5;
float p = 0.4;
int N = 10;
double X[10];


/* find one sample from a normal distribution with mean 0 and var = sigma^2 / (1 - phi^2) to find x_1 */
/* with a for loop for i = 2 ... N, sample from a normal distribution with mean phi * x_(i - 1) and variance sigma^2 = 1 and find all the other x_i 
   put all this values in a vector X */

int main()
{

  std::normal_distribution<double> normalDist(0,sigmasq / (1 - phi * phi));
  std::mt19937 generator(time(NULL));

  X[1] = normalDist(generator);
  
  for (int i = 2; i < N + 1; ++i){
  std::normal_distribution<double> normalDist(phi * X[i - 1], sigmasq);
    std::mt19937 generator(time(NULL));
    X[i] = normalDist(generator);
  }

printf("x_1 is %f", X[1]);
printf("x_3 is %f", X[3]);
printf("x_10 is %f", X[10]);

  return 0;
}

/* create N empty vectors R_(k_i) for k_i = 1, ..., N */
/* for loop for i = 1, ..., N sample from a geometric distribution with p = 0.4 and find the values t_i. */
/* If t_i > N - i do nothing. */
/* else calculate (k_i) = i + t_i and put the value i in the vector R_(k_i) */
/* create N empty vectors Z_(k_i) for (k_i) = 1, ..., N */
/* for each (k_i), take each element of R_(k_i) and put the corresponding value of x_i in Z_(k_i) */

/* calculate the expectation of the x_i */



/* this is the code for the method
   i here is the index for the current time that goes from 1 to N, j is the index for the particles that goes from 1 to n. I choose n to be 1000 */



/* create j vectors S_j of N 0s. These are the vectors that will list if it has been (1) or not (0) an observations at time i for the particle j*/
/* create j empty vectors Y_j for the sampled events */
/* create j empty vectors w_j for the unnormalised weights */
/* create j empty vectors W_j for the normalised weights */


/* for i = 1
   /* for j = 1, ..., n */
      /* sample j values for y_(1_j) and put each value as first element of the corresponding vector Y_j. 
         We define the fist observations to be 0 for every j (we don't need to change anything in the vecotrs S_j) */
      /* set j unnormalised importance weights w_(1_j) = 1, and save this as first elements of all w_j */
      /* set j normalised importance weight W_(1_j) = 1/N, and save this as first elements of all W_j */
/* for i = 2, ..., N */
   /* for j = 1, ..., n */
      /* sample from a normal distribution with mean = phi * y_(i - 1) and var = sigma^2 and find the values for y_(i_j), 
         put these as ith elements of the vectors Y_j */
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
