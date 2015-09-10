#include <iostream>
#include <string.h>						// string modification functions
#include <time.h>						// time functions
#include <math.h>						// math functions
#include <stdio.h>						// output functions
#include <stdlib.h>						// standard
#include <limits.h>
#include <gsl/gsl_rng.h>				// random number generator functions
#include <gsl/gsl_blas.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_randist.h>			// random number distributions
#include <gsl/gsl_sort_vector.h>		// vector sorting functions
#include <gsl/gsl_odeiv.h>				// differential equation solver
#include <gsl/gsl_errno.h>


using namespace std;

struct data{
	gsl_vector* alter;
	int size;
};

void variation(struct data dataArr[], int l, struct data temp);

int main()
{	
	int l = 5, i;
	gsl_vector* alter = gsl_vector_calloc(l);
	data temp = {alter,0};
	
	
	data dataArr[l];
	for(i = 0 ;i<l;i++)
	{
	  dataArr[i] = temp;
	}

	
	variation(dataArr,5, temp);
	
	cout << "alter von 3 ist " << gsl_vector_get(dataArr[3].alter,4) <<endl;
	cout << "size von 3 ist " << dataArr[3].size <<endl;
	
	
	return 0;
}


void variation(struct data dataArr[], int l, struct data temp)
{
  	for (int j = 0; j<l;j++)
	{	

 		gsl_vector* altertemp = gsl_vector_calloc(l);
// 		data tempData = {alter,0};
	
		//tempData.dA = [0,55,67,34];
		gsl_vector_set(altertemp,j,(2*j));
		temp.alter = altertemp;
		//cout << "test"<< endl;
		cout << "j ist " << j << " \talter " << gsl_vector_get(temp.alter,j) << endl;

		temp.size= 6+3*j;
		dataArr[j]=temp;
	}
}