#ifndef   	GET_INPUT_H_
# define   	GET_INPUT_H_
#include <stdio.h>
#include "../Classes/String_Matrices.h"
#include "../Classes/Double_Matrices_cont.h"
#include "../Classes/Int_Matrices_cont.h"
#include "./rand.h"
#include "./manip_dates.h"
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_blas.h>
#include <vector>
using namespace std;


void get_indicator_matrices(Int_Matrices_cont Is_SI_SI,
			    Int_Matrices_cont Is_SI_R,
			    Int_Matrices_cont Is_SI_M,
			    Int_Matrices_cont Is_SI_D,
			    Int_Matrices_cont Is_R_R,
			    Int_Matrices_cont Is_R_M,
			    Int_Matrices_cont Is_R_D,
			    gsl_vector *vect_DOB,
			    gsl_vector *vect_DODiag,
			    gsl_vector *vect_DOD,
			    Double_Matrices_cont mat_specs,
			    gsl_vector * vect_time_line);

void get_summary_from_indicator(vector < vector < int > > &summary,
				Int_Matrices_cont mat_indicator);

void display_summary(vector < vector < int > > &summary);


#endif 	    /* !GET_INPUT_H_ */
