#include <stdio.h>
#include "../Classes/String_Matrices.h"
#include "../Classes/Double_Matrices_cont.h"
#include "./rand.h"
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_blas.h>
#include <vector>

using namespace std; 

gsl_vector* get_Time_Line(gsl_vector *vect_DOB,
			  gsl_vector *vect_DODiag,
			  gsl_vector *vect_DOD);

double get_YOB(double DOB);

double get_MOB(double DOB,
	       double YOB);
double get_DayOB(double DOB,
		 double YOB,
		 double MOB);
double get_age(double DOB,
	       double date);
gsl_vector* get_vect_age_at_recr(gsl_vector *vect_DOB,
				 gsl_vector *vect_DORecr);
void get_matrix_age(Double_Matrices_cont mat_age,
		    gsl_vector * vect_time_line,
		    gsl_vector * vect_DOB);


void get_matrix_age_cont(Double_Matrices_cont mat_age,
			 gsl_vector * vect_time_line,
			 gsl_vector * vect_DOB);
