#include "manip_dates.h"

#define DEBUG 0

using namespace std;

gsl_vector* get_Time_Line(gsl_vector *vect_DOB,
			  gsl_vector *vect_DODiag,
			  gsl_vector *vect_DOD)
{
  double min_DOB=gsl_vector_min(vect_DOB);
  double max_DOB=gsl_vector_max(vect_DOB);
  double min_DODiag=gsl_vector_min(vect_DODiag);
  double max_DODiag=gsl_vector_max(vect_DODiag);
  double min_DOD=gsl_vector_min(vect_DOD);
  double max_DOD=gsl_vector_max(vect_DOD);
  
  cout << "Min DOB " << min_DOB
       << " -- Max DOB " << max_DOB
       << endl
       << "Min DODiag " << min_DODiag
       << " -- Max DODiag " << max_DODiag
       << " -- Min DOD " << min_DOD
       << " -- Max DOD " << max_DOD
        << endl;
  
  double min_YOB=trunc(min_DOB/10000);
  double max_YOB=trunc(max_DOB/10000);
  double min_YODiag=trunc(min_DODiag/10000);
  double max_YODiag=trunc(max_DODiag/10000);
  double min_YOD=trunc(min_DOD/10000);
  double max_YOD=trunc(max_DOD/10000);

  double latest_year=max(max_YOD,max_YODiag);
  unsigned int n_years=latest_year-min_YOB+1;

  cout << "Min YOB " << min_YOB
       << " -- Max YOB " << max_YOB
       << endl
       << "Min YODiag " << min_YODiag
       << " -- Max YODiag " << max_YODiag
       << " -- Min YOD " << min_YOD
       << " -- Max YOD " << max_YOD
       << " -- n_years " << n_years
       << endl;
  gsl_vector *vect_time_line=gsl_vector_alloc(n_years);
  for(unsigned int year=0;year<n_years;year++){
    vect_time_line->data[year]=year+min_YOB;
  }
  return vect_time_line;
  
}

gsl_vector* get_vect_age_at_recr(gsl_vector *vect_DOB,
				 gsl_vector *vect_DORecr)
{
  unsigned int n_ind=vect_DOB->size;
  gsl_vector *vect_age_at_recr=gsl_vector_alloc(n_ind);
  for(unsigned int ind=0;ind<n_ind;ind++){
    double current_DOB=vect_DOB->data[ind];
    double current_DORecr=vect_DORecr->data[ind];
    vect_age_at_recr->data[ind]=get_age(current_DOB,
					current_DORecr);
  }
  return vect_age_at_recr;
  
}

double get_YOB(double DOB)
{
  double YOB=trunc(DOB/10000);
  return YOB;
}
double get_MOB(double DOB,
	       double YOB)
{
  double MOB=trunc((DOB-YOB*10000)/100);
  return MOB;
}

double get_DayOB(double DOB,
		 double YOB,
		 double MOB)
{
  double DayOB=DOB-(YOB*10000+MOB*100);
  return DayOB;
}


double get_age(double DOB,
	       double date)
{
  double YOB=get_YOB(DOB);
  double MOB=get_MOB(DOB,
		     YOB);
  
  double DayOB=get_DayOB(DOB,
			 YOB,
			 MOB);

  double year=get_YOB(date);

  double month=get_MOB(date,
		       year);

  double day=get_DayOB(date,
		       year,
		       month);

  double age=year-YOB+((month-MOB)/12)+((day-DayOB)/365.25);

  return age;
}

void get_matrix_age(Double_Matrices_cont mat_age,
		    gsl_vector * vect_time_line,
		    gsl_vector * vect_DOB)
{

  unsigned int n_ind=mat_age.nb_rows;
  unsigned int n_years=mat_age.nb_columns;

  for(unsigned int ind=0;ind<n_ind;ind++){
    double current_DOB=vect_DOB->data[ind];
    double YOB=get_YOB(current_DOB);
    double MOB=get_MOB(current_DOB,
		       YOB);

    double DayOB=get_DayOB(current_DOB,
			   YOB,
			   MOB);
    double init_age=vect_time_line->data[0]-YOB+(1-MOB)/12+(1-DayOB)/365.25;
    //     cout << "Ind " << ind
    // 	 << " -- DOB " << current_DOB
    // 	 << " -- YOB " << YOB
    // 	 << " -- MOB " << MOB
    // 	 << " -- DayOB " << DayOB
    // 	 << " -- Age_init " << init_age
    // 	 << endl;
    if(init_age>0){
      mat_age.matrix[ind*mat_age.nb_columns]=init_age;
    }
    double current_age=init_age;
    for(unsigned int year=1;year<n_years;year++){
      current_age+=1.0;
      if(current_age>0.0){
	mat_age.matrix[ind*mat_age.nb_columns+year]=current_age;
      }
    }
  }


}


void get_matrix_age_cont(Double_Matrices_cont mat_age,
			 gsl_vector * vect_time_line,
			 gsl_vector * vect_DOB)
{

  unsigned int n_ind=mat_age.nb_rows;
  unsigned int n_years=mat_age.nb_columns;

  for(unsigned int ind=0;ind<n_ind;ind++){
    double current_DOB=vect_DOB->data[ind];
    double YOB=get_YOB(current_DOB);
    double MOB=get_MOB(current_DOB,
		       YOB);

    double DayOB=get_DayOB(current_DOB,
			   YOB,
			   MOB);
    double init_age=vect_time_line->data[0]-YOB+(1-MOB)/12+(1-DayOB)/365.25;
    //     cout << "Ind " << ind
    // 	 << " -- DOB " << current_DOB
    // 	 << " -- YOB " << YOB
    // 	 << " -- MOB " << MOB
    // 	 << " -- DayOB " << DayOB
    // 	 << " -- Age_init " << init_age
    // 	 << endl;
    if(init_age>0){
      mat_age.matrix[ind*mat_age.nb_columns]=init_age;
    }
    double current_age=init_age;
    for(unsigned int year=1;year<n_years;year++){
      current_age+=1.0;
      if(current_age>0.0){
	mat_age.matrix[ind*mat_age.nb_columns+year]=current_age;
      }
    }
  }


}


