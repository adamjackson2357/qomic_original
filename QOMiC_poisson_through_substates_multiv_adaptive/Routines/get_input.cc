#include "get_input.h"

#define DEBUG 0


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
			    gsl_vector * vect_time_line)
{
  unsigned int n_ind=Is_SI_SI.nb_rows;
  unsigned int n_years_fup=Is_SI_SI.nb_columns;
  double first_year=vect_time_line->data[0];
  
  for(unsigned int ind=0;ind<n_ind;ind++){
    unsigned int n_cols=mat_specs.nb_columns;
    unsigned int  current_CaCo=mat_specs.matrix[ind*n_cols+1];
    unsigned int current_smoke_status=mat_specs.matrix[ind*n_cols+2];
    unsigned int current_vit_status=mat_specs.matrix[ind*n_cols+5];
    double current_DOB=vect_DOB->data[ind];
    double current_DODiag=vect_DODiag->data[ind];
    double current_DOD=vect_DOD->data[ind];
    double current_YOB=get_YOB(current_DOB);
    double current_YOD=get_YOB(current_DOD);
    double current_YODiag=get_YOB(current_DODiag);
    unsigned int pos_DOB=current_YOB-first_year;
    unsigned int pos_last_year_SISI=pos_DOB;
    unsigned int pos_last_year_RR=0;



    //Simple data check
    if(current_CaCo==0 && current_YODiag!=1900.00){
      cout << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << endl
	   << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << endl;
      cout << "Ind: " << ind 
	   << " -- Control (" << current_CaCo
	   << ") with DODiag (" << current_YODiag 
	   << ")"
	   << endl
	   << "Run Stopped" << endl;
      cout << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << endl
	   << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << endl;
      exit(1);
    }

    if(current_CaCo==1 && current_YODiag==1900.00){
      cout << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << endl
	   << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << endl;
       cout << "Ind: " << ind 
	   << " -- Case (" << current_CaCo
	   << ") without DODiag (" << current_YODiag 
	   << ")"
	   << endl
	   << "Run Stopped" << endl;
     cout << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << endl
	   << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << endl;
       exit(1);
    }
    if(current_vit_status==2 && current_YOD==1900.00){
      cout << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << endl
	   << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << endl;
      cout << "Ind: " << ind 
	   << " -- Dead (" << current_vit_status
	   << ") without DOD (" << current_YOD 
	   << ")"
	   << endl
	   << "Run Stopped" << endl;
      cout << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << endl
	   << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << endl;
      exit(1);
    }
    if(current_vit_status==1 && current_YOD!=1900.00){
      cout << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << endl
	   << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << endl;
      cout << "Ind: " << ind 
	   << " -- Alive (" << current_vit_status
	   << ") with a DOD (" << current_YOD 
	   << ")"
	   << endl
	   << "Run Stopped" << endl;
      cout << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << endl
	   << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << endl;
      exit(1);
    }
    if(current_CaCo==1 && current_vit_status==2 &&current_YODiag>current_YOD){
      cout << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << endl
	   << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << endl;
      cout << "Ind: " << ind 
	   << " -- Case (" << current_CaCo
	   << ") DODiag (" << current_YODiag 
	   << ") > DODiag (" << current_YOD
	   << ")"
	   << endl
	   << "Run Stopped" << endl;
      cout << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << endl
	   << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << endl;
      exit(1);
    }

    if(current_CaCo==0.0){
      //Case 1: Controls
      if(current_vit_status==1){
	//Case 1.1: ALIVE SI_SI=1 all along
	//SI_M=0
	pos_last_year_SISI=n_years_fup;
      }
      else{
	//Case 1.2: DEAD SI_SI=1 until DOD
	pos_last_year_SISI=current_YOD-first_year;
	Is_SI_M.matrix[ind*Is_SI_M.nb_columns+pos_last_year_SISI]=1;
      }
    }//end of if controls
    else{
      //Case 2: Cases SISI=1 until YODiag-1
      pos_last_year_SISI=current_YODiag-first_year;
      Is_SI_R.matrix[ind*Is_SI_R.nb_columns+pos_last_year_SISI]=1;
      if(current_vit_status==1){
	//Case 2.1 ALIVE: RR=1 from DODiag+1 to today
	pos_last_year_RR=n_years_fup;
      }
      else{
	//Case 2.2 ALIVE: RR=1 from DODiag+1 to DOD-1
	pos_last_year_RR=current_YOD-first_year;
	Is_R_D.matrix[ind*Is_R_D.nb_columns+pos_last_year_RR]=1;
      }
    }
    for(unsigned int year=pos_DOB;year<pos_last_year_SISI;year++){
      Is_SI_SI.matrix[ind*Is_SI_SI.nb_columns+year]=1;
    }
    if(pos_last_year_SISI+1<=n_years_fup){
      for(unsigned int year=pos_last_year_SISI+1;year<pos_last_year_RR;year++){
	Is_R_R.matrix[ind*Is_R_R.nb_columns+year]=1;
      }
    }
    if(DEBUG==1){
      cout << "ind " << ind
	   << " -- CaCo " << current_CaCo
	   << " -- Smok_stat " << current_smoke_status
	   << " -- Vit_stat " << current_vit_status
	   << " -- YOB " << current_YOB
	   << " -- YODiag " << current_YODiag
	   << " -- YOD " << current_YOD
	   << " -- pos_last_year_SISI " << pos_last_year_SISI
	   << " -- pos_last_year_SISI+1 " << pos_last_year_SISI+1
	   << " -- SISI+1<=n_fup " << ((pos_last_year_SISI+1)<=n_years_fup)
	   << " -- pos_last_year_RR " << pos_last_year_RR
	   <<  " -- pos_last_RD " << pos_last_year_RR
	
	   << endl;
      for(unsigned int year=0;year<n_years_fup;year++){
	cout << "\tyear " << year 
	     << " -- date " << vect_time_line->data[year]
	     << " -- SISI " << Is_SI_SI.matrix[ind*Is_SI_SI.nb_columns+year]
	     << " -- SIR " << Is_SI_R.matrix[ind*Is_SI_R.nb_columns+year]
	     << " -- SIM " << Is_SI_M.matrix[ind*Is_SI_M.nb_columns+year]
	     << " -- SID " << Is_SI_D.matrix[ind*Is_SI_D.nb_columns+year]
	     << " -- RR " << Is_R_R.matrix[ind*Is_R_R.nb_columns+year]
	     << " -- RM " << Is_R_M.matrix[ind*Is_R_M.nb_columns+year]
	     << " -- RD " << Is_R_D.matrix[ind*Is_R_D.nb_columns+year]
	     << endl;
      }
    }
    
  }
  
  
  
  
}

void get_summary_from_indicator(vector < vector < int > > &summary,
				Int_Matrices_cont mat_indicator)
{
  unsigned int n_ind=summary.size();
  unsigned int n_years=mat_indicator.nb_columns;
  for(unsigned int ind=0;ind<n_ind;ind++){
    summary[ind].push_back(-1);
    for(unsigned int year=0;year<n_years;year++){
      if(mat_indicator.matrix[ind*mat_indicator.nb_columns+year]==1){
	summary[ind].push_back(year);
      }
    }
  }
}

void display_summary(vector < vector < int > > &summary)
{

  unsigned int n_ind=summary.size();
  for(unsigned int ind=0;ind<n_ind;ind++){
    unsigned int n_elements=summary[ind].size();
    for(unsigned int element=0;element<n_elements;element++){
      cout << summary[ind][element] << " ";
    }
    cout << endl;
  }
}

