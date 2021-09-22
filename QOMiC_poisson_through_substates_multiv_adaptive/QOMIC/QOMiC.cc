#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <ctype.h>
#include <iostream>
#include "xml_file_write.h"
#include "xml_file_read.h"
#include <vector>
#include <algorithm>
#include <memory>
#include <fstream>
#include <../Routines/struc.h>
#include <../Routines/manip_dates.h>
#include <../Routines/dyn_name.h>
#include <../Routines/matrix_handling.h>
#include <../Routines/rand.h>
#include "../Routines/get_input.h"
#include "../Routines/get_pbties.h"
#include <../Routines/xml_file_read.h>
#include <../Classes/String_Matrices.h>
#include <../Classes/Int_Matrices_cont.h>
#include <../Classes/Double_Matrices_cont.h>
#include <gsl/gsl_fit.h>
#include <gsl/gsl_multifit.h>
#include <gsl/gsl_cdf.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_sort.h>
#include <gsl/gsl_sort_vector.h>
#include <gsl/gsl_permutation.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_blas.h>

#define DEBUG 0

using namespace std;


// File Names
char nameIN[256]; // input file called by switch -i
char nameOUT[256]; // output  file called by switch -o


//******************************************************************
//*main
//******************************************************************


int main(int argc, char *  argv[])
{
    int na=0;
    char filename_in_IDs[1000];
    char filename_in_dates[1000];
    char filename_in_specs[1000];
    char filename_in_exposure[1000];
    char filename_in_mortalityF[1000];
    char filename_in_mortalityM[1000];
    char fic_out_MCMC[1000];
    
    long MY_SEED=0;
    
    
    // Initialisation of parameters (default values)
    
    unsigned int Nb_states_in_I=0;
    
    unsigned int n_iter=0;
    unsigned int burn_in=0;
    
    unsigned int adaptive=1;
    unsigned int n_iter_adaptive=100;
    unsigned int update_frequency=100;
    double adaptive_upperbound=50;
    double adaptive_lowerbound=40;
    
    double mu_init=0.0; // intercept
    double sigma_mu=0.0;
    double lambda1_init=0.0; // cumulative exposure
    double sigma_lambda1=0.0;
    double lambda2_init=0.0; // current age
    double sigma_lambda2=0.0;
    double lambda3_init=0.0; // age at starting
    double sigma_lambda3=0.0;
    double lambda4_init=0.0;
    double sigma_lambda4=0.0;
    double gamma_init=0.0; // transitions through sub-states of I
    double sigma_gamma=0.0;
    na++;
    
    
    // Reading parameters from command line
    
    while(na < argc){
        if ( 0 == strcmp(argv[na],"-ID") ){
            strcpy(filename_in_IDs,argv[++na]);
            if (na+1==argc) break;
            na++;
        }
        if ( 0 == strcmp(argv[na],"-dates") ){
            strcpy(filename_in_dates,argv[++na]);
            if (na+1==argc) break;
            na++;
        }
        else if ( 0 == strcmp(argv[na],"-specs") ){
            strcpy(filename_in_specs,argv[++na]);
            if (na+1==argc) break;
            na++;
        }
        else if ( 0 == strcmp(argv[na],"-exp") ){
            strcpy(filename_in_exposure,argv[++na]);
            if (na+1==argc) break;
            na++;
        }
        else if ( 0 == strcmp(argv[na],"-mort") ){
            strcpy(filename_in_mortalityF,argv[++na]);
            strcpy(filename_in_mortalityM,argv[++na]);
            if (na+1==argc) break;
            na++;
        }
        else if ( 0 == strcmp(argv[na],"-seed") ){
            MY_SEED=(long)((atoi(argv[++na])));
            if (na+1==argc) break;
            na++;
        }
        else if ( 0 == strcmp(argv[na],"-k") ){
            Nb_states_in_I=atoi(argv[++na]);
            if (na+1==argc) break;
            na++;
        }
        else if ( 0 == strcmp(argv[na],"-iter") ){
            n_iter=atoi(argv[++na]);
            if (na+1==argc) break;
            na++;
        }
        else if ( 0 == strcmp(argv[na],"-burn_in") ){
            burn_in=atoi(argv[++na]);
            if (na+1==argc) break;
            na++;
        }
        else if ( 0 == strcmp(argv[na],"-adaptive") ){
            adaptive=atoi(argv[++na]);
            if (na+1==argc) break;
            na++;
        }
        else if ( 0 == strcmp(argv[na],"-iter_adaptive") ){
            n_iter_adaptive=atoi(argv[++na]);
            if (na+1==argc) break;
            na++;
        }
        else if ( 0 == strcmp(argv[na],"-update_frequency") ){
            update_frequency=atoi(argv[++na]);
            if (na+1==argc) break;
            na++;
        }
        else if ( 0 == strcmp(argv[na],"-adaptive_upperbound") ){
            adaptive_upperbound=(double)(atof(argv[++na]));
            if (na+1==argc) break;
            na++;
        }
        else if ( 0 == strcmp(argv[na],"-adaptive_lowerbound") ){
            adaptive_lowerbound=(double)(atof(argv[++na]));
            if (na+1==argc) break;
            na++;
        }
        else if ( 0 == strcmp(argv[na],"-mu") ){
            mu_init=(double)(atof(argv[++na]));
            sigma_mu=(double)(atof(argv[++na]));
            if (na+1==argc) break;
            na++;
        }
        else if ( 0 == strcmp(argv[na],"-lambda1") ){
            lambda1_init=(double)(atof(argv[++na]));
            sigma_lambda1=(double)(atof(argv[++na]));
            if (na+1==argc) break;
            na++;
        }
        else if ( 0 == strcmp(argv[na],"-lambda2") ){
            lambda2_init=(double)(atof(argv[++na]));
            sigma_lambda2=(double)(atof(argv[++na]));
            if (na+1==argc) break;
            na++;
        }
        else if ( 0 == strcmp(argv[na],"-lambda3") ){
            lambda3_init=(double)(atof(argv[++na]));
            sigma_lambda3=(double)(atof(argv[++na]));
            if (na+1==argc) break;
            na++;
        }
        else if ( 0 == strcmp(argv[na],"-lambda4") ){
            lambda4_init=(double)(atof(argv[++na]));
            sigma_lambda4=(double)(atof(argv[++na]));
            if (na+1==argc) break;
            na++;
        }
        else if ( 0 == strcmp(argv[na],"-gamma") ){
            gamma_init=(double)(atof(argv[++na]));
            sigma_gamma=(double)(atof(argv[++na]));
            if (na+1==argc) break;
            na++;
        }
        else if ( 0 == strcmp(argv[na],"-out") ){
            strcpy(fic_out_MCMC,argv[++na]);
            if (na+1==argc) break;
            na++;
        }
        else{
            cout << "Unknown option: " << argv[na] << endl;
            exit(1);
        }
    }


    // Dynamic writing of output files

    // History file
    string Extension_out=".txt";
    string Name_number1="History_iter";
    string Name_number2="k";
    string OutputName_MCMC=Write_dyn_nname_2params(fic_out_MCMC,
                                                   Name_number1,
                                                   n_iter,
                                                   Name_number2,
                                                   Nb_states_in_I,
                                                   Extension_out);
    ofstream f_out;
    f_out.open(OutputName_MCMC.c_str(),ios::out);
    if(f_out.fail()){
        cout << "Invalid Path and/or permission rights for " << OutputName_MCMC << " -- run stopped." << endl;
        exit(1);
    }
    else{
        f_out << "Iter\tMu\tLambda1\tLambda2\tLambda3\tLambda4\tGamma\tL\tsigma_mu\tsigma_12\tsigma_13\tsigma_21\tsigma_lambda1\tsigma_23\tsigma_31\tsigma_32\tsigma_gamma"<<endl;
    }

    // Posterior likelihood file
    string Name_number3="Posterior_L_iter";
    string OutputName_post_L=Write_dyn_nname_2params(fic_out_MCMC,
                                                     Name_number3,
                                                     n_iter,
                                                     Name_number2,
                                                     Nb_states_in_I,
                                                     Extension_out);

    ofstream f_outL;
    f_outL.open(OutputName_post_L.c_str(),ios::out);
    if(f_outL.fail()){
        cout << "Invalid Path and/or permission rights for " << OutputName_post_L << " -- run stopped." << endl;
        exit(1);
    }
    else{
        f_outL << "BurnIn\tMu\tpost_Lambda1\tpost_Lambda2\tpost_Lambda3\tpost_Lambda4\tpost_Gamma\tpost_L"<<endl;
    }
    smyrand((long)(MY_SEED));


    // Summary of the adaptive parameters:

    if (adaptive==1){
        cout << "Adaptive algorithm: " << endl;
        cout << "adaptive_lowerbound: " << adaptive_lowerbound << endl;
        cout << "adaptive_upperbound: " << adaptive_upperbound << endl;
        cout << "Number of iterations for adaptive: " << n_iter_adaptive << endl;
        cout << "Update frequency: " << update_frequency << endl;
    }

    // Reading dates file

    //Col 1: DOB
    //Col 2: D_Recr
    //Col 3: DOD
    //Col 4: D_lung

    Double_Matrices_cont mat_dates;
    mat_dates.Read_from_file(filename_in_dates);
    //mat_dates.Display_matrix_header();
    gsl_matrix *mat_dates_work=Double_matrices_cont_2_gsl_matrix(mat_dates);
    unsigned int col_DOB=0;
    gsl_vector *vect_DOB=get_one_vect(col_DOB,
                                      mat_dates_work);

    unsigned int col_DOR=1;
    gsl_vector *vect_DOR=get_one_vect(col_DOR,
                                      mat_dates_work);
    
    
    unsigned int col_DOD=2;
    gsl_vector *vect_DOD=get_one_vect(col_DOD,
                                      mat_dates_work);
    unsigned int col_DODiag=3;
    gsl_vector *vect_DODiag=get_one_vect(col_DODiag,
                                         mat_dates_work);
    mat_dates.Free_double_matrix_cont();
    gsl_matrix_free(mat_dates_work);
    
    
    // Reading Specs file
    
    //Col 1: Gender
    //Col 2: Ca/Co
    //Col 3: Smok_status
    //Col 4: Country
    //Col 5: Center
    //Col 6: Vital_status
    //Col 7: Age_at_start
    //Col 8: Age_at_quitting
    
    Double_Matrices_cont mat_specs;
    mat_specs.Read_from_file(filename_in_specs);
    gsl_matrix *mat_specs_work=Double_matrices_cont_2_gsl_matrix(mat_specs);
    //mat_specs.Display_matrix_header();
    unsigned int col_AaS=6;
    gsl_vector *vect_AaS=get_one_vect(col_AaS,
                                      mat_specs_work);
    unsigned int col_AaQ=7;
    gsl_vector *vect_AaQ=get_one_vect(col_AaQ,
                                      mat_specs_work);
    gsl_matrix_free(mat_specs_work);
    
    
    // Reading exposure file
    
    Double_Matrices_cont mat_exposure;
    mat_exposure.Read_from_file(filename_in_exposure);
    
    //    cout << "Matrix of exposure:" << endl;
    //    mat_exposure.Display_matrix_header();
    
    gsl_vector *vect_time_line=get_Time_Line(vect_DOB,
                                             vect_DODiag,
                                             vect_DOD);
    display_gsl_vector(vect_time_line);
    
    unsigned int n_years_follow_up=vect_time_line->size;
    unsigned int n_ind=mat_exposure.nb_rows;
    cout << "N_year_follow-up= " << n_years_follow_up
    << " -- n_ind " << n_ind
    << endl;
    
    Double_Matrices_cont mat_mort_female;
    mat_mort_female.Read_from_file(filename_in_mortalityF);
    
    Double_Matrices_cont mat_mort_male;
    mat_mort_male.Read_from_file(filename_in_mortalityM);
    
    
    // Calculating state matrices from input
    
    Int_Matrices_cont Is_SI_SI;
    Is_SI_SI.Alloc_int_matrix_cont(n_ind,n_years_follow_up); //
    //    Is_SI_SI.Display_matrix_header();
    Int_Matrices_cont Is_SI_R;
    Is_SI_R.Alloc_int_matrix_cont(n_ind,n_years_follow_up);
    Int_Matrices_cont Is_SI_M;
    Is_SI_M.Alloc_int_matrix_cont(n_ind,n_years_follow_up);
    Int_Matrices_cont Is_SI_D;
    Is_SI_D.Alloc_int_matrix_cont(n_ind,n_years_follow_up);
    Int_Matrices_cont Is_R_R;
    Is_R_R.Alloc_int_matrix_cont(n_ind,n_years_follow_up);
    Int_Matrices_cont Is_R_M;
    Is_R_M.Alloc_int_matrix_cont(n_ind,n_years_follow_up);
    Int_Matrices_cont Is_R_D;
    Is_R_D.Alloc_int_matrix_cont(n_ind,n_years_follow_up);
    
    get_indicator_matrices(Is_SI_SI,
                           Is_SI_R,
                           Is_SI_M,
                           Is_SI_D,
                           Is_R_R,
                           Is_R_M,
                           Is_R_D,
                           vect_DOB,
                           vect_DODiag,
                           vect_DOD,
                           mat_specs,
                           vect_time_line);
    
    //    cout << "SI-SI transitions:" << endl;
    //    Is_SI_SI.Display_matrix_header();
    //
    //    cout << "SI-R transitions:" << endl;
    //    Is_SI_R.Display_matrix_header();
    
    //Int_Matrices_Var_Dim summary_SI_SI;
    //summary_SI_SI.Alloc_matrix(n_ind);
    vector < vector < int > > summary_SI_SI;
    vector < vector < int > > summary_SI_R;
    vector < vector < int > > summary_SI_M;
    vector < vector < int > > summary_R_R;
    
    summary_SI_SI.resize(n_ind);
    summary_SI_R.resize(n_ind);
    summary_SI_M.resize(n_ind);
    summary_R_R.resize(n_ind);
    
    get_summary_from_indicator(summary_SI_SI,
                               Is_SI_SI);
    get_summary_from_indicator(summary_SI_R,
                               Is_SI_R);
    get_summary_from_indicator(summary_SI_M,
                               Is_SI_M);
    get_summary_from_indicator(summary_R_R,
                               Is_R_R);
    
    //    cout << summary_SI_R[2][1] << endl;
    unsigned int ind=2;
    cout << "First year SI-SI: " << summary_SI_SI[ind][1] << endl;
    cout << "Second year SI-SI: " << summary_SI_SI[ind][2] << endl;
    cout << "Second year SI-SI: " << summary_SI_SI[ind][3] << endl;
    unsigned int n_trans_SI_SI=summary_SI_SI[ind].size();
    cout << "Number of SI-SI transitions: " << n_trans_SI_SI << endl;
    cout << summary_SI_SI[ind][n_trans_SI_SI-2] << endl;
    cout << summary_SI_SI[ind][n_trans_SI_SI-1] << endl;
    cout << summary_SI_SI[ind][n_trans_SI_SI] << endl;
    cout << "First year SI-R: " << summary_SI_R[ind][1] << endl;
    
    Is_SI_SI.Free_int_matrix_cont();
    Is_SI_R.Free_int_matrix_cont();
    Is_SI_M.Free_int_matrix_cont();
    Is_SI_D.Free_int_matrix_cont();
    Is_R_R.Free_int_matrix_cont();
    Is_R_M.Free_int_matrix_cont();
    Is_R_D.Free_int_matrix_cont();
    
    Double_Matrices_cont mat_age;
    mat_age.Alloc_double_matrix_cont(n_ind,
                                     n_years_follow_up);
    get_matrix_age_cont(mat_age,
                        vect_time_line,
                        vect_DOB);
    
    //    cout << "Matrix of age:" << endl;
    //    mat_age.Display_matrix_header();

    
    // Initialisation of probabilities
    
    // Step 1: Calculating P_SI_M=P_S_M=P_I_M other cause mortality (fixed, based on input mortality data)
    
    Double_Matrices_cont mat_P_SI_M;
    mat_P_SI_M.Alloc_double_matrix_cont(n_ind,
                                        n_years_follow_up);
    
    // mat_age.Display_matrix();
    get_P_SI_M_VBT(mat_P_SI_M,
                   mat_mort_female,
                   mat_mort_male,
                   mat_specs,
                   mat_age,
                   vect_time_line,
                   vect_DOR,
                   vect_AaS,
                   vect_AaQ);
    
    //    cout << "Matrix of probability of death:" << endl;
    //    mat_P_SI_M.Display_matrix_header();
    
    mat_mort_female.Free_double_matrix_cont();
    mat_mort_male.Free_double_matrix_cont();
    
    
    // Step 2: Initialisation of other transition probabilities (filled with zeroes for now)
    
    Double_Matrices_cont mat_P_S_I;
    mat_P_S_I.Alloc_double_matrix_cont(n_ind,
                                       n_years_follow_up);
    
    Double_Matrices_cont mat_P_S_S;
    mat_P_S_S.Alloc_double_matrix_cont(n_ind,
                                       n_years_follow_up);
    
    Double_Matrices_cont mat_P_Ii_Ij_baseline;
    mat_P_Ii_Ij_baseline.Alloc_double_matrix_cont(Nb_states_in_I+2,
                                                  Nb_states_in_I+2);
    
    Double_Matrices_cont mat_P_Ii_Ij;
    mat_P_Ii_Ij.Alloc_double_matrix_cont(Nb_states_in_I+2,
                                         Nb_states_in_I+2);
    
    Double_Matrices_cont mat_P_SI_R;
    mat_P_SI_R.Alloc_double_matrix_cont(n_ind,
                                        n_years_follow_up);
    
    Double_Matrices_cont mat_P_SI_SI;
    mat_P_SI_SI.Alloc_double_matrix_cont(n_ind,
                                         n_years_follow_up);
    
    
    // Step 3: Initialisation of probability of being in S or in one I sub-state (filled with zeroes for now)
    
    Double_Matrices_cont mat_H_k_per_ind;
    mat_H_k_per_ind.Alloc_double_matrix_cont(Nb_states_in_I+1,
                                             n_years_follow_up);
    
    //    cout << "Matrix of probability of being in S or one of the I sub-states:" << endl;
    //    mat_H_k_per_ind.Display_matrix_header();
    
    
    // Initialisation of current parameters

    // Number of parameters
    int p=2;
    int q=1;
    int N=p+q;
    int ncol=1+6+1+N*N;
    cout << "Number of columns in history matrix: "<< ncol << endl;
    
    // Matrix in which all visited parameters values and likelihoods will be stored
    Double_Matrices_cont mat_history;
    mat_history.Alloc_double_matrix_cont(n_iter,ncol);
    
    double current_mu=mu_init;
    double current_lambda1=lambda1_init;
    double current_lambda2=lambda2_init;
    double current_lambda3=lambda3_init;
    double current_lambda4=lambda4_init;
    double current_gamma=gamma_init;

    
    // Creating the covariance matrix Sigma between parameters in theta
    
    gsl_matrix* Sigma=gsl_matrix_calloc(N,N);
    Sigma->data[0 * Sigma->tda + 0]=sigma_mu;
    Sigma->data[1 * Sigma->tda + 1]=sigma_lambda1;
    Sigma->data[2 * Sigma->tda + 2]=sigma_gamma;
    
    cout<<"Current covariance Sigma:"<<endl;
    display_gsl_matrix(Sigma);
    
    
    // Cholesky decomposition of the covariance Sigma
    
    gsl_matrix* L=gsl_matrix_calloc(N,N);
    L->data[0 * L->tda + 0]=sigma_mu;
    L->data[1 * L->tda + 1]=sigma_lambda1;
    //    Sigma->data[2 * Sigma->tda + 2]=sigma_lambda2;
    //    Sigma->data[3 * Sigma->tda + 3]=sigma_lambda3;
    //    Sigma->data[4 * Sigma->tda + 4]=sigma_lambda4;
    L->data[2 * L->tda + 2]=sigma_gamma;
    
    cout<<"After Cholesky decomposition:"<<endl;
    gsl_linalg_cholesky_decomp(L);
    display_gsl_matrix(L);
    
    double cd=pow(2.4,2)/N;
    cout << "Scaling factor: " << cd << endl;
    
    
    
    // Initialisation of the MCMC
    
    // Computing the likelihood with initial parameters
    
    double current_L=wrapped_L_calculation_5_Params_model5(mat_P_S_I,
                                                           mat_P_S_S,
                                                           mat_exposure,
                                                           mat_P_Ii_Ij_baseline,
                                                           mat_P_Ii_Ij,
                                                           mat_P_SI_M,
                                                           mat_P_SI_R,
                                                           mat_P_SI_SI,
                                                           mat_H_k_per_ind,
                                                           mat_age,
                                                           mat_specs,
                                                           vect_time_line,
                                                           vect_DOR,
                                                           vect_AaS,
                                                           vect_AaQ,
                                                           vect_DOB,
                                                           vect_DOD,
                                                           vect_DODiag,
                                                           summary_SI_SI,
                                                           summary_SI_R,
                                                           summary_SI_M,
                                                           current_mu,
                                                           current_lambda1,
                                                           current_lambda2,
                                                           current_lambda3,
                                                           current_lambda4,
                                                           current_gamma);
    
    
    if(DEBUG==1){
        cout << "Mu " << current_mu
        << "Lambda1 " << current_lambda1
        << " -- Lambda2 " << current_lambda2
        << " -- Lambda3 " << current_lambda3
        << " -- Lambda4 " << current_lambda4
        << " -- gamma " << current_gamma
        << " -- initial L " << current_L
        << endl;
    }
    
    
    // Storing parameters values at current iteration
    
    mat_history.matrix[0*mat_history.nb_columns]=0;
    
    mat_history.matrix[0*mat_history.nb_columns+1]=current_mu;
    mat_history.matrix[0*mat_history.nb_columns+2]=current_lambda1;
    mat_history.matrix[0*mat_history.nb_columns+3]=current_lambda2;
    mat_history.matrix[0*mat_history.nb_columns+4]=current_lambda3;
    mat_history.matrix[0*mat_history.nb_columns+5]=current_lambda4;
    mat_history.matrix[0*mat_history.nb_columns+6]=current_gamma;
    
    mat_history.matrix[0*mat_history.nb_columns+7]=current_L;
    
    cout << "Filling covariance values:" << endl;
    for(int k=0;k<N*N;k++){
        mat_history.matrix[0*mat_history.nb_columns+8+k]=Sigma->data[k];
        cout << Sigma->data[k] << endl;
    }
    
//    mat_history.matrix[0*mat_history.nb_columns+13]=sigma_mu;
//    mat_history.matrix[0*mat_history.nb_columns+14]=sigma_lambda1;
//    mat_history.matrix[0*mat_history.nb_columns+15]=sigma_lambda2;
//    mat_history.matrix[0*mat_history.nb_columns+16]=sigma_lambda3;
//    mat_history.matrix[0*mat_history.nb_columns+17]=sigma_lambda4;
//    mat_history.matrix[0*mat_history.nb_columns+18]=sigma_gamma;
    
    
    // Writing parameters values in dynamic output file
    
    for(unsigned int col=0;col<mat_history.nb_columns;col++){
        f_out << mat_history.matrix[0*mat_history.nb_columns+col]
        << "\t";
    }
    f_out << endl;
    
    
    // Initialisation of acceptance counts for each parameter (for computation of acceptance probabilities at the end of the run)
    
    unsigned int nb_accepted_theta=0;
    unsigned int nb_accepted_theta_batch=0;
    unsigned int nb_accepted_theta_after_BI=0;
    
    unsigned int nb_accepted_mu=0;
    unsigned int nb_accepted_mu_batch=0;
    unsigned int nb_accepted_mu_after_BI=0;
    unsigned int nb_accepted_lambda1=0;
    unsigned int nb_accepted_lambda1_batch=0;
    unsigned int nb_accepted_lambda1_after_BI=0;
    unsigned int nb_accepted_lambda2=0;
    unsigned int nb_accepted_lambda2_batch=0;
    unsigned int nb_accepted_lambda2_after_BI=0;
    unsigned int nb_accepted_lambda3=0;
    unsigned int nb_accepted_lambda3_batch=0;
    unsigned int nb_accepted_lambda3_after_BI=0;
    unsigned int nb_accepted_lambda4=0;
    unsigned int nb_accepted_lambda4_batch=0;
    unsigned int nb_accepted_lambda4_after_BI=0;
    unsigned int nb_accepted_gamma=0;
    unsigned int nb_accepted_gamma_batch=0;
    unsigned int nb_accepted_gamma_after_BI=0;
    //    unsigned int accept_candidate=0;
    
    
    // Initialisation of adaptive iterations (batches of 50 iterations)
    unsigned int iter_adaptive=0;
    
    // Starting MCMC algorithm
    
    for(unsigned int iter=1;iter<n_iter;iter++){
        // Writing iteration ID
        
        mat_history.matrix[iter*mat_history.nb_columns]=iter;
        
        if(DEBUG==1){
            cout << "**************" << endl
            << "Iter= " << iter << endl
            << "**************" << endl;
        }
        else{
            if(iter%100==0){
                cout << "**************" << endl
                << "Iter= " << iter << endl
                << "**************" << endl;
            }
        }
        
        
        // Storing previous parameters values
        
        //        double previous_mu=mat_history.matrix[(iter-1)*mat_history.nb_columns+1];
        //        double previous_lambda1=mat_history.matrix[(iter-1)*mat_history.nb_columns+2];
        //        double previous_lambda2=mat_history.matrix[(iter-1)*mat_history.nb_columns+3];
        //        double previous_lambda3=mat_history.matrix[(iter-1)*mat_history.nb_columns+4];
        //        double previous_lambda4=mat_history.matrix[(iter-1)*mat_history.nb_columns+5];
        //        double previous_gamma=mat_history.matrix[(iter-1)*mat_history.nb_columns+6];
        
        // Storing parameters in the vector theta
        gsl_vector *theta_previous=gsl_vector_alloc(N);
        theta_previous->data[0]=mat_history.matrix[(iter-1)*mat_history.nb_columns+1];
        theta_previous->data[1]=mat_history.matrix[(iter-1)*mat_history.nb_columns+2];
        theta_previous->data[2]=log(mat_history.matrix[(iter-1)*mat_history.nb_columns+6]);
        
        
        // Metropolis step: sampling new candidate, computing its likelihood and accepting/rejecting it
        
        unsigned int accept_candidate=metropolis_step_multiv(mat_P_S_I,
                                                             mat_P_S_S,
                                                             mat_exposure,
                                                             mat_P_Ii_Ij_baseline,
                                                             mat_P_Ii_Ij,
                                                             mat_P_SI_M,
                                                             mat_P_SI_R,
                                                             mat_P_SI_SI,
                                                             mat_H_k_per_ind,
                                                             mat_age,
                                                             mat_specs,
                                                             vect_time_line,
                                                             vect_DOR,
                                                             vect_AaS,
                                                             vect_AaQ,
                                                             vect_DOB,
                                                             vect_DOD,
                                                             vect_DODiag,
                                                             summary_SI_SI,
                                                             summary_SI_R,
                                                             summary_SI_M,
                                                             mat_history,
                                                             iter,
                                                             theta_previous,
                                                             L);
        
        
        // Updating acceptance count for this parameter
        
        if(accept_candidate==1){
            nb_accepted_theta++;
            nb_accepted_theta_batch++;
            if(iter>burn_in){
                nb_accepted_theta_after_BI++;
            }
        }
        
        
        // Storing the covariance values
        
        for(int k=0;k<N*N;k++){
            mat_history.matrix[iter*mat_history.nb_columns+8+k]=Sigma->data[k];
        }
        
        
        // Updating scaling of the proposal every n_iter_adaptive iterations
        
        if (adaptive==1){
            if(iter%update_frequency==0){
                gsl_matrix* hat_sigma=gsl_matrix_alloc(3,3);
                sigma_mu=adaptive_step_multiv(iter,
                                              n_iter_adaptive,
                                              mat_history,
                                              hat_sigma);
                
                // Scaling factor
                
                for (int i=0;i<(N*N);i++){
                    L->data[i]=cd*hat_sigma->data[i];
                    Sigma->data[i]=cd*hat_sigma->data[i];
                }

                // Adding a constant to ensure the matrix is positive definite

                for (int i=0;i<N;i++){
                    L->data[i * L->tda + i]+=0.00001;
                    Sigma->data[i * Sigma->tda + i]+=0.00001;
                }
                
                if (DEBUG==0){
                    cout << "Updated covariance:" << endl;
                    display_gsl_matrix(Sigma);
                }
                 
                gsl_linalg_cholesky_decomp(L);
                
                if (DEBUG==0){
                    cout<<"After Cholesky decomposition:"<<endl;
                    display_gsl_matrix(L);
                }
                
                nb_accepted_mu_batch=0;
            }
        }
        
        for(unsigned int col=0;col<mat_history.nb_columns;col++){
            f_out << mat_history.matrix[iter*mat_history.nb_columns+col]
            << "\t";
        }
        
        
        // Close line in dynamic history file
        f_out << endl;
        
        if(DEBUG==1){
            cout << endl;
            cout << "nb_accepted_mu " <<  nb_accepted_mu
            << " -- nb_accepted_lambda1 " << nb_accepted_lambda1
            << " -- nb_accepted_lambda2 " << nb_accepted_lambda2
            << " -- nb_accepted_lambda3 " << nb_accepted_lambda3
            << " -- nb_accepted_lambda4 " << nb_accepted_lambda4
            << " -- nb_accepted_gamma " << nb_accepted_gamma
            << endl;
            
            cout << "nb_accepted_mu_after_BI " <<  nb_accepted_mu_after_BI
            << " -- nb_accepted_lambda1_after_BI " << nb_accepted_lambda1_after_BI
            << " -- nb_accepted_lambda2_after_BI " << nb_accepted_lambda2_after_BI
            << " -- nb_accepted_lambda3_after_BI " << nb_accepted_lambda3_after_BI
            << " -- nb_accepted_lambda4_after_BI " << nb_accepted_lambda4_after_BI
            << " -- nb_accepted_gamma_after_BI " << nb_accepted_gamma_after_BI
            << endl
            << endl
            << endl;
        }
    }//end of for iter
    
    cout << "//////////////////////////////////////////////////////////////////" << endl
    << "//////////////////////////////////////////////////////////////////" << endl
    << "                      MCMC estimation done" << endl
    << "\tAcceptance rate: " << 100.0*(double)(nb_accepted_theta)/(double)(n_iter-1) << endl
    << "\tAcceptance rate (after BI): " << 100.0*(double)(nb_accepted_theta_after_BI)/(double)(n_iter-burn_in-1) << endl
    << "\tAdapted sigma_mu: " << sigma_mu << endl
    << "\tAdapted sigma_lambda1: " << sigma_lambda1 << endl
    << "\tAdapted sigma_lambda2: " << sigma_lambda2 << endl
    << "\tAdapted sigma_lambda3: " << sigma_lambda3 << endl
    << "\tAdapted sigma_lambda4: " << sigma_lambda4 << endl
    << "\tAdapted sigma_gamma: " << sigma_gamma << endl;
    
    cout << "//////////////////////////////////////////////////////////////////" << endl
    << "//////////////////////////////////////////////////////////////////" << endl;
    
    // Close dynamic history file
    f_out.close();
    
    cout << "//////////////////////////////////////////////////////////////////" << endl
    << "//////////////////////////////////////////////////////////////////" << endl
    << "                      Post-processing" << endl;
    
    
    // Post-processing: getting the posterior means after different numbers of iterations and computing corresponding likelihoods
    
    bool stop_iter=false;
    unsigned int current_iter=0;
    unsigned int max_BI=100000;
    
    while(!stop_iter){
        // Writing in log file
        cout << "\titer" << current_iter << endl;
        
        
        // Starting at current_BI iterations
        unsigned int current_BI=1000*current_iter;
        
        
        if(current_BI>=n_iter || current_BI>max_BI){
            stop_iter=true;
        } else {
            // Initialisation of the posterior means
            
            double post_mean_mu=0.0;
            double post_mean_lambda1=0.0;
            double post_mean_lambda2=0.0;
            double post_mean_lambda3=0.0;
            double post_mean_lambda4=0.0;
            double post_mean_gamma=0.0;
            
            
            // Initialisation of the iterations count
            
            double count=0.0;
            
            
            // Computing the posterior means
            
            // Sum of the visited values
            for(unsigned int loop=current_BI;loop<n_iter;loop++){
                count+=1.0;
                post_mean_mu+=mat_history.matrix[loop*mat_history.nb_columns+1];
                post_mean_lambda1+=mat_history.matrix[loop*mat_history.nb_columns+2];
                post_mean_lambda2+=mat_history.matrix[loop*mat_history.nb_columns+3];
                post_mean_lambda3+=mat_history.matrix[loop*mat_history.nb_columns+4];
                post_mean_lambda4+=mat_history.matrix[loop*mat_history.nb_columns+5];
                post_mean_gamma+=mat_history.matrix[loop*mat_history.nb_columns+6];
            }
            
            // Division by number of iterations
            post_mean_mu/=count;
            post_mean_lambda1/=count;
            post_mean_lambda2/=count;
            post_mean_lambda3/=count;
            post_mean_lambda4/=count;
            post_mean_gamma/=count;
            
            
            // Computing the corresponding likelihood
            
            double post_L=wrapped_L_calculation_5_Params_model5(mat_P_S_I,
                                                                mat_P_S_S,
                                                                mat_exposure,
                                                                mat_P_Ii_Ij_baseline,
                                                                mat_P_Ii_Ij,
                                                                mat_P_SI_M,
                                                                mat_P_SI_R,
                                                                mat_P_SI_SI,
                                                                mat_H_k_per_ind,
                                                                mat_age,
                                                                mat_specs,
                                                                vect_time_line,
                                                                vect_DOR,
                                                                vect_AaS,
                                                                vect_AaQ,
                                                                vect_DOB,
                                                                vect_DOD,
                                                                vect_DODiag,
                                                                summary_SI_SI,
                                                                summary_SI_R,
                                                                summary_SI_M,
                                                                post_mean_mu,
                                                                post_mean_lambda1,
                                                                post_mean_lambda2,
                                                                post_mean_lambda3,
                                                                post_mean_lambda4,
                                                                post_mean_gamma);
            
            
            // Writing in the dynamic posterior file
            
            f_outL << setprecision(20)
            << current_BI <<"\t"
            << post_mean_mu <<"\t"
            << post_mean_lambda1 <<"\t"
            << post_mean_lambda2 <<"\t"
            << post_mean_lambda3 <<"\t"
            << post_mean_lambda4 <<"\t"
            << post_mean_gamma <<"\t"
            << post_L << endl;
        }
        // Incrementing iteration
        current_iter++;
    }
    
    // Writing in the log file
    cout << "//////////////////////////////////////////////////////////////////" << endl
    << "//////////////////////////////////////////////////////////////////" << endl;
    
    
    // Closing dynamic posterior file
    f_outL.close();
    
    
    // Removing all objects
    
    mat_history.Free_double_matrix_cont();
    mat_P_SI_M.Free_double_matrix_cont();
    mat_P_S_I.Free_double_matrix_cont();
    mat_P_S_S.Free_double_matrix_cont();
    mat_P_Ii_Ij_baseline.Free_double_matrix_cont();
    mat_P_Ii_Ij.Free_double_matrix_cont();
    mat_H_k_per_ind.Free_double_matrix_cont();
    mat_P_SI_R.Free_double_matrix_cont();
    mat_P_SI_SI.Free_double_matrix_cont();
    mat_exposure.Free_double_matrix_cont();
    gsl_vector_free(vect_time_line);
    gsl_vector_free(vect_DOB);
    gsl_vector_free(vect_DOR);
    gsl_vector_free(vect_DOD);
    gsl_vector_free(vect_AaS);
    gsl_vector_free(vect_AaQ);
    
    mat_specs.Free_double_matrix_cont();
    gsl_vector_free(vect_DODiag);
    
    for(unsigned int ind=0;ind<n_ind;ind++){
        summary_SI_SI[ind].clear();
        summary_SI_R[ind].clear();
        summary_SI_M[ind].clear();
        summary_R_R[ind].clear();
    }
    
    summary_SI_SI.clear();
    summary_SI_R.clear();
    summary_SI_M.clear();
    summary_R_R.clear();
    mat_age.Free_double_matrix_cont();
    
    return 0;
}

