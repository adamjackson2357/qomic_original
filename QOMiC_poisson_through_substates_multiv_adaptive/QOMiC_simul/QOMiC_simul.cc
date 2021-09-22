/* This file is part of QOMiC.
 *      Copyright (c) Marc Chadeau-Hyam (m.chadeau@imperial.ac.uk)
 *      2013
 *
 * QOMiC is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * QOMiC is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with QOMiC.  If not, see <http://www.gnu.org/licenses/>.
 */
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <ctype.h>
#include <iostream>
#include <vector>
#include <algorithm>
#include <memory>
#include <fstream>
#include <../Routines/manip_dates.h>
#include <../Routines/dyn_name.h>
#include <../Routines/matrix_handling.h>
#include <../Routines/rand.h>
#include "../Routines/get_input.h"
#include "../Routines/get_pbties.h"
#include <../Classes/Double_Matrices_cont.h>
#include <../Classes/Int_Matrices_cont.h>
#include <gsl/gsl_fit.h>
#include <gsl/gsl_multifit.h>
#include <gsl/gsl_cdf.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_sort.h>
#include <gsl/gsl_sort_vector.h>
#include <gsl/gsl_permutation.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_blas.h>

#define DEBUG 1

#define first_year_mort 1980

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
    char filename_in_dates[1000];
    char filename_in_specs[1000];
    char filename_in_exposure[1000];
    char filename_in_mortalityF[1000];
    char filename_in_mortalityM[1000];
    char filename_in_MCMC[1000];
    char filename_scenario[1000];
    
    //   char path_name_out[1000];
    
    long MY_SEED=0;
    
    
    // Initialisation of parameters (default values)
    
    unsigned int n_iter=0;
    
    unsigned int Nb_states_in_I=0;
    unsigned int nb_iter=0;
    unsigned int burn_in=0;
    
    na++;
    while(na < argc){
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
        else if ( 0 == strcmp(argv[na],"-iter") ){ // Number of simulations
            n_iter=atoi(argv[++na]);
            if (na+1==argc) break;
            na++;
        }
        else if ( 0 == strcmp(argv[na],"-MCMC") ){
            strcpy(filename_in_MCMC,argv[++na]);
            if (na+1==argc) break;
            na++;
        }
        else if (0 == strcmp(argv[na],"-scenario")){
            strcpy(filename_scenario,argv[++na]);
            if (na+1==argc) break;
            na++;
        }
        else if ( 0 == strcmp(argv[na],"-k") ){
            Nb_states_in_I=atoi(argv[++na]);
            if (na+1==argc) break;
            na++;
        }
        else if ( 0 == strcmp(argv[na],"-sweeps") ){ // Number of iterations in QOMiC output
            nb_iter=atoi(argv[++na]);
            if (na+1==argc) break;
            na++;
        }
        else if ( 0 == strcmp(argv[na],"-burn_in") ){ // Burn-in of QOMiC output
            burn_in=atoi(argv[++na]);
            if (na+1==argc) break;
            na++;
        }
        else{
            cout << "Unknown option: " << argv[na] << endl;
            exit(1);
        }
    }
    
    
    // Dynamic writing of output files
    
    // Summary file (with numbers of transitions)
    string Extension=".txt";
    string OutputName_Summary=Get_scen_name("./Results/Summary_simul",
                                            filename_scenario,
                                            Extension);
    
    ofstream f_out_summary;
    f_out_summary.open(OutputName_Summary.c_str(),ios::out);
    if(f_out_summary.fail()){
        cout << "Invalid Path and/or permission rights for " << OutputName_Summary << " -- run stopped." << endl;
        exit(1);
    }
    else{
        f_out_summary << "Ind\tCaCo\tSmok\tGender\tnSI\tpSI\tnIR\tpIR\tDelta_t_cum\tmean_IP"<<endl;
    }
    
    
    // File with mean probability of S-I transition for all individuals and years
    string OutputName_PSI=Get_scen_name("./Results/Mean_PSI",
                                        filename_scenario,
                                        Extension);
    ofstream f_out_PSI;
    f_out_PSI.open(OutputName_PSI.c_str(),ios::out);
    if(f_out_PSI.fail()){
        cout << "Invalid Path and/or permission rights for " << OutputName_PSI << " -- run stopped." << endl;
        exit(1);
    }
    else{
        cout <<  OutputName_PSI << endl;
    }
    
    smyrand((long)(MY_SEED));
    
    
    // Reading dates matrix:
    
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
    
    
    // Reading Specs matrix:
    
    //Col 1: Gender
    //Col 2: Ca/Co
    //Col 3: Smok_status
    //Col 4: Vital_status
    //Col 5: Age_at_start
    //Col 6: Age_at_quitting
    Double_Matrices_cont mat_specs;
    mat_specs.Read_from_file(filename_in_specs);
    gsl_matrix *mat_specs_work=Double_matrices_cont_2_gsl_matrix(mat_specs);
    unsigned int col_AaS=4;
    gsl_vector *vect_AaS=get_one_vect(col_AaS,
                                      mat_specs_work);
    unsigned int col_AaQ=5;
    gsl_vector *vect_AaQ=get_one_vect(col_AaQ,
                                      mat_specs_work);
    gsl_matrix_free(mat_specs_work);
    
    Double_Matrices_cont mat_exposure;
    mat_exposure.Read_from_file(filename_in_exposure);
    
    
    gsl_vector *vect_time_line=get_Time_Line(vect_DOB,
                                             vect_DODiag,
                                             vect_DOD);
    unsigned int n_years_follow_up=vect_time_line->size;
    unsigned int n_ind=mat_exposure.nb_rows;
    cout << "N_year_follow-up= " << n_years_follow_up
    << " -- n_ind " << n_ind
    << endl;
    
    Double_Matrices_cont mat_mort_female;
    mat_mort_female.Read_from_file(filename_in_mortalityF);
    
    Double_Matrices_cont mat_mort_male;
    mat_mort_male.Read_from_file(filename_in_mortalityM);
    
    Double_Matrices_cont mat_MCMC;
    mat_MCMC.Read_from_file(filename_in_MCMC);
    if(DEBUG==1){
        cout << "MCMC matrix" << endl;
        mat_MCMC.Display_matrix();
    }
    
    Double_Matrices_cont mat_age;
    mat_age.Alloc_double_matrix_cont(n_ind,
                                     n_years_follow_up);
    get_matrix_age_cont(mat_age,
                        vect_time_line,
                        vect_DOB);
    
    
    // Initialisation of probabilities
    
    // Step 1: Calculating P_SI_M=P_S_M=P_I_M other cause mortality (fixed, based on input mortality data)
    
    Double_Matrices_cont mat_P_SI_M;
    mat_P_SI_M.Alloc_double_matrix_cont(n_ind,
                                        n_years_follow_up);
    
    get_P_SI_M_VBT(mat_P_SI_M,
                   mat_mort_female,
                   mat_mort_male,
                   mat_specs,
                   mat_age,
                   vect_time_line,
                   vect_DOR,
                   vect_AaS,
                   vect_AaQ);
    
    mat_mort_female.Free_double_matrix_cont();
    mat_mort_male.Free_double_matrix_cont();
    
    
    // Step 2: Initialisation of other transition probabilities
    
    Double_Matrices_cont mat_P_S_I;
    mat_P_S_I.Alloc_double_matrix_cont(n_ind,
                                       n_years_follow_up);
    
    // Matrix storing all P_S_I for computation of the mean P_S_I
    Double_Matrices_cont mat_P_S_I_cum;
    mat_P_S_I_cum.Alloc_double_matrix_cont(n_ind,
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
    
    
    // Initialisation of the vectors of transition counts
    
    vector < unsigned > vect_nSI;
    vect_nSI.resize(n_ind);
    vector < unsigned > vect_nIR;
    vect_nIR.resize(n_ind);
    vector < unsigned > vect_nDeltaT_cum;
    vect_nDeltaT_cum.resize(n_ind);
    
    for(unsigned int iter=0;iter<n_iter;iter++){
        // Sample a random simulation ID from estimation file (after burn-in)
        
        unsigned int current_sim=(unsigned int)(myrandRange(burn_in,nb_iter-1));
        
        
        // Extract the parameters of the sampled simulation
        
        double current_mu=mat_MCMC.matrix[current_sim*mat_MCMC.nb_columns+1];
        double current_lambda0=mat_MCMC.matrix[current_sim*mat_MCMC.nb_columns+2];
        double current_lambda1=mat_MCMC.matrix[current_sim*mat_MCMC.nb_columns+3];
        double current_lambda2=mat_MCMC.matrix[current_sim*mat_MCMC.nb_columns+4];
        double current_lambda3=mat_MCMC.matrix[current_sim*mat_MCMC.nb_columns+5];
        double current_gamma=mat_MCMC.matrix[current_sim*mat_MCMC.nb_columns+6];
        
        
        // Computing the probability of S-I transition given these parameters
        
        get_P_S_I_5_param_model5(mat_P_S_I,
                                 mat_exposure,
                                 mat_P_SI_M,
                                 mat_age,
                                 mat_specs,
                                 vect_time_line,
                                 vect_DOR,
                                 vect_AaS,
                                 vect_AaQ,
                                 current_mu,
                                 current_lambda0,
                                 current_lambda1,
                                 current_lambda2,
                                 current_lambda3);
        
        
        // Storing the P_S_I values for all individuals and time points
        
        for(unsigned int MyRow=0;MyRow<mat_P_S_I.nb_rows;MyRow++){
            for(unsigned int MyCol=0;MyCol<mat_P_S_I.nb_columns;MyCol++){
                mat_P_S_I_cum.matrix[MyRow*mat_P_S_I_cum.nb_columns+MyCol]+=mat_P_S_I.matrix[MyRow*mat_P_S_I_cum.nb_columns+MyCol];
            }
        }
        
        
        // Computing the probability of S-S transition
        
        get_P_S_S(mat_P_S_S,
                  mat_P_SI_M,
                  mat_P_S_I);
        
        
        // Computing the transitions through sub-states of I
        
        get_P_Ii_Ij_baseline(mat_P_Ii_Ij_baseline,
                             current_gamma);
        
        if(DEBUG==1){
            cout << "gamma: " << current_gamma << endl;
        }
        
        
        // Simulate S to I transitions for each individual and time point
        
        for(unsigned int ind=0;ind<n_ind;ind++){
            if(DEBUG==1){
                cout << "Individual: " << ind << endl;
            }
            
            // Initialisation
            
            unsigned int year=0;
            unsigned int Stop_simul=0;
            vector < double > vect_S_trans;
            vector < double > vect_I_trans;
            unsigned int year_SI=0;
            
            
            // Loop over the years while staying in S
            while(Stop_simul==0 && year<n_years_follow_up){
                if(DEBUG==1){
                    cout << "Year: " << year << endl;
                }
                
                // Extract the current probabilities of transition
                
                double current_P_S_S=mat_P_S_S.matrix[ind*mat_P_S_S.nb_columns+year];
                double current_P_S_I=mat_P_S_I.matrix[ind*mat_P_S_I.nb_columns+year];
                double current_P_S_M=mat_P_SI_M.matrix[ind*mat_P_SI_M.nb_columns+year];
                
                if(DEBUG==1){
                    cout << "Proability of transition S-S: " << current_P_S_S << endl;
                    cout << "Proability of transition S-I: " << current_P_S_I << endl;
                    cout << "Proability of transition S-M: " << current_P_S_M << endl;
                }
                
                
                // Storing all three probabilities in vect_S_trans
                
                vect_S_trans.push_back(current_P_S_S);
                vect_S_trans.push_back(current_P_S_I);
                vect_S_trans.push_back(current_P_S_M);
                
                
                // Sampling the state based on transition probabilities
                
                int sampled_state=SampleFromDiscrete_non_cum(vect_S_trans);
                
                if(DEBUG==1){
                    cout << "Sampled state: " << sampled_state << endl;
                }
                
                if(sampled_state==0){ // Staying in S
                    year++;
                } else { // Leaving S
                    // Stoping the current while loop
                    Stop_simul=1;
                    
                    if(sampled_state==1){
                        // Storing the year of S-I transition for this individual
                        year_SI=year;
                        
                        // Counting the number of S-I transitions for this individual
                        vect_nSI[ind]+=1.0;
                    }
                }
                
                // Emptying the vector vect_S_trans
                vect_S_trans.clear();
            }//end of while year
            
            
            // If transition S-I has happened for this individual
            if(year_SI!=0){
                // Get year of S-I transition
                unsigned int current_year=year_SI;
                
                // Get current sub-state of I (0 means S)
                unsigned int current_substate=0;
                
                while(current_substate<mat_P_Ii_Ij.nb_columns-2 && current_year<n_years_follow_up){
                    if(DEBUG==1){
                        cout << "Year: " << current_year << endl;
                    }
                    
                    // Computing the probability of transitions through the sub-states of I
                    
                    update_P_Ii_Ij(mat_P_Ii_Ij,
                                   mat_P_Ii_Ij_baseline,
                                   mat_P_SI_M,
                                   ind,
                                   current_year-1);
                    
                    
                    // Storing the probabilities of transitions across sub-states of I
                    
                    if(DEBUG==1){
                        cout << "Maximum number of sub-states: " << mat_P_Ii_Ij.nb_columns << endl;
                    }
                    
                    for(unsigned int substate=0;substate<mat_P_Ii_Ij.nb_columns;substate++){
                        if(DEBUG==0){
                            cout << "Substate: " << substate << endl;
                            cout << "Corresponding probability of transition: " << mat_P_Ii_Ij.matrix[current_substate*mat_P_Ii_Ij.nb_columns+substate] << endl;
                        }
                    
                        vect_I_trans.push_back(mat_P_Ii_Ij.matrix[current_substate*mat_P_Ii_Ij.nb_columns+substate]);
                    }
                    
                    
                    // Sampling the state based on transition probabilities
                    
                    current_substate=SampleFromDiscrete_non_cum(vect_I_trans);
                    
                    if(DEBUG==1){
                        cout << "Sampled state: " << current_substate << endl;
                    }
                    
                    if(current_substate < mat_P_Ii_Ij.nb_columns-2){ // Staying in one of the I sub-states
                        current_year++;
                    } else {
                        if(current_substate==mat_P_Ii_Ij.nb_columns-2){ // Transition to R
                            // Count of I-R transitions
                            vect_nIR[ind]+=1;
                            
                            // Time spent in I sub-states (i.e. time to diagnosis)
                            vect_nDeltaT_cum[ind]+=(current_year-year_SI);
                        }
                        // Otherwise: transition to M (not stored)
                    }
                    vect_I_trans.clear();
                }
            }//end of if S_I!=0
        }//end of for ind
    }
    
    
    // Computing the average probability of S-I transition for each individual and time point
    
    for(unsigned int MyRow=0;MyRow<mat_P_S_I.nb_rows;MyRow++){ // For loop over individuals
        for(unsigned int MyCol=0;MyCol<mat_P_S_I.nb_columns;MyCol++){ // For loop over time points
            mat_P_S_I_cum.matrix[MyRow*mat_P_S_I_cum.nb_columns+MyCol]/=n_iter;
            f_out_PSI << mat_P_S_I_cum.matrix[MyRow*mat_P_S_I_cum.nb_columns+MyCol]
            << " ";
        }
        f_out_PSI << endl;
    }
    
    for(unsigned int ind=0;ind<n_ind;ind++){
        // Storing some individual characteristics from the input data
        
        unsigned int current_CaCo=(unsigned int)(mat_specs.matrix[ind*mat_specs.nb_columns+1]);
        unsigned int current_smok=(unsigned int)(mat_specs.matrix[ind*mat_specs.nb_columns+2]);
        unsigned int current_gender=(unsigned int)(mat_specs.matrix[ind*mat_specs.nb_columns]);
        
        
        // Storing the S-I and I-R transitions counts and proportions, cumulative and average time-to-diagnosis over simulations
        
        f_out_summary << ind+1 << "\t"
        << current_CaCo << "\t"
        << current_smok << "\t"
        << current_gender << "\t"
        << vect_nSI[ind] << "\t"
        << (double)(vect_nSI[ind])/(double)(n_iter) << "\t"
        << vect_nIR[ind] << "\t"
        << (double)(vect_nIR[ind])/(double)(n_iter) << "\t"
        << vect_nDeltaT_cum[ind] << "\t"
        << (double)(vect_nDeltaT_cum[ind])/(double)(vect_nIR[ind])
        << endl;
    }
    
    
    // Closing the connections
    
    f_out_summary <<endl;
    f_out_summary.close();
    f_out_PSI.close();
    
    
    // Clearing all objects
    
    vect_nSI.clear();
    vect_nDeltaT_cum.clear();
    vect_nIR.clear();
    
    mat_P_S_I.Free_double_matrix_cont();
    mat_P_S_I_cum.Free_double_matrix_cont();
    mat_P_S_S.Free_double_matrix_cont();
    mat_P_Ii_Ij_baseline.Free_double_matrix_cont();
    mat_P_Ii_Ij.Free_double_matrix_cont();
    mat_exposure.Free_double_matrix_cont();
    mat_MCMC.Free_double_matrix_cont();
    gsl_vector_free(vect_time_line);
    gsl_vector_free(vect_DOB);
    gsl_vector_free(vect_DOD);
    
    mat_specs.Free_double_matrix_cont();
    gsl_vector_free(vect_DODiag);
    mat_age.Free_double_matrix_cont();
    
    return 0;
}

