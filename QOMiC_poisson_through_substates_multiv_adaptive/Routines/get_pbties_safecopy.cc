#include "get_pbties.h"

#define DEBUG 0


void get_P_SI_M_VBT(Double_Matrices_cont mat_PSI_M,
                    Double_Matrices_cont mat_mort_female,
                    Double_Matrices_cont mat_mort_male,
                    Double_Matrices_cont mat_specs,
                    Double_Matrices_cont mat_age,
                    gsl_vector * vect_time_line,
                    gsl_vector * vect_DOR,
                    gsl_vector * vect_AaS,
                    gsl_vector * vect_AaQ)
{
    unsigned int n_ind=mat_PSI_M.nb_rows;
    unsigned int n_years= mat_PSI_M.nb_columns;
    
    for(unsigned int ind=0;ind<n_ind;ind++){
        unsigned int current_gender=(unsigned int)(mat_specs.matrix[ind*mat_specs.nb_columns]);
        unsigned int current_smok=(unsigned int)(mat_specs.matrix[ind*mat_specs.nb_columns+2]);
        unsigned int current_YOR=(unsigned int)(vect_DOR->data[ind]/10000);
        unsigned int current_AaS=(unsigned int)(vect_AaS->data[ind]);
        unsigned int current_AaQ=(unsigned int)(vect_AaQ->data[ind]);
        if(DEBUG==1){
            cout << "Ind " << ind
            << " -- Gender " << current_gender
            << " -- Smoking " << current_smok
            << " -- YOR " << current_YOR
            << " -- AaS " << current_AaS
            << " -- AaQ " << current_AaQ
            << endl;
        }
        for(unsigned int year=current_YOR-vect_time_line->data[0];year<n_years;year++){
            double current_age=mat_age.matrix[ind*mat_age.nb_columns+year];
            double current_year=vect_time_line->data[year];
            unsigned int col_pos=0;
            unsigned int line_pos=(unsigned int)(current_age);
            if(current_smok==1){
                col_pos=0;
            }
            else{
                if(current_smok==3 && current_age>=current_AaS){
                    col_pos=1;
                }
                if(current_smok==2 && (current_age>=current_AaS && current_age<=current_AaQ)){
                    col_pos=1;
                }
            }
            if(current_gender==2){//Gender=2; Female
                mat_PSI_M.matrix[ind*mat_PSI_M.nb_columns+year]=mat_mort_female.matrix[line_pos*mat_mort_female.nb_columns+col_pos];
                
            }
            else{
                mat_PSI_M.matrix[ind*mat_PSI_M.nb_columns+year]=mat_mort_male.matrix[line_pos*mat_mort_male.nb_columns+col_pos];
            }
            
            if(DEBUG==1){
                cout << "\tyear " << year
                << " -- date " << current_year
                << " -- col_pos " << col_pos
                << " -- age " << current_age
                << " -- line_pos " << line_pos
                << " -- mort " << mat_PSI_M.matrix[ind*mat_PSI_M.nb_columns+year]
                << endl;
            }
        }//end of for year
    }//end of for ind
}


void get_P_S_S(Double_Matrices_cont mat_P_S_S,
               Double_Matrices_cont mat_P_SI_M,
               Double_Matrices_cont mat_P_S_I)
{
    unsigned int n_ind=mat_P_S_I.nb_rows;
    unsigned int n_years= mat_P_S_I.nb_columns;
    
    for(unsigned int ind=0;ind<n_ind;ind++){
        if(DEBUG==1){
            cout << "Ind " << ind
            << endl;
        }
        
        for(unsigned int year=0;year<n_years;year++){
            mat_P_S_S.matrix[ind*mat_P_S_S.nb_columns+year]=1.0-(mat_P_SI_M.matrix[ind*mat_P_SI_M.nb_columns+year]+
                                                                 mat_P_S_I.matrix[ind*mat_P_S_I.nb_columns+year]);
            if(mat_P_S_S.matrix[ind*mat_P_S_S.nb_columns+year]<0.0){
                cout << "WARNING: PEM+PEI>1.0 run stopped" << endl;
                exit(1);
            }
        }
    }
}



//void get_P_Ii_Ij_baseline(Double_Matrices_cont mat_Ii_Ij_baseline,
//                          double gamma)
//{
//    unsigned int n_states=mat_Ii_Ij_baseline.nb_rows;
//    //   cout << "Ind " << ind
//    //        << " -- year " << year
//    //        << " -- tmp mort " << tmp_mort
//    //        << endl;
//
//    double tmp_pdf_comp=gsl_cdf_gamma_Q(1,1,gamma);
//    for(unsigned int state1=0;state1<n_states-2;state1++){
//        double tmp_sum=0.0;
//        mat_Ii_Ij_baseline.matrix[state1*n_states+(n_states-1)]=0.0;
//        for(unsigned int state2=state1+1;state2<n_states-1;state2++){
//            double tmp_pdf=gsl_cdf_gamma_P(1,state2-state1,gamma);
//            //double tmp_pdf=gsl_ran_poisson_pdf(state2-state1,gamma);
//            mat_Ii_Ij_baseline.matrix[state1*n_states+state2]=tmp_pdf;
//            tmp_sum+=tmp_pdf;
//        }
//        //cout << "tmp_sum" << tmp_sum << endl;
//        for(unsigned int state2=state1+1;state2<n_states-1;state2++){
//            mat_Ii_Ij_baseline.matrix[state1*n_states+state2]*=((1- tmp_pdf_comp)/tmp_sum);
//        }
//        mat_Ii_Ij_baseline.matrix[state1*n_states+state1]=tmp_pdf_comp;
//        mat_Ii_Ij_baseline.matrix[n_states*n_states-1]=1.0;
//        mat_Ii_Ij_baseline.matrix[(n_states-1)*n_states-2]=1.0;
//    }
//}


void get_P_Ii_Ij_baseline(Double_Matrices_cont mat_P_Ii_Ij_baseline,
                          double gamma)
{
    unsigned int n_states=mat_P_Ii_Ij_baseline.nb_rows;
    
    for(unsigned int state1=0;state1<n_states-2;state1++){
        double tmp_sum=0.0;
        mat_P_Ii_Ij_baseline.matrix[state1*n_states+(n_states-1)]=0.0; // Mortality is included later
        
        // Probability of having exactly (state2-state1) transitions in Poisson process with intensity 1/gamma
        for(unsigned int state2=state1;state2<n_states-1;state2++){ // For I_1, state2 is I_1, I_2, ..., I_K or R
            double tmp_poisson=exp(-1/gamma)*pow(1/gamma, state2-state1)/factorial(state2-state1);
            mat_P_Ii_Ij_baseline.matrix[state1*n_states+state2]=tmp_poisson;
            tmp_sum+=tmp_poisson;
            
            if (DEBUG==1){
                cout << "Proba: " << tmp_poisson << endl;
                cout << "Cumulative proba: " << tmp_sum << endl;
            }
        }
        
        // Condition the transition probability on the maximum number of transitions (K) using cumulative probability
        for(unsigned int state2=state1;state2<n_states-1;state2++){ // For I_1, state2 is I_1, I_2, ..., I_K or R
            mat_P_Ii_Ij_baseline.matrix[state1*n_states+state2]*=1/tmp_sum;
        }
    }
    
    // R and M are absorbing states
    mat_P_Ii_Ij_baseline.matrix[n_states*n_states-1]=1.0;
    mat_P_Ii_Ij_baseline.matrix[(n_states-1)*n_states-2]=1.0;
    
    if (DEBUG==1){
        mat_P_Ii_Ij_baseline.Display_matrix_header();
    }
}


void update_P_Ii_Ij(Double_Matrices_cont mat_Ii_Ij,
                    Double_Matrices_cont mat_Ii_Ij_baseline,
                    Double_Matrices_cont mat_P_SI_M,
                    unsigned int ind,
                    unsigned int year)
{
    mat_Ii_Ij.Reset_double_matrix_cont();
    unsigned int n_states=mat_Ii_Ij.nb_rows;
    double tmp_mort=mat_P_SI_M.matrix[ind*mat_P_SI_M.nb_columns+year];
    if(DEBUG==1){
        cout << "tmp_mort " << tmp_mort << endl;
    }
    for(unsigned int state1=0;state1<n_states-2;state1++){
        mat_Ii_Ij.matrix[state1*n_states+(n_states-1)]=tmp_mort;
        for(unsigned int state2=state1;state2<n_states-1;state2++){
            mat_Ii_Ij.matrix[state1*n_states+state2]=mat_Ii_Ij_baseline.matrix[state1*n_states+state2]*(1.0-tmp_mort);
            
        }
    }
    mat_Ii_Ij.matrix[n_states*n_states-1]= mat_Ii_Ij_baseline.matrix[n_states*n_states-1];
    mat_Ii_Ij.matrix[(n_states-1)*n_states-2]= mat_Ii_Ij_baseline.matrix[n_states*n_states-1];
}


void get_P_SI_R_and_H_k(Double_Matrices_cont mat_P_SI_R,
                        Double_Matrices_cont mat_P_SI_SI,
                        Double_Matrices_cont mat_H_k_per_ind,
                        Double_Matrices_cont mat_P_S_S,
                        Double_Matrices_cont mat_P_S_I,
                        Double_Matrices_cont mat_P_SI_M,
                        Double_Matrices_cont mat_P_Ii_Ij,
                        Double_Matrices_cont mat_P_Ii_Ij_baseline,
                        gsl_vector * vect_time_line,
                        gsl_vector *vect_DOB,
                        gsl_vector *vect_DOR,
                        gsl_vector *vect_DOD,
                        gsl_vector *vect_DODiag,
                        Double_Matrices_cont mat_specs)
{
    unsigned int n_states=mat_H_k_per_ind.nb_rows;
    unsigned int n_ind=mat_P_S_S.nb_rows;
    unsigned int n_years=mat_H_k_per_ind.nb_columns;
    
    for(unsigned int ind=0;ind<n_ind;ind++){
        // Initialisation of mat_H_k_per_ind
        mat_H_k_per_ind.Reset_double_matrix_cont();
        
        unsigned int first_year_SI_SI=0;
//        unsigned int n_yearsbis=summary_SI_SI[ind].size();
//        unsigned int first_year_SI_SI=summary_SI_SI[ind][1];
//        unsigned int last_year_SI_SI=summary_SI_SI[ind][n_yearsbis-1];
        unsigned int current_CaCo=(unsigned int)(mat_specs.matrix[ind*mat_specs.nb_columns+1]);
        unsigned int current_Vit_stat=(unsigned int)(mat_specs.matrix[ind*mat_specs.nb_columns+5]);
//
//        if(DEBUG==1){
//            cout << "Ind " << ind
//            << " -- CaCo " << current_CaCo
//            << " -- Vit_stat " << current_Vit_stat
//            << " -- n_yearsbis " << n_yearsbis
//            << " -- first_year_SI_SI " << first_year_SI_SI
//            << " -- DOR " <<  trunc(vect_DOR->data[ind]/10000)
//            << " -- DOB " <<  vect_time_line->data[first_year_SI_SI]
//            << " -- DOBInput " <<  trunc(vect_DOB->data[ind]/10000)
//            << " -- last_year_SI_SI " << last_year_SI_SI;
//            if(current_CaCo==0 && current_Vit_stat==2){
//                cout << " -- DOD " <<  vect_time_line->data[last_year_SI_SI+1]
//                << " -- DODInput " <<  trunc(vect_DOD->data[ind]/10000);
//            }
//            else if(current_CaCo==1){
//                cout << " -- DODiag " <<  vect_time_line->data[last_year_SI_SI+1]
//                << " -- DODiagInput " <<  trunc(vect_DODiag->data[ind]/10000);
//            }
//            else{
//                cout << " -- End in SI: " << vect_time_line->data[last_year_SI_SI];
//            }
//            cout << endl;
//        }
//
        
        // Filling H_S at starting year
        
        // H_S(first_year)=1.0
        mat_H_k_per_ind.matrix[first_year_SI_SI]=1;
        
        // Filling H_k at starting year
        
        // H_k(first_year)=0.0
        for(unsigned int k=1;k<n_states;k++){
            mat_H_k_per_ind.matrix[k*mat_H_k_per_ind.nb_columns+first_year_SI_SI]=0.0;
        }
        
        if(DEBUG==1){
            cout << "\tyear " <<  first_year_SI_SI
            << " -- date " << vect_time_line->data[first_year_SI_SI];
            
            for(unsigned int k=0;k<n_states;k++){
                //cout << " -- H_" << k << "(t-1) " << mat_H_k_per_ind.matrix[k* mat_H_k_per_ind.nb_columns+first_year_SI_SI-1]
                cout << " -- H_" << k << "(t) " << mat_H_k_per_ind.matrix[k* mat_H_k_per_ind.nb_columns+first_year_SI_SI];
            }
            cout << endl;
            cout << "P_S_S" << endl;
            for(unsigned int year=0;year<n_years;year++){
                cout << mat_P_S_S.matrix[ind*mat_P_S_S.nb_columns+year]
                << " ";
            }
            cout << endl;
            cout << "P_S_I" << endl;
            for(unsigned int year=0;year<n_years;year++){
                cout << mat_P_S_I.matrix[ind*mat_P_S_S.nb_columns+year]
                << " ";
            }
            cout << endl;
            cout << "P_SI_M" << endl;
            for(unsigned int year=0;year<n_years;year++){
                cout << mat_P_SI_M.matrix[ind*mat_P_S_S.nb_columns+year]
                << " ";
            }
            cout << endl;
        }
        
        
        // Probability of remaining in SI at (first year + 1) is probability of not dying
        mat_P_SI_SI.matrix[ind*mat_P_SI_R.nb_columns+first_year_SI_SI]=1.0-mat_P_SI_M.matrix[ind*mat_P_SI_M.nb_columns+first_year_SI_SI];
        
        
        // Impossible to be in R at (first year + 1) because need to enter in I sub-states first
        mat_P_SI_R.matrix[ind*mat_P_SI_R.nb_columns+first_year_SI_SI]=0.0;
        
        
        double tmp_prod=1.0;
        for(unsigned int year=first_year_SI_SI+1;year<n_years;year++){
            if(DEBUG==1){
                cout << "***********************************************" << endl;
                cout << "year " <<  year
                << " -- date " << vect_time_line->data[year]
                << endl;
            }
            
            
            // Computing H_S: first line
            
            tmp_prod*=mat_P_S_S.matrix[ind*mat_P_S_S.nb_columns+year-1];
            mat_H_k_per_ind.matrix[year]=tmp_prod;
            
            if(DEBUG==1){
                cout << "P_S_S " << mat_P_S_S.matrix[ind*mat_P_S_S.nb_columns+year]
                << " -- previous H_S " << mat_H_k_per_ind.matrix[year-1]
                << " -- H_S " << mat_H_k_per_ind.matrix[year]
                << endl;
                
                cout << endl
                << "previous PII baseline" << endl;
                mat_P_Ii_Ij_baseline.Display_matrix_header();
                
                cout << endl;
                cout << "previous P_SI_M " << mat_P_SI_M.matrix[ind*mat_P_SI_M.nb_columns+year-1] << endl;
            }
            
            
            // Computing H_1 (H_I1): second line
            
            update_P_Ii_Ij(mat_P_Ii_Ij,
                           mat_P_Ii_Ij_baseline,
                           mat_P_SI_M,
                           ind,
                           year-1);
            
            if(DEBUG==1){
                cout << endl
                << "previous PII updated" << endl;
                mat_P_Ii_Ij.Display_matrix_header();
            }
            
            double previous_P_S_I=mat_P_S_I.matrix[ind*mat_P_S_I.nb_columns+year-1];
            double previous_P_1_1=mat_P_Ii_Ij.matrix[0];
            double previous_H_S=mat_H_k_per_ind.matrix[year-1];
            double previous_H_1=mat_H_k_per_ind.matrix[mat_H_k_per_ind.nb_columns+year-1];
            
            mat_H_k_per_ind.matrix[mat_H_k_per_ind.nb_columns+year]=(previous_H_S*previous_P_S_I+
                                                                     previous_H_1*previous_P_1_1);
            
            if(DEBUG==1){
                cout << "previous P_S_I " << previous_P_S_I
                << " -- previous P_1_1 " << previous_P_1_1
                << " -- previous_H_S " << previous_H_S
                << " -- previous_H_1 " << previous_H_1
                << " -- current H_1 " << mat_H_k_per_ind.matrix[mat_H_k_per_ind.nb_columns+year]
                << endl;
            }
            
            
            // Computing H_k for k between 2 and number of sub-states (H_Ik)
            
            for(unsigned int k1=2;k1<n_states;k1++){
                double tmp_cum=0.0;
                
                if(DEBUG==1){
                    cout << "Filling line " << k1
                    << " -- tmp_cum " << tmp_cum
                    << endl;
                }
                
                
                // Product over all previous sub-states
                
                for(unsigned int k2=1;k2<=k1;k2++){
                    double previous_Hk=mat_H_k_per_ind.matrix[k2*mat_H_k_per_ind.nb_columns+year-1];
                    double previous_P=mat_P_Ii_Ij.matrix[(k2-1)*mat_P_Ii_Ij.nb_columns+(k1-1)];
                    
                    tmp_cum+=previous_Hk*previous_P;
                    
                    if(DEBUG==1){
                        cout << "\tI" << k2
                        << "-I" << k1
                        << " transitions"
                        << " -- previous H " << previous_Hk
                        << " -- previous P " << previous_P
                        << " -- tmp_cum " << tmp_cum
                        << endl;
                    }
                }
                
                mat_H_k_per_ind.matrix[k1*mat_H_k_per_ind.nb_columns+year]=tmp_cum;
                
                if(DEBUG==1){
                    cout << "current H_" << k1
                    << " = " << mat_H_k_per_ind.matrix[k1*mat_H_k_per_ind.nb_columns+year]
                    << endl;
                }
            }//end of for k1
            
            
            // Getting P_Ii_Ij at time year
            
            update_P_Ii_Ij(mat_P_Ii_Ij,
                           mat_P_Ii_Ij_baseline,
                           mat_P_SI_M,
                           ind,
                           year);
            
            
            if(DEBUG==1){
                cout << endl
                << "current PII updated" << endl;
                mat_P_Ii_Ij.Display_matrix_header();
            }
            
            
            // Computing P_SI_R for a given year/ind
            
            double tmp_prod_cum=0.0;
            for(unsigned int k=1;k<n_states;k++){
                double tmp_prod = (mat_H_k_per_ind.matrix[k*mat_H_k_per_ind.nb_columns+year]*
                                   mat_P_Ii_Ij.matrix[(k-1)*mat_P_Ii_Ij.nb_columns+mat_P_Ii_Ij.nb_columns-2]);
                tmp_prod_cum+=tmp_prod;
                
                if(DEBUG==1){
                    cout << " -- current H" << k << " " << mat_H_k_per_ind.matrix[k*mat_H_k_per_ind.nb_columns+year]
                    << " -- current P_Ik_R " <<  mat_P_Ii_Ij.matrix[(k-1)*mat_P_Ii_Ij.nb_columns+mat_P_Ii_Ij.nb_columns-2]
                    << " -- tmp_prod " << tmp_prod
                    << " -- tmp_cum " << tmp_prod_cum
                    << endl;
                    cout << endl;
                }
            }
            //            cout << "Probability of being in SI inter R:" << endl;
            //            cout << tmp_prod_cum << endl;
            
            double current_mort=mat_P_SI_M.matrix[ind*mat_P_SI_M.nb_columns+year];
            
            double tmp=0;
            //            cout << tmp << endl;
            for(unsigned int k=0;k<n_states;k++){
                tmp+=mat_H_k_per_ind.matrix[k*mat_H_k_per_ind.nb_columns+year];
                //                cout << tmp << endl;
            }
            //            cout << "Probability of being in S or in an I sub-state:" << endl;
            //            cout << tmp << endl;
            
            mat_P_SI_R.matrix[ind*mat_P_SI_R.nb_columns+year]=tmp_prod_cum/tmp; // revised formula
            //            cout << "Probability of transition SI-R:" << endl;
            //            cout << mat_P_SI_R.matrix[ind*mat_P_SI_R.nb_columns+year] << endl;
            
            if(DEBUG==1){
                cout << "current P_SI_M: " << current_mort << endl;
                cout << "current P_SI_inter_R: " << tmp_prod_cum << endl;
                cout << "current P(SI): " << tmp << endl;
                cout << "current P_SI_R: " << mat_P_SI_R.matrix[ind*mat_P_SI_R.nb_columns+year] << endl;
            }
            
            mat_P_SI_SI.matrix[ind*mat_P_SI_SI.nb_columns+year]=1.0-(mat_P_SI_R.matrix[ind*mat_P_SI_R.nb_columns+year]+current_mort);
            
//            // Precision issues
//            if(mat_P_SI_SI.matrix[ind*mat_P_SI_R.nb_columns+year]<0.0){
//                mat_P_SI_SI.matrix[ind*mat_P_SI_R.nb_columns+year]=0.0;
//            }
            
            if(mat_P_SI_R.matrix[ind*mat_P_SI_R.nb_columns+year]>1 || mat_P_SI_R.matrix[ind*mat_P_SI_R.nb_columns+year]<0.0){
                cout << "Problem P_SI_R not in [0-1]: " << tmp_prod_cum
                << endl
                << "Run Stopped"
                << endl;
                exit(1);
            }
            
            if(mat_P_SI_SI.matrix[ind*mat_P_SI_R.nb_columns+year]<0.0 || mat_P_SI_SI.matrix[ind*mat_P_SI_R.nb_columns+year]>1.0){
                cout << "Problem P_SI_SI not in [0-1]: " << mat_P_SI_SI.matrix[ind*mat_P_SI_R.nb_columns+year]
                << endl
                << "Run Stopped"
                << endl;
                exit(1);
            }
        }//end of for year
        if(DEBUG==1){
            mat_H_k_per_ind.Display_matrix_header();
        }
    }//end of for ind
}


double get_likelihood(Double_Matrices_cont mat_P_SI_R,
                      Double_Matrices_cont mat_P_SI_SI,
                      Double_Matrices_cont mat_P_SI_M,
                      vector < vector < int > > &summary_SI_SI,
                      vector < vector < int > > &summary_SI_R,
                      vector < vector < int > > &summary_SI_M)
{
    unsigned int n_ind=mat_P_SI_SI.nb_rows;
    double tmp_cum=0.0;
    
    
    for(unsigned int ind=0;ind<n_ind;ind++){
        //SI_SI transitions
        unsigned int n_trans_SI_SI=summary_SI_SI[ind].size();
        unsigned int n_trans_SI_R=summary_SI_R[ind].size();
        unsigned int n_trans_SI_M=summary_SI_M[ind].size();
        
        if(DEBUG==1){
            cout << "Ind " << ind
            << " -- n_SI_SI " << n_trans_SI_SI
            << " -- n_SI_R " << n_trans_SI_R
            << " -- n_SI_M " << n_trans_SI_M
            << endl;
        }
        if(n_trans_SI_SI>1){
            for(unsigned int trans=1;trans<n_trans_SI_SI;trans++){
                unsigned int year_trans=(unsigned int)(summary_SI_SI[ind][trans]);
                double logP=log(mat_P_SI_SI.matrix[ind*mat_P_SI_SI.nb_columns+year_trans]);
                tmp_cum+=logP;
                if(DEBUG==1){
                    cout << "\tYear " << year_trans
                    << " -- P_SI_SI " << mat_P_SI_SI.matrix[ind*mat_P_SI_SI.nb_columns+year_trans]
                    << " -- P_SI_R " << mat_P_SI_R.matrix[ind*mat_P_SI_SI.nb_columns+year_trans]
                    << " -- P_SI_M " << mat_P_SI_M.matrix[ind*mat_P_SI_SI.nb_columns+year_trans]
                    << " -- log " << logP
                    << " -- cum " << tmp_cum
                    << endl;
                }
            }
        }
        if(n_trans_SI_R>1){
            unsigned int year_trans_SI_R=(unsigned int)(summary_SI_R[ind][1]);
            double logP=log(mat_P_SI_R.matrix[ind*mat_P_SI_R.nb_columns+year_trans_SI_R]);
            tmp_cum+=logP;
            if(DEBUG==1){
                cout << "\tYear " << year_trans_SI_R
                << " -- P_SI_R " << mat_P_SI_R.matrix[ind*mat_P_SI_R.nb_columns+year_trans_SI_R]
                << " -- log " << logP
                << " -- cum " << tmp_cum
                << endl;
            }
        }
        if(n_trans_SI_M>1){
            unsigned int year_trans_SI_M=(unsigned int)(summary_SI_M[ind][1]);
            double logP=log(mat_P_SI_M.matrix[ind*mat_P_SI_M.nb_columns+year_trans_SI_M]);
            tmp_cum+=logP;
            if(DEBUG==1){
                cout << "\tYear " << year_trans_SI_M
                << " -- P_SI_M " << mat_P_SI_M.matrix[ind*mat_P_SI_M.nb_columns+year_trans_SI_M]
                << " -- log " << logP
                << " -- cum " << tmp_cum
                << endl;
            }
            
        }
    }
    return tmp_cum;
}


void get_P_S_I_5_param_model5(Double_Matrices_cont mat_P_S_I,
                              Double_Matrices_cont mat_exposure,
                              Double_Matrices_cont mat_P_SI_M,
                              Double_Matrices_cont mat_age,
                              Double_Matrices_cont mat_specs,
                              gsl_vector * vect_time_line,
                              gsl_vector * vect_DOR,
                              gsl_vector * vect_AaS,
                              gsl_vector * vect_AaQ,
                              double mu,
                              double lambda1,
                              double lambda2,
                              double lambda3,
                              double lambda4)
{
    unsigned int n_ind=mat_P_S_I.nb_rows;
    unsigned int n_years= mat_P_S_I.nb_columns;
    
    for(unsigned int ind=0;ind<n_ind;ind++){
        unsigned int current_smok=(unsigned int)(mat_specs.matrix[ind*mat_specs.nb_columns+2]);
        unsigned int current_YOR=(unsigned int)(vect_DOR->data[ind]/10000);
        unsigned int current_AaS=(unsigned int)(vect_AaS->data[ind]);
        unsigned int current_AaQ=vect_AaQ->data[ind];
        unsigned int first_year=current_YOR-vect_time_line->data[0];
        
        if(DEBUG==1){
            cout << "Ind " << ind
            << " -- Smoking " << current_smok
            << " -- AaS " << current_AaS
            << " -- YOR " << current_YOR
            << " -- first year " << first_year
            << " -- test " << vect_time_line->data[first_year]
            << endl;
        }
        
        for(unsigned int year=first_year;year<n_years;year++){
            // Get current exposure, age and probability of dying
            
            double tmp_exp=mat_exposure.matrix[ind*mat_exposure.nb_columns+year];
            double current_age=mat_age.matrix[ind*mat_age.nb_columns+year];
            double tmp_mort=mat_P_SI_M.matrix[ind*mat_P_SI_M.nb_columns+year];
            
            
            // Compute time since quiting smoking
            
            double current_T_since_Q=0.0;
            if(current_AaQ<120 && current_age>current_AaQ){
                current_T_since_Q=(double)(current_age)-(double)(current_AaQ);
            }
            
            if(tmp_exp>0.0){
                // Compute probability of S-I transition
                
                double arg1=lambda1*log(tmp_exp);
                double arg2=lambda2*(double)(current_age);
                
                double arg3=0.0;
                if(current_smok>1){
                    arg3=lambda3*(double)(current_AaS);
                }
                
                double argTsQ=current_T_since_Q*lambda4;
                
                double arg4=mu;
                
                double arg_tot=arg1+arg2+arg3+arg4+argTsQ;
                
                double pbty=0.0;
                if(isinf(exp(arg_tot))==0){
                    pbty=exp(arg_tot)/(1.0+exp(arg_tot));
                }
                else{
                    pbty=1.0;
                }
                
                mat_P_S_I.matrix[ind*mat_P_S_I.nb_columns+year]=pbty*(1.0-tmp_mort);
                
                if(DEBUG==1){
                    cout << "\tyear " << year
                    << " -- age " << current_age
                    << " -- exp " << tmp_exp
                    << " -- mu " << mu
                    << " -- lambda1 " << lambda1
                    << " -- tmp_exp " << tmp_exp
                    << " -- arg1 " << arg1
                    << " -- lambda2 " << lambda2
                    << " -- arg2 " << arg2
                    << " -- lambda3 " << lambda3
                    << " -- arg3 " << arg3
                    << " -- arg4 " << arg4
                    << " -- lambda4 " << lambda4
                    << " -- argTsQ " << argTsQ
                    << " -- arg_tot " << arg_tot
                    << " -- pbty " << pbty
                    << " -- 1-tmp_mort " << 1.0-tmp_mort
                    << " -- final_pbty " << mat_P_S_I.matrix[ind*mat_P_S_I.nb_columns+year]
                    << endl;
                }
                
            }
            else{
                cout << "!!!!!!!Null Exposure Run Stopped" << endl;
                cout << "Year " << year << " -- date " << vect_time_line->data[year] << endl;
                exit(0);
                mat_P_S_I.matrix[ind*mat_P_S_I.nb_columns+year]=0.0;
            }
        }
    }
}


double wrapped_L_calculation_5_Params_model5(Double_Matrices_cont mat_P_S_I,
                                             Double_Matrices_cont mat_P_S_S,
                                             Double_Matrices_cont mat_exposure,
                                             Double_Matrices_cont mat_P_Ii_Ij_baseline,
                                             Double_Matrices_cont mat_P_Ii_Ij,
                                             Double_Matrices_cont mat_P_SI_M,
                                             Double_Matrices_cont mat_P_SI_R,
                                             Double_Matrices_cont mat_P_SI_SI,
                                             Double_Matrices_cont mat_H_k_per_ind,
                                             Double_Matrices_cont mat_age,
                                             Double_Matrices_cont mat_specs,
                                             gsl_vector * vect_time_line,
                                             gsl_vector * vect_DOR,
                                             gsl_vector * vect_AaS,
                                             gsl_vector * vect_AaQ,
                                             gsl_vector *vect_DOB,
                                             gsl_vector *vect_DOD,
                                             gsl_vector *vect_DODiag,
                                             vector < vector < int > > &summary_SI_SI,
                                             vector < vector < int > > &summary_SI_R,
                                             vector < vector < int > > &summary_SI_M,
                                             double current_mu,
                                             double current_lambda1,
                                             double current_lambda2,
                                             double current_lambda3,
                                             double current_lambda4,
                                             double current_gamma)
{
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
                             current_lambda1,
                             current_lambda2,
                             current_lambda3,
                             current_lambda4);
    
    get_P_S_S(mat_P_S_S,
              mat_P_SI_M,
              mat_P_S_I);
    
    get_P_Ii_Ij_baseline(mat_P_Ii_Ij_baseline,
                         current_gamma);
    
//    cout << "P_Ii_Ij:" << endl;
//    mat_P_Ii_Ij_baseline.Display_matrix_header();
    
    get_P_SI_R_and_H_k(mat_P_SI_R,
                       mat_P_SI_SI,
                       mat_H_k_per_ind,
                       mat_P_S_S,
                       mat_P_S_I,
                       mat_P_SI_M,
                       mat_P_Ii_Ij,
                       mat_P_Ii_Ij_baseline,
//                       summary_SI_SI,
                       vect_time_line,
                       vect_DOB,
                       vect_DOR,
                       vect_DOD,
                       vect_DODiag,
                       mat_specs);
    
    double current_L=get_likelihood(mat_P_SI_R,
                                    mat_P_SI_SI,
                                    mat_P_SI_M,
                                    summary_SI_SI,
                                    summary_SI_R,
                                    summary_SI_M);
    
    return current_L;
}


unsigned int candidate_acceptance(double current_L,
                                  double previous_L)
{
    unsigned int accept=0;
    if(isinf(current_L)==0){
        if(current_L>previous_L){
            if(DEBUG==1){
                cout << "\t\tHigher L, Candidate accepted" << endl;
            }
            accept=1;
        }
        else{
            double acc_ratio=exp(current_L-previous_L);
            double rand_test=myrand();
            if(DEBUG==1){
                cout << "\t\tLower L, likelihood ratio " << acc_ratio << " -- test_sample " << rand_test << endl;
            }
            if(rand_test<=acc_ratio){
                accept=1;
                if(DEBUG==1){
                    cout << "\t\tCandidate accepted" << endl;
                }
            }//end of if (test<ratio)
            else{
                if(DEBUG==1){
                    cout << "\t\tCandidate rejected" << endl;
                }
            }//end of else (test>ratio)
        }//end of else (ratio)
    }
    return accept;
}


unsigned int candidate_acceptance_log(double current_L,
                                      double previous_L,
                                      double current_param,
                                      double previous_param)
{
    unsigned int accept=0;
    if(isinf(current_L)==0){
        //        // Probability g(previous|proposed) where g is the proposal
        //        double proposal_proba_previous=1/(2.50663*current_sigma*previous_param)*exp(-0.5*pow((log(previous_param)-log(current_param)), 2)/pow(current_sigma, 2));
        //        cout << "proposal_proba_previous, i.e. g(previous|proposed): " << proposal_proba_previous << endl;
        //
        //        // Probability g(proposed|previous)
        //        double proposal_proba_current=1/(2.50663*current_sigma*current_param)*exp(-0.5*pow((log(current_param)-log(previous_param)), 2)/pow(current_sigma, 2));
        //        cout << "proposal_proba_current, i.e. g(proposed|previous): " << proposal_proba_current << endl;
        //
        //        double acc_ratio=exp(current_L-previous_L)*(current_param/previous_param);
        //        cout << "acc_ratio: " << acc_ratio << endl;
        
        double acc_ratio=exp(current_L-previous_L)*(current_param/previous_param);
        
        if (DEBUG==1){
            cout << "acc_ratio: " << acc_ratio << endl;
        }
        
        // Formula from original code:
        //        if(current_L>previous_L){
        if(acc_ratio>1){
            if(DEBUG==1){
                cout << "\t\tHigher L, Candidate accepted" << endl;
            }
            accept=1;
        }
        else{
            double rand_test=myrand();
            if(DEBUG==1){
                cout << "\t\tLower L, likelihood ratio " << acc_ratio << " -- test_sample " << rand_test << endl;
            }
            if(rand_test<=acc_ratio){
                accept=1;
                if(DEBUG==1){
                    cout << "\t\tCandidate accepted" << endl;
                }
            }//end of if (test<ratio)
            else{
                if(DEBUG==1){
                    cout << "\t\tCandidate rejected" << endl;
                }
            }//end of else (test>ratio)
        }//end of else (ratio)
    }
    return accept;
}


unsigned int metropolis_step(Double_Matrices_cont mat_P_S_I,
                             Double_Matrices_cont mat_P_S_S,
                             Double_Matrices_cont mat_exposure,
                             Double_Matrices_cont mat_P_Ii_Ij_baseline,
                             Double_Matrices_cont mat_P_Ii_Ij,
                             Double_Matrices_cont mat_P_SI_M,
                             Double_Matrices_cont mat_P_SI_R,
                             Double_Matrices_cont mat_P_SI_SI,
                             Double_Matrices_cont mat_H_k_per_ind,
                             Double_Matrices_cont mat_age,
                             Double_Matrices_cont mat_specs,
                             gsl_vector * vect_time_line,
                             gsl_vector * vect_DOR,
                             gsl_vector * vect_AaS,
                             gsl_vector * vect_AaQ,
                             gsl_vector *vect_DOB,
                             gsl_vector *vect_DOD,
                             gsl_vector *vect_DODiag,
                             vector < vector < int > > &summary_SI_SI,
                             vector < vector < int > > &summary_SI_R,
                             vector < vector < int > > &summary_SI_M,
                             Double_Matrices_cont mat_history,
                             unsigned int iter,
                             unsigned int param_id,
                             double sigma)
{
    // Keeping previous values for other parameters
    
    double current_mu=mat_history.matrix[(iter-1)*mat_history.nb_columns+1];
    double current_lambda1=mat_history.matrix[(iter-1)*mat_history.nb_columns+2];
    double current_lambda2=mat_history.matrix[(iter-1)*mat_history.nb_columns+3];
    double current_lambda3=mat_history.matrix[(iter-1)*mat_history.nb_columns+4];
    double current_lambda4=mat_history.matrix[(iter-1)*mat_history.nb_columns+5];
    double current_gamma=mat_history.matrix[(iter-1)*mat_history.nb_columns+6];
    
    if (param_id>0){
        current_mu=mat_history.matrix[(iter)*mat_history.nb_columns+1];
    }
    
    if (param_id>1){
        current_lambda1=mat_history.matrix[(iter)*mat_history.nb_columns+2];
    }
    
    if (param_id>2){
        current_lambda2=mat_history.matrix[(iter)*mat_history.nb_columns+3];
    }
    
    if (param_id>3){
        current_lambda3=mat_history.matrix[(iter)*mat_history.nb_columns+4];
    }
    
    if (param_id>4){
        current_lambda4=mat_history.matrix[(iter)*mat_history.nb_columns+5];
    }
    
    
    // Step 1.1: Sampling candidate for parameter: random walk
    
    double previous_value=mat_history.matrix[(iter-1)*mat_history.nb_columns+(param_id+1)];
    
    double current_value=0;
    if (param_id==5){
        //        current_value=gennor(log(previous_value),sigma);
        //        current_value=exp(current_value);
        
        double log_gamma=gennor(0,sigma);
        log_gamma+=log(previous_value);
        current_value=exp(log_gamma);
        
        if(DEBUG==1){
            cout << "Previous gamma: " << previous_value << endl;
            cout << "Previous log(gamma): " << log(previous_value) << endl;
            cout << "New gamma: " << current_value << endl;
        }
    } else {
        current_value=gennor(previous_value,sigma);
    }
    
    if (param_id==0){
        current_mu=current_value;
    }
    else if (param_id==1){
        current_lambda1=current_value;
    }
    else if (param_id==2){
        current_lambda2=current_value;
    }
    else if (param_id==3){
        current_lambda3=current_value;
    }
    else if (param_id==4){
        current_lambda4=current_value;
    }
    else if (param_id==5){
        current_gamma=current_value;
    }
    
    
    // Step 1.2: Computing the likelihood for the candidate mu
    
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
        cout << "\tcurrent mu " <<  current_mu
        << " -- current lambda1 " <<  current_lambda1
        << " -- current lambda2 " <<  current_lambda2
        << " -- current lambda3 " <<  current_lambda3
        << " -- current lambda4 " <<  current_lambda4
        << " -- current gamma " <<  current_gamma
        << " -- current L " << current_L
        << endl;
    }
    
    
    // Step 1.3: Accepting/Rejecting the candidate
    
    double previous_L=mat_history.matrix[(iter-1)*mat_history.nb_columns+12];
    if (param_id>0){
        previous_L=mat_history.matrix[(iter)*mat_history.nb_columns+7+param_id-1];
    }
    
    unsigned int accept_candidate=0;
    if (param_id<5){
        accept_candidate=candidate_acceptance(current_L,
                                              previous_L);
    } else {
        if (DEBUG==1){
            cout << "Iteration: " << iter << endl;
            cout << "current_L:" << current_L << endl;
            cout << "previous_L:" << previous_L << endl;
            cout << "current_value:" << current_value << endl;
            cout << "previous_value:" << previous_value << endl;
        }
        accept_candidate=candidate_acceptance_log(current_L,
                                                  previous_L,
                                                  current_value,
                                                  previous_value);
    }
    
    if(DEBUG==1){
        cout << "\taccept_candidate " <<  accept_candidate
        << endl;
    }
    
    
    // Storing values in history matrix
    
    if(accept_candidate==1){
        mat_history.matrix[iter*mat_history.nb_columns+param_id+1]=current_value;
        mat_history.matrix[iter*mat_history.nb_columns+param_id+7]=current_L;
    } else {
        mat_history.matrix[iter*mat_history.nb_columns+param_id+1]=previous_value;
        mat_history.matrix[iter*mat_history.nb_columns+param_id+7]=previous_L;
    }
    mat_history.matrix[iter*mat_history.nb_columns+param_id+13]=sigma;
    
    return accept_candidate;
}


unsigned int metropolis_step_multiv(Double_Matrices_cont mat_P_S_I,
                                    Double_Matrices_cont mat_P_S_S,
                                    Double_Matrices_cont mat_exposure,
                                    Double_Matrices_cont mat_P_Ii_Ij_baseline,
                                    Double_Matrices_cont mat_P_Ii_Ij,
                                    Double_Matrices_cont mat_P_SI_M,
                                    Double_Matrices_cont mat_P_SI_R,
                                    Double_Matrices_cont mat_P_SI_SI,
                                    Double_Matrices_cont mat_H_k_per_ind,
                                    Double_Matrices_cont mat_age,
                                    Double_Matrices_cont mat_specs,
                                    gsl_vector * vect_time_line,
                                    gsl_vector * vect_DOR,
                                    gsl_vector * vect_AaS,
                                    gsl_vector * vect_AaQ,
                                    gsl_vector *vect_DOB,
                                    gsl_vector *vect_DOD,
                                    gsl_vector *vect_DODiag,
                                    vector < vector < int > > &summary_SI_SI,
                                    vector < vector < int > > &summary_SI_R,
                                    vector < vector < int > > &summary_SI_M,
                                    Double_Matrices_cont mat_history,
                                    unsigned int iter,
                                    gsl_vector * theta_previous,
                                    gsl_matrix * Sigma)
{
    // Number of parameters
    
    int N=theta_previous->size;
    
    if(DEBUG==1){
        cout<<"Previous theta:"<<endl;
        display_gsl_vector(theta_previous);
    }
    
    if(DEBUG==1){
        cout<<"Cholesky of current Sigma:"<<endl;
        display_gsl_matrix(Sigma);
    }
    
    
    // Step 1.1: Sampling candidate for parameter: random walk
    
    static gsl_rng *RandomNumberGenerator = gsl_rng_alloc( gsl_rng_taus );
    gsl_vector *theta=gsl_vector_alloc(N);
    gsl_ran_multivariate_gaussian(RandomNumberGenerator, theta_previous, Sigma, theta);
    
    if(DEBUG==1){
        cout<<"Generated sample:"<<endl;
        display_gsl_vector(theta);
    }
    
    
    // Step 1.2: Computing the likelihood for the candidate mu
    
    double current_mu=theta->data[0];
    double current_lambda1=theta->data[1];
//    double current_lambda2=theta->data[2];
//    double current_lambda3=theta->data[3];
//    double current_lambda4=theta->data[4];
//    double current_gamma=exp(theta->data[5]);
    double current_lambda2=0;
    double current_lambda3=0;
    double current_lambda4=0;
    double current_gamma=exp(theta->data[2]);

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
        cout << "\tcurrent mu " <<  current_mu
        << " -- current lambda1 " <<  current_lambda1
        << " -- current lambda2 " <<  current_lambda2
        << " -- current lambda3 " <<  current_lambda3
        << " -- current lambda4 " <<  current_lambda4
        << " -- current gamma " <<  current_gamma
        << " -- current L " << current_L
        << endl;
    }


    // Step 1.3: Accepting/Rejecting the candidate

    double previous_L=mat_history.matrix[(iter-1)*mat_history.nb_columns+7];

    unsigned int accept_candidate=0;
    double previous_gamma=exp(theta_previous->data[2]);
    accept_candidate=candidate_acceptance_log(current_L,
                                              previous_L,
                                              current_gamma,
                                              previous_gamma);

    if(DEBUG==1){
        cout << "Previous L: " << previous_L << endl;
        cout << "\taccept_candidate " <<  accept_candidate
        << endl;
    }


    // Storing values in history matrix

    if(accept_candidate==1){
        mat_history.matrix[iter*mat_history.nb_columns+1]=current_mu;
        mat_history.matrix[iter*mat_history.nb_columns+2]=current_lambda1;
        mat_history.matrix[iter*mat_history.nb_columns+3]=current_lambda2;
        mat_history.matrix[iter*mat_history.nb_columns+4]=current_lambda3;
        mat_history.matrix[iter*mat_history.nb_columns+5]=current_lambda4;
        mat_history.matrix[iter*mat_history.nb_columns+6]=current_gamma;
        
        mat_history.matrix[iter*mat_history.nb_columns+7]=current_L;
    } else {
        mat_history.matrix[iter*mat_history.nb_columns+1]=theta_previous->data[0];
        mat_history.matrix[iter*mat_history.nb_columns+2]=theta_previous->data[1];
//        mat_history.matrix[iter*mat_history.nb_columns+3]=theta_previous->data[2];
//        mat_history.matrix[iter*mat_history.nb_columns+4]=theta_previous->data[3];
//        mat_history.matrix[iter*mat_history.nb_columns+5]=theta_previous->data[4];
//        mat_history.matrix[iter*mat_history.nb_columns+6]=exp(theta_previous->data[5]);
        mat_history.matrix[iter*mat_history.nb_columns+6]=exp(theta_previous->data[2]);
        
        mat_history.matrix[iter*mat_history.nb_columns+7]=previous_L;
    }
    
//    mat_history.matrix[iter*mat_history.nb_columns+param_id+13]=sigma;

    return accept_candidate;
}


double adaptive_step(unsigned int nb_accepted,
                     unsigned int iter,
                     unsigned int iter_adaptive,
                     double sigma,
                     double adaptive_lowerbound,
                     double adaptive_upperbound)
{
    // Computing acceptance rate at this stage (over completed iterations)
    double acc_rate=100.0*(double)(nb_accepted)/(double)(iter);
    
    // Computing increment to add/remove to log(sigma), making sure its amount is decreasing over iterations
    double inc=pow((double)(iter_adaptive), (-0.5));
    
    if (inc>0.01){
        inc=0.01;
    }
    
    if(DEBUG==1){
        cout << "***" << endl;
        cout << "Adaptive iteration: " << iter_adaptive << endl;
        cout << "Number of iterations in batch (always the same): " << iter << endl;
        cout << "Number of accepted values: " << nb_accepted << endl;
        cout << "Acceptance rate: " << acc_rate << endl;
        cout << "Previous sigma: " << sigma << endl;
        cout << "Increment of sigma: " << inc << endl;
    }
    
    if(acc_rate>adaptive_upperbound){
        if(DEBUG==1){
            cout << "Increasing sigma" << endl;
        }
        double tmp=log(sigma)+inc;
        sigma=exp(tmp);
    }
    if(acc_rate<adaptive_lowerbound) {
        if(DEBUG==1){
            cout << "Decreasing sigma" << endl;
        }
        double tmp=log(sigma)-inc;
        sigma=exp(tmp);
    }
    
    if(DEBUG==1){
        cout << "Updated sigma: " << sigma << endl;
        cout << "***" << endl;
    }
    
    return sigma;
}


double adaptive_step_multiv(unsigned int iter,
                            unsigned int n_iter_adaptive,
                            Double_Matrices_cont mat_history,
                            gsl_matrix* hat_sigma)
{
    cout << "Updating covariance of the proposal -- iter " << iter << endl;
    int N=3;
    
    
    // Storing values of theta over "n_iter_adaptive" previous iterations in X
    
    gsl_matrix* X=gsl_matrix_alloc(n_iter_adaptive,N);
    for (unsigned int i=0;i<n_iter_adaptive;i++){
        X->data[i * X->tda + 0]=mat_history.matrix[(iter-i)*mat_history.nb_columns+1];
        X->data[i * X->tda + 1]=mat_history.matrix[(iter-i)*mat_history.nb_columns+2];
        X->data[i * X->tda + 2]=log(mat_history.matrix[(iter-i)*mat_history.nb_columns+6]);
    }
    
    if(DEBUG==1){
        cout << "Matrix X:" << endl;
        display_gsl_matrix(X);
    }

    // Computing the covariance
    
    gsl_ran_multivariate_gaussian_vcov(X, hat_sigma);

    if(DEBUG==0){
        cout << "Empirical covariance:" << endl;
        display_gsl_matrix(hat_sigma);
    }

//    // Scaling the matrix X to compute the empirical correlation over previous iterations
//
//    // Estimating the mean
//    gsl_vector* hat_mean=gsl_vector_alloc(N);
//    gsl_ran_multivariate_gaussian_mean(X, hat_mean);
//
//    if(DEBUG==1){
//        cout<<"Empirical mean:"<<endl;
//        display_gsl_vector(hat_mean);
//    }
//
//    // X = X - empirical mean
//    for (unsigned int i=0;i<n_iter_adaptive;i++){
//        X->data[i * X->tda + 0]-=hat_mean->data[0];
//        X->data[i * X->tda + 1]-=hat_mean->data[1];
//        X->data[i * X->tda + 2]-=hat_mean->data[2];
//    }
//
//    if(DEBUG==1){
//        cout << "Scaled matrix X:" << endl;
//        display_gsl_matrix(X);
//    }
//
//
//    // Computing the covariance
//
//    gsl_ran_multivariate_gaussian_vcov(X, hat_sigma);
//
//    if(DEBUG==0){
//        cout << "Empirical covariance2:" << endl;
//        display_gsl_matrix(hat_sigma);
//    }
//
//    // Computing increment to add/remove to log(sigma), making sure its amount is decreasing over iterations
//    double inc=pow((double)(iter_adaptive), (-0.5));
//
//    if (inc>0.01){
//        inc=0.01;
//    }
//
//    if(DEBUG==1){
//        cout << "***" << endl;
//        cout << "Adaptive iteration: " << iter_adaptive << endl;
//        cout << "Number of iterations in batch (always the same): " << iter << endl;
//        cout << "Number of accepted values: " << nb_accepted << endl;
//        cout << "Acceptance rate: " << acc_rate << endl;
//        cout << "Previous sigma: " << sigma << endl;
//        cout << "Increment of sigma: " << inc << endl;
//    }
//
//    if(acc_rate>adaptive_upperbound){
//        if(DEBUG==1){
//            cout << "Increasing sigma" << endl;
//        }
//        double tmp=log(sigma)+inc;
//        sigma=exp(tmp);
//    }
//    if(acc_rate<adaptive_lowerbound) {
//        if(DEBUG==1){
//            cout << "Decreasing sigma" << endl;
//        }
//        double tmp=log(sigma)-inc;
//        sigma=exp(tmp);
//    }
//
//    if(DEBUG==1){
//        cout << "Updated sigma: " << sigma << endl;
//        cout << "***" << endl;
//    }
    
//    return hat_sigma;
}


int factorial(int n)
{
    if(n > 1)
        return n * factorial(n - 1);
    else
        return 1;
}

