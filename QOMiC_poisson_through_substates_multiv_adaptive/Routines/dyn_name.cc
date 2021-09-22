#include "dyn_name.h"
#include "struc.h"
#include <string>
#include <iostream>
#include <sstream> 



//using namespace std; 

string Get_dyn_name(string File_name,
		    int Number,
		    string Extension)
{
  string Result;
  ostringstream Int_to_String;

  Int_to_String << Number;

  Result= File_name+Int_to_String.str()+Extension;

  return Result;
}


string Get_stddzed_name(string File_name,
			int Number,
			string Name_number,
			string Extension)
{
  string Result;
  ostringstream Int_to_String;
  string separator="_";
  Int_to_String << Number;

  Result= File_name+separator+Int_to_String.str()+separator+Name_number+Extension;

  return Result;
}


string Get_simple_name(string File_name,
		       string Extension)
{
  string Result;
  ostringstream Int_to_String;

  Result= File_name+Extension;

  return Result;
}
string Get_scen_name(string File_name,
		     char *scenario,
		     string Extension)
{
  string Result;
  ostringstream Int_to_String;
  string separator="_";
  string temp_scen=scenario;

  Result= File_name+separator+temp_scen+Extension;
  
  return Result;
}

string Write_dyn_nname_2params(char * File_name,
			     string Name_number1 ,
			     int Number1,
			     string Name_number2 ,
			     int Number2,
			     string Extension)
{
  string separator="_";
  string temp_res1=File_name;
  ostringstream ostr1, ostr2;
  
  ostr1 << Number1;
  ostr2 << Number2;
  
  
  string Result=temp_res1+separator+Name_number1+separator+ostr1.str()+separator+Name_number2+separator+ostr2.str()+Extension;
  return Result;  
}

void Write_dyn_double(data_double * object,
		      char * File_name,
		      string Name_number1 ,
		      int Number1,
		      string Name_number2 ,
		      int Number2,
		      string Extension)
{
  char *temp;
  
  string temp_res1=File_name;
  ostringstream ostr1, ostr2;
  
  ostr1 << Number1;
  ostr2 << Number2;
  
  
  string temp_res2=temp_res1+Name_number1+ostr1.str()+Name_number2+ostr2.str()+Extension;
  (temp)=(char*)(temp_res2.c_str());
  
  
  Write_matrix_double(object,temp);
  
}

void Write_dyn_double_bis(data_double * object,
			  char * File_name,
			  string Name_number1 ,
			  int Number1,
			  string Name_number2 ,
			  int Number2,
			  string Name_number3 ,
			  int Number3,
			  string Extension)
{
  char *temp;
  
  string temp_res1=File_name;
  ostringstream ostr1, ostr2, ostr3;
  
  ostr1 << Number1;
  ostr2 << Number2;
  ostr3 << Number3;
  
  
  string temp_res2=temp_res1+Name_number1+ostr1.str()+Name_number2+ostr2.str()+Name_number3+ostr3.str()+Extension;
  (temp)=(char*)(temp_res2.c_str());
  
  
  Write_matrix_double(object,temp);
  
}



void Write_matrix_double(data_double *M, char *file_name)
{
  int i;
  int j;
  FILE *flot;
  
  /* flot est un pointeur vers un fichier de type FILE*/
  
  flot=fopen(file_name,"w");
  fprintf(flot, "%d \n",(*M).nb_rows);
  fprintf(flot, "%d \n",(*M).nb_columns);
  
  for(i=0;i<(*M).nb_rows;i++) {
    for(j=0;j<(*M).nb_columns;j++) {
      fprintf(flot,"%f ",(double)((*M).matrix[i][j]));
    }
    fprintf(flot,"\n");
  }
  fclose(flot);
}
