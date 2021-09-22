#ifndef DYN_NAME_H
#define DYN_NAME_H


#include <stdlib.h>
#include <stddef.h>
#include <stdarg.h>
#include <string.h>
#include <stdio.h>
#include <cstring>
#include <iostream>
#include <sstream> 
#include <cstdarg>
#include "struc.h"
using namespace std; 

string Get_dyn_name(string File_name,
		    int Number,
		    string Extension);

void Write_dyn_double(data_double * object,
		      char * File_name,
		      string Name_number1 ,
		      int Number1,
		      string Name_number2 ,
		      int Number2,
		      string Extension);


string Get_scen_name(string File_name,
		     char *scenario,
		     string Extension);

void Write_dyn_double_bis(data_double * object,
			  char * File_name,
			  string Name_number1 ,
			  int Number1,
			  string Name_number2 ,
			  int Number2,
			  string Name_number3 ,
			  int Number3,
			  string Extension);

string Get_stddzed_name(string File_name,
			int Number,
			string Name_number,
			string Extension);
string Get_simple_name(string File_name,
		       string Extension);
void Write_matrix_double(data_double *M, char *file_name);
string Write_dyn_nname_2params(char * File_name,
			     string Name_number1 ,
			     int Number1,
			     string Name_number2 ,
			     int Number2,
			       string Extension);
#endif /* DYN_NAME_H */
