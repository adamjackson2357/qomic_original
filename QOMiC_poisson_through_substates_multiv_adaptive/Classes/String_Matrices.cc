#include "String_Matrices.h"
#define DEBUG 0

using namespace std;


String_Matrices::String_Matrices()
{
  nb_rows=0;
  nb_columns=0;
}

void String_Matrices::Alloc_string_matrix(int Rows,
					  int Columns)
{
  
  nb_rows=Rows;
  nb_columns=Columns;

  matrix=new string*[Rows];
  for(int row=0;row<Rows;row++){
    matrix[row]=new string[Columns];
  }
}


void String_Matrices::Free_string_matrix()
{
  for(int row=0;row<nb_rows;row++){
    delete[] matrix[row];
  }
  delete[] matrix;
}


void String_Matrices::Read_from_file(char *filename)
{

  ifstream INFILE;
  INFILE.open(filename, ios::in);
  //Checking the file path

  if(INFILE.fail()){
    cout << "Invalid Path and/or permission rights for " << filename << " -- run stopped." << endl;
    exit(1);
  }

  //  Reading the first element: nb_rows
  INFILE >> nb_rows;
  //  Reading the second element: nb_columns
  INFILE >> nb_columns;
  cout << "Reading file : " << filename
       << " -- #rows: " << nb_rows
       << " -- #cols: " << nb_columns << endl;

  //Allocating memory for the Matrix
  Alloc_string_matrix(nb_rows,
		      nb_columns);

  //Storing the map file in the matrix.
  for(int current_line=0;current_line<nb_rows;current_line++){
    for(int current_column=0;current_column<nb_columns;current_column++){
      if(!INFILE.eof()){
	if(DEBUG){
	  cout << " ligne " << current_line
	       << " -- colonne " << current_column
	       << " --Test " << INFILE.eof() << endl;
	}
	  INFILE>>matrix[current_line][current_column];
      }
      else{
	cout << "Missing element at line " << current_line
	     << " and column " << current_column
	     << " in file " << filename << "." << endl
	     << "!!!!Run Stopped!!!!" << endl;
	exit(1);
      }
    }
  } 
}


void String_Matrices::Display_matrix()
{
  for(int current_line=0;current_line<nb_rows;current_line++){
    for(int current_column=0;current_column<nb_columns;current_column++){
      cout << matrix[current_line][current_column] << " ";
    }
    cout << endl;
  }
  cout << endl;
}


void String_Matrices::Display_matrix_header()
{

  cout << nb_rows << endl;
  cout << nb_columns << endl;

  for(int current_line=0;current_line<nb_rows;current_line++){
    for(int current_column=0;current_column<nb_columns;current_column++){
      cout << matrix[current_line][current_column] << " ";
    }
    cout << endl;
  }
  cout << endl;

}


void String_Matrices::Write_to_file(char *filename)
{

  ofstream OUTFILE;
  OUTFILE.open(filename, ios::out);
  //Checking the file path

  if(OUTFILE.fail()){
    cout << "Invalid Path and/or permission rights for " << filename << " -- run stopped." << endl;
    exit(1);
  }

  //Writing the first element: nb_rows
  OUTFILE << nb_rows << endl;
  //Writing the second element: nb_columns
  OUTFILE << nb_columns << endl;

  //Writing the core of the matrix.
  for(int current_line=0;current_line<nb_rows;current_line++){
    for(int current_column=0;current_column<nb_columns;current_column++){
      OUTFILE << matrix[current_line][current_column] << " ";
    }
    OUTFILE << endl;
  } 
  OUTFILE.close();
}
