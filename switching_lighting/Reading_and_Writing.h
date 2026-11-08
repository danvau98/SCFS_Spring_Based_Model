#pragma once
#include <sstream>
//Generates a matrix into a datafile
void Writing_to_File_Mat(double** Mat, int row, int col, std::string Filename);
//Generates an array into a datafile
void Writing_to_File_Array(double* Array, int size, std::string Filename);
//Reads a matrix from a datafile in a specific format
double** Reading_from_File_Mat(double** Mat, int row, int col, std::string Filename);
//Reads an array from a datafile in a specific format
double* Reading_from_File_Array(double* Array, int size, std::string Filename);