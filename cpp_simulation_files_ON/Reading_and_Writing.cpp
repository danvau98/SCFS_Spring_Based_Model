#include <iostream>
#include <cmath>
#include <sstream>
#include <fstream>
#include <cstring>
#include <cctype>
#include "Reading_and_Writing.h"
//Generates a matrix into a datafile
void Writing_to_File_Mat(double** Mat, int row, int col, std::string Filename)
{
	//Function to write data to a textfile
	std::ofstream myfile;
	myfile.open(Filename);
    for (int i = 0; i < row; i++) 
    {
        for (int j = 0; j < col; j++) 
        {
            myfile << Mat[i][j] << " ";
        }
        myfile << "\n";
    }
}
//Generates an array into a datafile
void Writing_to_File_Array(double* Array, int size, std::string Filename)
{
    std::ofstream myfile;
    myfile.open(Filename);
    for (int i = 0; i < size; i++)
    {
        myfile << Array[i] << "\n";
    }
}
//Reads a matrix from a datafile in a specific format
double** Reading_from_File_Mat(double** Mat, int row, int col, std::string Filename)
{
    //Function to read data from a file and produce a matrix. 
    std::ifstream infile(Filename);
    for (int i = 0; i < row; i++) 
    {
        for (int j = 0; j < col; j++)
        {
            infile >> Mat[i][j];
        }
    }
    infile.close();
    return Mat;
}
//Reads an array from a datafile in a specific format
double* Reading_from_File_Array(double* Array, int size, std::string Filename)
{
    //Function to read data from a file and produce a matrix. 
    std::ifstream infile(Filename); 
    for (int i = 0; i < size; i++)
    {
        infile >> Array[i];
    }
    infile.close();
    return Array;
}
