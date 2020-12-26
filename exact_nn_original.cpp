#include <iostream>
#include <string.h>
#include <string>
#include <sstream>
#include <fstream>
#include <iterator>


#include "Point_Table.hpp"
#include "Point.hpp"
#include "utilities.hpp"
#include "NN_Functions.hpp"


int main(){

	/******************************************
	 * Reading input and query datasets.
	*****************************************/

    string filename = "train-images-idx3-ubyte";
	string filename2 = "t10k-images-idx3-ubyte";
	
	int option;
	cout << "Press 1 you want to run programme with default search and output data file. Press 2 if you want to choose other dataset and queryset. Any other option will exit programme" << endl;
	cin >> option;
	
	if (option == 1){
		cout << "You have chosen the default inputfiles!" << endl;
	}
	else if(option == 2){
		cout << "Insert dataset filename: "; 
		cin >> filename;
		cout << "Insert queryset filename: "; 
		cin >> filename2;
	}
	else{
		cout << "Programme is exiting!" << endl;
		exit(1);
	}
	cout << "Search will be done using dataset with name: " << filename << endl;
	cout << "Search will be done using queryset with name: " << filename2 << endl;


	int input_count = NumberOfPoints(filename);   // number of input points
	int dimension = DimensionofPoint(filename);	

	int queries_count = NumberOfPoints(filename2); // number of query points

	cout << "Number of points is : " << input_count <<endl;
	cout << "Number of queries is : " << queries_count <<endl;

	Point_Array input(input_count, dimension);
	Point_Array queries(queries_count, dimension);
	
	if(input.FillPoints(filename) == 0) cout << "Filling input points successful"<<endl;
	else exit(-1);
	
	if(queries.FillPoints(filename2) == 0) cout << "Filling query points successful"<<endl;
	else exit(-1);


	cout << endl << "Dimension = "<< dimension <<endl;
	
	/******************************************
	 * Exporting results to output file with the required format
	*****************************************/

    ofstream outfile;
    outfile.open("exact_results_original.txt", ios::out | ios::trunc );

    Exact_NN(input, queries, input_count, queries_count, outfile);

    outfile.close();
}