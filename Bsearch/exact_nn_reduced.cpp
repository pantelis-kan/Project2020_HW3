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

    string filename = "out_dataset.bin";
	string filename2 = "out_queryset.bin";

	int input_count = NumberOfPoints(filename);   // number of input points
	int dimension = DimensionofPoint(filename);	

	int queries_count = NumberOfPoints(filename2); // number of query points

	cout << "Number of points is : " << input_count <<endl;
	cout << "Number of queries is : " << queries_count <<endl;

	Point_Array input(input_count, dimension);
	Point_Array queries(queries_count, dimension);
	
	if(input.FillPoints_reduced(filename) == 0) cout << "Filling input points successful"<<endl;
	else exit(-1);
	
	if(queries.FillPoints_reduced(filename2) == 0) cout << "Filling query points successful"<<endl;
	else exit(-1);

	cout << endl << "Dimension = "<< dimension <<endl;

	/******************************************
	 * Exporting results to output file with the required format
	*****************************************/

    ofstream outfile;
    outfile.open("exact_results_reduced.txt", ios::out | ios::trunc );

    Exact_NN(input, queries, input_count, queries_count, outfile);

    // int time_passed;
    // Exact_NN(input, queries, input_count, queries_count, outfile, &time_passed);

    // outfile << time_passed << endl;
    outfile.close();
}