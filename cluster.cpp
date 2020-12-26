#include <iostream>
#include <string.h>
#include <string>
#include <algorithm>
#include <sstream>
#include <fstream>
#include <iterator>
#include <limits>
#include <cmath>
#include <vector>
#include <chrono>
#include <random>
#include <vector>

#include "Point_Table.hpp"
#include "Point.hpp"
#include "NN_Functions.hpp"
#include "Hash_Table.hpp"
#include "utilities.hpp"
#include "Cluster.hpp"
#include "Cluster_Functions.hpp"
//#include "Hypercube.hpp"

using namespace std;

int k = 10; // k = number of clusters

string filename = "train-images-idx3-ubyte";
string filename2 = "t10k-images-idx3-ubyte";
string configuration_filename = "cluster.conf";
string method = "Classic";
string outputfile = "cluster_output.txt";
string inputfile_reduced = "out_queryset.bin";
string class_file = "class_file.txt";

bool complete = false;

int main(int argc, char* argv[]){


	for (int i = 1; i < argc; i+=2){
		string arg = argv[i];
		if (arg == "-d"){
			filename = argv[i+1];
		}
		else if (arg == "-c"){
			configuration_filename = argv[i+1];
		}
		else if (arg == "-i"){
			inputfile_reduced = argv[i+1];
		}
		else if (arg == "-o"){
			outputfile = argv[i+1];
		}
		else if (arg == "-n"){
			class_file = argv[i+1];
		}
	}

    int N = NumberOfPoints(filename2);   // number of input points
	int TableSize = N/8;
	cout << "Number of points is : " << N <<endl;
	cout << "TableSize = " <<TableSize <<endl;

	int dimension_original = DimensionofPoint(filename2);
	Point_Array input(N,dimension_original);

	int dimension_new = DimensionofPoint(inputfile_reduced);
	Point_Array input_new(N,dimension_new);
	
	//cout << dimension_original << "  " << dimension_new << endl;
	if(input.FillPoints(filename2) == 0) cout << "Filling input points successful"<<endl;
	else exit(-1);
	

	//int dimension = input.get_dimension();
	//cout << endl << "Dimension = "<< dimension <<endl;

	if(input_new.FillPoints_reduced(inputfile_reduced) == 0) cout << "Filling input points successful"<<endl;
	else{ cout << "Error at FillPoints()" << endl; exit(-1);}

	Configuration_File(configuration_filename,&k);
	
	

	// ----------- ORIGINAL SPACE -----------

    Cluster* clusters = new Cluster[k];  

	for(int i = 0; i < k; i++)
		clusters[i].Inititalize_Centroid_Points(dimension_original);

	cout << "Starting k-means++ Initialization" << endl;
	Initialize_Centroids(input,clusters,k); // k-means++
	cout << "Initialization complete"<<endl;

	auto t1 = std::chrono::high_resolution_clock::now();


	//Lloyds
	cout << endl << "Performing Lloyd's assignment" <<endl <<endl;
	Loyds_Clusters(input,clusters,k);


	auto t2 = std::chrono::high_resolution_clock::now();		
	auto duration = std::chrono::duration_cast<std::chrono::seconds>( t2 - t1 ).count();
	cout << "Time taken : " << duration << " seconds" <<endl;

	double time = (double)duration;

	double* s = new double[input.get_ArraySize()];
	cout << endl << "Starting silhouette " <<endl;

	double s_total;

	Silhouette(input,clusters,k,s,&s_total);

	cout << endl << "Printing output to file " << outputfile << endl;
	string original_header = "ORIGINAL SPACE";
	Output_Results(input,clusters,k,s,outputfile,time,s_total,original_header);

	delete[] clusters;

	// ---------------------- NEW SPACE ----------------------------

	string outputfile_new_dimension = "cluster_output_new_dimension.txt";
	clusters = new Cluster[k];  

	for(int i = 0; i < k; i++)
		clusters[i].Inititalize_Centroid_Points(dimension_new);

	cout << "Starting k-means++ Initialization" << endl;
	Initialize_Centroids(input_new,clusters,k); // k-means++
	cout << "Initialization complete"<<endl;

	auto t1_new = std::chrono::high_resolution_clock::now();


	//Lloyds
	cout << endl << "Performing Lloyd's assignment" <<endl <<endl;
	Loyds_Clusters(input_new,clusters,k);


	auto t2_new = std::chrono::high_resolution_clock::now();		
	auto duration_new = std::chrono::duration_cast<std::chrono::seconds>( t2_new - t1_new ).count();
	cout << "Time taken : " << duration_new << " seconds" <<endl;


	double* s1 = new double[input_new.get_ArraySize()];
	cout << endl << "Starting silhouette " <<endl;

	Silhouette(input_new,clusters,k,s1,&s_total);

	double time_new = (double)duration_new;

	cout << endl << "Printing output to file " << outputfile << endl;
	string new_header = "NEW SPACE";
	Output_Results(input_new,clusters,k,s1,outputfile_new_dimension,time_new,s_total,new_header);
	

	string class_output = "class_file_output.txt";
	Class_File_Output(input,class_file,class_output);

	//cout << Distance(input.Retrieve(3),input.Retrieve(10),1) << endl;

	delete[] clusters;
	delete[] s;
	delete[] s1;
}


