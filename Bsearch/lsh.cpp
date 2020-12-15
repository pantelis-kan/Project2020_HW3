#include <iostream>
#include <string.h>
#include <string>
#include <sstream>
#include <fstream>
#include <iterator>
#include <limits>
#include <cmath>
#include <vector>
#include <chrono>
#include <random>

#include "Point_Table.hpp"
#include "Point.hpp"
#include "utilities.hpp"
#include "Hash_Table.hpp"
#include "NN_Functions.hpp"

using namespace std;

const long long int m = 4294967291; // 2^32  - 5

int N = 1;
double R = 10000.0;
int k = 4, L = 5;
string inputfile_original = "1train-images-idx3-ubyte";
string inputfile_reduced = "out.bin";
string queryfile_original = "t10k-images-idx3-ubyte";
string queryfile_reduced;
string outputfile = "lsh_results.txt";

double w = 30000.0;

std::default_random_engine rand_generator(time(NULL));

int main(int argc, char* argv[]){

	/******************************************
	 * Parse execution arguments on runtime.
	*******************************************/

	for (int i = 1; i < argc; i+=2){
		string arg = argv[i];
		if (arg == "-d"){
			inputfile_original = argv[i+1];
		}
		else if (arg == "-i"){
			inputfile_reduced = argv[i+1];
		}
		else if (arg == "-q"){
			queryfile_original = argv[i+1];
		}
		else if (arg == "-s"){
			queryfile_reduced = argv[i+1];
		}
		else if (arg == "-k"){
			k = atoi(argv[i+1]);
			if (k <= 0){
				cout << "Wrong -k parameter! Please try again!" << endl;
				exit(1);
			}
		}
		else if (arg == "-L"){
			L = atoi(argv[i+1]);
			if (L <= 0){
				cout << "Wrong -L parameter! Please try again!" << endl;
				exit(1);
			}				
		}
		else if (arg == "-o"){
			outputfile = argv[i+1];
		}
		else{
			cout << "Wrong parameteres! Please try again!" << endl;
			exit(1);
		}

	}

	int M = pow(2,32/k);
	
	cout << "Programme will run with inputfile_original: " << inputfile_original << " for input_original data" << endl; 




	/******************************************
	 * Reading ORIGINAL input dataset.
	*******************************************/


	int input_count_original = NumberOfPoints(inputfile_original);   // number of input_original points
	int TableSize_original = input_count_original/8;

	int dimension_original = DimensionofPoint(inputfile_original);	

	cout << "Number of points is : " << input_count_original <<endl;
	cout << "TableSize_original = " << TableSize_original <<endl;
	cout << endl << "Dimension ORIGINAL= "<< dimension_original <<endl;

	Point_Array input_original(input_count_original, dimension_original);
	
	if(input_original.FillPoints(inputfile_original) == 0) cout << "Filling input_original points successful"<<endl;
	else exit(-1);


	//int dimension_original = input_original.get_dimension();
	cout << endl << "Dimension ORIGINAL= "<< dimension_original <<endl;

	/******************************************
	 * Reading REDUCED input dataset.
	*******************************************/

	int input_count_reduced = NumberOfPoints(inputfile_reduced);   // number of input points
	int TableSize_reduced = input_count_reduced/8;

	int dimension_reduced = DimensionofPoint(inputfile_reduced);	

	cout << "Number of points in REDUCED is : " << input_count_reduced <<endl;
	cout << "TableSize REDUCED = " << TableSize_reduced <<endl;

	Point_Array input_reduced(input_count_reduced, dimension_reduced);
	
	if(input_reduced.FillPoints(inputfile_reduced) == 0) cout << "Filling input_reduced points successful"<<endl;
	else exit(-1);
	

	//int dimension_reduced = input_reduced.get_dimension();
	cout << endl << "Dimension REDUCED= "<< dimension_reduced <<endl;

	exit(-1);

	/******************************************
	 * Building si parameters needed for amplification - ORIGINAL dataset
	*******************************************/

	// every h (h1, h2, ..., hk) has its own parameters for the amplification
	// definition of si parameter with i = 0,1,...,d-1
	double**  s_params_original = new double*[L*k];
	
	for (int i = 0; i < L*k; i++){
		s_params_original[i] = new double[dimension_original];
	} 


	//Limiting rand function to take values from 0.0 to w
	std::uniform_real_distribution<double> distribution(0.0,w);

	// create s parameters for each h(i). Since there are L hash tables, there are L*k rows in the table
	for(int i = 0; i < L*k; i++){
		
		for(int j = 0; j < dimension_original; j++) {
			double rand = distribution(rand_generator);
			s_params_original[i][j] = rand;
		}
	}


	/******************************************
	 * Building si parameters needed for amplification - REDUCED dataset
	*******************************************/

	// every h (h1, h2, ..., hk) has its own parameters for the amplification
	// definition of si parameter with i = 0,1,...,d-1
	double**  s_params_reduced = new double*[L*k];
	
	for (int i = 0; i < L*k; i++){
		s_params_reduced[i] = new double[dimension_reduced];
	} 


	//Limiting rand function to take values from 0.0 to w
	//std::uniform_real_distribution<double> distribution(0.0,w);

	// create s parameters for each h(i). Since there are L hash tables, there are L*k rows in the table
	for(int i = 0; i < L*k; i++){
		
		for(int j = 0; j < dimension_reduced; j++) {
			double rand = distribution(rand_generator);
			s_params_reduced[i][j] = rand;
		}
	}



	/******************************************
	 * Building LSHashtable and storing ORIGINAL input data
	 * Finding nearest neighbors with approximate and Range methods.
	*******************************************/	

	cout << endl << "Stage 1: Preprocessing stage for ORIGINAL dataset... "<<endl;

	//create L hash tables
	Hash_Table**  H_Tables_original = new Hash_Table*[L];
	
	for (int i = 0; i < L; i++){	
		H_Tables_original[i] = new Hash_Table(TableSize_original);
	}

	auto t1 = std::chrono::high_resolution_clock::now();		
	Preprocessing(H_Tables_original, input_original, input_count_original, TableSize_original, s_params_original, L, k, M, m, w);
	auto t2 = std::chrono::high_resolution_clock::now();

	auto duration = std::chrono::duration_cast<std::chrono::seconds>( t2 - t1 ).count();

	cout << "Stage 1 completed in " << duration << " seconds" << endl;


	/******************************************
	 * Building LSHashtable and storing REDUCED input data
	 * Finding nearest neighbors with approximate method.
	*******************************************/	

	cout << endl << "Stage 1: Preprocessing stage for REDUCED dataset... "<<endl;

	//create L hash tables
	Hash_Table**  H_Tables_reduced = new Hash_Table*[L];
	
	for (int i = 0; i < L; i++){
		H_Tables_reduced[i] = new Hash_Table(TableSize_reduced);
	}

	t1 = std::chrono::high_resolution_clock::now();		
	Preprocessing(H_Tables_reduced, input_reduced, input_count_reduced, TableSize_reduced, s_params_reduced, L, k, M, m, w);
	t2 = std::chrono::high_resolution_clock::now();

	duration = std::chrono::duration_cast<std::chrono::seconds>( t2 - t1 ).count();

	cout << "Stage 1 completed in " << duration << " seconds" << endl;


	/******************************************
	 * Reading query set and output filenames from user
	 * Creating necessary structures for the queries
	 * Running until user terminates the programme
	*******************************************/	
	int option;
	do{
		cout << "Press 1 you want to run programme with default search and output data file. Press 2 if you want to choose other inputfile_original. Any other option will exit programme" << endl;
		cin >> option;
		if (option == 1){
			cout << "You have chosen the default inputfile_original!" << endl;
		}
		else if(option == 2){
			cout << "Insert queries_original filename: "; 
			cin >> queryfile_original;
			cout << "Insert output filename: "; 
			cin >> outputfile;
		}
		else{
			cout << "Programme is exiting!" << endl;
			exit(1);
		}

		cout << "Search will be done using queryfile_original with name: " << queryfile_original << endl;
		cout << "Search will be done using queryfile_reduced with name: " << queryfile_reduced << endl;
		cout << "Output will be exported to file with name: " << outputfile << endl;


		/******************************************
		 * Reading queries for ORIGINAL and REDUCED.
		 * Checking that the files have the same number of queries
		*******************************************/

		int queries_count_original = NumberOfPoints(queryfile_original); // number of query points

		cout << "Number of queries_original is : " << queries_count_original <<endl;
		Point_Array queries_original(queries_count_original, dimension_original);
		if(queries_original.FillPoints(queryfile_original) == 0) cout << "Filling query points successful"<<endl;
		else exit(-1);		

		int queries_count_reduced = NumberOfPoints(queryfile_reduced); // number of query points
		cout << "Number of queries_reduced is : " << queries_count_reduced <<endl;
		Point_Array queries_reduced(queries_count_reduced, dimension_original);
		if(queries_reduced.FillPoints(queryfile_reduced) == 0) cout << "Filling query points successful"<<endl;
		else exit(-1);

		if (queries_count_original != queries_count_reduced){
			cout << "Error: queries_count_original is not the same as queries_count_reduced" << endl;
			exit(-1);
		}
		

		/******************************************
		 * Finding nearest neighbors with approximate method.
		*******************************************/

		cout << "Stage 2: Finding Nearest Neighbors" <<endl;
		
		Results results_original[queries_count_original];
		LSH_Nearest_Neighbors(results_original, H_Tables_original, input_original, queries_original, queries_count_original, TableSize_original, s_params_original, L, k, M, m, w, N);

		Results results_reduced[queries_count_reduced];
		LSH_Nearest_Neighbors(results_reduced, H_Tables_reduced, input_reduced, queries_reduced, queries_count_reduced, TableSize_reduced, s_params_reduced, L, k, M, m, w, N);

		cout << "Stage 2 completed!" << endl;

		/******************************************
		 * Exporting results to output file with the required format
		*******************************************/

		cout << "Stage 3: Exporting results to file" << endl;
		string exact_NN_fp = "exact_results.txt";
		Exact_NN_readonly(results_original, queries_count_original, N, exact_NN_fp);
//here probably needs to be added the reading of exact results for the NeuralNet


		ofstream final_results;
		final_results.open(outputfile, ios::out | ios::trunc);

		for (int i = 0; i < queries_count_original; i++){
			final_results << "Query: " << results_original[i].get_query_id() << endl;
			
			vector <int> temp_N_nearest_id_original = results_original[i].get_N_nearest_id();
			vector <double> temp_N_nearest_distance_original = results_original[i].get_N_nearest_distance();
			vector <double> temp_exact_N_nearest_original = results_original[i].get_exact_N_nearest();


			int counter = 1;
			auto it_distance_original = temp_N_nearest_distance_original.cbegin();
			auto it_exact_distance_original = temp_exact_N_nearest_original.cbegin();

			for(auto it_id = temp_N_nearest_id_original.cbegin(); it_id != temp_N_nearest_id_original.cend(); ++it_id){
				final_results << "Nearest neighbor-" << counter << ": " << *it_id << endl;
				final_results << "distanceLSH: " << *it_distance_original << endl;
				final_results << "distanceTrue: " << *it_exact_distance_original << endl;

				it_distance_original++;
				it_exact_distance_original++;
				counter++;
			}

			final_results << "tLSH: " << results_original[i].get_t_NN() << endl; 
			final_results << "tTrue: " << results_original[i].get_tTrue() << endl;


		}

		final_results.close();

		cout << "Stage 3 - Completed!" << endl;



		cout << "If you want to repeat search press 1! Any other answer will terminate the programme!" << endl; 
		cin >> option; 

	}while (option == 1);

	for(int i = 0; i < L;  i++) delete H_Tables_original[i]; 
	delete[] H_Tables_original;
	
	for (int i = 0; i <  L * k; i++){
		delete[] s_params_original[i];
	}
	delete[] s_params_original;


	return 0;
}


