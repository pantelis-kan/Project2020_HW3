#ifndef NN_FUNCTIONS_HPP
#define NN_FUNCTIONS_HPP


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


void Preprocessing(Hash_Table** H_Tables,Point_Array& input, int input_count, int TableSize, double** s_params, int L,
					int k, int M, long long int m, double w);

void LSH_Nearest_Neighbors(Results* results, Hash_Table** H_Tables, Point_Array& input, Point_Array& queries,
				int queries_count, int TableSize, double** s_params,int L, int k, int M, long long int m, double w, int N);

void LSH_Range_Search(Results* results, Hash_Table** H_Tables, Point_Array& input, Point_Array& queries,
				int queries_count, int TableSize, double** s_params,int L, int k, int M, long long int m, double w, double R);

void Exact_NN_readonly(Results* results, int queries_count, int N, string& input_fp);

void Exact_NN(Point_Array& input, Point_Array& queries,int N,int N_q,ofstream& outfile,int* time_passed);
double compute_w(Point_Array& input, int N);


#endif
