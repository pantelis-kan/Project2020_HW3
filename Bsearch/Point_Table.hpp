
// a collection of class Point. Fixed size, which is chosen at the construction

#ifndef POINT_TABLE_HPP
#define POINT_TABLE_HPP

#include <iostream>
#include <random>

#include "Point.hpp"

class Hypercube;
using namespace std;

class Point_Array{

	public:
		Point_Array(int num);
		~Point_Array();
		int get_dimension(); // returns dimension of an arbitrary point (first), because all points have the same dimension

		int get_ArraySize();

		// compute h(x) for x = points[position] 
		//long long int LSH_Manhattan( int M, const long long int m, double w, int position); 
		int FillPoints(string &input_fp);

		void PrintPoint(int position);
		bool check_centroid(int pos);
		void set_centroid(int pos,bool val);
		void AddtoPoint(int point_position,int coordinate_pos,int val);

		// compute g(x) , by concatenation of h1(x)...hk(x)
		int Compute_g(int position,int k,int M, const long long int m, double w,int tablesize, double** s_params,int l); 

		// used by the Hypercube structure. Computes f1(x)...fd'(x)
		string Compute_f(int position, int k, int M, const long long int m, 
					double w, double** s_params, Hypercube* hcube);

		Point& Retrieve(int position);
		Point* Retrieve_ptr(int position);

	private:
		Point* points;
		int ArraySize = 0;

};

#endif
