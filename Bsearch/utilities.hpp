#ifndef UTILITIES_HPP
#define UTILITIES_HPP


#include <iostream>
#include <cmath>
#include <chrono>
#include <random>
#include <string.h>
#include <string>
#include <sstream>
#include <fstream>
#include <vector>
#include <algorithm>

using namespace std;


class Results{
    public:
        Results();
        ~Results();

        void set_query_id(int id);
        void insert_N_nearest(int nearest_neibhor_id, double distance);
        void insert_exact_N_nearest(int nearest_neighbor_id);
        void insert_Range_nearest(int nearest_neibhor_id); 
        void insert_tTrue(int tTrue);     
        void insert_t_NN(int t_NN);  
        
        int get_query_id();
        vector <int> get_N_nearest_id();
        vector <double> get_N_nearest_distance();
        vector <double> get_exact_N_nearest();
        vector <int> get_Range_nearest();
        int get_tTrue();
        int get_t_NN();

    private:
        int query_id;
        vector <int> N_nearest_id;
        vector <double> N_nearest_distance;
        vector <double> exact_N_nearest;
        vector <int> Range_nearest;
        int exact_NN_time;
        int NN_time;
};



/******************************************
 * Other various Utilities
 *****************************************/


int modular(int base,unsigned int exp , unsigned int mod);
int powxy(long long int x, int y,int M);
int string_mod(string num, int a);
int NumberOfPoints(string&);
double FRandomGen(double low, double high, std::default_random_engine generator);
int reverseInteger (int);

int power(long long int x, int y, int p);
 int mod( int k,  int n);
int bigMod(int  a,int  b,int c);

int moduloMultiplication(long long int a, 
                            int b, 
                            int mod);

int partition (int* arr, int low, int high);
void quickSort(int* arr, int low, int high);
void swap(int* a, int* b);
int binary_search_index(double *A,double key,int left, int right);



#endif
