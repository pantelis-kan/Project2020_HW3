/* Methods required globally in order to perform different computations*/

#include "utilities.hpp"

using namespace std;

/******************************************
 * Functions related to class Results
 *****************************************/

Results::Results(){

}

Results::~Results(){

}

void Results::set_query_id(int id){
    query_id = id;
}

void Results::insert_N_nearest(int nearest_neighbor_id, double distance){
	N_nearest_id.emplace_back(nearest_neighbor_id);
    N_nearest_distance.emplace_back(distance);
}

void Results::insert_exact_N_nearest(int nearest_neighbor_id){
    exact_N_nearest.emplace_back(nearest_neighbor_id);
}

void Results::insert_tTrue(int tTrue){
    exact_NN_time = tTrue;
}     

void Results::insert_t_NN(int t_NN){
    NN_time = t_NN;
} 

void Results::insert_Range_nearest(int nearest_neighbor_id){
	Range_nearest.emplace_back(nearest_neighbor_id);
}

int Results::get_query_id(){
    return query_id;
}

vector <int> Results::get_N_nearest_id(){
    return N_nearest_id;
}

vector <double> Results::get_N_nearest_distance(){
    return N_nearest_distance;
}

vector <double> Results::get_exact_N_nearest(){
    return exact_N_nearest;

}

/*
 * Removes possible dublicates
 * Returns the vector with all the range nearests
 */
vector <int> Results::get_Range_nearest(){
    sort(Range_nearest.begin(), Range_nearest.end()); 
    auto last = std::unique(Range_nearest.begin(), Range_nearest.end());
    Range_nearest.erase(last, Range_nearest.end());
    
    return Range_nearest;
}

int Results::get_tTrue(){
    return exact_NN_time;
}

int Results::get_t_NN(){
    return NN_time;
}


/******************************************
 * Other various Utilities
 *****************************************/

/* Reversing an integer */
int reverseInteger (int i){
    unsigned char c1, c2, c3, c4;

    c1 = i & 255;
    c2 = (i >> 8) & 255;
    c3 = (i >> 16) & 255;
    c4 = (i >> 24) & 255;

    return ((int)c1 << 24) + ((int)c2 << 16) + ((int)c3 << 8) + c4;
}

/* Finding the modulo of a number given in string format */
int string_mod(string num, int a){ 
    // Initialize result 
    int res = 0; 
  
    // One by one process all digits of 'num' 
    for (int i = 0; i < num.length(); i++) 
         res = (res*10 + (int)num[i] - '0') %a; 
  
    return res; 
}


int mod( int k,  int n) {
    return ((k %= n) < 0) ? k+n : k;
}


/* x^y mod M */
int powxy(long long int x, int y,int M){

	//cout << "Called with arguments : " << x << "   " << y <<endl;
	if(y==0) return 1;
	if(y%2 == 1) return moduloMultiplication(x, powxy(x,y-1,M), M);
	long long int t = powxy(x,y/2,M);
	//cout << "t = "<< t <<endl;
	return moduloMultiplication(t, t, M);

}


/*Modular exponentiation*/
int power(long long int x, int y,  int p){  
    int res = 1;     // Initialize result  
  
    x = x % p; // Update x if it is more than or  
                // equal to p 
   
    if (x == 0) return 0; // In case x is divisible by p; 
  
    while (y > 0)  
    {  
        // If y is odd, multiply x with result  
        if (y & 1)  
            res = (res*x) % p;  
  
        // y must be even now  
        y = y>>1; // y = y/2  
        x = (x*x) % p;  
    }  

	//cout << "Modular exponentiation returning " << res <<endl;
    return res;  
}  


int moduloMultiplication(long long int a, 
                            int b, 
                            int mod) 
{ 
    int res = 0; // Initialize result 
  
    // Update a if it is more than 
    // or equal to mod 
    a %= mod; 
  
    while (b) 
    { 
        // If b is odd, add a with result 
        if (b & 1) 
            res = (res + a) % mod; 
  
        // Here we assume that doing 2*a 
        // doesn't cause overflow 
        a = (2 * a) % mod; 
  
        b >>= 1; // b = b / 2 
    } 
  
    return res; 
} 

int bigMod(int a,int  b,int c) {
    if (a == 0 || b == 0) {
        return 0;
    }
    if (a == 1) {
        return b;
    }
    if (b == 1) {
        return a;
    } 

    // Returns: (a * b/2) mod c
    long a2 = bigMod(a, b / 2, c);

    // Even factor
    if ((b & 1) == 0) {
        // [((a * b/2) mod c) + ((a * b/2) mod c)] mod c
        return (a2 + a2) % c;
    } else {
        // Odd exponent
        // [(a mod c) + ((a * b/2) mod c) + ((a * b/2) mod c)] mod c
        return ((a % c) + (a2 + a2)) % c;
    }
}



int modular(int base,unsigned int exp , unsigned int mod){

	int x = 1;
	int i;
	int power = base % mod;

	for(i = 0; i< sizeof(int)*8; i++){

		int least_sig_bit = 0x00000001 & (exp >> i);
		
		if(least_sig_bit)
			x = (x*power) % mod;

		power = (power*power) % mod;

	}

	return x;

}


/* Counting the total number of images available in the input file */
int NumberOfPoints(string& input_fp){
	ifstream myfile;
	
	myfile.open(input_fp, ios::out | ios::binary);

	if (!myfile.is_open()) {						/* probably exception handling here */
	    cout << "Cannot open file!" << endl;
		return 1;									/* what return value is proper? NULL? */
	}


	/* Initialize all with zero */
	int magic_number = 0;								/* Read and ignore it! */
	int number_of_images = 0;
	int num_rows = 0;
	int num_cols = 0;
	
	myfile.read((char*)&magic_number, sizeof(magic_number)); 
	magic_number = reverseInteger(magic_number);
	
	myfile.read((char*)&number_of_images, sizeof(number_of_images));
	number_of_images = reverseInteger(number_of_images);

	return number_of_images;

	myfile.close();

}


double FRandomGen(double low, double high, std::default_random_engine generator){

	//Create random (uniform distribution) float number	
	
	uniform_real_distribution<double> distribution (low,high);

	double k = distribution(generator);
	
	return k;
}


int partition (int* arr, int low, int high) 
{ 
    int pivot = arr[high];    // pivot 
    int i = (low - 1);  // Index of smaller element 
  
    for (int j = low; j <= high- 1; j++) 
    { 
        // If current element is smaller than or 
        // equal to pivot 
        if (arr[j] <= pivot) 
        { 
            i++;    // increment index of smaller element 
            swap(&arr[i], &arr[j]); 
        } 
    } 
    swap(&arr[i + 1], &arr[high]); 
    return (i + 1); 
} 
  

void quickSort(int* arr, int low, int high) 
{ 
    if (low < high) 
    { 
        /* pi is partitioning index, arr[p] is now 
           at right place */
        int pi = partition(arr, low, high); 
  
        // Separately sort elements before 
        // partition and after partition 
        quickSort(arr, low, pi - 1); 
        quickSort(arr, pi + 1, high); 
    } 
} 


/*
Returns the right index of key
example: 34.0 < key < 40.0
The function will return the index of 40.0 
*/
int binary_search_index(double *A,double key,int left, int right)
{

  while (left <= right) {
        int middle = left + (right - left) / 2;
        if (A[middle] < key)
            left = middle+1;
        else if(A[middle] > key)
            right = middle-1;
        else
            return middle;
    }
    return (right + 1);
}



void swap(int* a, int* b) 
{ 
    int t = *a; 
    *a = *b; 
    *b = t; 
} 
