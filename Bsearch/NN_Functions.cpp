
#include "NN_Functions.hpp"

#include <algorithm>
#include <chrono>
#include <map>

using namespace std;
using namespace std::chrono;


/* 	Assigns all images to buckets
	int input_count : is the number of total images in the input 
	int TableSize : the size of the hashtable
	double** s_params: it's the s(i) necessary for the aplification of LSH
*/
void Preprocessing(Hash_Table** H_Tables, Point_Array& input, int input_count, int TableSize, double** s_params, int L, 
					int k, int M, long long int m, double w){


	int bucket_position;

	// for each point x in the input
	for(int j = 0; j < input_count; j++ ){

		auto t1 = std::chrono::high_resolution_clock::now();

		// Compute g1(x) ... gL(x)
		for (int l = 0; l < L; l++){
		
			bucket_position = input.Compute_g(j, k, M, m, w, TableSize, s_params, l);
			

			//cout << "Function g(x) computed in " << duration <<endl;

			H_Tables[l]->InsertToBucket(bucket_position, j+1); // Insert to bucket point with id j+1	

			//cout << " point id " << j+1 << " has bucket " << bucket_position <<endl;
		}

		auto t2 = std::chrono::high_resolution_clock::now();		
		auto duration = std::chrono::duration_cast<std::chrono::microseconds>( t2 - t1 ).count();

		//cout << "Point " << j+1 << " inserted in " << duration << " microseconds " << endl;
		//cout << endl;

	}


}




/* 	Searches for up to N nearest neghbors
	Stores results to temporary structure
	Returns results through class Results
*/
void LSH_Nearest_Neighbors(Results* results, Hash_Table** H_Tables,Point_Array& input, Point_Array& queries,
				int queries_count,int TableSize, double** s_params,int L,int k,int M, long long int m, double w, int N){

	int query_bucket_position, nearest_neighbor_id;
	double min_distance;


	// for each query
	for (int q = 0; q < queries_count; q++){

		min_distance  = std::numeric_limits<double>::max();
		bool found_nn = false;
		multimap<double, int> all_NN_storage; //temporary storage of results
		double min_distance_previous = 0;	//to avoid taking everytime the minimum -> we set the minimum to the latest found

		auto t1 = std::chrono::high_resolution_clock::now();		

		//for up to N neighbors (we may stop earlier)
		for(int c = 0; c < N; c++){
			// for each hash table
			for (int l = 0; l < L; l++){

				// find the position of the query in the lth hash table
				query_bucket_position = queries.Compute_g(q, k, M, m, w, TableSize, s_params, l);

				// find how many elements the bucket has
				int size_of_bucket = H_Tables[l]->SizeofBucket(query_bucket_position);

				// if the query fell on an empty bucket, ignore
				if (size_of_bucket == 0) continue; 

				// for each element in the bucket
				for(int i = 0; i < size_of_bucket; i++ ){

					if (i > 10*L) break; // 10*L are enough points

					// pop id from the query's bucket
					int id = H_Tables[l]->Pop_ID(query_bucket_position,i); 

					Point& query_point = queries.Retrieve(q);
					Point& input_point = input.Retrieve(id-1);

					// compute Manhattan Distance for the query and the popped id
					double distance = Distance(query_point,input_point,1); 

					//check if the input we retrieved is a suitable nearest neighbor
					if (distance < min_distance && distance > min_distance_previous){
						min_distance = distance;
						nearest_neighbor_id = id;
						found_nn = true;
					}

				}
			}	

			//if a nearest neighbor was found then store it and reset values to continue to next one
			if(found_nn == true){
				//cout << "Approximate NN for query " << q+1 << " = " << nearest_neighbor_id << " with distance " << min_distance <<endl;
				//make a multimap and insert there. then in the end results take N from there
  				all_NN_storage.insert( pair<double, int>(min_distance, nearest_neighbor_id) );
				min_distance_previous = min_distance;
				min_distance  = numeric_limits<double>::max();
			}
			else{
				//cout << "Could not find any other approximate nearest neighbor for query " << q+1 <<endl;
				break;
			}
		}

		auto t2 = std::chrono::high_resolution_clock::now();
		auto duration = std::chrono::duration_cast<std::chrono::microseconds>( t2 - t1 ).count();

		//store findings to results, after completing the search.
		results[q].insert_t_NN(duration);

		results[q].set_query_id(q+1);
		
		for (auto it = all_NN_storage.cbegin(); it != all_NN_storage.cend(); it++){
			results[q].insert_N_nearest((*it).second, (*it).first);
		}

		//check if empty, means that no nearest found
		if (all_NN_storage.empty()){
			cout << "No Approximate nearest found for query image: " << q+1 << endl; 
		}
		

	}


}


/* 	Searches for all images that are up to R range far from the query
	Returns results through class Results
*/
void LSH_Range_Search(Results* results, Hash_Table** H_Tables, Point_Array& input, Point_Array& queries, int queries_count, int TableSize, double** s_params,int L, int k, int M, long long int m, double w, double R){
	int query_bucket_position, nearest_neighbor_id;
	double min_distance;

	// for each query
	for (int q = 0; q < queries_count; q++){

		bool found_nn;
		multimap<double, int> all_NN_storage;

		// for each hash table
		for (int l = 0; l < L; l++){
			min_distance  = R;
			// find the position of the query in the lth hash table
			query_bucket_position = queries.Compute_g(q, k, M, m, w, TableSize, s_params, l);

			// find how many elements the bucket has
			int size_of_bucket = H_Tables[l]->SizeofBucket(query_bucket_position);

			// if the query fell on an empty bucket, ignore
			if (size_of_bucket == 0) continue; 

			// for each element in the bucket
			for(int i = 0; i < size_of_bucket; i++ ){

				if (i > 10*L) break; // 10*L are enough points

				// pop id from the query's bucket
				int id = H_Tables[l]->Pop_ID(query_bucket_position,i); 

				Point& query_point = queries.Retrieve(q);
				Point& input_point = input.Retrieve(id-1);

				// compute Manhattan Distance for the query and the popped id
				double distance = Distance(query_point,input_point,1); 
				if (distance < min_distance){
					nearest_neighbor_id = id;
					results[q].insert_Range_nearest(nearest_neighbor_id);
				}

			}

		}	

	}
}

/*	Exact_NN is prerun (because it takes long time) and results are stored in a .txt file
	This function reads from the file the N required number of Exact NN for every query. 
	N can be up to 50 neighbors. For more -1 appears next to the result
*/
void Exact_NN_readonly(Results* results, int queries_count, int N, string& input_fp){
	string line;
	int nearest_neighbor_id;
	int distance;
	int tTrue;

	ifstream myfile;
	myfile.open(input_fp, ios::out);

	if (!myfile.is_open()) {
	    cout << "Cannot open file with Exact NN results!" << endl;
		exit(1);
	}

	for (int i = 0; i < queries_count; i++){
		for (int j = 0; j < 1; j++){
			getline(myfile, line);
			istringstream iss(line);
        	iss >> nearest_neighbor_id >> distance >> tTrue;
			results[i].insert_exact_N_nearest(nearest_neighbor_id, distance);
		}

		results[i].insert_tTrue(tTrue);
		// if (N > 50){
		// 	int remaining = N - 50;
		// 	for (int q = 0; q < remaining; q++){
		// 		results[i].insert_exact_N_nearest(-1);
		// 	}
		// }
	}

	
	myfile.close();

}

/* 	Finding 1 exact nearest neghbor for every query
	Saves the tTrue time
*/
// void Exact_NN(Point_Array& input, Point_Array& queries, int input_count, int queries_count, ofstream& outfile,int* time_passed){
void Exact_NN(Point_Array& input, Point_Array& queries, int input_count, int queries_count, ofstream& outfile){

	int nearest_neighbor_id;
	double min_distance;
	double distance;

	long double dist_sum = 0.0;
	
    multimap <double, int> distances; // empty multimap container - <distance,ID>

	
	for(int j = 0; j < queries_count; j++){
	// for(int j = 0; j < 600; j++){
		auto t1 = std::chrono::high_resolution_clock::now();
		Point& query_point = queries.Retrieve(j);

		for(int i = 0; i < input_count; i++){
		// for(int i = 0; i < 100; i++){
	
			Point& input_point = input.Retrieve(i);
			distance = Distance(query_point,input_point,1);
			
			// distances[i] = distance;
			distances.insert(pair<double,int> (distance, i+1));

	 	}

		auto t2 = std::chrono::high_resolution_clock::now();		
		auto duration = std::chrono::duration_cast<std::chrono::microseconds>( t2 - t1 ).count();	 
	
		multimap <double, int> :: iterator itr;
		itr = distances.begin();
		for (int i = 0; i < 1; i++){
		// for (int i = 0; i < 10; i++){
			
			if(i == 0) outfile << itr->second << " " << itr->first << " " << duration << endl;
			else outfile << itr->second << " " << itr->first << endl;

			++itr;
			if(itr == distances.end())
				break;

		}

		distances.clear();

	}

}


/* 	Computes w in a brute force way
	Takes a lot of time so it has not been finally used
*/
double compute_w(Point_Array& input, int input_count){
	double distance_sum = 0.0;
	double distance_a_b = 0.0;
	int count_of_sums = 0;

	Point* image_a;
	Point* image_b;

	cout << "Calculating started now!!" << endl;

	for (int i = 0; i < input_count; i += 2){
		image_a = input.Retrieve_ptr(i);
		image_b = input.Retrieve_ptr(i+1);
		distance_a_b = Distance(*image_a, *image_b, 1);
		distance_sum += distance_a_b;
		count_of_sums++;
		cout << "we are at " << i << endl;
	}


	cout << "average distance is:" << distance_sum/count_of_sums << endl;
	cout << "count of sums:" << count_of_sums << endl;

	return distance_sum/count_of_sums;

}
