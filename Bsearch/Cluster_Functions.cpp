
#include "Cluster_Functions.hpp"
#include <random>
#include <set>
#include <chrono>
#include "utilities.hpp"

using namespace std;

// Initilize uniform random generator
std::default_random_engine rand_generator(time(NULL));


int k_lsh = 4,L = 5;
double w = 30000.0;
int k_hypercube = 3;
int probes = 2;
int Max_elements = 10;

//const int M = 4294967291; 3.43597e+10

int M = pow(2,32/k_lsh);

const long long int m = 4294967291; // 2^32  - 5


std::default_random_engine param_rand_generator(time(NULL));

/* Min distance from a point to all the available centroids
    size = number of already assigned centroids
    Returns the nearest cluster index (0-k) 
*/

int Min_Centroid_From_Point(Point& point, Cluster* clusters,int size){

	double min_distance = std::numeric_limits<double>::max();
    double distance;

    int nearest_centroid_id = 0;

    for(int i = 0; i < size; i++){
        Point& centr = *(clusters[i].get_centroid());

        distance = Distance(point,centr,1);

        if (distance < min_distance){
			min_distance = distance;
			nearest_centroid_id = i;
		}	 
    }

    return nearest_centroid_id;
}


double Min_Centroid_Distance(Point& point, Cluster* clusters,int size){

	double min_distance = std::numeric_limits<double>::max();
    double distance;

    int nearest_centroid_id = 0;

    for(int i = 0; i < size; i++){
        Point& centr = *(clusters[i].get_centroid());

        distance = Distance(point,centr,1);

        if (distance < min_distance) min_distance = distance;
    }

    return min_distance;
}




void Loyds_Clusters(Point_Array& input,Cluster* clusters,int k){

	bool not_converged = true;
	int loops = 0;

	int old_assigned;
	// in the first loop, all the points will be assigned
	int new_assigned = input.get_ArraySize();

	while(not_converged == true){

		++loops;
		if(loops > 40) break;

		cout << endl << "Loop " << loops << endl;
		old_assigned = new_assigned;
		bool changed = Loyds_Assign(input,clusters,k,&new_assigned);

		if(loops > 1){
			if(abs(new_assigned - old_assigned) <= 10){
				not_converged = false;
				break;
			}
		}

		if(changed == false){
			not_converged = false;
			break;
		}

		Update(input,clusters,k);
	}

	if(not_converged == true) cout << "Clustering failed"<<endl;
	else cout << "Clustering successful" << endl;
}


void Update(Point_Array& input,Cluster* clusters,int k){

	cout << endl;
	for(int j = 0; j < k; j++){
		cout << "Cluster " << j << " has : " << clusters[j].Cluster_Size() << " points " << " and objective function value = " << clusters[j].Get_Objective() <<endl;
	}

	// for each cluster
	for(int i = 0; i < k; i++){

		//long double old_obj = clusters[i].Get_Objective();

		//cout << "Cluster " << i << " has obj function before update : " << old_obj <<endl;
		
		clusters[i].Compute_New_Centroid(input,i);
		//cout << "Cluster " << i << " updated" << endl;
		//cout << "Cluster " << i << " has obj function after update : " << clusters[i].Get_Objective()<<endl;
	
	}

}


bool Loyds_Assign(Point_Array& input,Cluster* clusters,int k,int* assigned){

	int input_points = input.get_ArraySize();
	int nearest_cluster = 0;
	double min_distance;
    double distance;

	bool changed = true;

	int points_not_changed = 0;

	// for each point
	for(int i = 0; i < input_points; i++){
		
		min_distance = std::numeric_limits<double>::max();

		Point& point = input.Retrieve(i);
		//if(point.check_centroid() == true) continue;

		// for each cluster
		for(int j = 0; j < k; j++){

			// compute the distance from a point to all the clusters
			Point& centr = *(clusters[j].get_centroid());
			distance = Distance(point,centr,1);
			
			if (distance < min_distance){
				min_distance = distance;
				nearest_cluster = j;

			}
			
		}

		int point_nearest_cluster = point.Nearest_Cluster_id();
		
		// check if point needs to change cluster
		if( point_nearest_cluster != nearest_cluster){

			// it's not the first time the point has been assigned to a cluster
			// remove the point from the old cluster
			if(point_nearest_cluster != -1){
				clusters[point_nearest_cluster].Remove_Point(i);

				Point& old_centr = *(clusters[point_nearest_cluster].get_centroid());
				int old_distance = Distance(point,old_centr,1);
				

				 // decrease the obj function of the old cluster
				clusters[point_nearest_cluster].Add_Objective( - old_distance);
			}

			//cout << "Nearest cluster for point " << i << " is " << nearest_cluster <<endl;
			
			clusters[nearest_cluster].Assign_Point(i);
			
			clusters[nearest_cluster].Add_Objective(min_distance);
			point.Assign_Cluster(nearest_cluster);
			
		}
		else{
			++points_not_changed;
		}
	}

	// There is no change in assignment
	if(points_not_changed == input_points) changed = false;

	*assigned = input_points - points_not_changed;
	cout << "New points assigned = " << input_points - points_not_changed <<endl;
	return changed;
}



void Initialize_Centroids(Point_Array& input,Cluster* clusters,int k){

    int input_points = input.get_ArraySize();

    // Pick the first centroid uniformly at random
    std::uniform_real_distribution<double> first_centroid(0.0,(double)input_points - 1); 
    double first = first_centroid(rand_generator);
	

    int first_centr = (int)first;

    input.set_centroid(first_centr,true); // mark the first centroid

    Point* point = input.Retrieve_ptr(first_centr);
    clusters[0].Assign_Centroid(point,first_centr); // assign first centroid

    cout << "1st centroid chosen at index " << first_centr + 1 << endl;

    // conserve memory space 
    double prob[input_points];
    prob[0] = 0;

    for(int t = 1; t < k; t++ ){


        /* Compute distances and probabilities.
            They must be re-calculated each time a new
            centroid is selected.
        */

        int position = 0;
        for(int i = 0; i < input_points; i++){

            
            // if the point is a centroid, ignore
            if(input.check_centroid(i) == true) continue;
            
            Point& pt = input.Retrieve(i);

            // P[0] = 0
            // P[i] = P[i-1] + (D[i])^2
            if(position > 0){
                prob[position] = prob[position-1] + pow(Min_Centroid_Distance(pt,clusters,t),2);
            }

            ++position;
        }

        // the maximum element is always the last element
        double max_P = prob[input_points -t-1];

        // Set up a uniform random distribution from (0,P(r)]
        std::uniform_real_distribution<double> distribution(1.0,max_P); 
        double x = distribution(rand_generator);

        int index_of_x = binary_search_index(prob, x, 0, input_points - t);
        bool is_centr = true;

        /*
            If the index of x is an already assigned centroid,
            increment the index until you find a point that isn't assigned as centroid
        */
        while(is_centr == true){

            Point* next_centroid = input.Retrieve_ptr(index_of_x);
            
            if(next_centroid->check_centroid() == true) index_of_x++;
            else{
                // centroid found. Assign to cluster
                cout << t+1 << " centroid chosen at index "<< index_of_x << endl;
                clusters[t].Assign_Centroid(next_centroid,index_of_x);
                is_centr = false;
            }
        }


    }

}


void Distance_From_Centroids(Cluster* clusters, int k){

    	
	for(int i = 0; i < k; i++){
		Point& centr1 = *(clusters[i].get_centroid());

		for(int j = 0; j < k; j++){
			if(i==j) continue;
			Point& centr2 = *(clusters[j].get_centroid());

			cout << "Distance from centroid " << i << " to centroid " << j << " : " << Distance(centr1,centr2,1) <<endl;
		}
	}
	
}


void Reverse_Assignment(Point_Array& input,Cluster* clusters,int k,bool lsh){

	bool not_converged = true;
	int loops = 0;

	int old_assigned;
	// in the first loop, all the points will be assigned
	int new_assigned = input.get_ArraySize();

	int TableSize = input.get_ArraySize() / 8;
	int dimension = input.get_dimension();

	Hash_Table**  H_Tables;
	Hypercube* cube;

	double**  s_params_cube = new double*[k_lsh];
	for (int i = 0; i < k_lsh; i++){
		s_params_cube[i] = new double[dimension];
	} 

	std::uniform_real_distribution<double> distribution(0.0,w);

	// create s parameters for each hi. Since there is only one hash table, there are k rows in the table
	for(int i = 0; i < k_lsh; i++){
		for(int j = 0; j < dimension; j++) {
			double rand = distribution(rand_generator);
			s_params_cube[i][j] = rand;
		}
	}


	// Both the lsh and the hypercube use the same parameters
	double**  s_params = new double*[L*k];
	
	for (int i = 0; i < L*k_lsh; i++){
		s_params[i] = new double[dimension];
	} 


	// create s parameters for each h(i). Since there are L hash tables, there are L*k rows in the table
	for(int i = 0; i < L*k_lsh; i++){
		
		for(int j = 0; j < dimension; j++) {
			double rand = distribution(param_rand_generator);
			s_params[i][j] = rand;
		}
	}



	// start LSH routine
	if(lsh == true){
	
		//create L hash tables
		H_Tables = new Hash_Table*[L];

		for (int i = 0; i < L; i++){	
			H_Tables[i] = new Hash_Table(TableSize);
		}

		cout << "Hashing input to buckets " <<endl;

		auto t1 = std::chrono::high_resolution_clock::now();	
		// hashing the input to buckets
		Preprocessing(H_Tables, input, input.get_ArraySize(), TableSize, s_params, L, k_lsh, M, m, w);
		
		auto t2 = std::chrono::high_resolution_clock::now();		
		auto duration = std::chrono::duration_cast<std::chrono::seconds>( t2 - t1 ).count();
		cout << "Time taken : " << duration << " seconds" << endl;
		
	}	 
	// start HyperCube routine
	else{

		// create a Hypercube instance
		cube = new Hypercube(k_hypercube);
		
		//Map points to Hypercube 
		cube->Map_images(input, input.get_ArraySize(), k_hypercube, s_params, M, m, w, cube);
	}

	while(not_converged == true){
		++loops;
		if(loops > 50) break;

		cout << endl << "Loop " << loops << endl;
		old_assigned = new_assigned;
		bool changed;
		
		if(lsh == true) changed = LSH_Reverse_Assignment(input,clusters,k,H_Tables,s_params,&new_assigned, TableSize);
		else changed = Hypercube_Reverse_Assignment(input,clusters,k,s_params_cube,cube,&new_assigned,probes,Max_elements);


		if(loops > 1){
			if(abs(new_assigned - old_assigned) <= 10){
				not_converged = false;
				break;
			}
		}

		if(changed == false){
			not_converged = false;
			break;
		}

		Update(input,clusters,k);
		cout << "Update() finished" <<endl;
	}
	if(not_converged == true) cout << "Clustering failed"<<endl;
	else cout << "Clustering successful" << endl;

	
	if(lsh == true){

		for(int i = 0; i < L;  i++) delete H_Tables[i]; 
		for(int i = 0; i < L*k_lsh; i++) delete[] s_params[i]; 
    	
    	//Free the array of pointers
    	delete[] H_Tables;
		delete[] s_params;
	}
	else{
	
		for (int i = 0; i < k_lsh; i++){
			delete s_params_cube[i];
		} 
		delete s_params_cube;
		delete cube;		
	}
	

}


bool LSH_Reverse_Assignment(Point_Array& input,Cluster* clusters,int k,
				Hash_Table** H_Tables,double** s_params,int* assigned,int TableSize){
	
	int dimension = input.get_dimension();
	int query_bucket_position;
	int points_not_changed = 0;
	int input_points = input.get_ArraySize();
	bool changed = true;
	int points_already_assigned = 0;

	// Insert clusters to a Point_Array class
	Point_Array cluster_table(k);

	// stores distinct point ids that are assigned to a cluster
	std::set<int> already_assigned;

	// for each cluster
	for(int i = 0; i < k; i++){		
		Point* centr = clusters[i].get_centroid();
		Point& centroid_point = *centr;

		// Fill the cluster Point_Array class with the centroids.
		// They must be re-filled at each assignment, since
		// the centroids change coordinates in each Update() call
		for(int j = 0; j < dimension; j++){

			int val = centr->get_coordinate(j);
			cluster_table.AddtoPoint(i,j,val);
		}

		for (int l = 0; l < L; l++){
			

			// find the position of the query in the lth hash table
			query_bucket_position = cluster_table.Compute_g(i,k_lsh,M,m,w,TableSize,s_params,l);

			// find how many elements the bucket has
			int size_of_bucket = H_Tables[l]->SizeofBucket(query_bucket_position);

			// if the query fell on an empty bucket, ignore
			if (size_of_bucket == 0) continue; 

			// for each element in the bucket
			for(int n = 0; n < size_of_bucket; n++ ){

				int id = H_Tables[l]->Pop_ID(query_bucket_position,n);

				Point& bucket_point = input.Retrieve(id-1);
				int nearest_cluster = bucket_point.Nearest_Cluster_id();

				// It's the first assignment for point
				if(nearest_cluster == -1){

					bucket_point.Assign_Cluster(i);
					clusters[i].Assign_Point(id-1);
					already_assigned.insert(id-1);
				}
				else{
					// Ignore if it's the same cluster
					if(nearest_cluster == i){ 
						++points_not_changed;
						continue;
					}

					// The point was already assigned to a different cluster
					// Calculate distance from point to the two clusters
					Point& old_centroid_point = *(clusters[nearest_cluster].get_centroid());
					
					int old_cluster_distance = Distance(old_centroid_point,bucket_point,1);
					int new_cluster_distance = Distance(centroid_point,bucket_point,1);

					// If the new distance is better, assign the point to the new
					// cluster, and dissasociate the old cluster from the point
					if(new_cluster_distance < old_cluster_distance){
						
						bucket_point.Assign_Cluster(i);
						clusters[i].Assign_Point(id-1);	
						clusters[nearest_cluster].Remove_Point(id-1);
						already_assigned.insert(id-1);
					}
					else{ 
						
						++points_not_changed;
					}	
				}

			}

		}	

	
	}


	int not_visited = 0;
	int not_assigned = 0;

	// for the points that were not assigned, compute loyd's distance
	for(int i = 0; i < input_points; i++ ){

		std::set<int>::iterator it = already_assigned.find(i);
		// if i doesn't belong in the set
		if(it == already_assigned.end()){
			++not_visited;
			Point& point = input.Retrieve(i);
			int old_point_cluster = point.Nearest_Cluster_id();

			int nearest_cluster_id = Min_Centroid_From_Point(point,clusters,k);
			if(old_point_cluster != nearest_cluster_id){
				
				++not_assigned;
				//cout << "Point id " << i << " found a better cluster: " << nearest_cluster_id << " .Old cluster id = " << old_point_cluster <<endl;
				point.Assign_Cluster(nearest_cluster_id);
				clusters[nearest_cluster_id].Assign_Point(i);

				if(old_point_cluster != -1) clusters[old_point_cluster].Remove_Point(i);
			}
			
		}

	}

	cout << "Set has " << already_assigned.size() << " points " <<endl;
	cout << "Non visited points : " << not_visited <<endl;
	cout << "Non assigned points : " << not_assigned <<endl;

	points_already_assigned = already_assigned.size() + not_assigned;

	if(points_already_assigned == 0) changed = false;

	*assigned = points_already_assigned;
	cout << "New points assigned = " << *assigned <<endl;
	return changed;

}


bool Hypercube_Reverse_Assignment(Point_Array& input,Cluster* clusters,int k,double** s_params_cube,
									Hypercube* hcube,int* assigned,int probes, int Max_elements){

	vector<int> *bucket_records;

	// Insert clusters to a Point_Array class
	Point_Array cluster_table(k);
	int dimension = input.get_dimension();
	int input_points = input.get_ArraySize();
	int remaining;
	int probes_count;

	// stores distinct point ids that are assigned to a cluster
	std::set<int> already_assigned;

	// for each cluster
	for(int i = 0; i < k; i++){	

		string query_label;	
		Point* centr = clusters[i].get_centroid();
		Point& centroid_point = *centr;

		// Fill the cluster Point_Array class with the centroids.
		// They must be re-filled at each assignment, since
		// the centroids change coordinates in each Update() call
		for(int j = 0; j < dimension; j++){

			int val = centr->get_coordinate(j);
			cluster_table.AddtoPoint(i,j,val);
		}

		// find the position of the query in the cube table
		query_label = cluster_table.Compute_f(i, k_hypercube, M, m, w, s_params_cube, hcube);

		//Initialize Hamming class needed for the probes. Make and delete for every query
		Hamming* hamming = new Hamming(query_label, probes);

		remaining = Max_elements;

		
		do{
			//retrieve pointer to a Vertex which is the actual bucket corresponding to the query_label
			bucket_records = hcube->retrieve_records_vector(query_label);

			while(bucket_records == NULL){
				cout << "Records vector is empty!" << endl;

				if (hamming->get_usedprobes() == probes)
					break;		//Thresold reached: we cannot go further so searching has to stop
				
				//move_to_next: should actually check next in map, change the current_in_use and increase used_probes
				//Returns the new label of the bucket we move to
				query_label = hamming->move_to_next();
				cout << "New query label after probing is: " << query_label << endl;
				//Change bucket to the next one to be checked
				bucket_records = hcube->retrieve_records_vector(query_label);
			
			}

			
			int bucket_size = bucket_records->size();
			cout << "Cluster " << i << " fell in vertex " << query_label << " with size " << bucket_size <<endl;
			
			if(bucket_size == 0) continue;
			

			// for each element in the vertex
			for(int j = 0; j < bucket_size; j++ ){

				
				if(remaining <= 0) break;
				// pop id from the query's bucket
				int id = bucket_records->at(j); 
				Point& bucket_point = input.Retrieve(id-1);

				int nearest_cluster = bucket_point.Nearest_Cluster_id();

				// It's the first assignment for point
				if(nearest_cluster == -1){

					--remaining;
					bucket_point.Assign_Cluster(i);
					clusters[i].Assign_Point(id-1);
					already_assigned.insert(id-1);
				}
				else{
					// Ignore if it's the same cluster
					if(nearest_cluster == i){ 
						//already_assigned.insert(id-1);
						continue;
					}

					// The point was already assigned to a different cluster
					// Calculate distance from point to the two clusters
					Point& old_centroid_point = *(clusters[nearest_cluster].get_centroid());
					
					int old_cluster_distance = Distance(old_centroid_point,bucket_point,1);
					int new_cluster_distance = Distance(centroid_point,bucket_point,1);

					// If the new distance is better, assign the point to the new
					// cluster, and dissasociate the old cluster from the point
					if(new_cluster_distance < old_cluster_distance){
						
						--remaining;
						bucket_point.Assign_Cluster(i);
						clusters[i].Assign_Point(id-1);	
						clusters[nearest_cluster].Remove_Point(id-1);
						already_assigned.insert(id-1);
					}

				}
				
			}		
	
			query_label = hamming->move_to_next();
			
		}while(hamming->get_usedprobes() < probes);

		remaining = Max_elements;
		
		delete hamming;
	}

	
	int not_visited = 0;
	int not_assigned = 0;

	// for the points that were not assigned, compute loyd's distance
	for(int i = 0; i < input_points; i++ ){

		std::set<int>::iterator it = already_assigned.find(i);
		// if i doesn't belong in the set
		if(it == already_assigned.end()){
			++not_visited;
			Point& point = input.Retrieve(i);
			int old_point_cluster = point.Nearest_Cluster_id();

			int nearest_cluster_id = Min_Centroid_From_Point(point,clusters,k);
			if(old_point_cluster != nearest_cluster_id){
				
				++not_assigned;
				//cout << "Point id " << i << " found a better cluster: " << nearest_cluster_id << " .Old cluster id = " << old_point_cluster <<endl;
				point.Assign_Cluster(nearest_cluster_id);
				clusters[nearest_cluster_id].Assign_Point(i);

				if(old_point_cluster != -1) clusters[old_point_cluster].Remove_Point(i);
			}
			
		}

	}

	cout << "Set has " << already_assigned.size() << " points " <<endl;
	cout << "Non visited points : " << not_visited <<endl;
	cout << "Non assigned points : " << not_assigned <<endl;

	int points_already_assigned = already_assigned.size() + not_assigned;

	bool changed = true;
	if(points_already_assigned == 0) changed = false;

	*assigned = points_already_assigned;
	cout << "New points assigned = " << *assigned <<endl;
	return changed;
}


void Silhouette(Point_Array& input,Cluster* clusters,int k,double* s,double* s_total){

	int point_size = input.get_ArraySize();
	cout <<"Input size "<< point_size <<endl;
	int second_nearest = 0;
	double min_distance; 
	int nearest;
	double distance;
	

	// for each point 
	for(int i = 0; i < point_size; i++){

		Point& pt = input.Retrieve(i);

		nearest = pt.Nearest_Cluster_id();
		//cout << i << "nearest cluster = " << nearest <<endl; 
		min_distance = std::numeric_limits<double>::max();


		for(int j = 0; j < k; j++){
			if(j == nearest) continue;

			Point& centr = *(clusters[j].get_centroid());
			distance = Distance(pt,centr,1);

			if (distance < min_distance){
				min_distance = distance;
				second_nearest = j;
			}
			//cout << min_distance << " " << second_nearest << endl;
		}

		//cout << nearest << "  " << second_nearest <<endl;
		
		double average_nearest = Average_Distance(input,pt,clusters[nearest]);
		double average_second_nearest = Average_Distance(input,pt,clusters[second_nearest]);
		//cout << i << " after average "<<endl;

		s[i] = (average_nearest - average_second_nearest) / (std::max(average_nearest,average_second_nearest)) ;

		//if(s[i] <= 1 && s[i] >= 0) cout << "Point "<< i << " was correctly assigned to cluster, s[i] = " << s[i] << endl;
		//else cout << "Point "<< i << " was INcorrectly assigned to cluster, s[i] = " << s[i] << endl;
			
	}

	double average_total = 0.0;

	for(int i = 0; i < point_size; i++){
		average_total += s[i];
	}	 

	average_total = average_total/point_size;
	cout << "Silhouette total average : " << average_total << endl;
	*s_total = average_total;

}

double Average_Distance(Point_Array& input,Point& pt,Cluster& cluster){

	int cluster_size = cluster.Cluster_Size();
	double average = 0.0;

	for(int j = 0; j < cluster_size; j++){
		
		int id = cluster.Retrieve_ID(j);
		Point& point = input.Retrieve(id);

		average += Distance(pt,point,1);
	}

	return (average/cluster_size);
}

void Configuration_File(string filename,int* K){

	string line;

    ifstream infile;
    infile.open(filename, ios::in );
    int line_count =  1;

    string description,value;

    while(getline(infile, line)){

        istringstream iss(line);
        iss >> description >> value;
        cout << description << " " << value << endl;

        switch(line_count){

            case 1:
                *K = stoi(value);
            break;

            case 2:
                L = stoi(value);
            break;

            case 3:
                k_lsh = stoi(value);
            break;

            case 4:
                Max_elements = stoi(value);
            break;

            case 5:
                k_hypercube = stoi(value);
            break;

            case 6:
                probes = stoi(value);
            break;

            default:
            break;
        }

        ++line_count;
    }

   // cout << k << " " << L << " " << k_lsh << " " << M_cube << " " << k_cube << " " << probes << endl;

    infile.close();
}

void Output_Results(Point_Array& input,Cluster* clusters, int k ,double *s, string outputfile,string method,
					double time,double s_total,bool complete){

	int input_points = input.get_ArraySize();
	ofstream outfile;
    outfile.open(outputfile, ios::out | ios::trunc );

	outfile << "Algorithm: ";

	if(method == "Classic") outfile << "Lloyds" <<endl;
	else if(method == "LSH") outfile << "Range Search LSH" <<endl;
	else outfile << "Range Search Hypercube" <<endl;

	int clustersize;

	for(int i = 0; i < k; i++){

		clustersize = clusters[i].Cluster_Size();
		outfile << "CLUSTER-" << i+1 << " {size: " << clustersize << ", centroid: ";

		Point* centr = clusters[i].get_centroid();
		int dimension = centr->get_dimension();

		for(int j = 0; j < dimension; j++){
			outfile << centr->get_coordinate(j) << " ";
		}
		outfile << "}" << endl;
	}

	outfile << "clustering_time: " << time << endl;
	outfile << "Silhouette: [" << endl;

	for(int i = 0; i < input_points; i++){
		outfile << s[i] << ",";
	}
	outfile << " " << s_total << "]" << endl;


	if(complete == true){
		
		for(int i = 0; i < k; i++){

			clustersize = clusters[i].Cluster_Size();
			outfile << "CLUSTER-" << i+1 << "{centroid, ";

			for(int j = 0; j < clustersize; j++){
				outfile << clusters[i].Retrieve_ID(j) << ", ";
			}
			outfile << "}" << endl;
		}
	}

	outfile.close();
}