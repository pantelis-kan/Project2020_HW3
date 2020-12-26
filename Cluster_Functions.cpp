
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
	//point->PrintPoint();
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

void Output_Results(Point_Array& input,Cluster* clusters, int k ,double *s, string outputfile,
					double time,double s_total,string message){

	int input_points = input.get_ArraySize();
	ofstream outfile;
    outfile.open(outputfile, ios::out | ios::trunc );

	outfile << message << endl;

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
	outfile << " " << s_total << "]" << endl << endl;

	
	for(int i = 0; i < k; i++)
		outfile << "CLUSTER-"<<i+1 << " value of Objective Function: " <<  clusters[i].Get_Objective() << endl;
	
	outfile << endl << "Silhouette Total Average: " << s_total << endl;

	outfile.close();
}


void Class_File_Output(Point_Array& input,string class_file,string class_output){

	ofstream outfile;
    outfile.open(class_output, ios::out | ios::trunc );

	outfile << "CLASSES AS CLUSTERS" << endl;

	ifstream infile;
    infile.open(class_file);

    string line;

    int cluster_id = 0;
    int** images_id;
    
	int* cluster_sizes = new int[10];
	int* centroids = new int[10];

    images_id = new int*[10];

    for(int i = 0; i < 10; i++)
        images_id[i] = NULL;

    while(getline(infile, line)){
       // cout << line.size() << endl;
        //line_c = (char*)malloc(line.size()*sizeof(char));
        char line_c[(int)line.size()];
        strcpy(line_c,line.c_str());
        char* token = strtok(line_c, " ");
        char cluster_size_str[] = "size:"; 
        bool found_size = false;
        int position = 0;

        while(token != NULL){

            if(images_id[cluster_id] != NULL){
                images_id[cluster_id][position] = atoi(token);
                ++position;
            }

            if(found_size == true){
                int size = atoi(token);
				cluster_sizes[cluster_id] = size;
                //cout << "Size: " << size << endl;
                found_size = false;
                images_id[cluster_id] = new int[size];
            }

            if(strcmp(token,cluster_size_str) == 0){
                found_size = true;
            }

            token = strtok(NULL," ");
            
        }
        ++cluster_id;
        if(cluster_id == 10) break;
    }

	cout << endl << "Writing to " << class_output << "..." << endl;

	for(int k = 0; k < 10; k++){

		double global_min_distance = std::numeric_limits<double>::max();
		double dist = 0.0;
		int min_pos = 0;

		for(int i = 0; i < cluster_sizes[k]; i++){
			Point& pt1 = input.Retrieve(images_id[k][i]);

			for(int j = 0; j < cluster_sizes[k]; j++){
				if(i != j){
					Point& pt2 = input.Retrieve(images_id[k][j]);
					dist += Distance(pt1,pt2,1);
					//outfile << Distance(pt1,pt2,1) << endl;

				}
			}
			
			if(dist < global_min_distance){
				global_min_distance = dist;
				min_pos = i;
			}

			dist = 0.0;
		}

		centroids[k] = min_pos;

		//cout << "The point with the min distance from every other point is : " <<images_id[k][min_pos] << endl;

		Point& centr = input.Retrieve(images_id[k][min_pos]);
		//cout << endl << "Objective function value for cluster " << k+1 << " : " << global_min_distance << endl;
		outfile << "CLUSTER-"<<k+1 << " value of Objective Function: " <<  global_min_distance << endl;
	}

	cout << "Finding Silhouette for class file" << endl;
	outfile << "Silhouette: [" << endl;

	int point_size = input.get_ArraySize();
	int nearest,second_nearest;
	int centroid_distances[10];
	int centroid_distances_copy[10];
	double dist;
	double s[point_size];

	for(int i = 0; i < point_size; i++){

		Point& pt = input.Retrieve(i);

		for(int j = 0; j < 10; j++){
			Point& centr = input.Retrieve(centroids[j]);
			dist = Distance(pt,centr,1);
			centroid_distances[j] = dist;
			centroid_distances_copy[j] = dist;

		}

		quickSort(centroid_distances,0,9);
		for(int j = 0; j < 10; j++){
			if(centroid_distances[0] == centroid_distances_copy[j])
				nearest = j;
			else if (centroid_distances[1] == centroid_distances_copy[j])
				second_nearest = j;
		}

		double average_nearest;
		double average_second_nearest;
		
		double average = 0.0;
		for(int k = 0; k < cluster_sizes[nearest]; k++){
			Point& pt1 = input.Retrieve(images_id[nearest][k]);
			average += Distance(pt,pt1,1);

		}
		average_nearest = average/cluster_sizes[nearest];

		average = 0.0;
		for(int k = 0; k < cluster_sizes[second_nearest]; k++){
			Point& pt1 = input.Retrieve(images_id[second_nearest][k]);
			average += Distance(pt,pt1,1);

		}
		average_second_nearest = average/cluster_sizes[second_nearest];

		s[i] = (average_nearest - average_second_nearest) / (std::max(average_nearest,average_second_nearest)) ;

		outfile << s[i] << ",";
	
	}	

	double s_total = 0.0;

	for(int i = 0; i < point_size; i++){
		s_total += s[i];
	}

	outfile << " " << s_total << "]" << endl << endl;
	outfile << endl << "Silhouette Total Average: " << s_total << endl;


	for(int i = 0; i < 10; i++)
		delete[] images_id[i];

	delete[] images_id;
	delete[] cluster_sizes;
	delete[] centroids;
	infile.close();
	outfile.close();
}