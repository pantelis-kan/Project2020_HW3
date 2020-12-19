

#include <iostream>
#include <string>
#include <iterator>
#include <random>
#include <algorithm>
#include <limits>
#include <vector>

using namespace std;

#include "Cluster.hpp"
#include "Point_Table.hpp"
#include "Point.hpp"
#include "utilities.hpp"



Cluster::Cluster(){
    centroid = new Point;
    objective_function_value = 0.0;
}

void Cluster::Add_Objective(double val){
    //cout << "Adding to objective function : " << val<<endl;
    objective_function_value += val;
}

long double Cluster::Get_Objective(){
    return objective_function_value;
}


Cluster::~Cluster(){
    delete centroid;
}

void Cluster::Assign_Centroid(Point* c,int id){
   centroid_id = id; 
   Copy_Points(c,centroid);
}

Point* Cluster::get_centroid(){
    return centroid;
}


int Cluster::Cluster_Size(){

    return points.size();
}

void Cluster::Assign_Point(int id){
   
    points.push_back(id);
   // cout << "push back successful current size = " << points.size() << endl;
}

void Cluster::Compute_New_Centroid(Point_Array& input,int cluster_num){

    int dimension = input.get_dimension();

    int points_in_cluster = points.size();
    int nth_position = ceil(points_in_cluster/2); 

    int coordinates[points_in_cluster]; 

	// for each dimension
	for(int i = 0; i < dimension; i++){

		//for each id that belongs in the cluster
		for (int j = 0; j < points_in_cluster; j++){
			int id = points.at(j);

			Point* p = input.Retrieve_ptr(id);
			coordinates[j] = p->get_coordinate(i);
		}

        // sort the coordinate array and put the n/2-th maximum value in the centroid
        // the function sorts in ascending order 
		quickSort(coordinates, 0, points_in_cluster - 1);
		centroid->set_coordinate(i,coordinates[nth_position]);
	}
    
}


void Cluster::Compute_Objective(Point_Array& input){

    long double sum = 0.0;
    int points_in_cluster = points.size();

    Point& centr = *centroid;

    for (int j = 0; j < points_in_cluster; j++){
        int id = points.at(j);

        Point& p = input.Retrieve(id);
        sum += Distance(centr,p,1);
    }

    //cout << endl << "Old objective function value : " << objective_function_value << endl;
   // cout << "New objective function value : " << sum << endl;
    objective_function_value = sum;

}


void Cluster::Remove_Point(int id){

	//points.remove(id);

    auto it = std::find(points.begin(), points.end(), id);

    points.erase(it);
}

int Cluster::Retrieve_ID(int position){

    return points.at(position);
}