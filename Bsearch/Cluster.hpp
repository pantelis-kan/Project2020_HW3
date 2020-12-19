#ifndef CLUSTER_HPP
#define CLUSTER_HPP

#include <iostream>
#include <vector>

#include "Point.hpp"
#include "Point_Table.hpp"

using namespace std;



class Cluster{

    public:
        Cluster();
        ~Cluster();
        void Assign_Point(int id); // assign the point id to the cluster
        void Assign_Centroid(Point* c,int id);
        void Compute_New_Centroid(Point_Array& input,int cluster_num);
        Point* get_centroid();
        void Remove_Point(int id);
        int Cluster_Size();
        int Retrieve_ID(int position);
        void Add_Objective(double val);
        long double Get_Objective();
        void Compute_Objective(Point_Array& input);

    private:

        Point* centroid; // centroid is a Point
        vector<int> points; // the point ids that belong to the cluster
        int centroid_id; // the id of the centroid

        long double objective_function_value;


};

#endif