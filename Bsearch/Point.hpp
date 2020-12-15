#ifndef POINT_HPP
#define POINT_HPP

#include <iostream>
#include <string>
#include <vector>
#include <iterator>
#include <random>

using namespace std;

class Point{


	public:
		Point();
		void Point_init(int dim);
		~Point();
		int get_dimension(); // getter
		int PrintPoint();  // prints the coordinates of the image
//		void AddtoPoint(int pos,int val);  // adds the integer 'val' to the coordinate in position 'pos'
void AddtoPoint(int pos, int* tempdata);


		bool check_centroid();   // checks wether this point is a centroid
		void set_centroid(bool val);  // sets/unsets the point as a centroid 
		void Assign_Cluster(int c);  // assigns the point to a cluster pointed by 'c'
		int Nearest_Cluster_id();    // returns the id of the point's current cluster 
		int get_coordinate(int pos); // returns the coordinate in position 'pos'
		void set_coordinate(int pos,int val); // sets 'val' as a coordinate in position 'pos'
		
		// computes h_currentk(x) hash value for the point x. 0 <= currentk <= k
		int LSH_Manhattan(int M, const long long int m, double w, double** s_params,int current_k);

		// returns the distance between two points with a specific metric
		friend double Distance(Point& p1, Point& p2, int metric); 
		friend void Copy_Points(Point* from, Point* to);

	private:

		int * point;
		//int point[784]; // a collection of integers for each point
		int dimension;
		bool is_centroid; // wether the point is randomly chosen as a centroid in k-means++
		int nearest_cluster;  // 
};

#endif
