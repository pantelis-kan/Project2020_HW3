
#include "Point.hpp"
#include "utilities.hpp"


using namespace std;

Point::Point(){
	dimension = 784;
	is_centroid = false;
	nearest_cluster = -1;
}

int Point::Nearest_Cluster_id(){

	return nearest_cluster;
}

void Point::Assign_Cluster(int c){
	nearest_cluster = c;
}

Point::~Point(){

}


int Point::get_dimension(){
	return dimension;
}


void Point::PrintPoint(){
		for(int j = 0; j < dimension;  j++)
			cout << point[j] << ' ';

		cout << endl << endl;
}

// add val to point vector
void Point::AddtoPoint(int pos, int val){ 

	//cout << "Pushing " << val <<endl;
	point[pos] = val;

}



 // computes h(x) 
 int Point::LSH_Manhattan(int M, const long long int m, double w, double** s_params, int current_k){ 	

	int Dimension = 784;

	int* coeff = new int[Dimension];
	int hash = 0;
	int result;

	int k = Dimension-1;

	for(int i = 0; i < Dimension; i++){


		coeff[i] = floor( (  point[i] - s_params[current_k][i]) / w );

		//cout << m << "  " << k << "  " << M << endl;
		int power_m = power(m,k,M);
		//cout << "powxy returned " << power_m <<endl;

		int coeff_param = mod(coeff[i],M);
		//cout << "coeff[i] " << coeff[i]<< " mod M = " << coeff_param <<endl;


		//result = moduloMultiplication(power_m, coeff[j]%M, M); // (m^k mod M * aj mod M ) mod M
		result = mod(bigMod(power_m , coeff_param,M) , M);
		
		//cout << "Adding " << result <<endl;
		hash += (result % M); 
		--k;
		
	}

	delete[] coeff;
	//coeff = NULL;
	return (hash % M);

}

bool Point::check_centroid(){
	return is_centroid;
}

void Point::set_centroid(bool val){

	is_centroid = val;
}


double Distance(Point& v1, Point& v2, int metric){

	double dist = 0.0;
	double sum = 0.0;

	int Dimension = v1.get_dimension();
	

	for(int i = 0; i < Dimension; i++){
		
		sum = abs(v1.point[i] - v2.point[i]);
		if(metric > 1) sum = pow(sum,metric);


		//cout << "Iteration "<< i << ". Vector numbers : " << (*it1)<<"  "<<(*it2)<< ". Adding  "<<sum<<endl;
		dist += sum;
		//cout << "Dist in while : " << dist <<endl;

	}

	
	//cout << "dist before root : "<<dist <<endl;
	if(metric > 1) dist = pow(dist,1.0/metric);

	return dist;

}

void Copy_Points(Point* from, Point* to){

	int dimension = from->get_dimension();

	for(int i = 0; i < dimension; i++){

		to->point[i] = from->point[i];
	}
}

int Point::get_coordinate(int pos){
	return point[pos];
}

void Point::set_coordinate(int pos,int val){
	point[pos] = val;
}