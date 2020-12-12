#include "Bucket.hpp"
#include <iostream>
#include <vector>

using namespace std;

Bucket::Bucket(){
}

Bucket::~Bucket(){
}

void Bucket::Insert(int id){
	records.push_back(id);
}

int Bucket::SizeofBucket(){
	return records.size();
}

int Bucket::Pop_ID(int position){
	int point_id = records.at(position);
	return point_id;
}


