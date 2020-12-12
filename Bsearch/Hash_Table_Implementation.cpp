#include "Hash_Table.hpp"
#include "Bucket.hpp"

using namespace std;


Hash_Table::Hash_Table(int tablesize){
	buckets = new Bucket[tablesize]; // array of Buckets

}

Hash_Table::~Hash_Table(){
	delete[] buckets;
}

void Hash_Table::InsertToBucket(int position,int id){

	buckets[position].Insert(id);	

}

int Hash_Table::SizeofBucket(int bucket_position){

	return buckets[bucket_position].SizeofBucket();
}

int Hash_Table::Pop_ID(int bucket_position,int element){

	return buckets[bucket_position].Pop_ID(element);
}


