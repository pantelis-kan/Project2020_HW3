#ifndef HASH_TABLE_HPP
#define HASH_TABLE_HPP

#include "Bucket.hpp"

using namespace std;

class Hash_Table{

	public:
		Hash_Table(int tablesize);
		~Hash_Table();
		void InsertToBucket(int position,int id);
		int SizeofBucket(int bucket_position); // return the ammount of elements from a specific bucket
		int Pop_ID(int bucket_position,int element); // pop the id from a specific bucket 

	private:
		Bucket* buckets;
		
};


#endif
