#ifndef BUCKET_HPP
#define BUCKET_HPP

#include <vector>

using namespace std;

class Bucket{

	public:

		Bucket();
		~Bucket();
		void Insert(int id);
		int SizeofBucket();
		int Pop_ID(int position);

	private:

		// store the id of each point
		// no information repetition, since the id is only store once in the bucket

		vector<int> records;
		

};

#endif
