CC=g++
DEPS = Bucket.hpp Hash_Table.hpp NN_Functions.hpp Point_Table.hpp Point.hpp utilities.hpp Hypercube.hpp Vertex.hpp
OBJ = lsh.o Bucket_Implementation.o Hash_Table_Implementation.o NN_Functions.o Point_Table_Implementation.o Point_Implementation.o utilities.o Hypercube_Implementation.o Vertex_Implementation.o
CFLAGS  = -g -Wall

all: search cluster exact_reduced exact_original

exact_original: exact_nn_original.o Point_Table_Implementation.o Point_Implementation.o utilities.o NN_Functions.o Hash_Table_Implementation.o Bucket_Implementation.o
	$(CC) $(CFLAGS) -o exact_original exact_nn_reduced.o Point_Table_Implementation.o Point_Implementation.o utilities.o NN_Functions.o Hash_Table_Implementation.o Bucket_Implementation.o 

exact_reduced: exact_nn_reduced.o Point_Table_Implementation.o Point_Implementation.o utilities.o NN_Functions.o Hash_Table_Implementation.o Bucket_Implementation.o
	$(CC) $(CFLAGS) -o exact_reduced exact_nn_reduced.o Point_Table_Implementation.o Point_Implementation.o utilities.o NN_Functions.o Hash_Table_Implementation.o Bucket_Implementation.o

search: lsh.o Bucket_Implementation.o Hash_Table_Implementation.o NN_Functions.o Point_Table_Implementation.o Point_Implementation.o utilities.o 
	$(CC) $(CFLAGS) -o search lsh.o Bucket_Implementation.o Hash_Table_Implementation.o NN_Functions.o Point_Table_Implementation.o Point_Implementation.o utilities.o

cluster: cluster.o Cluster_Implementation.o Point_Implementation.o Point_Table_Implementation.o utilities.o Hash_Table_Implementation.o NN_Functions.o Bucket_Implementation.o Cluster_Functions.o 
	$(CC) $(CFLAGS) -o cluster cluster.o Cluster_Implementation.o Point_Implementation.o Point_Table_Implementation.o utilities.o Hash_Table_Implementation.o NN_Functions.o Bucket_Implementation.o Cluster_Functions.o 

cluster.o: cluster.cpp 
	$(CC) -c cluster.cpp

Cluster_Implementation.o: Cluster_Implementation.cpp
	$(CC) -c Cluster_Implementation.cpp

Cluster_Functions.o: Cluster_Functions.cpp
	$(CC) -c Cluster_Functions.cpp

lsh.o: lsh.cpp
	$(CC) -c lsh.cpp

Bucket_Implementation.o: Bucket_Implementation.cpp Bucket.hpp
	$(CC) -c Bucket_Implementation.cpp

Hash_Table_Implementation.o: Hash_Table_Implementation.cpp Hash_Table.hpp
	$(CC) -c Hash_Table_Implementation.cpp

NN_Functions.o: NN_Functions.cpp NN_Functions.hpp
	$(CC) -c NN_Functions.cpp

Point_Implementation.o: Point_Implementation.cpp Point.hpp
	$(CC) -c Point_Implementation.cpp

Point_Table_Implementation.o: Point_Table_Implementation.cpp Point_Table.hpp
	$(CC) -c Point_Table_Implementation.cpp

utilities.o: utilities.cpp utilities.hpp
	$(CC) -c utilities.cpp

.PHONY: clean //necessary in case file with name clean exists
clean:
	rm *.o search cluster exact_original exact_reduced