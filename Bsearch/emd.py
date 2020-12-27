import sys
import matplotlib.pyplot as plt
import numpy as np
import tensorflow as tf
import keras
import time
import math
from keras.optimizers import Adam,RMSprop
from keras.layers import Dense, Dropout, Activation, Flatten
from keras.layers import UpSampling2D,Conv2D, MaxPooling2D, BatchNormalization
from sklearn.model_selection import train_test_split
from keras.utils import to_categorical
import scipy
import warnings


def read_labels(filename, num_images):
    f = open(filename, "rb")
    f.read(4)

    num_images = int.from_bytes(f.read(4), 'big')
    print("Total images in labels file: ", num_images)

    buf = f.read(1 * num_images)
    labels = np.frombuffer(buf, dtype=np.uint8).astype(np.int64)
    return labels


def read_input(filename, num_images):

    f = open(filename, "rb")
    image_size = 28

    #ignore header
    f.read(4)

    num_images = int.from_bytes(f.read(4), 'big')
#    print("Total images in file: ", num_images)

    f.read(8)

    buf = f.read(image_size * image_size * num_images)
    #data = np.frombuffer(buf, dtype=np.uint8).astype(np.float32) / 255
    data = np.frombuffer(buf, dtype=np.uint8).astype(np.float32) / 255
    data = data.reshape(num_images, image_size, image_size,1)
    return data

def show_image(img):
    image = np.asarray(img).squeeze()
    img = plt.imshow(image, cmap='gray')
    #plt.savefig("/content/drive/MyDrive/ML/data/image.png")
    #plt.show()
    img.set_cmap('gray')
    plt.axis('off')
    #plt.savefig("/content/drive/MyDrive/ML/data/test.png", bbox_inches='tight')


def split(array, nrows, ncols):
    #Split a matrix into sub-matrices.

    r, h,z = array.shape
    return (array.reshape(h//nrows, nrows, -1, ncols)
                 .swapaxes(1, 2)
                 .reshape(-1, nrows, ncols))    



def compute_signature(image,cluster_size):
  
  clusters = split(image,cluster_size,cluster_size)

  clusters = clusters.reshape(clusters.shape[0],clusters.shape[1],clusters.shape[2],1)
  #print(clusters.shape)

  total_clusters = clusters.shape[0]

  #show_image(clusters[1])

  total_weights = 0.0

  for i in range(total_clusters):
    total_weights += np.sum(clusters[i])

  #print(total_weights)

  cluster_weights = []
  for i in range(total_clusters):
    cluster_weights.append(np.sum(clusters[i]) / total_weights)
    #cluster_weights.append(np.sum(clusters[i]))

  #print(cluster_weights)


  cluster_centroids = []
  cluster_centroids.append((np.ceil(cluster_size/2)-1,np.ceil(cluster_size/2)-1))

  for i in range(total_clusters-1):

    previous_y = cluster_centroids[i][1]
    previous_x = cluster_centroids[i][0]

    # If the column index exceeds the size of the array, 
    # reposition the y coordinate at the start and change the x coordinate
    if previous_y + cluster_size > image_dimension-1:
      cluster_centroids.append((previous_x + cluster_size , np.ceil(cluster_size/2)-1))
    else:
      cluster_centroids.append((previous_x , previous_y + cluster_size))

  #print(cluster_centroids)

  signature = (list(zip(cluster_centroids,cluster_weights)))
  return signature    

def compute_distances(input_img,query_img,input_signature,query_signature):
  distances = []

  query_length = len(query_signature) 
  input_length = len(input_signature)

  for i in range(input_length):
    input_centroid = np.array(input_signature[i][0])
    for j in range(query_length):
      query_centroid = np.array(query_signature[j][0])
      distances.append(np.linalg.norm(query_centroid - input_centroid))

  return distances    



def compute_EMD(input_img,query_img,input_signature,query_signature,aub,distances):

  #query_length = len(query_signatures[0]) 
  #input_length = len(input_signatures[0])

  query_length = len(query_signature) 
  input_length = len(input_signature)

  #print(query_length)
  #print(input_length)
  #print(len(distances))

  b = []
  query_weight_sum = 0.0

  for i in range(input_length):
    b.append(input_signature[i][1])

  for i in range(query_length):
    query_weight_sum += query_signature[i][1]
    b.append(query_signature[i][1])
    
  #print(len(aub))
  #print(len(b))
  res = scipy.optimize.linprog(distances,A_eq = aub, b_eq = b,options={'cholesky':False,'sym_pos':False} )
  return res.fun




# Handle command line arguments
if len(sys.argv) < 10:
    exit("Please give ALL the necessairy arguments. Exiting...")

args = 0

if "-d" in sys.argv:
    input_filename = sys.argv[sys.argv.index("-d") + 1]
    args += 1
if "-q" in sys.argv: 
    test_filename = sys.argv[sys.argv.index("-t") + 1]
    args += 1
if "-l1" in sys.argv:
    input_label_filename = sys.argv[sys.argv.index("-dl") + 1]
    args += 1
if "-l2" in sys.argv: 
    test_label_filename = sys.argv[sys.argv.index("-tl") + 1]
    args += 1     
if "-o" in sys.argv: 
    outputfile = sys.argv[sys.argv.index("-o") + 1]
    args += 1 

if args != 5:
    exit("Please give ALL the necessairy arguments. Exiting...")

# SciPy keeps showing a warning that doesn't affect the 
# solution, so it is better to ignore this warning.
warnings.filterwarnings('ignore')

input_num_images = 0
test_num_images = 0

train_data = read_input(input_filename,input_num_images)
test_data = read_input(test_filename,test_num_images)

train_labels = read_labels(input_label_filename,input_num_images)
test_labels = read_labels(test_label_filename, test_num_images)

print(train_data.shape)
print(train_data.dtype)

image_dimension = train_data.shape[1]


# define the cluster size in pixels
# example : cluster_size = 7 means 7x7 sized cluster
# 7x7 = 16 clusters for one image
cluster_size = 7
input_signatures = []
query_signatures = []

input_num_images = train_data.shape[0]
query_num_images = test_data.shape[0]

# Partition the 600 input images into clusters and compute the centroid and weight
for i in range(input_num_images // 100):
  sig = compute_signature(train_data[i],cluster_size)
  input_signatures.append(sig)

# Partition 100 of the query images into clusters and compute the centroid and weight
for i in range(query_num_images // 100):
  sig = compute_signature(test_data[i],cluster_size)
  query_signatures.append(sig)

#print(input_signatures[0])
#print(len(query_signatures))



input_length = len(input_signatures[0])
query_length = len(query_signatures[0])

# All distances are the same for all clusters
distance_matrix = compute_distances(train_data[0],test_data[0],input_signatures[0],query_signatures[0])

print(distance_matrix)
print(len(distance_matrix))

  # for a 7x7 cluster, dimensions are 16^2 = 256
emd_dimensions = len(query_signatures[0])*len(input_signatures[0])
aub = []

# this is for the first constraint of sums
for i in range(input_length):
  arr = np.zeros(emd_dimensions,dtype=int)
  arr[input_length*i:input_length*i + input_length ] = 1
  aub.append(list(arr))  

# this is for the second constraint of sums
for i in range(input_length):
  arr = np.zeros(emd_dimensions,dtype=int)
  for j in np.arange(0,emd_dimensions,input_length):
    arr[i+j] = 1
    #print(i+j)
  aub.append(list(arr))


#print(len(aub))
#print(aub)  

#_=scipy.seterr(all='ignore')

total_score = 0

for i in range(len(query_signatures)):
  obj_functions = {}
  start_time = time.time()
  score = 0
  print("Query number %d" %test_labels[i])
  #show_image(test_data[i])

  for j in range(len(input_signatures)):
    obj_value = compute_EMD(train_data[j],test_data[i],input_signatures[j],query_signatures[i],aub,distance_matrix)
    obj_functions[j] = obj_value

  #print(obj_functions)
  elapsed_time = time.time() - start_time
  print("Query finished in %d seconds" %elapsed_time)  
  sorted_obj_functions = {l: v for l, v in sorted(obj_functions.items(), key=lambda item: item[1])}
  #print(sorted_obj_functions)
  for k in range(10):
    key = list(sorted_obj_functions.keys())
    neighbor = train_labels[key[k]]

    if neighbor == test_labels[i]:
      score+=1
      print("Neighbor number %d for Query number %d is number %d - SUCCESS" %(k+1,test_labels[i],neighbor))
    else:
      print("Neighbor number %d for Query number %d is number %d - FAILED" %(k+1,test_labels[i],neighbor)) 
    #show_image(train_data[key[k]])
  print("Scored %d / 10 \n" %score)
  total_score += score

print("Total score : %d" %total_score)

f = open("exact_results_reduced_10Ν.txt", "r")

Lines = f.readlines()

query_count = 0
counter = 0
score = 0
total_score_NN = 0

for i in range(0,len(Lines),10):

  query_label = test_labels[query_count]

  for j in range(10):
    value_list = Lines[i+j].split()
    neighbor = train_labels[int(value_list[0])-1]

    if neighbor == query_label:
      score+=1
      print("Neighbor number %d for Query number %d is number %d - SUCCESS" %(j+1,query_label,neighbor))
    else:
      print("Neighbor number %d for Query number %d is number %d - FAILED" %(j+1,query_label,neighbor)) 

  print("Scored %d / 10 \n" %score)
  total_score_NN += score

  score = 0
  query_count += 1

  if query_count == 100:
    break

print("Total score : %d" %total_score_NN)
f.close()

f = open("output.txt", "w")
f.write("Average Correct Search Results EMD: %d / 1000 \n" %total_score)
f.write("Average Correct Search Results MANHATTAN: %d / 1000 " %total_score_NN)
f.close()


