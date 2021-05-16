#!/usr/bin/env python
# coding: utf-8

# In[60]:


import pandas as pd
import numpy as np
from math import radians, cos, sin, asin, sqrt
import matplotlib.pyplot as mtp 
import os


# In[2]:


#Indian cities dataset consist cities and their latitude and longitude
csv_path='Indian Cities Database.csv'
df=pd.read_csv(csv_path)


# In[3]:


df.index


# In[4]:


data=[]
lat_all=[]
long_all=[]
for ind in df.index:
    city=df['City'][ind]
    lat=df['Lat'][ind]
    long=df['Long'][ind]
    lat_all.append(lat)
    long_all.append(long)
    temp_tuple=(city,(lat,long))
    data.append(temp_tuple)
data


# In[5]:


fig=mtp.figure()
mtp.scatter(long_all, lat_all)  
mtp.title('latitude and longitude of cities')  
mtp.xlabel('longitude')  
mtp.ylabel('latitude')  
mtp.show() 
fig.clear()
mtp.close(fig)


# In[6]:


def earth_dist(lat1,long1,lat2,long2):
    #converting decimal degree to radians
    lat1,long1,lat2,long2=map(radians,[lat1,long1,lat2,long2])
    
    # Using Haversine formula for distance calculation
    dlat=lat2-lat1
    dlong=long2-long1
    a=sin(dlat/2)**2 + cos(lat1) * cos(lat2) * sin(dlong/2)**2
    c=2 * asin(sqrt(a))
    r=6371 #radius of earth in km
    return c*r


# In[7]:


dist_list=[]
i=0
j=0
for i in range(len(data)):
    tup1=data[i][1]
    dist_all=[]
    for j in range(len(data)):
        tup2=data[j][1]
        dist=earth_dist(tup1[0],tup1[1],tup2[0],tup2[1])
        dist_all.append(dist)
        j+=1
    dist_list.append(dist_all)
    i+=1
dist_matrix=np.array(dist_list)
dist_matrix


# In[8]:


import kmedoids


# In[9]:


max_k=dist_matrix.shape[0]//3
def WCSS(dist):
    wcss_list=[]
    for n in range(1,max_k):
        M, C = kmedoids.kMedoids(dist, n)
        wcss=0
        i=0
        for i in range(n):
            cluster_wcss=0
            for j in C[i]:
                cluster_wcss+=(dist[j][M[i]])**2
            wcss+=cluster_wcss
            i+=1
        wcss_list.append(wcss)
        n+=1
    return wcss_list
wcss_list=WCSS(dist_matrix)
print(wcss_list)


# In[10]:


fig=mtp.figure()
mtp.plot(range(1, max_k), wcss_list)  
mtp.title('The Elobw Method Graph')  
mtp.xlabel('Number of clusters(k)')  
mtp.ylabel('wcss_list')  
mtp.show()
#fig.clear()
#mtp.close(fig)


# In[11]:


M, C = kmedoids.kMedoids(dist_matrix, 10)
print(M)
print(C)


# In[12]:


label=[0]*len(data)
long_cluster=[]
lat_cluster=[]
cluster_data=[]
cluster_name=[]
z=1
for key,item in C.items():
    long_cluster.append(list(map(lambda i:long_all[i],item)))
    lat_cluster.append(list(map(lambda i:lat_all[i],item)))
    k=1
    c_name="cluster"+str(z)
    z+=1
    cluster_name.append(c_name)
    single_cluster=[]
    for i,j in zip(long_cluster[key],lat_cluster[key]):
        single_cluster.append((str(k), str(i), str(j)))
        k+=1
    cluster_data.append(single_cluster)
    for i in item:
        label[i]=key
colors=['blue','orange','green','red','purple','brown','grey','pink','olive','cyan']
print(label)


# In[14]:


fig, ((ax0, ax1, ax2, ax3, ax4), (ax5, ax6, ax7, ax8, ax9)) = mtp.subplots(2, 5)
#fig.suptitle('Individual Clusters')
ax_list=[ax0,ax1,ax2,ax3,ax4,ax5,ax6,ax7,ax8,ax9]
i=0
for ax in fig.get_axes():
    ax_list[i].scatter(long_cluster[i], lat_cluster[i],c=colors[i])
    ax_list[i].set_title('Cluster'+str(i+1))
    ax.label_outer()
    i+=1


# In[15]:


fig.clear()
mtp.close(fig)
fig=mtp.figure()
mtp.scatter(long_all, lat_all, c=list(map(lambda i:colors[i],label)))  
mtp.title('latitude and longitude of cities')  
mtp.xlabel('longitude')  
mtp.ylabel('latitude')  
mtp.show()
fig.clear()
mtp.close(fig)


# In[16]:


idc=list(range(1,len(M)+1))
idc=list(map(lambda i:str(i),idc))
long_center=list((map(lambda i:str(long_all[i]),M)))
lat_center=list((map(lambda i:str(lat_all[i]),M)))
coord_center=list(zip(idc,long_center,lat_center))
print(coord_center)
cluster_name.append("CentroidCluster")
cluster_data.append(coord_center)


# In[17]:


#del df
#del wcss_list
#del M
import gc
gc.collect()
from julia import Main as jl
jl.include("Firefly.jl")


# In[18]:


solution=jl.main(cluster_data,cluster_name)


# In[19]:


print(solution)


# In[73]:


city_path_id=list(map(lambda i:solution[i][0],range(0,len(solution))))
new_list=[]
for i in range(0,len(city_path_id)):
    temp_list=[]
    for j in range(0,len(city_path_id[i])):
        temp_list.append(int(city_path_id[i][j])-1)
    new_list.append(temp_list)
city_path_id=new_list
city_paths=[]
n=0
for path_id in city_path_id:
    if(n==len(M)):
        break
    city_id=list(map(lambda i:C[n][i],path_id))
    city_names=list(map(lambda i:data[i][0],city_id))
    city_paths.append(city_names)
    n+=1
cluster_id=city_path_id[-1]
city_id=list(map(lambda i:M[i],cluster_id))
centroid_cities_path=list(map(lambda i:data[i][0],city_id))


# In[72]:


final_path=[]
final_path_latitude=[]
final_path_longitude=[]
label_final_path=[]
cluster_path_lat=[]
cluster_path_long=[]
for i in city_path_id[-1]:
    cluster_path_lat.append(solution[i][2])
    cluster_path_long.append(solution[i][1])
    for t in range(0,len(solution[i][0])):
        sol_long_list=solution[i][1]
        sol_lat_list=solution[i][2]
        final_path_longitude.append(sol_long_list[t])
        final_path_latitude.append(sol_lat_list[t])
        label_final_path.append(i)
    for k in city_paths[i]:
        final_path.append(k)
print(final_path)


# In[71]:


from IPython.display import Image
graphs_path=os.getcwd()+"\\graphs"
path_images=[]
iteration_images=[]
for filename in os.listdir(graphs_path):
    image_name=os.path.join(graphs_path, filename)
    if filename.find("path") !=-1:
        path_images.append(Image(image_name))
    elif filename.find("iteration") !=-1:
        iteration_images.append(Image(image_name))
print("Individual cluster path images:\n\n")
for img in path_images:
    display(img)
print("Iteration cost images:\n\n")
for img in iteration_images:
    display(img)


# In[59]:


for i in range(len(cluster_path_lat)):
    if(i!=len(cluster_path_lat)-1):
            x=[cluster_path_long[i][-1], cluster_path_long[i+1][0]]
            y=[cluster_path_lat[i][-1], cluster_path_lat[i+1][0]]
            mtp.plot(x, y,c='black')
    mtp.plot(cluster_path_long[i], cluster_path_lat[i],c=colors[i])
    
mtp.title('path between cities')  
mtp.xlabel('longitude')  
mtp.ylabel('latitude')  

