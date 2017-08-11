# -*- coding: utf-8 -*-
"""
Created on Tue Jul 11 12:58:14 2017

@author: Charlotte
"""

import numpy as np
import networkx as nx
import random
import scipy.sparse as sp
import matplotlib.pyplot as plt
import matplotlib as mpl
import scipy.io as sio
import pandas as pd
import h5py
import json
#.mat file is a adjacency matrix with the actual weights given for each element
#for matlab files saved in version 7.3, use the h5py library to read the file
f=h5py.File('DisKNN.mat')
#the adjacency matrix is stored under the variable PATH
x=f["PATH"]
adjmat_disease=x

#if the file containing the matlab matrix is not saved in version 7.3, use sio.loadmat
# # mat_contents=sio.loadmat('TFgraph.mat')
# # adjmat_disease=mat_contents['PATH']

#example graph to test layout on
#provides a scale free graph
G=nx.erdos_renyi_graph(50,0.1)


# send in graph with number of iterations to run, which
#method you want to use for repulsion, the strength of repulsion you want,
#and the strength of each centrality(betweeness or degree centrality) you want to
#count for the type of repulsion used
#for example, if centrality fraction is set to 0.001, that means that 99.9% of the effect
#of repulsion is due to betweenness centrality of each node

def layout(G,iter,repulstype,repulsfactor,centralityfraction):
    #given a graph object or adjacency matrix convert to adj mat
    adjmat=nx.adjacency_matrix(G)
    adjmat=adjmat.todense()
    a,b=adjmat.shape
    adjmat=np.asarray((adjmat),dtype='float64')
    iterations=iter
    #set initial positions of nodes randomly
    #use a seed value to get repeatable results
    pos_array=np.random.RandomState(42)
    pos=pos_array.rand(a,2)
    pos=pos.astype(adjmat.dtype)
    #k is the spring constant 
    k=np.sqrt(1.0/a)
    t = max(max(pos.T[0]) - min(pos.T[0]), max(pos.T[1]) - min(pos.T[1]))*0.1
#if we want the 'temperature' of the system to linearly decrease, use commented dt
#    dt=t/float(iterations+1)
    dt=0.9   
    #delta will become the x and y distance between every node and every other node,
   # giving a node x node x 2 size matrix
        delta = np.zeros((pos.shape[0],pos.shape[0],pos.shape[1]),dtype=adjmat.dtype)

    #betweenness and degree are arrays giving the specific centrality measure for each node
    #these are made into matrices, where the centrality factor between any two nodes
    # is equal to some linear equation of the product of their centrality values
    #this means that nodes that both have high values of centrality will repulse each other more
    #the centrality fraction variable provides the factor of combination for the two centrality products
    #to get the final relational matrix in centralitymat
    degree1=nx.degree_centrality(G)
    betweenness=nx.betweenness_centrality(G)
    betweenness=np.array(betweenness.values(),dtype='float64')*(1-centralityfraction)
    degree=np.array(degree1.values(),dtype='float64')*centralityfraction
    degree2=np.array(degree1.values(),dtype='float64')
    degreemat2=[[i*j for i in degree2]for j in degree2]
    degreemat=[[i*j/centralityfraction for i in degree]for j in degree]
    betweenmat=[[i*j/(1-centralityfraction) for i in betweenness] for j in betweenness]
    centralitymat=[[x+y for x,y in zip(degreemat[i],betweenmat[i])]for i in range(len(betweenmat))]


#the energy after each iteration will be tracked in an array to see how the energy changes over time
#this way we can adjust the temperature of the system based on how much the energy has changed over
#multiple iterations
    energy=[]
    progress=0
    for iteration in range(iterations):
        #find the displacement matrix between everynode and everyother node
        for i in range(pos.shape[1]):
            #each 2d array is for one coordinate
            #in order to get an array of displacements between one node and every other node
            delta[:,:,i]=pos[:,i,None]-pos[:,i]
        delta=np.array(delta,dtype=np.float64)
        #turn the displacement matrix into a distance matrix 
        #if the distance between two nodes is less than 0.01, reset the min value
        #to 0.01 so 
        distance=np.sqrt((delta**2).sum(axis=-1))
        distance=np.where(distance<0.01,0.01,distance)


        #repulsion function, only based on node to node interaction
        if repulstype =='edgespring1':
            edgeforce=degree + betweenness
            edgeforce=edgeforce/distance**2
            springrepuls=k*k/distance**2
            edgerepuls=edgeforce/(repulsfactor)
            force1=springrepuls+edgerepuls
            force1=force1/repulsfactor

        if repulstype == 'edgespring2':           
            force_edge=centralitymat
            force_edge=force_edge/(distance**2)
            force1=force_edge/repulsfactor/10+(k*k/distance**2)/repulsfactor
       
        #this layout has been the most successful
        if repulstype == 'trial':
            force_edge=centralitymat
            force_edge=force_edge/(distance**2)
            force1=force_edge/repulsfactor

        #follows linlog documentation for edge repulsion
        if repulstype == 'linlog':
            force_edge=centralitymat
            force_edge=force_edge*np.log(distance)
            force1=-1*force_edge

        #basic fructerman reihngold repulsion
        if repulstype == 'noedge':
            force1=k*k/distance**2
            force1*=1/repulsfactor

        # basic repulsion from spring layout
        if repulstype == 'spring':
            force1=k*k/distance
            force1*=1/repulsfactor


        #attraction function
        #basic fruchterman reihngold formula, only affects node-node interaction that have a nonzero 
        #relationship
        force2=adjmat*distance/k

        #gives the total force contribution from each node to every other node
        forces=force1-force2

        #displacement1 gives a 2column array for every node, applies the force along
        #the radial direction of the other nodes
        displacement1=np.transpose(np.transpose(delta)*(forces))
        #the forces from each node are summed together to get one final force vector
        #for each node
        displacement=displacement1.sum(axis=1)
        #length is the magnitude of the displacement vector
        length=np.sqrt((displacement**2).sum(axis=1))

        # following eq in http://yifanhu.net/PUB/graph_draw_small.pdf page 5 for energy
        energy.append(length.sum(axis=0))
        length=np.where(length<0.01,0.01,length)
        #gives a unit vector for each node in the position that they should move
        #this unit vector is multiplied by the temperature of the system
        #to get the final displacement for every node
        delta_pos=np.transpose(np.transpose(displacement)*t/length)
        pos+=delta_pos

        #calculate the energy of the system, if the energy has continuously decreased 5 times
        #increase the temperature of the system so that more movement is possible
        #do this once every 5 times, if energy increases instead of decreases reset the progress counter
        if (iteration>0 and energy[iteration]<energy[iteration-1]):
            progress+=1
            if progress>=5:
                progress=0
                t=t/dt
        else: 
            progress=0
            t*=dt                            
#if temperature should decrease linearly and unimpacted by stem energy, use t-=dt
#        t-=dt

#if there are no fixed nodes in the layout, rescale to fit in [0,1]
        pos=_rescale_layout(pos)
        #if delta_pos<tolerance
    

    plt.plot(energy[5:])
    return pos
    # return  coords dict with id original as key 
        
def _rescale_layout(pos,scale=1.):
    maxlim=0
    for i in range(pos.shape[1]):
        pos[:,i]-= pos[:,i].min()
        maxlim=max(maxlim,pos[:,i].max())
    if maxlim>0:
        for i in range(pos.shape[1]):
            pos[:,i]*=scale/maxlim
    return pos

#call to method
L2=layout(graph_disease,50,'trial',1,.001 )
#take the position array and input into format similar to cyjs
jsondata={'elements':{'nodes':[]}}
counter=1
for rows in L2:
    jsondata['elements']['nodes'].append({'data':{'name':counter},'position':{'x':rows[0],'y':rows[1]}})
    counter+=1
json1=json.dumps(jsondata)
with open('DisKnn(0.001).txt','w') as outfile:
    json.dump(jsondata,outfile)


#display graph with given position array
plt.figure('TF')
nx.draw_networkx(graph_disease,pos=L2)

# #
plt.show()  