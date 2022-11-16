import numpy as np
import operations as op
import cubic_lattice as cl

'''
SUMMARY:
This code takes a list of nodes that represent points on a 3D lattice. 

First part of this script identifies and labels separate networks. Nodes in a network
are all separated by SOME degree, i.e from one nucleus, there needs edges that can be traced 
to any of the network members. What determines if there is an edge between two nodes? 
In this model (Edwards-Anderson spin glass model), there is an edge if two nodes are only one point away from 
each other in any direction -- because this is 3D, each node can have up to 6 edges.  

Second part calculates whether a network percolates (has a path from one face to the opposite 
face). It then visualizes the network in a 3D lattice, as well as the path. 
@author: mrahm32
'''
N = 1000
''' do not change N if running with test data'''

class networkInfo:
    def __init__(self, spinIndices, N):
        self.spinIndices = spinIndices
        self.N = N
        self.L = np.round(N**(1/3.0))
        self.spinsNum = len(spinIndices)
        self.str_coordMatrix = []
        self.build_matrix(spinIndices,N)      
        self.label_networks()  
        self.is_network_spanning()
        
    def build_matrix(self, spinIndices, N):
        ''' Parameters: N = int, system size. spinIndices = 1D array of integers
        takes the data (list of N 1D indices) and converts 
        into array of string(3D coordinates)''' 
        
        for count, spinIndex in enumerate(spinIndices):
            self.coord = op.index_conversion_1D_3D(N, spinIndex)
            coord_array = np.array([int(self.coord[0]), int(self.coord[1]), int(self.coord[2])])
            #self.indexMatrix[count] = spinIndex
            self.str_coordMatrix.append(str(coord_array))
       

    def generate_neighbors(self, coord):
        '''Parameters: coord : (i,j,k) is a string. 
        Function: Uses a "neighborhood kernel" to calculate all 6 possible neighboring
        nodes connected by just 1 edge.'''
                
        r_Neighbors = np.identity(3)
        l_Neighbors = (-1)*r_Neighbors
        allNeighbors = np.vstack((r_Neighbors, l_Neighbors))
        int_coord = np.array(op.str_2_int(coord))
        localNeighbors = np.tile(int_coord, (6,1)) + allNeighbors
                
        localNeighbors[localNeighbors == -1]= self.L-1
        localNeighbors[localNeighbors == self.L] = 0

        localNeighbors_str = op.int_2_str(localNeighbors)
        unoccupiedNeighbors = [] 
        
        '''Neighbors that are not part of the real dataset or are already cached
        are deleted.'''
        
        for neighbor in localNeighbors_str:
            if neighbor in self.cached_coords or neighbor not in self.str_coordMatrix:
                unoccupiedNeighbors.append(neighbor)
        for neighbor in unoccupiedNeighbors:
            localNeighbors_str.remove(neighbor)

        if len(localNeighbors_str) > 0:
            self.neighbors[coord]= self.neighbors[coord] + localNeighbors_str
            self.node_queue.append(coord)

    def add_node_edge(self, label, coord, node):
        '''node is added to a graph structure '''
        
        try:
            self.networks_graph[label][coord].append(node)
        except:
            self.networks_graph[label][coord] = []
            self.networks_graph[label][coord].append(node)

    def add_coord(self, label, node, coord):
        '''node is added to network and also list of 3D coords'''
        
        
        vec_diff= np.array(op.str_2_int(node)) - np.array(op.str_2_int(coord))        
        node1DIndex = op.index_conversion_3D_1D(self.N, op.str_2_int(node))
        coord1DIndex = op.index_conversion_3D_1D(self.N, op.str_2_int(coord))
        
        ''' nodes are ignored if they are related through periodic boundary 
        condition --- this is because eventually this network will be used
        to find a path from boundary to boundary'''
        
        if np.abs(sum(vec_diff)) == 1:
            self.networks_3Dcoord[label].append(op.str_2_int(node))
           
            self.add_node_edge(label, coord1DIndex, node1DIndex)
            self.add_node_edge(label, node1DIndex, coord1DIndex)
        
        try:
            self.str_coordMatrix.remove(node)
            self.cached_coords.append(node)
        except:
               pass
           
            
    def is_network_spanning(self):
        ''' lays out a structure for how to evaluate whether a networks spans from
        boundary to boundary by iterating through all the cluster labels, evaluating via
        the "find path" function and terminating if any path is found'''
        
        self.perc = 0
         
        for network_label in self.networks_graph.keys():
            network_coords = np.array(self.networks_3Dcoord[network_label])
            if len(network_coords) > 1:
                for axs_indx in range(3):
                    self.find_path(network_coords, network_label, axs_indx)
       
        if self.perc > 0:
            self.path = []
            for i in self.visited:
                coord = op.index_conversion_1D_3D(self.N, i)
                self.path.append(coord)
           
            cl.make_vis(self.networks_3Dcoord, self.L, self.path)
            

    def label_networks(self):
        ''' * SIMILAR TO BURNING ALGORITHM *
            Nodes in the same network can be connected by any degree of separation. So first pick
            a random coord, look for whether it has any neirest neighbors, and add to network, and for 
            THOSE neighbors, look for THEIR nearest neighbors and add to the same network... terminate
            until there are no more nearest neighbors. And then start with a disconnected node ...and repeat.
        
            Objects below:
            networks_3Dcoord = a dictionary with network labels, and 3D coordinates the member
            to preserve geometric information in the lattice
            networks_graph = dictionary representation of the graph structure where each key is a node(vertex), 
            and dict values are other nodes separated by 1 edge
            neighbors = nearest neighbors of a given node 
            node_queue = list of nodes to be assessed 
            cached_coords = coordinates that have been added to the graph, and we do not want to revist
            or else we would loop forever... 
            '''
            
            
        self.networks_3Dcoord = {0:[]}
        self.networks_graph = {}
        self.neighbors ={}
        self.node_queue=[]
        self.cached_coords = []
        
        while len(self.str_coordMatrix) != 0:
            for coord in self.str_coordMatrix:
                self.neighbors[coord] = []
                self.generate_neighbors(coord)
                
                if len(self.neighbors[coord]) !=0:
                    label = np.max(np.array(list(self.networks_3Dcoord.keys()))) + 1
                    self.networks_3Dcoord[label] = []
                    self.networks_graph[label] = {}
                    self.add_coord(label, coord, coord)

                    while len(set(self.node_queue)) != 0:
                        for point in self.node_queue:                        
                            all_neighbors = self.neighbors[point]                            
                            if len(all_neighbors)> 0:
                                for count, neighbor in enumerate(all_neighbors):                               
                                    self.neighbors[neighbor] = []
                                    self.add_coord(label, neighbor, point)
                                    self.generate_neighbors(neighbor)                                                                      
                                    if count == len(all_neighbors) - 1:
                                        self.node_queue.remove(point)
                            else:
                                self.node_queue.remove(point)
                else:
                    
                    self.cached_coords.append(coord)
                    self.str_coordMatrix.remove(coord)
        
        
    def find_path(self, network_coords, label, axs_indx):
        
        ''' Parameters: cluster_coords == 3D coordinates of a network
        label = label of the network
        axs_indx = since we want boundary to boundary spanning, we can 
        evaluate each plane to see if there are points on opposite sides
        by transposing coordinate matrix, and looking at i,j,k (0,1,2)
        
        
        Check to see whether any network percolates -- I.E Does it span completely
        from boundary to boundary? If I take one point on the surface, can I trace an
        exact path of points to the bottom? 
        
        Since there can be a lot of points, the first technique is to eliminate many of
        them, by ONLY looking for points that are on opposite walls -- that is where 
        networks_3D comes in handy!!!
        
        Then, for each pair of points on the opposite sides, check if there is a path
        between them. '''
           
        plane = network_coords.T[axs_indx] # transposing coords to get information along one axis     
        if np.max(plane) - np.min(plane) >= self.L-1:
            if len(np.unique(plane))>= self.L:           
                ## look for all points that are touching opposite walls ##
                all_initialINDX = np.where(plane == 0)[0]
                all_finalINDX = np.where(plane == self.L-1)[0]                
                for initialINDX in all_initialINDX:
                    for finalINDX in all_finalINDX:
                        if self.perc == 0:
                            initCoord = network_coords[initialINDX]
                            finalCoord = network_coords[finalINDX]                        
                            ## convert coordinates to 1D indices to navigate self.networks_graph##
                            init1D = op.index_conversion_3D_1D(self.N, initCoord)
                            final1D = op.index_conversion_3D_1D(self.N, finalCoord)
                                
                            ''' elimination is done -- now the real work begins.
                            create a branching queue to store nodes we have yet to 
                            branch from. from the initial node, explore one connected
                            node, from that connected node...explore another until
                            final point is reached. Having a branch quue will make
                            backtracking possible if there is a better route'''
                            
                            self.visited = []
                            branch_queue= []                           
                            branch_queue.append(init1D)  
                                                   
                            while branch_queue:
                                branch_queue = list(set(branch_queue))
                                
                                init1D = branch_queue.pop(0)
                                self.visited.append(init1D)
                               
                                if init1D == final1D:
                                    self.perc = 1  
                                    # terminate once a path is found. that's all we need.
                                    break 
                                else:
                                    for point in self.networks_graph[label][init1D]:                                                                     
                                        if point not in self.visited:                                             
                                            branch_queue.append(point)
                                                                       
    


data = np.loadtxt('test_data.txt')
network = networkInfo(data, N)
        

