""" Implementation of the methods described in https://kevgildea.github.io/blog/Kinematic-Chain-Mapping/"""

import numpy as np
import matplotlib.pyplot as plt
import mpl_toolkits.mplot3d as plt3d
import copy
import networkx as nx

# functions
# ----------------------------------------------------------------------------------------------------------------
def jnt_path(graph, start, end, path=[]):
    path = path + [start]
    if start == end:
        return path
    if not start in graph:
        return None
    for node in graph[start]:
        if node not in path:
            newpath = jnt_path(graph, node, end, path)
            if newpath: return newpath
    return None

def takeSecond(elem):
    return elem[1]

def parent_child(graph):
    parent_child_pairs = []
    for i in range(max(graph[max(graph)])+1):
        for j in range(max(graph[max(graph)])+1):
            pair = jnt_path(dir_graph, i, j)
            if pair != None: 
                if len(pair)==2:
                    parent_child_pairs.append(pair)
    parent_child_pairs.sort(key=takeSecond)
    return parent_child_pairs

def Vector_mapping_Euler_Axis_Space(vec1, vec2): 
    """ Calculate all rotation matrices that map vector a to align with vector b"""
    a, b = (vec1 / np.linalg.norm(vec1)), (vec2 / np.linalg.norm(vec2))
    p1 = np.array([0,0,0])
    p2 = b / np.linalg.norm(b) + (a / np.linalg.norm(a))
    p3 = np.cross(a, b) / np.linalg.norm(np.cross(a, b))
    n=[]
    # create a list of candidate Euler axes (discritised to 1 degree)
    for i in range(0,360,1):
        Φ=np.radians(i)
        x_Φ=p1[0]+np.cos(Φ)*(p2[0]-p1[0])+np.sin(Φ)*(p3[0]-p1[0])
        y_Φ=p1[1]+np.cos(Φ)*(p2[1]-p1[1])+np.sin(Φ)*(p3[1]-p1[1])
        z_Φ=p1[2]+np.cos(Φ)*(p2[2]-p1[2])+np.sin(Φ)*(p3[2]-p1[2])
        n_Φ=[x_Φ,y_Φ,z_Φ]
        n.append(n_Φ / np.linalg.norm(n_Φ))
    # project vectors to form a cone around the Euler axis, and determine required angle for mapping
    rotation_matrices=[]
    Euler_axes=[]
    Euler_angles=[]
    for i in range(len(n)):
        Euler_axes.append(n[i])
        a_Φ, b_Φ = (np.cross(a,n[i]) / np.linalg.norm(np.cross(a,n[i]))), (np.cross(b,n[i]) / np.linalg.norm(np.cross(b,n[i])))
        θ = np.arccos(np.dot(a_Φ,b_Φ))
        θ = θ*np.sign(np.dot(n[i], np.cross(a_Φ,b_Φ)))
        Euler_angles.append(θ)
        if θ != θ: # if θ is NaN
            rotation_matrices.append(np.array([[ 1, 0, 0],
                                               [ 0, 1, 0],
                                               [ 0, 0, 1]]))
        else:
            rotation_matrices.append(np.array([[n[i][0]**2+(n[i][1]**2+n[i][2]**2)*(np.cos(θ)),n[i][0]*n[i][1]*(1-np.cos(θ))-n[i][2]*np.sin(θ),n[i][0]*n[i][2]*(1-np.cos(θ))+n[i][1]*np.sin(θ)],
                                        [n[i][0]*n[i][1]*(1-np.cos(θ))+n[i][2]*np.sin(θ),n[i][1]**2+(n[i][0]**2+n[i][2]**2)*(np.cos(θ)),n[i][1]*n[i][2]*(1-np.cos(θ))-n[i][0]*np.sin(θ)],
                                        [n[i][0]*n[i][2]*(1-np.cos(θ))-n[i][1]*np.sin(θ),n[i][1]*n[i][2]*(1-np.cos(θ))+n[i][0]*np.sin(θ),n[i][2]**2+(n[i][0]**2+n[i][1]**2)*(np.cos(θ))]]))
    
    return Euler_axes, Euler_angles, rotation_matrices

def FK_MDH(chain,dir_graph): 
    """ perform forward kinematics to to convert joint orientations and position vectors into the global coordinate system"""
    oris=[]
    poss=[]
    for i in range (len(chain)):
        path = jnt_path(dir_graph, 0, i)
        T=np.array([[ 1, 0, 0, 0],
                    [ 0, 1, 0, 0],
                    [ 0, 0, 1, 0],
                    [ 0, 0, 0, 1]])
        for j in path:
            T= T @ chain[j][1]
            pos=np.array([T[0][3],T[1][3],T[2][3]])
            ori=np.array([[T[0][0],T[0][1],T[0][2]],[T[1][0],T[1][1],T[1][2]],[T[2][0],T[2][1],T[2][2]]])
        oris.append(ori)
        poss.append(pos)
    return oris, poss

def Serial_open_chain_mapping(chain_a,chain_b,dir_graph):
    """ perform inverse kinematics to map a serial kinematic chain a (an open chain without forking) to chain b, which must have the same directed graph"""
    ori_a, pos_a = FK_MDH(chain_a,dir_graph)
    _, pos_b = FK_MDH(chain_b,dir_graph)

    pos_a= pos_a - np.array([chain_a[0][1][0][3],chain_a[0][1][1][3],chain_a[0][1][2][3]])
    pos_b= pos_b - np.array([chain_b[0][1][0][3],chain_b[0][1][1][3],chain_b[0][1][2][3]])
    
    #convert directed graph for use in nx
    dir_graph = nx.DiGraph(dir_graph)

    repos_a_EAS=[]
    reori_a_EAS=[]
    Euler_axis=[]
    Euler_axes=[]

    #I think I will need to create another function here for looping every solution in upchain joints to output the solution space for downchain joints
    for i in range(1,360,10): # discretised further here
        reori_a=copy.deepcopy(ori_a)
        repos_a=copy.deepcopy(pos_a)
        Euler_axis=[]
        for j in range(len(chain_a)):
            downchain = list(nx.nodes(nx.dfs_tree(dir_graph, j)))
            if j < len(chain_a)-1:
                ROTM= Vector_mapping_Euler_Axis_Space(repos_a[j+1]-repos_a[j], pos_b[j+1]-pos_b[j])[2][i]
                Euler_axis.append(Vector_mapping_Euler_Axis_Space(repos_a[j+1]-repos_a[j], pos_b[j+1]-pos_b[j])[0][i])
            for k in range(len(downchain)):
                if downchain[0] != downchain[-1]:
                    reori_a[downchain[k]]= ROTM @ reori_a[downchain[k]]
                if downchain[k] != downchain[-1]:
                    repos_a[downchain[k+1]]= repos_a[j] + ROTM @ (repos_a[downchain[k+1]]-repos_a[j])
            repos_a_EAS.append(repos_a)
            reori_a_EAS.append(reori_a)
            Euler_axes.append(Euler_axis)
    
    for i in range(len(repos_a_EAS)):
        repos_a_EAS[i]=repos_a_EAS[i]+np.array([chain_a[0][1][0][3],chain_a[0][1][1][3],chain_a[0][1][2][3]])


    return reori_a_EAS, repos_a_EAS, Euler_axes

def IK_complex_open_chain(chain_a, chain_b,dir_graph): 
    """ perform inverse kinematics to map a complex kinematic chain a (an open chain with forking) to chain b, which must have the same directed graph"""

    #TO DO!

    return 
# ----------------------------------------------------------------------------------------------------------------





# choose vectors a and b
# ----------------------------------------------------------------------------------------------------------------
a=[1,3,4]
b=[3,1,2]
# ----------------------------------------------------------------------------------------------------------------



# example choosing an initial local coordinate system with Y axis alligned with vector a
# ----------------------------------------------------------------------------------------------------------------
fig = plt.figure(0)
ax = fig.add_subplot(111, projection='3d')
fig.suptitle('Example choosing an initial local coordinate system with Y axis alligned with vector a', fontsize=12)

Euler_axes, Euler_angles, rotation_matrices = Vector_mapping_Euler_Axis_Space(a, b)

ax.quiver(0,0,0,a[0],a[1],a[2],color='black')
ax.text(a[0],a[1],a[2],  '%s' % (r'$\vec a$'), size=12, zorder=1,color='black')
ax.quiver(0,0,0,b[0],b[1],b[2],color='grey')
ax.text(b[0],b[1],b[2],  '%s' % (r'$\vec b$'), size=12, zorder=1,color='grey')

dirxglobal=[1,0,0]
diryglobal=[0,1,0]
dirzglobal=[0,0,1]

ax.quiver(-5,-5,-5,dirxglobal[0],dirxglobal[1],dirxglobal[2],color='r')
ax.quiver(-5,-5,-5,diryglobal[0],diryglobal[1],diryglobal[2],color='g')
ax.quiver(-5,-5,-5,dirzglobal[0],dirzglobal[1],dirzglobal[2],color='b')

_, _, Localcoordsys = Vector_mapping_Euler_Axis_Space(np.array([0,1,0]), a)
Localcoordsys=Localcoordsys[100]

dirxlocal=Localcoordsys @ [1,0,0]
dirylocal=Localcoordsys @ [0,1,0]
dirzlocal=Localcoordsys @ [0,0,1]

ax.quiver(0,0,0,dirxlocal[0],dirxlocal[1],dirxlocal[2],color='r')
ax.quiver(0,0,0,dirylocal[0],dirylocal[1],dirylocal[2],color='g')
ax.quiver(0,0,0,dirzlocal[0],dirzlocal[1],dirzlocal[2],color='b')

for i in range(0,len(rotation_matrices),2):#(0,len(rotation_matrices),2):
    a_mapped = rotation_matrices[i] @ a
    ax.quiver(0,0,0,a_mapped[0],a_mapped[1],a_mapped[2],color='black',linestyle='--',alpha=0.3)
    ax.text(a_mapped[0],a_mapped[1],a_mapped[2],  '%s' % (r'$\vec a_{mapped}$'), size=12, zorder=1,color='black',alpha=0.3)
    θ = Euler_angles[i]
    n = Euler_axes[i]

    ax.quiver(0,0,0,n[0],n[1],n[2],color='orange',linestyle='--', alpha=0.3)

    dirxlocal_new=rotation_matrices[i] @ dirxlocal
    dirylocal_new=rotation_matrices[i] @ dirylocal
    dirzlocal_new=rotation_matrices[i] @ dirzlocal

    ax.quiver(0,0,0,dirxlocal_new[0],dirxlocal_new[1],dirxlocal_new[2],color='r',linestyle='--', alpha=0.1)
    ax.quiver(0,0,0,dirylocal_new[0],dirylocal_new[1],dirylocal_new[2],color='g',linestyle='--', alpha=0.1)
    ax.quiver(0,0,0,dirzlocal_new[0],dirzlocal_new[1],dirzlocal_new[2],color='b',linestyle='--', alpha=0.1)

ax.set_xlim3d(-5,5)
ax.set_ylim3d(-5,5)
ax.set_zlim3d(-5,5)

ax.set_xlabel('Global X')
ax.set_ylabel('Global Y')
ax.set_zlabel('Global Z')

plt.show()
# ----------------------------------------------------------------------------------------------------------------



# Define serial chains a and b
# ----------------------------------------------------------------------------------------------------------------
# The chans are expressed using local orientations (relative to parent) as 3x3 rotation matrices, and local position offset vectors (position of child relative to parent)
chain_a=[]
chain_b=[]

chain_a.append(['jnt_a0', np.array([[ 1, 0, 0, 0],
                                    [ 0, 1, 0, 0.2],
                                    [ 0, 0, 1, 0],
                                    [ 0, 0, 0, 1]])])
chain_a.append(['jnt_a1', np.array([[ 0, 1, 0, 0],
                                    [ 0, 0, 1, 0.2],
                                    [ 1, 0, 0, 0],
                                    [ 0, 0, 0, 1]])])               
chain_a.append(['jnt_a2', np.array([[ 1, 0, 0, 0.5],
                                    [ 0, 0, 1, 0],
                                    [ 0, 1, 0, 0],
                                    [ 0, 0, 0, 1]])])
chain_a.append(['jnt_a3', np.array([[ 0, 0, 1, 0],
                                    [ 1, 0, 0, 0.2],
                                    [ 0, 1, 0, 0],
                                    [ 0, 0, 0, 1]])])
chain_a.append(['jnt_a4', np.array([[ 1, 0, 0, 0],
                                    [ 0, 0, 1, 0],
                                    [ 0, 1, 0, 0.2],
                                    [ 0, 0, 0, 1]])])

chain_b.append(['jnt_b0', np.array([[ 1, 0, 0, 0],
                                    [ 0, 1, 0, -0.2],
                                    [ 0, 0, 1, 0],
                                    [ 0, 0, 0, 1]])])
chain_b.append(['jnt_b1', np.array([[ 1, 0, 0, 0],
                                    [ 0, 1, 0, 0],
                                    [ 0, 0, 1, -0.2],
                                    [ 0, 0, 0, 1]])])             
chain_b.append(['jnt_b2', np.array([[ 1, 0, 0, 0],
                                    [ 0, 1, 0, -0.2],
                                    [ 0, 0, 1, -0.2],
                                    [ 0, 0, 0, 1]])])
chain_b.append(['jnt_b3', np.array([[ 1, 0, 0, 0.5],
                                    [ 0, 1, 0, 0],
                                    [ 0, 0, 1, 0],
                                    [ 0, 0, 0, 1]])])
chain_b.append(['jnt_b4', np.array([[ 1, 0, 0, 0],
                                    [ 0, 1, 0, -0.5],
                                    [ 0, 0, 1, 0],
                                    [ 0, 0, 0, 1]])])


# a directed graph with parent/child relationship knowledge is needed for constructing the kinematic chain
dir_graph = {0: [1],
             1: [2],
             2: [3],
             3: [4]}

#print(nx.incidence_matrix(nx.DiGraph(dir_graph), nodelist=None, edgelist=None, oriented=False, weight=None).toarray())

# find parent/child pairs from the directed graph
parent_child_pairs=parent_child(dir_graph)#[[0,1],[1,2],[2,3],[0,4],[4,5],[5,6],[0,7],[7,8],[8,9],[7,10],[9,10]] # REMOVE
# ----------------------------------------------------------------------------------------------------------------


# Calculate global orientations and positions for sinple chains a and b
# ----------------------------------------------------------------------------------------------------------------

ori_a, pos_a = FK_MDH(copy.deepcopy(chain_a),dir_graph)

ori_b, pos_b = FK_MDH(copy.deepcopy(chain_b),dir_graph)

# ----------------------------------------------------------------------------------------------------------------



# Plot serial chain a
# ----------------------------------------------------------------------------------------------------------------
fig = plt.figure(1)
ax = fig.add_subplot(111, projection='3d')
fig.suptitle('Serial open kinematic chain a, plotted using forward kinematics', fontsize=12)

dirxglobal=[1,0,0]
diryglobal=[0,1,0]
dirzglobal=[0,0,1]

ax.quiver(0,0,0,dirxglobal[0],dirxglobal[1],dirxglobal[2],length=0.05,color='r')
ax.quiver(0,0,0,diryglobal[0],diryglobal[1],diryglobal[2],length=0.05,color='g')
ax.quiver(0,0,0,dirzglobal[0],dirzglobal[1],dirzglobal[2],length=0.05,color='b')
#ax.text(-0.5,-0.5,-0.5,  '%s' % ('global'), size=8, zorder=1,color='dimgrey')

ax.text(0,0,0, 'ref. space', size=8, zorder=1,color='black')

x=[]
y=[]
z=[]
# plot chain a
for i in range(len(pos_a)):
    x.append(pos_a[i][0])
    y.append(pos_a[i][1])
    z.append(pos_a[i][2])
    
    dirxlocal=ori_a[i]@[1,0,0]
    dirylocal=ori_a[i]@[0,1,0]
    dirzlocal=ori_a[i]@[0,0,1]

    ax.quiver(x[i],y[i],z[i],dirxlocal[0],dirxlocal[1],dirxlocal[2],length=0.05,color='r')
    ax.quiver(x[i],y[i],z[i],dirylocal[0],dirylocal[1],dirylocal[2],length=0.05,color='g')
    ax.quiver(x[i],y[i],z[i],dirzlocal[0],dirzlocal[1],dirzlocal[2],length=0.05,color='b')
ax.scatter(x, y, z, c=['gray']*len(chain_a), marker='o')

for i in range(len(chain_a)):
    ax.scatter(x[i], y[i], z[i], c='black', marker='.')
    ax.text(x[i], y[i], z[i],  '%s' % (str(chain_a[i][0])), size=8, zorder=1,color='black')

for j in range(len(parent_child_pairs)):
    xs = [x[parent_child_pairs[j][0]],x[parent_child_pairs[j][1]]]
    ys = [y[parent_child_pairs[j][0]],y[parent_child_pairs[j][1]]]
    zs = [z[parent_child_pairs[j][0]],z[parent_child_pairs[j][1]]]
    
    line = plt3d.art3d.Line3D(xs, ys, zs, c='black')
    ax.add_line(line)

#ax.quiver(0,0,0,pos_a[0][0],pos_a[0][1],pos_a[0][2],color='gray',alpha=0.6,linestyle='--')
ax.add_line(plt3d.art3d.Line3D([0,pos_a[0][0]], [0,pos_a[0][1]], [0,pos_a[0][2]] ,c='gray',alpha=0.6,linestyle='--'))

ax.set_xlim3d(-2,2)
ax.set_ylim3d(-2,2)
ax.set_zlim3d(-2,2)

ax.set_xlabel('Global X')
ax.set_ylabel('Global Y')
ax.set_zlabel('Global Z')

plt.show()
# ----------------------------------------------------------------------------------------------------------------



# Plot serial chains a and b  --> Plot b without local joint orientations
# ----------------------------------------------------------------------------------------------------------------
fig = plt.figure(2)
ax = fig.add_subplot(111, projection='3d')
fig.suptitle('Serial open kinematic chains a & b', fontsize=12)

dirxglobal=[1,0,0]
diryglobal=[0,1,0]
dirzglobal=[0,0,1]

ax.quiver(0,0,0,dirxglobal[0],dirxglobal[1],dirxglobal[2],length=0.05,color='r')
ax.quiver(0,0,0,diryglobal[0],diryglobal[1],diryglobal[2],length=0.05,color='g')
ax.quiver(0,0,0,dirzglobal[0],dirzglobal[1],dirzglobal[2],length=0.05,color='b')
#ax.text(-0.5,-0.5,-0.5,  '%s' % ('global'), size=8, zorder=1,color='dimgrey')

ax.text(0,0,0, 'ref. space', size=8, zorder=1,color='black')

x=[]
y=[]
z=[]
# plot chain a
for i in range(len(pos_a)):
    x.append(pos_a[i][0])
    y.append(pos_a[i][1])
    z.append(pos_a[i][2])
    
    dirxlocal=ori_a[i]@[1,0,0]
    dirylocal=ori_a[i]@[0,1,0]
    dirzlocal=ori_a[i]@[0,0,1]

    ax.quiver(x[i],y[i],z[i],dirxlocal[0],dirxlocal[1],dirxlocal[2],length=0.05,color='r')
    ax.quiver(x[i],y[i],z[i],dirylocal[0],dirylocal[1],dirylocal[2],length=0.05,color='g')
    ax.quiver(x[i],y[i],z[i],dirzlocal[0],dirzlocal[1],dirzlocal[2],length=0.05,color='b')
ax.scatter(x, y, z, c=['gray']*len(chain_a), marker='o')

for i in range(len(chain_a)):
    ax.scatter(x[i], y[i], z[i], c='black', marker='.')
    ax.text(x[i], y[i], z[i],  '%s' % (str(chain_a[i][0])), size=8, zorder=1,color='black')

for j in range(len(parent_child_pairs)):
    xs = [x[parent_child_pairs[j][0]],x[parent_child_pairs[j][1]]]
    ys = [y[parent_child_pairs[j][0]],y[parent_child_pairs[j][1]]]
    zs = [z[parent_child_pairs[j][0]],z[parent_child_pairs[j][1]]]
    
    line = plt3d.art3d.Line3D(xs, ys, zs, c='black')
    ax.add_line(line)

#ax.quiver(0,0,0,pos_a[0][0],pos_a[0][1],pos_a[0][2],color='gray',alpha=0.6,linestyle='--')
ax.add_line(plt3d.art3d.Line3D([0,pos_a[0][0]], [0,pos_a[0][1]], [0,pos_a[0][2]] ,c='gray',alpha=0.6,linestyle='--'))

x=[]
y=[]
z=[]
# plot chain b
for i in range(len(pos_b)):
    x.append(pos_b[i][0])
    y.append(pos_b[i][1])
    z.append(pos_b[i][2])
    
    #dirxlocal=ori_b[i]@[1,0,0]
    #dirylocal=ori_b[i]@[0,1,0]
    #dirzlocal=ori_b[i]@[0,0,1]

    #ax.quiver(x[i],y[i],z[i],dirxlocal[0],dirxlocal[1],dirxlocal[2],length=0.05,color='r')
    #ax.quiver(x[i],y[i],z[i],dirylocal[0],dirylocal[1],dirylocal[2],length=0.05,color='g')
    #ax.quiver(x[i],y[i],z[i],dirzlocal[0],dirzlocal[1],dirzlocal[2],length=0.05,color='b')
ax.scatter(x, y, z, c=['black']*len(chain_b), marker='o')

for i in range(len(chain_b)):
    ax.scatter(x[i], y[i], z[i], c='black', marker='.')
    ax.text(x[i], y[i], z[i],  '%s' % (str(chain_b[i][0])), size=8, zorder=1,color='black')

for j in range(len(parent_child_pairs)):
    xs = [x[parent_child_pairs[j][0]],x[parent_child_pairs[j][1]]]
    ys = [y[parent_child_pairs[j][0]],y[parent_child_pairs[j][1]]]
    zs = [z[parent_child_pairs[j][0]],z[parent_child_pairs[j][1]]]
    
    line = plt3d.art3d.Line3D(xs, ys, zs, c='gray')
    ax.add_line(line)

#ax.quiver(0,0,0,pos_b[0][0],pos_b[0][1],pos_b[0][2],color='gray',alpha=0.6,linestyle='--')
ax.add_line(plt3d.art3d.Line3D([0,pos_b[0][0]], [0,pos_b[0][1]], [0,pos_b[0][2]] ,c='gray',alpha=0.6,linestyle='--'))

ax.set_xlim3d(-2,2)
ax.set_ylim3d(-2,2)
ax.set_zlim3d(-2,2)

ax.set_xlabel('Global X')
ax.set_ylabel('Global Y')
ax.set_zlabel('Global Z')

plt.show()
# ----------------------------------------------------------------------------------------------------------------



# Apply inverse kinematics to map serial chain a to chain b
# ----------------------------------------------------------------------------------------------------------------
reori_a, repos_a, Euler_axes = Serial_open_chain_mapping(copy.deepcopy(chain_a),copy.deepcopy(chain_b),dir_graph)
# ----------------------------------------------------------------------------------------------------------------



# Plot serial chain a repositioned
# ----------------------------------------------------------------------------------------------------------------
fig = plt.figure(3)
ax = fig.add_subplot(111, projection='3d')
fig.suptitle('Solution space for mapping kinematic chain a to chain b', fontsize=12)

dirxglobal=[1,0,0]
diryglobal=[0,1,0]
dirzglobal=[0,0,1]

ax.quiver(0,0,0,dirxglobal[0],dirxglobal[1],dirxglobal[2],length=0.05,color='r')
ax.quiver(0,0,0,diryglobal[0],diryglobal[1],diryglobal[2],length=0.05,color='g')
ax.quiver(0,0,0,dirzglobal[0],dirzglobal[1],dirzglobal[2],length=0.05,color='b')

ax.text(0,0,0, 'ref. space', size=8, zorder=1,color='black')

x=[]
y=[]
z=[]

# plot chain a for the EAS of the chain
for s in range(len(reori_a)):
    x=[]
    y=[]
    z=[]
    for i in range(len(reori_a[s])):
        x.append(repos_a[s][i][0])
        y.append(repos_a[s][i][1])
        z.append(repos_a[s][i][2])
        
        dirxlocal=reori_a[s][i]@[1,0,0]
        dirylocal=reori_a[s][i]@[0,1,0]
        dirzlocal=reori_a[s][i]@[0,0,1]

        ax.quiver(x[i],y[i],z[i],dirxlocal[0],dirxlocal[1],dirxlocal[2],length=0.05,color='r',linestyle='--', alpha=0.3)
        ax.quiver(x[i],y[i],z[i],dirylocal[0],dirylocal[1],dirylocal[2],length=0.05,color='g',linestyle='--', alpha=0.3)
        ax.quiver(x[i],y[i],z[i],dirzlocal[0],dirzlocal[1],dirzlocal[2],length=0.05,color='b',linestyle='--', alpha=0.3)

        # plot the Euler axis
        if i<len(Euler_axes[s]):
            n= Euler_axes[s][i]
            ax.quiver(x[i],y[i],z[i],n[0],n[1],n[2],length=0.05,color='orange')
        
ax.scatter(x, y, z, c=['gray']*len(chain_a), marker='o')

for i in range(len(chain_a)):
    ax.scatter(x[i], y[i], z[i], c='black', marker='.')
    ax.text(x[i], y[i], z[i],  '%s' % (str(chain_a[i][0])), size=8, zorder=1,color='black')

for j in range(len(parent_child_pairs)):
    xs = [x[parent_child_pairs[j][0]],x[parent_child_pairs[j][1]]]
    ys = [y[parent_child_pairs[j][0]],y[parent_child_pairs[j][1]]]
    zs = [z[parent_child_pairs[j][0]],z[parent_child_pairs[j][1]]]
    
    line = plt3d.art3d.Line3D(xs, ys, zs, c='black')
    ax.add_line(line)

#ax.quiver(0,0,0,repos_a[0][0],repos_a[0][1],repos_a[0][2],color='gray',alpha=0.6,linestyle='--')
ax.add_line(plt3d.art3d.Line3D([0,pos_a[0][0]], [0,pos_a[0][1]], [0,pos_a[0][2]] ,c='gray',alpha=0.6,linestyle='--'))

x=[]
y=[]
z=[]
# plot chain b
for i in range(len(pos_b)):
    x.append(pos_b[i][0])
    y.append(pos_b[i][1])
    z.append(pos_b[i][2])
    
    #dirxlocal=ori_b[i]@[1,0,0]
    #dirylocal=ori_b[i]@[0,1,0]
    #dirzlocal=ori_b[i]@[0,0,1]

    #ax.quiver(x[i],y[i],z[i],dirxlocal[0],dirxlocal[1],dirxlocal[2],length=0.05,color='r')
    #ax.quiver(x[i],y[i],z[i],dirylocal[0],dirylocal[1],dirylocal[2],length=0.05,color='g')
    #ax.quiver(x[i],y[i],z[i],dirzlocal[0],dirzlocal[1],dirzlocal[2],length=0.05,color='b')
ax.scatter(x, y, z, c=['black']*len(chain_b), marker='o')

for i in range(len(chain_a)):
    ax.scatter(x[i], y[i], z[i], c='black', marker='.')
    ax.text(x[i], y[i], z[i],  '%s' % (str(chain_b[i][0])), size=8, zorder=1,color='black')

for j in range(len(parent_child_pairs)):
    xs = [x[parent_child_pairs[j][0]],x[parent_child_pairs[j][1]]]
    ys = [y[parent_child_pairs[j][0]],y[parent_child_pairs[j][1]]]
    zs = [z[parent_child_pairs[j][0]],z[parent_child_pairs[j][1]]]
    
    line = plt3d.art3d.Line3D(xs, ys, zs, c='gray')
    ax.add_line(line)

#ax.quiver(0,0,0,pos_b[0][0],pos_b[0][1],pos_b[0][2],color='gray',alpha=0.6,linestyle='--')
ax.add_line(plt3d.art3d.Line3D([0,pos_b[0][0]], [0,pos_b[0][1]], [0,pos_b[0][2]] ,c='gray',alpha=0.6,linestyle='--'))

ax.set_xlim3d(-2,2)
ax.set_ylim3d(-2,2)
ax.set_zlim3d(-2,2)

ax.set_xlabel('Global X')
ax.set_ylabel('Global Y')
ax.set_zlabel('Global Z')

plt.show()
# ----------------------------------------------------------------------------------------------------------------



# Define complex chain a
# ----------------------------------------------------------------------------------------------------------------
# The chans are expressed using local orientations (relative to parent) as 3x3 rotation matrices, and local position offset vectors (position of child relative to parent)
chain_a=[]
chain_a_repos=[]


chain_a.append(['jnt_a0', np.array([[ 0.1466, 0.7779, -0.6111, 0],
                                    [ -0.04568, -0.6118, -0.7897, 0.2],
                                    [ -0.9881, 0.1437, -0.05418, 0],
                                    [ 0, 0, 0, 1]])])
chain_a.append(['jnt_a1', np.array([[ -0.531, 0.7242, -0.4399, -0.1],
                                    [ -0.267, -0.6357, -0.7243, -0.1],
                                    [ -0.8042, -0.2672, 0.5309, -0.1],
                                    [ 0, 0, 0, 1]])])               
chain_a.append(['jnt_a2', np.array([[ -0.1898, 0.7854, -0.5892, 0.2],
                                    [ -0.395, -0.6104, -0.6865, -0.2],
                                    [ -0.8989, 0.1024, 0.4261, 0],
                                    [ 0, 0, 0, 1]])])
chain_a.append(['jnt_a3', np.array([[ 0.8989, -0.1024, -0.4261, 0.1],
                                    [ -0.395, -0.6104, -0.6865, 0],
                                    [ -0.1898, 0.7854, -0.5892, 0.1],
                                    [ 0, 0, 0, 1]])])
chain_a.append(['jnt_a4', np.array([[ -0.6784, -0.2447, 0.6928, 0],
                                    [ -0.4144, -0.6512, -0.6358, 0.1],
                                    [ 0.6067, -0.7184, 0.3404, 0.1],
                                    [ 0, 0, 0, 1]])])
chain_a.append(['jnt_a5', np.array([[ 0.3584, 0.5254, -0.7717, 0.1],
                                    [ 0.4144, 0.6512, 0.6358, 0.2],
                                    [ 0.8366, -0.5476, 0.01571, 0],
                                    [ 0, 0, 0, 1]])])
chain_a.append(['jnt_a6', np.array([[ -0.3819, 0.8207, -0.425, 0],
                                    [ -0.2589, -0.5365, -0.8032, 0],
                                    [ -0.8872, -0.1967, 0.4174, -0.1],
                                    [ 0, 0, 0, 1]])])
chain_a.append(['jnt_a7', np.array([[ -0.3384, 0.8237, -0.455, 0],
                                    [ -0.688, -0.5464, -0.4776, 0.2],
                                    [ -0.642, 0.1514, 0.7516, -0.1],
                                    [ 0, 0, 0, 1]])])
chain_a.append(['jnt_a8', np.array([[ 0.1306, 0.5971, 0.7915, 0.1],
                                    [ -0.4515, 0.7465, -0.4888, 0.1],
                                    [ -0.8827, -0.2935, 0.367, -0.1],
                                    [ 0, 0, 0, 1]])])
chain_a.append(['jnt_a9', np.array([[ -0.4515, 0.7465, -0.4888, -0.1],
                                    [ -0.1306, -0.5971, -0.7915, -0.1],
                                    [ -0.8827, -0.2935, 0.367, 0],
                                    [ 0, 0, 0, 1]])])
chain_a.append(['jnt_a10', np.array([[ -0.02267, 0.9877, -0.1547, 0],
                                    [ 0.9976,  0.03246, 0.06109, 0],
                                    [ 0.06537,  -0.153, -0.9861, -0.1],
                                    [ 0, 0, 0, 1]])])

# a directed graph with parent/child relationship knowledge is needed for constructing the kinematic chain
dir_graph = {0: [1],
             1: [2, 5],
             2: [3, 4],
             5: [6, 7],
             7: [8, 9, 10]}

# find parent/child pairs from the directed graph
parent_child_pairs=parent_child(dir_graph)
# ----------------------------------------------------------------------------------------------------------------



# Calculate global orientations and positions for complex chain a
# ----------------------------------------------------------------------------------------------------------------

ori_a, pos_a = FK_MDH(copy.deepcopy(chain_a),dir_graph)
# ----------------------------------------------------------------------------------------------------------------



# Plot complex chain a
# ----------------------------------------------------------------------------------------------------------------
fig = plt.figure(4)
ax = fig.add_subplot(111, projection='3d')
fig.suptitle('Complex open kinematic chain a, with branching, plotted using forward kinematics', fontsize=12)

dirxglobal=[1,0,0]
diryglobal=[0,1,0]
dirzglobal=[0,0,1]

ax.quiver(0,0,0,dirxglobal[0],dirxglobal[1],dirxglobal[2],length=0.05,color='r')
ax.quiver(0,0,0,diryglobal[0],diryglobal[1],diryglobal[2],length=0.05,color='g')
ax.quiver(0,0,0,dirzglobal[0],dirzglobal[1],dirzglobal[2],length=0.05,color='b')

ax.text(0,0,0, 'ref. space', size=8, zorder=1,color='black')

# plot chain a
x=[]
y=[]
z=[]
for i in range(len(pos_a)):
    x.append(pos_a[i][0])
    y.append(pos_a[i][1])
    z.append(pos_a[i][2])
    
    dirxlocal=ori_a[i]@[1,0,0]
    dirylocal=ori_a[i]@[0,1,0]
    dirzlocal=ori_a[i]@[0,0,1]

    ax.quiver(x[i],y[i],z[i],dirxlocal[0],dirxlocal[1],dirxlocal[2],length=0.05,color='r')
    ax.quiver(x[i],y[i],z[i],dirylocal[0],dirylocal[1],dirylocal[2],length=0.05,color='g')
    ax.quiver(x[i],y[i],z[i],dirzlocal[0],dirzlocal[1],dirzlocal[2],length=0.05,color='b')
ax.scatter(x, y, z, c=['gray']*len(chain_a), marker='o')

for i in range(len(chain_a)):
    ax.scatter(x[i], y[i], z[i], c='black', marker='.')
    ax.text(x[i], y[i], z[i],  '%s' % (str(chain_a[i][0])), size=8, zorder=1,color='black')

for j in range(len(parent_child_pairs)):
    xs = [x[parent_child_pairs[j][0]],x[parent_child_pairs[j][1]]]
    ys = [y[parent_child_pairs[j][0]],y[parent_child_pairs[j][1]]]
    zs = [z[parent_child_pairs[j][0]],z[parent_child_pairs[j][1]]]
    
    line = plt3d.art3d.Line3D(xs, ys, zs, c='black')
    ax.add_line(line)

#ax.quiver(0,0,0,pos_a[0][0],pos_a[0][1],pos_a[0][2],color='gray',alpha=0.6,linestyle='--')
ax.add_line(plt3d.art3d.Line3D([0,pos_a[0][0]], [0,pos_a[0][1]], [0,pos_a[0][2]] ,c='gray',alpha=0.6,linestyle='--'))

ax.set_xlim3d(-2,2)
ax.set_ylim3d(-2,2)
ax.set_zlim3d(-2,2)

ax.set_xlabel('Global X')
ax.set_ylabel('Global Y')
ax.set_zlabel('Global Z')

plt.show()
# ----------------------------------------------------------------------------------------------------------------


# Forward Kinematics repositioning example
# Define complex chain a, and the repositioned 
# ----------------------------------------------------------------------------------------------------------------
# The chans are expressed using local orientations (relative to parent) as 3x3 rotation matrices, and local position offset vectors (position of child relative to parent)
chain_a=[]

chain_a.append(['jnt_a0', np.array([[ 0.1466, 0.7779, -0.6111, 0],
                                    [ -0.04568, -0.6118, -0.7897, 0.2],
                                    [ -0.9881, 0.1437, -0.05418, 0],
                                    [ 0, 0, 0, 1]])])
chain_a.append(['jnt_a1', np.array([[ -0.531, 0.7242, -0.4399, -0.1],
                                    [ -0.267, -0.6357, -0.7243, -0.1],
                                    [ -0.8042, -0.2672, 0.5309, -0.1],
                                    [ 0, 0, 0, 1]])])               
chain_a.append(['jnt_a2', np.array([[ -0.1898, 0.7854, -0.5892, 0.2],
                                    [ -0.395, -0.6104, -0.6865, -0.2],
                                    [ -0.8989, 0.1024, 0.4261, 0],
                                    [ 0, 0, 0, 1]])])
chain_a.append(['jnt_a3', np.array([[ 0.8989, -0.1024, -0.4261, 0.1],
                                    [ -0.395, -0.6104, -0.6865, 0],
                                    [ -0.1898, 0.7854, -0.5892, 0.1],
                                    [ 0, 0, 0, 1]])])
chain_a.append(['jnt_a4', np.array([[ -0.6784, -0.2447, 0.6928, 0],
                                    [ -0.4144, -0.6512, -0.6358, 0.1],
                                    [ 0.6067, -0.7184, 0.3404, 0.1],
                                    [ 0, 0, 0, 1]])])
chain_a.append(['jnt_a5', np.array([[ 0.3584, 0.5254, -0.7717, 0.1],
                                    [ 0.4144, 0.6512, 0.6358, 0.2],
                                    [ 0.8366, -0.5476, 0.01571, 0],
                                    [ 0, 0, 0, 1]])])
chain_a.append(['jnt_a6', np.array([[ -0.3819, 0.8207, -0.425, 0],
                                    [ -0.2589, -0.5365, -0.8032, 0],
                                    [ -0.8872, -0.1967, 0.4174, -0.1],
                                    [ 0, 0, 0, 1]])])
chain_a.append(['jnt_a7', np.array([[ -0.3384, 0.8237, -0.455, 0],
                                    [ -0.688, -0.5464, -0.4776, 0.2],
                                    [ -0.642, 0.1514, 0.7516, -0.1],
                                    [ 0, 0, 0, 1]])])
chain_a.append(['jnt_a8', np.array([[ 0.1306, 0.5971, 0.7915, 0.1],
                                    [ -0.4515, 0.7465, -0.4888, 0.1],
                                    [ -0.8827, -0.2935, 0.367, -0.1],
                                    [ 0, 0, 0, 1]])])
chain_a.append(['jnt_a9', np.array([[ -0.4515, 0.7465, -0.4888, -0.1],
                                    [ -0.1306, -0.5971, -0.7915, -0.1],
                                    [ -0.8827, -0.2935, 0.367, 0],
                                    [ 0, 0, 0, 1]])])
chain_a.append(['jnt_a10', np.array([[ -0.02267, 0.9877, -0.1547, 0],
                                    [ 0.9976,  0.03246, 0.06109, 0],
                                    [ 0.06537,  -0.153, -0.9861, -0.1],
                                    [ 0, 0, 0, 1]])])


#reorient a joint by 45 degrees about the local coordinate system
chain_a_reori=copy.deepcopy(chain_a)
chain_a_reori[5][1]= chain_a_reori[5][1] @ np.array([[ 1, 0, 0, 0],
                                                      [ 0, 0.7071068, -0.7071068, 0],
                                                      [ 0, 0.7071068, 0.7071068, 0],
                                                      [ 0, 0, 0, 1]])

# a directed graph with parent/child relationship knowledge is needed for constructing the kinematic chain
dir_graph = {0: [1],
             1: [2, 5],
             2: [3, 4],
             5: [6, 7],
             7: [8, 9, 10]}

# find parent/child pairs from the directed graph
parent_child_pairs=parent_child(dir_graph)
# ----------------------------------------------------------------------------------------------------------------



# Calculate global orientations and positions for complex chain a
# ----------------------------------------------------------------------------------------------------------------

ori_a, pos_a = FK_MDH(copy.deepcopy(chain_a),dir_graph)
reori_a, repos_a = FK_MDH(copy.deepcopy(chain_a_reori),dir_graph)
# ----------------------------------------------------------------------------------------------------------------



# Plot complex chain a
# ----------------------------------------------------------------------------------------------------------------
fig = plt.figure(5)
ax = fig.add_subplot(111, projection='3d')
fig.suptitle('Reoriented complex open kinematic chain a, with branching, plotted using forward kinematics', fontsize=12)

dirxglobal=[1,0,0]
diryglobal=[0,1,0]
dirzglobal=[0,0,1]

ax.quiver(0,0,0,dirxglobal[0],dirxglobal[1],dirxglobal[2],length=0.05,color='r')
ax.quiver(0,0,0,diryglobal[0],diryglobal[1],diryglobal[2],length=0.05,color='g')
ax.quiver(0,0,0,dirzglobal[0],dirzglobal[1],dirzglobal[2],length=0.05,color='b')

ax.text(0,0,0, 'ref. space', size=8, zorder=1,color='black')

# plot chain a
x=[]
y=[]
z=[]
for i in range(len(pos_a)):
    x.append(pos_a[i][0])
    y.append(pos_a[i][1])
    z.append(pos_a[i][2])
    
    dirxlocal=ori_a[i]@[1,0,0]
    dirylocal=ori_a[i]@[0,1,0]
    dirzlocal=ori_a[i]@[0,0,1]

    ax.quiver(x[i],y[i],z[i],dirxlocal[0],dirxlocal[1],dirxlocal[2],length=0.05,linestyle='--', alpha=0.3,color='r')
    ax.quiver(x[i],y[i],z[i],dirylocal[0],dirylocal[1],dirylocal[2],length=0.05,linestyle='--', alpha=0.3,color='g')
    ax.quiver(x[i],y[i],z[i],dirzlocal[0],dirzlocal[1],dirzlocal[2],length=0.05,linestyle='--', alpha=0.3,color='b')
ax.scatter(x, y, z, c=['gray']*len(chain_a), marker='o')

#for i in range(len(chain_a)):
#    ax.scatter(x[i], y[i], z[i], c='black', marker='.')
#    ax.text(x[i], y[i], z[i],  '%s' % (str(chain_a[i][0])), size=8, zorder=1,color='black')

for j in range(len(parent_child_pairs)):
    xs = [x[parent_child_pairs[j][0]],x[parent_child_pairs[j][1]]]
    ys = [y[parent_child_pairs[j][0]],y[parent_child_pairs[j][1]]]
    zs = [z[parent_child_pairs[j][0]],z[parent_child_pairs[j][1]]]
    
    line = plt3d.art3d.Line3D(xs, ys, zs, linestyle='--', alpha=0.3, c='black')
    ax.add_line(line)

#ax.quiver(0,0,0,pos_a[0][0],pos_a[0][1],pos_a[0][2],color='gray',alpha=0.6,linestyle='--')
ax.add_line(plt3d.art3d.Line3D([0,pos_a[0][0]], [0,pos_a[0][1]], [0,pos_a[0][2]] ,c='gray',alpha=0.6,linestyle='--'))

# plot chain a with reorientation
x=[]
y=[]
z=[]
for i in range(len(repos_a)):
    x.append(repos_a[i][0])
    y.append(repos_a[i][1])
    z.append(repos_a[i][2])
    
    dirxlocal=reori_a[i]@[1,0,0]
    dirylocal=reori_a[i]@[0,1,0]
    dirzlocal=reori_a[i]@[0,0,1]

    ax.quiver(x[i],y[i],z[i],dirxlocal[0],dirxlocal[1],dirxlocal[2],length=0.05,color='r')
    ax.quiver(x[i],y[i],z[i],dirylocal[0],dirylocal[1],dirylocal[2],length=0.05,color='g')
    ax.quiver(x[i],y[i],z[i],dirzlocal[0],dirzlocal[1],dirzlocal[2],length=0.05,color='b')
ax.scatter(x, y, z, c=['gray']*len(chain_a_reori), marker='o')

for i in range(len(chain_a_reori)):
    ax.scatter(x[i], y[i], z[i], c='black', marker='.')
    ax.text(x[i], y[i], z[i],  '%s' % (str(chain_a_reori[i][0])), size=8, zorder=1,color='black')

for j in range(len(parent_child_pairs)):
    xs = [x[parent_child_pairs[j][0]],x[parent_child_pairs[j][1]]]
    ys = [y[parent_child_pairs[j][0]],y[parent_child_pairs[j][1]]]
    zs = [z[parent_child_pairs[j][0]],z[parent_child_pairs[j][1]]]
    
    line = plt3d.art3d.Line3D(xs, ys, zs, c='black')
    ax.add_line(line)

ax.set_xlim3d(-2,2)
ax.set_ylim3d(-2,2)
ax.set_zlim3d(-2,2)

ax.set_xlabel('Global X')
ax.set_ylabel('Global Y')
ax.set_zlabel('Global Z')

plt.show()


# ----------------------------------------------------------------------------------------------------------------
