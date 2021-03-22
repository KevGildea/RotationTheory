""" Implementation of the methods described in https://kevgildea.github.io/blog/Euler-Axis-Vector-Mapping/"""

import numpy as np
import matplotlib.pyplot as plt
import mpl_toolkits.mplot3d as plt3d
import math


def Vector_mapping_cross_product(vec1, vec2):
    """ Calculate the rotation matrix that maps unit vector a to align with unit vector b about an axis aligned with the cross product"""
    a, b = vec1 / np.linalg.norm(vec1), vec2 / np.linalg.norm(vec2)
    n = np.cross(a, b) / np.linalg.norm(np.cross(a, b))
    rotation_matrix = np.array([[n[0]**2+(n[1]**2+n[2]**2)*(np.dot(a, b)),n[0]*n[1]*(1-np.dot(a, b))-n[2]*np.linalg.norm(np.cross(a, b)),n[0]*n[2]*(1-np.dot(a, b))+n[1]*np.linalg.norm(np.cross(a, b))],
                                [n[0]*n[1]*(1-np.dot(a, b))+n[2]*np.linalg.norm(np.cross(a, b)),n[1]**2+(n[0]**2+n[2]**2)*(np.dot(a, b)),n[1]*n[2]*(1-np.dot(a, b))-n[0]*np.linalg.norm(np.cross(a, b))],
                                [n[0]*n[2]*(1-np.dot(a, b))-n[1]*np.linalg.norm(np.cross(a, b)),n[1]*n[2]*(1-np.dot(a, b))+n[0]*np.linalg.norm(np.cross(a, b)),n[2]**2+(n[0]**2+n[1]**2)*(np.dot(a, b))]])
    return rotation_matrix


def Vector_mapping_bisect(vec1, vec2):
    """ Calculate the rotation matrix that maps unit vector a to align with unit vector b about an axis aligned with norm(a) + norm(b)"""
    a, b = vec1 / np.linalg.norm(vec1), vec2 / np.linalg.norm(vec2)
    n = (a+b) / np.linalg.norm(a+b)
    θ = math.pi
    rotation_matrix = np.array([[n[0]**2+(n[1]**2+n[2]**2)*(np.cos(θ)),n[0]*n[1]*(1-np.cos(θ))-n[2]*np.sin(θ),n[0]*n[2]*(1-np.cos(θ))+n[1]*np.sin(θ)],
                                [n[0]*n[1]*(1-np.cos(θ))+n[2]*np.sin(θ),n[1]**2+(n[0]**2+n[2]**2)*(np.cos(θ)),n[1]*n[2]*(1-np.cos(θ))-n[0]*np.sin(θ)],
                                [n[0]*n[2]*(1-np.cos(θ))-n[1]*np.sin(θ),n[1]*n[2]*(1-np.cos(θ))+n[0]*np.sin(θ),n[2]**2+(n[0]**2+n[1]**2)*(np.cos(θ))]])
    return rotation_matrix


def Vector_mapping_UDaxis(vec1, vec2, axis): 
    """ Calculate the rotation matrix that maps unit vector a to align with unit vector b along an user defined axis"""
    a, b = vec1 / np.linalg.norm(vec1), vec2 / np.linalg.norm(vec2)
    n = axis / np.linalg.norm(axis)
    # project vectors to form the base of a right cone around the Euler axis
    a, b = np.cross(a,n) / np.linalg.norm(np.cross(a,n)), np.cross(b,n) / np.linalg.norm(np.cross(b,n))
    θ = np.arccos(np.dot(a,b))*np.sign(np.dot(n, np.cross(a,b)))
    rotation_matrix = np.array([[n[0]**2+(n[1]**2+n[2]**2)*(np.cos(θ)),n[0]*n[1]*(1-np.cos(θ))-n[2]*np.sin(θ),n[0]*n[2]*(1-np.cos(θ))+n[1]*np.sin(θ)],
                                [n[0]*n[1]*(1-np.cos(θ))+n[2]*np.sin(θ),n[1]**2+(n[0]**2+n[2]**2)*(np.cos(θ)),n[1]*n[2]*(1-np.cos(θ))-n[0]*np.sin(θ)],
                                [n[0]*n[2]*(1-np.cos(θ))-n[1]*np.sin(θ),n[1]*n[2]*(1-np.cos(θ))+n[0]*np.sin(θ),n[2]**2+(n[0]**2+n[1]**2)*(np.cos(θ))]])
    return rotation_matrix


# Plot vectors a and b
fig = plt.figure(1)
ax = fig.add_subplot(111, projection='3d')
fig.suptitle('Vectors a and b', fontsize=12)

a=[0.1,0.3,0.4]
b=[0.3,0.1,0.2]

# addition of a scalar multiple of a to b (and vice-versa) gives the same cross product, however, the rotation angle approaches 0 (and A21 approaches the identity matrix) as the scalar multiple approaches infinity:
# print(np.cross(a, b))
# scalar1=-1.5
# scalar2=3.5
# a=np.array(a)+np.array(b)*scalar1
# b=np.array(b)+np.array(a)*scalar2
# print(np.cross(a, b))

axb=np.cross(a, b)

ax.quiver(0,0,0,a[0],a[1],a[2],color='black')
ax.text(a[0],a[1],a[2],  '%s' % (r'$\vec a$'), size=12, zorder=1,color='black')
ax.quiver(0,0,0,b[0],b[1],b[2],color='grey')
ax.text(b[0],b[1],b[2],  '%s' % (r'$\vec b$'), size=12, zorder=1,color='grey')
ax.quiver(0,0,0,axb[0],axb[1],axb[2],color='orange',linestyle='--')
ax.text(axb[0],axb[1],axb[2],  '%s' % (r'$\vec {axb}$'), size=12, zorder=1,color='orange')

dirxglobal=[1,0,0]
diryglobal=[0,1,0]
dirzglobal=[0,0,1]

ax.quiver(-0.5,-0.5,-0.5,dirxglobal[0],dirxglobal[1],dirxglobal[2],color='r')
ax.quiver(-0.5,-0.5,-0.5,diryglobal[0],diryglobal[1],diryglobal[2],color='g')
ax.quiver(-0.5,-0.5,-0.5,dirzglobal[0],dirzglobal[1],dirzglobal[2],color='b')

ax.set_xlim3d(-1,1)
ax.set_ylim3d(-1,1)
ax.set_zlim3d(-1,1)

ax.set_xlabel('Global X')
ax.set_ylabel('Global Y')
ax.set_zlabel('Global Z')

plt.show()


# Plot rotation to map vectors about axis aligned with the cross product
fig = plt.figure(2)
ax = fig.add_subplot(111, projection='3d')
fig.suptitle('Rotation to map vectors about axis aligned with the cross product', fontsize=12)

ax.quiver(0,0,0,a[0],a[1],a[2],color='black')
ax.text(a[0],a[1],a[2],  '%s' % (r'$\vec a$'), size=12, zorder=1,color='black')
ax.quiver(0,0,0,b[0],b[1],b[2],color='grey')
ax.text(b[0],b[1],b[2],  '%s' % (r'$\vec b$'), size=12, zorder=1,color='grey')
ax.quiver(0,0,0,axb[0],axb[1],axb[2],color='orange',linestyle='--')
ax.text(axb[0],axb[1],axb[2],  '%s' % (r'$\vec {axb}$'), size=12, zorder=1,color='orange')

a_mapped = Vector_mapping_cross_product(a, b) @ a
ax.quiver(0,0,0,a_mapped[0],a_mapped[1],a_mapped[2],color='black',linestyle='--')
ax.text(a_mapped[0],a_mapped[1],a_mapped[2],  '%s' % (r'$\vec a_{mapped}$'), size=12, zorder=1,color='black')

θ = np.arccos((np.trace(Vector_mapping_cross_product(a, b))-1)/2)

dirxglobal=[1,0,0]
diryglobal=[0,1,0]
dirzglobal=[0,0,1]

ax.quiver(-0.5,-0.5,-0.5,dirxglobal[0],dirxglobal[1],dirxglobal[2],color='r')
ax.quiver(-0.5,-0.5,-0.5,diryglobal[0],diryglobal[1],diryglobal[2],color='g')
ax.quiver(-0.5,-0.5,-0.5,dirzglobal[0],dirzglobal[1],dirzglobal[2],color='b')

n = (axb / np.linalg.norm(axb)).reshape(3)
ax.quiver(-0.5,-0.5,-0.5,n[0],n[1],n[2],color='orange',linestyle='--')
ax.text(n[0]-0.5,n[1]-0.5,n[2]-0.5,  '%s' % (r'$\vec {n}$'+',  θ='+ str(round(θ, 5))), size=12, zorder=1,color='orange')

dirxglobal=Vector_mapping_cross_product(a, b) @ [1,0,0]
diryglobal=Vector_mapping_cross_product(a, b) @ [0,1,0]
dirzglobal=Vector_mapping_cross_product(a, b) @ [0,0,1]

ax.quiver(-0.5,-0.5,-0.5,dirxglobal[0],dirxglobal[1],dirxglobal[2],color='r',linestyle='--')
ax.quiver(-0.5,-0.5,-0.5,diryglobal[0],diryglobal[1],diryglobal[2],color='g',linestyle='--')
ax.quiver(-0.5,-0.5,-0.5,dirzglobal[0],dirzglobal[1],dirzglobal[2],color='b',linestyle='--')

ax.set_xlim3d(-1,1)
ax.set_ylim3d(-1,1)
ax.set_zlim3d(-1,1)

ax.set_xlabel('Global X')
ax.set_ylabel('Global Y')
ax.set_zlabel('Global Z')

plt.show()


# Plot rotation to map vectors about axis bisecting vectors a and b

# choose a rotation axis in the plane bisecting a and b i.e. where a and b are symmetrical
x=1
axis = x* (b / np.linalg.norm(b) + (a / np.linalg.norm(a)))

fig = plt.figure(2)
ax = fig.add_subplot(111, projection='3d')
fig.suptitle('Rotation to map vectors about axis bisecting vectors a and b', fontsize=12)

ax.quiver(0,0,0,a[0],a[1],a[2],color='black')
ax.text(a[0],a[1],a[2],  '%s' % (r'$\vec a$'), size=12, zorder=1,color='black')
ax.quiver(0,0,0,b[0],b[1],b[2],color='grey')
ax.text(b[0],b[1],b[2],  '%s' % (r'$\vec b$'), size=12, zorder=1,color='grey')
ax.quiver(0,0,0,axis[0],axis[1],axis[2],color='orange',linestyle='--')
ax.text(axis[0],axis[1],axis[2],  '%s' % (r'$\vec {b+a}$'), size=12, zorder=1,color='orange')

a_mapped = Vector_mapping_bisect(a, b) @ a
ax.quiver(0,0,0,a_mapped[0],a_mapped[1],a_mapped[2],color='black',linestyle='--')
ax.text(a_mapped[0],a_mapped[1],a_mapped[2],  '%s' % (r'$\vec a_{mapped}$'), size=12, zorder=1,color='black')

θ = np.arccos((np.trace(Vector_mapping_bisect(a, b))-1)/2)

dirxglobal=[1,0,0]
diryglobal=[0,1,0]
dirzglobal=[0,0,1]

ax.quiver(-0.5,-0.5,-0.5,dirxglobal[0],dirxglobal[1],dirxglobal[2],color='r')
ax.quiver(-0.5,-0.5,-0.5,diryglobal[0],diryglobal[1],diryglobal[2],color='g')
ax.quiver(-0.5,-0.5,-0.5,dirzglobal[0],dirzglobal[1],dirzglobal[2],color='b')

n = (axis / np.linalg.norm(axis)).reshape(3)
ax.quiver(-0.5,-0.5,-0.5,n[0],n[1],n[2],color='orange',linestyle='--')
ax.text(n[0]-0.5,n[1]-0.5,n[2]-0.5,  '%s' % (r'$\vec {n}$'+',  θ='+ str(round(θ, 5))), size=12, zorder=1,color='orange')

dirxglobal=Vector_mapping_bisect(a, b) @ [1,0,0]
diryglobal=Vector_mapping_bisect(a, b) @ [0,1,0]
dirzglobal=Vector_mapping_bisect(a, b) @ [0,0,1]

ax.quiver(-0.5,-0.5,-0.5,dirxglobal[0],dirxglobal[1],dirxglobal[2],color='r',linestyle='--')
ax.quiver(-0.5,-0.5,-0.5,diryglobal[0],diryglobal[1],diryglobal[2],color='g',linestyle='--')
ax.quiver(-0.5,-0.5,-0.5,dirzglobal[0],dirzglobal[1],dirzglobal[2],color='b',linestyle='--')

ax.set_xlim3d(-1,1)
ax.set_ylim3d(-1,1)
ax.set_zlim3d(-1,1)

ax.set_xlabel('Global X')
ax.set_ylabel('Global Y')
ax.set_zlabel('Global Z')

plt.show()

# Plot possible rotation axis space - any vector on the bisecting symmetric plane for a and b

# i.e. create the plane which symmetrically bisects vectors a and b

fig = plt.figure(2)
ax = fig.add_subplot(111, projection='3d')
fig.suptitle('Possible rotation axis space - any vector on the bisecting symmetric plane for a and b', fontsize=12)

p1 = np.array([0,0,0])
p2 = b / np.linalg.norm(b) + (a / np.linalg.norm(a))
p3 = np.cross(a, b)

cp = np.cross(p2, p3)

#print('The equation is {0}x + {1}y + {2}z = {3}'.format(cp[0], cp[1], cp[2], np.dot(cp, p3)))
ax.text(0,0.7,0.7,  '{0}x + {1}y + {2}z = {3}'.format(round(cp[0],3), round(cp[1],3), round(cp[2],3), round(np.dot(cp, p3),3)), size=10, zorder=1,color='grey')

x = np.linspace(-4, 4, 200)
y = np.linspace(-4, 4, 200)
X, Y = np.meshgrid(x, y)

Z = (np.dot(cp, p3) - cp[0] * X - cp[1] * Y) / cp[2]

ax.plot(X.flatten(),
        Y.flatten(),
        Z.flatten(), 'bo', marker='.', alpha=0.3)

ax.plot(*zip(p1, p2, p3), color='r', linestyle=' ', marker='o')

ax.quiver(0,0,0,a[0],a[1],a[2],color='black')
ax.text(a[0],a[1],a[2],  '%s' % (r'$\vec a$'), size=12, zorder=1,color='black')
ax.quiver(0,0,0,b[0],b[1],b[2],color='grey')
ax.text(b[0],b[1],b[2],  '%s' % (r'$\vec b$'), size=12, zorder=1,color='grey')

dirxglobal=[1,0,0]
diryglobal=[0,1,0]
dirzglobal=[0,0,1]

ax.quiver(-0.5,-0.5,-0.5,dirxglobal[0],dirxglobal[1],dirxglobal[2],color='r')
ax.quiver(-0.5,-0.5,-0.5,diryglobal[0],diryglobal[1],diryglobal[2],color='g')
ax.quiver(-0.5,-0.5,-0.5,dirzglobal[0],dirzglobal[1],dirzglobal[2],color='b')

ax.set_xlim3d(-1,1)
ax.set_ylim3d(-1,1)
ax.set_zlim3d(-1,1)

ax.set_xlabel('Global X')
ax.set_ylabel('Global Y')
ax.set_zlabel('Global Z')

plt.show()

# Plot rotation to map vectors about axis bisecting vectors a and b

# p2 = b / np.linalg.norm(b) + (a / np.linalg.norm(a))
# p3 = np.cross(a, b) / np.linalg.norm(np.cross(a, b))
# s1 = 0.2
# s2 = -1
# axis = s1*p3 + s2*p2

# chose a vector on this plane for a new rotation axis choosing x and y
X = 0.5
Y = -0.7
Z = (np.dot(cp, p3) - cp[0] * X - cp[1] * Y) / cp[2]

axis = np.array([X,Y,Z])

fig = plt.figure(2)
ax = fig.add_subplot(111, projection='3d')
fig.suptitle('Rotation to map vectors about axis bisecting vectors a and b', fontsize=12)

ax.quiver(0,0,0,a[0],a[1],a[2],color='black')
ax.text(a[0],a[1],a[2],  '%s' % (r'$\vec a$'), size=12, zorder=1,color='black')
ax.quiver(0,0,0,b[0],b[1],b[2],color='grey')
ax.text(b[0],b[1],b[2],  '%s' % (r'$\vec b$'), size=12, zorder=1,color='grey')
ax.quiver(0,0,0,axis[0],axis[1],axis[2],color='orange',linestyle='--')
ax.text(axis[0],axis[1],axis[2],  '%s' % (r'$\vec {ax}$'), size=12, zorder=1,color='orange')

a_mapped = Vector_mapping_UDaxis(a, b, axis) @ a
ax.quiver(0,0,0,a_mapped[0],a_mapped[1],a_mapped[2],color='black',linestyle='--')
ax.text(a_mapped[0],a_mapped[1],a_mapped[2],  '%s' % (r'$\vec a_{mapped}$'), size=12, zorder=1,color='black')

θ = np.arccos((np.trace(Vector_mapping_UDaxis(a, b, axis))-1)/2)

dirxglobal=[1,0,0]
diryglobal=[0,1,0]
dirzglobal=[0,0,1]

ax.quiver(-0.5,-0.5,-0.5,dirxglobal[0],dirxglobal[1],dirxglobal[2],color='r')
ax.quiver(-0.5,-0.5,-0.5,diryglobal[0],diryglobal[1],diryglobal[2],color='g')
ax.quiver(-0.5,-0.5,-0.5,dirzglobal[0],dirzglobal[1],dirzglobal[2],color='b')

n = (axis / np.linalg.norm(axis)).reshape(3)
ax.quiver(-0.5,-0.5,-0.5,n[0],n[1],n[2],color='orange',linestyle='--')
ax.text(n[0]-0.5,n[1]-0.5,n[2]-0.5,  '%s' % (r'$\vec {n}$'+',  θ='+ str(round(θ, 5))), size=12, zorder=1,color='orange')

dirxglobal=Vector_mapping_UDaxis(a, b, axis) @ [1,0,0]
diryglobal=Vector_mapping_UDaxis(a, b, axis) @ [0,1,0]
dirzglobal=Vector_mapping_UDaxis(a, b, axis) @ [0,0,1]

ax.quiver(-0.5,-0.5,-0.5,dirxglobal[0],dirxglobal[1],dirxglobal[2],color='r',linestyle='--')
ax.quiver(-0.5,-0.5,-0.5,diryglobal[0],diryglobal[1],diryglobal[2],color='g',linestyle='--')
ax.quiver(-0.5,-0.5,-0.5,dirzglobal[0],dirzglobal[1],dirzglobal[2],color='b',linestyle='--')

ax.set_xlim3d(-1,1)
ax.set_ylim3d(-1,1)
ax.set_zlim3d(-1,1)

ax.set_xlabel('Global X')
ax.set_ylabel('Global Y')
ax.set_zlabel('Global Z')

plt.show()
