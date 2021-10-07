""" Implementation of the methods described in https://kevgildea.github.io/blog/Kinematic-differential-equations/"""

import matplotlib.pyplot as plt
import numpy as np
import mpl_toolkits.mplot3d as plt3d
import sys
from math import isclose

# Functions
# ----------------------------------------------------------------------------------------------------------------

def isRotationMatrix(M):
    """ Checks whether a rotation matrix is orthogonal, and has a determinant of 1 (tolerance of 1e-3)"""
    tag = False
    I = np.identity(M.shape[0])
    check1 = np.isclose((np.matmul(M, M.T)), I, atol=1e-3)
    check2 = isclose(np.linalg.det(M),1,abs_tol=1e-3)
    if check1.all()==True and check2==True: tag = True
    return tag

def SkewSym(a):
    """ Converts a 1x3 vector to a 3x3 skew symmetric matrix"""
    a_tilda = np.array([[0, -a[2], a[1]],
                        [a[2], 0, a[0]],
                        [-a[1], a[0], 0]])
    return a_tilda

def RevSkewSym(a_tilda):
    """ Converts a 3x3 skew symmetric matrix to a 1x3 vector"""
    if a_tilda[0][0] != 0 or a_tilda[1][1] != 0 or a_tilda[2][2] != 0:
        sys.exit('Error: Matix is not skew symmetric')
    a = np.array([a_tilda[1][2], a_tilda[0][2], a_tilda[1][0]])
    return a

def EulerInt(A0,ω,r0,v,t_step,t_end):
    """ Performs Euler Integration to determine position and velocity given initial conditions, 1st derivatives, time step, and end time"""
    At=A0
    Ats=[]
    rt=r0
    rts=[]
    ω_tilda = SkewSym(ω)
    for t in range(int(t_end/t_step)):
        ω = At @ RevSkewSym(ω_tilda)
        ω_tilda = SkewSym(ω)
        At = At + t_step*(-ω_tilda @ At)
        Ats.append(At)
        correction_matrix = (3*np.array([[1, 0, 0],[0, 1, 0],[0, 0, 1]]) - (At @ At.T))/2 # correction matrix to ensure that the rotation matrix is orthogonal (see http://www.euclideanspace.com/maths/algebra/matrix/orthogonal/ )
        At = correction_matrix @ At
        if isRotationMatrix(At)==False:
            sys.exit('Error: Chosen integration time step is too large - try a smaller value (generally a step of <=1e-2 is recommended)')
        rt = rt + t_step*(At@v)
        rts.append(rt)
    return Ats,rts
# ----------------------------------------------------------------------------------------------------------------
 
  
# Initial conditions
# ----------------------------------------------------------------------------------------------------------------
t_end = 2
t_step = .01
A0 = np.array([[1, 0, 0],
               [0, 1, 0],
               [0, 0, 1]])
ω = np.array([1, 4, 2]) # specified in the local coordinate system
r0 = np.array([0, 1, 0])
vel = np.array([10, 10, 10]) # specified in the local coordinate system
ω0 = ω
# ----------------------------------------------------------------------------------------------------------------


# Euler integration
# ----------------------------------------------------------------------------------------------------------------
At,rt = EulerInt(A0,ω,r0,vel,t_step,t_end)

# ----------------------------------------------------------------------------------------------------------------




# Plot the ball at every time step
# ----------------------------------------------------------------------------------------------------------------
fig = plt.figure(1)
ax = fig.add_subplot(111, projection='3d')
fig.suptitle('Kinematic differential equations', fontsize=12)

# Plot the initial coordinate system
dirxinitial=A0 @ [1,0,0]
diryinitial=A0 @ [0,1,0]
dirzinitial=A0 @ [0,0,1]

ax.quiver(r0[0],r0[1],r0[2],dirxinitial[0],dirxinitial[1],dirxinitial[2],color='r',linestyle='--')
ax.quiver(r0[0],r0[1],r0[2],diryinitial[0],diryinitial[1],diryinitial[2],color='g',linestyle='--')
ax.quiver(r0[0],r0[1],r0[2],dirzinitial[0],dirzinitial[1],dirzinitial[2],color='b',linestyle='--')
ax.text(r0[0],r0[1],r0[2], '  t = 0ms', size=12, zorder=1)

#draw sphere   
centre = [r0[0],r0[1],r0[2]]
radius = 1
u, v = np.mgrid[0:2*np.pi:10*1j, 0:np.pi:10*1j]
ax.plot_surface(centre[0] + radius * np.cos(u) * np.sin(v),centre[1] + radius * np.sin(u) * np.sin(v),centre[2] + radius * np.cos(v), color="gray", alpha=0.1)

# plot the sphere after each timestep up to the t_end
for i in range(len(At)):
    dirx= At[i] @ [1,0,0]
    diry= At[i] @ [0,1,0]
    dirz= At[i] @ [0,0,1]

    ax.quiver(rt[i][0],rt[i][1],rt[i][2],dirx[0],dirx[1],dirx[2],color='r',linestyle='--')
    ax.quiver(rt[i][0],rt[i][1],rt[i][2],diry[0],diry[1],diry[2],color='g',linestyle='--')
    ax.quiver(rt[i][0],rt[i][1],rt[i][2],dirz[0],dirz[1],dirz[2],color='b',linestyle='--')
    t = round((i+1)*t_step, 3)
    if i in range(0,len(At))[0::10]:
        ax.text(rt[i][0],rt[i][1],rt[i][2], '  t = '+str(t)+'ms', size=12, zorder=1)
    if i == len(At)-1:
        ax.text(rt[i][0],rt[i][1],rt[i][2], '  t = '+str(t)+'ms', size=12, zorder=1)

    #draw sphere   
    centre = [rt[i][0],rt[i][1],rt[i][2]]
    radius = 1
    u, v = np.mgrid[0:2*np.pi:10*1j, 0:np.pi:10*1j]
    ax.plot_surface(centre[0] + radius * np.cos(u) * np.sin(v),centre[1] + radius * np.sin(u) * np.sin(v),centre[2] + radius * np.cos(v), color="gray", alpha=0.3)

# Axis limits and lables and legend
ax.set_xlim3d(-100,100)
ax.set_ylim3d(-100,100)
ax.set_zlim3d(-100,100)

ax.set_xlabel('Global X')
ax.set_ylabel('Global Y')
ax.set_zlabel('Global Z')

plt.show()
# ----------------------------------------------------------------------------------------------------------------



# Plot the ball at each time step and output images ( I used ffmpeg to create a gif)
# ----------------------------------------------------------------------------------------------------------------
filenames = []
fig = plt.figure(2)
ax = fig.add_subplot(111, projection='3d')

# Axis limits and lables and legend
ax.set_xlim3d(-10,10)
ax.set_ylim3d(0,20)
ax.set_zlim3d(-10,10)

ax.set_xlabel('Global X')
ax.set_ylabel('Global Y')
ax.set_zlabel('Global Z')

ax.view_init(10, -60)

dirx= A0 @ [1,0,0]
diry= A0 @ [0,1,0]
dirz= A0 @ [0,0,1]

ax.quiver(r0[0],r0[1],r0[2],dirx[0],dirx[1],dirx[2],color='r',linestyle='--')
ax.quiver(r0[0],r0[1],r0[2],diry[0],diry[1],diry[2],color='g',linestyle='--')
ax.quiver(r0[0],r0[1],r0[2],dirz[0],dirz[1],dirz[2],color='b',linestyle='--')
#ax.text(r0[0]+dirx[0],r0[1]+dirx[1],r0[2]+dirx[2], 'local x', size=12, zorder=1, color='r')
#ax.text(r0[0]+diry[0],r0[1]+diry[1],r0[2]+diry[2], 'local y', size=12, zorder=1, color='g')
#ax.text(r0[0]+dirz[0],r0[1]+dirz[1],r0[2]+dirz[2], 'local z', size=12, zorder=1, color='b')

# plot time/title
fig.suptitle('t = 0.0ms', fontsize=12)

#draw sphere   
centre = [r0[0],r0[1],r0[2]]
radius = 1
u, v = np.mgrid[0:2*np.pi:100*1j, 0:np.pi:100*1j]
ax.plot_surface(centre[0] + radius * np.cos(u) * np.sin(v),centre[1] + radius * np.sin(u) * np.sin(v),centre[2] + radius * np.cos(v), color="gray", alpha=0.6)

# Plot the local velocity vector (vt)
vt= A0 @ vel
ax.quiver(r0[0],r0[1],r0[2],vel[0],vel[1],vel[2],color='black')
ax.text(r0[0]+vel[0],r0[1]+vel[1],r0[2]+vel[2], 'v(t)', size=12, zorder=1, color='black')

# Plot the angular velocity (ω)
ω = A0 @ ω0
ax.quiver(r0[0],r0[1],r0[2],ω[0],ω[1],ω[2],color='orange',linestyle='--')
ax.text(r0[0]+ω[0],r0[1]+ω[1],r0[2]+ω[2], 'ω(t)', size=12, zorder=1,color='orange')

# build file name and append to list of file names
filename = f'KDE images/frame_{1000}.png'
filenames.append(filename)

# save img
plt.savefig(filename)
plt.close()

for i in range(len(At)):
    fig = plt.figure(2)
    ax = fig.add_subplot(111, projection='3d')

    dirx= At[i] @ [1,0,0]
    diry= At[i] @ [0,1,0]
    dirz= At[i] @ [0,0,1]

    ax.quiver(rt[i][0],rt[i][1],rt[i][2],dirx[0],dirx[1],dirx[2],color='r',linestyle='--')
    ax.quiver(rt[i][0],rt[i][1],rt[i][2],diry[0],diry[1],diry[2],color='g',linestyle='--')
    ax.quiver(rt[i][0],rt[i][1],rt[i][2],dirz[0],dirz[1],dirz[2],color='b',linestyle='--')
    #ax.text(rt[i][0]+dirx[0],rt[i][1]+dirx[1],rt[i][2]+dirx[2], 'local x', size=12, zorder=1, color='r')
    #ax.text(rt[i][0]+diry[0],rt[i][1]+diry[1],rt[i][2]+diry[2], 'local y', size=12, zorder=1, color='g')
    #ax.text(rt[i][0]+dirz[0],rt[i][1]+dirz[1],rt[i][2]+dirz[2], 'local z', size=12, zorder=1, color='b')
    
    # plot time/title
    t = round((i+1)*t_step, 3)  
    fig.suptitle('t = '+str(t)+'ms', fontsize=12)
    
    # draw sphere   
    centre = [rt[i][0],rt[i][1],rt[i][2]]
    radius = 1
    u, v = np.mgrid[0:2*np.pi:100*1j, 0:np.pi:100*1j]
    ax.plot_surface(centre[0] + radius * np.cos(u) * np.sin(v),centre[1] + radius * np.sin(u) * np.sin(v),centre[2] + radius * np.cos(v), color="gray", alpha=0.6)

    # Plot the local velocity vector (vt)
    vt= At[i] @ vel
    ax.quiver(rt[i][0],rt[i][1],rt[i][2],vt[0],vt[1],vt[2],color='black')
    ax.text(rt[i][0]+vt[0],rt[i][1]+vt[1],rt[i][2]+vt[2], 'v(t)', size=12, zorder=1, color='black')

    # Plot the angular velocity (ω)
    ω = At[i] @ ω0
    ax.quiver(rt[i][0],rt[i][1],rt[i][2],ω[0],ω[1],ω[2],color='orange',linestyle='--')
    ax.text(rt[i][0]+ω[0],rt[i][1]+ω[1],rt[i][2]+ω[2], 'ω(t)', size=12, zorder=1,color='orange')

    # build file name and append to list of file names
    filename = f'KDE images/frame_{1000+i+1}.png'
    filenames.append(filename)
    # last frame of each viz stays longer
    if (i == len(At)):
        for i in range(5):
            filenames.append(filename)

    # plot position history
    xs=[]
    ys=[]
    zs=[]
    for j in range(i):
        xs.append(rt[j][0])
        ys.append(rt[j][1])
        zs.append(rt[j][2])

    ax.scatter3D(xs, ys, zs, c=zs, cmap='Greys',marker='.')
    
    # Axis limits and lables and legend
    ax.set_xlim3d(-10,10)
    ax.set_ylim3d(0,20)
    ax.set_zlim3d(-10,10)

    ax.set_xlabel('Global X')
    ax.set_ylabel('Global Y')
    ax.set_zlabel('Global Z')

    ax.view_init(10, i-60)

    # save img
    plt.savefig(filename)
    plt.close()
# ----------------------------------------------------------------------------------------------------------------
