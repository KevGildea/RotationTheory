import matplotlib.pyplot as plt
import numpy as np
import mpl_toolkits.mplot3d as plt3d


#### 4B17 Week 4 tutorial – Kinematic differential equations
# ----------------------------------------------------------------------------------------------------------------

# Function for performing Euler integration on the rotation matrix
def EulerInt_A(A0,ω,t_step,t_end):
    At=A0
    Ats=[]
    ω_skewsym = np.array([[0, -ω[2], ω[1]],
                         [ω[2], 0,-ω[0]],
                         [-ω[1], ω[0], 0]])
    for t in range(int(t_end/t_step)):
        At = At + t_step*(-ω_skewsym @ At)
        Ats.append(At)
    return Ats

############## Q1 ##############

# Initial conditions
t_end = 0.01
t_step = 0.0005
A0 = np.array([[1, 0, 0],
               [0, 1, 0],
               [0, 0, 1]])
ω = np.array([0.1, 4, 0.2])

# # Euler integration
At = EulerInt_A(A0,ω,t_step,t_end)

fig = plt.figure(2)
ax = fig.add_subplot(111, projection='3d')
fig.suptitle('4B17 Week 4 tutorial – Kinematic differential equations: Q1', fontsize=12)

# Plot the initial coordinate system
dirxinitial=A0 @ [1,0,0]
diryinitial=A0 @ [0,1,0]
dirzinitial=A0 @ [0,0,1]

ax.quiver(0,0,0,dirxinitial[0],dirxinitial[1],dirxinitial[2],color='r',linestyle='--')
ax.quiver(0,0,0,diryinitial[0],diryinitial[1],diryinitial[2],color='g',linestyle='--')
ax.quiver(0,0,0,dirzinitial[0],dirzinitial[1],dirzinitial[2],color='b',linestyle='--')
ax.text(dirxinitial[0],dirxinitial[1],dirxinitial[2], '  t = 0ms', size=12, zorder=1, color='r')
ax.text(diryinitial[0],diryinitial[1],diryinitial[2], '  t = 0ms', size=12, zorder=1, color='g')
ax.text(dirzinitial[0],dirzinitial[1],dirzinitial[2], '  t = 0ms', size=12, zorder=1, color='b')

# plot the orientation after each timestep up to the end
for i in range(0,len(At)):
    dirx= At[i] @ [1,0,0]
    diry= At[i] @ [0,1,0]
    dirz= At[i] @ [0,0,1]

    ax.quiver(0,0,0,dirx[0],dirx[1],dirx[2],color='r', alpha=0.3)
    ax.quiver(0,0,0,diry[0],diry[1],diry[2],color='g', alpha=0.3)
    ax.quiver(0,0,0,dirz[0],dirz[1],dirz[2],color='b', alpha=0.3)
    t = (i+1)*t_step
    ax.text(dirx[0],dirx[1],dirx[2], '  t = '+str(t)+'ms', size=12, zorder=1, color='r')
    ax.text(diry[0],diry[1],diry[2], '  t = '+str(t)+'ms', size=12, zorder=1, color='g')
    ax.text(dirz[0],dirz[1],dirz[2], '  t = '+str(t)+'ms', size=12, zorder=1, color='b')

# Axis limits and lables and legend
ax.set_xlim3d(-2,2)
ax.set_ylim3d(-2,2)
ax.set_zlim3d(-2,2)

ax.set_xlabel('Global X')
ax.set_ylabel('Global Y')
ax.set_zlabel('Global Z')

print('Q1 answer', At[-1])
plt.show()
# ----------------------------------------------------------------------------------------------------------------

############## Q2 ##############

# Initial conditions
t_end = 0.005
t_step = 0.005
A0 = np.array([[0.8847, 0.2557, -0.3898],
               [-0.2060, 0.9645, 0.1651],
              [0.4182, -0.0658, 0.9060]])
ω =  A0 @ np.array([1, 7, 2]) # Must convert to local

# # Euler integration
At = EulerInt_A(A0,ω,t_step,t_end)

fig = plt.figure(3)
ax = fig.add_subplot(111, projection='3d')
fig.suptitle('4B17 Week 4 tutorial – Kinematic differential equations: Q2', fontsize=12)

# Plot the initial coordinate system
dirxinitial=A0 @ [1,0,0]
diryinitial=A0 @ [0,1,0]
dirzinitial=A0 @ [0,0,1]

ax.quiver(0,0,0,dirxinitial[0],dirxinitial[1],dirxinitial[2],color='r',linestyle='--')
ax.quiver(0,0,0,diryinitial[0],diryinitial[1],diryinitial[2],color='g',linestyle='--')
ax.quiver(0,0,0,dirzinitial[0],dirzinitial[1],dirzinitial[2],color='b',linestyle='--')
ax.text(dirxinitial[0],dirxinitial[1],dirxinitial[2], '  t = 0ms', size=12, zorder=1, color='r')
ax.text(diryinitial[0],diryinitial[1],diryinitial[2], '  t = 0ms', size=12, zorder=1, color='g')
ax.text(dirzinitial[0],dirzinitial[1],dirzinitial[2], '  t = 0ms', size=12, zorder=1, color='b')

# plot the orientation after each timestep up to the end
for i in range(0,len(At)):
    dirx= At[i] @ [1,0,0]
    diry= At[i] @ [0,1,0]
    dirz= At[i] @ [0,0,1]

    ax.quiver(0,0,0,dirx[0],dirx[1],dirx[2],color='r', alpha=0.3)
    ax.quiver(0,0,0,diry[0],diry[1],diry[2],color='g', alpha=0.3)
    ax.quiver(0,0,0,dirz[0],dirz[1],dirz[2],color='b', alpha=0.3)
    ax.text(dirx[0],dirx[1],dirx[2], '  i = '+str(i+1), size=12, zorder=1, color='r')
    ax.text(diry[0],diry[1],diry[2], '  i = '+str(i+1), size=12, zorder=1, color='g')
    ax.text(dirz[0],dirz[1],dirz[2], '  i = '+str(i+1), size=12, zorder=1, color='b')
    t = (i+1)*t_step
    ax.text(dirx[0],dirx[1],dirx[2], '  t = '+str(t)+'ms', size=12, zorder=1, color='r')
    ax.text(diry[0],diry[1],diry[2], '  t = '+str(t)+'ms', size=12, zorder=1, color='g')
    ax.text(dirz[0],dirz[1],dirz[2], '  t = '+str(t)+'ms', size=12, zorder=1, color='b')

# Axis limits and lables and legend
ax.set_xlim3d(-2,2)
ax.set_ylim3d(-2,2)
ax.set_zlim3d(-2,2)

ax.set_xlabel('Global X')
ax.set_ylabel('Global Y')
ax.set_zlabel('Global Z')

print('Q2 answer', At[-1])
plt.show()

# ----------------------------------------------------------------------------------------------------------------

############## Q3 ##############
print('#### Question 3 ####')

# Initial conditions
t_end = 0.01
t_step = 0.005
A0 = np.array([[1, 0, 0],
               [0, 1, 0],
               [0, 0, 1]])
ω = np.array([0.5,-5,1])

# # Euler integration
At = EulerInt_A(A0,ω,t_step,t_end)

fig = plt.figure(4)
ax = fig.add_subplot(111, projection='3d')
fig.suptitle('4B17 Week 4 tutorial – Kinematic differential equations: Q3', fontsize=12)

# Plot the initial coordinate system
dirxinitial=A0 @ [1,0,0]
diryinitial=A0 @ [0,1,0]
dirzinitial=A0 @ [0,0,1]

ax.quiver(0,0,0,dirxinitial[0],dirxinitial[1],dirxinitial[2],color='r',linestyle='--')
ax.quiver(0,0,0,diryinitial[0],diryinitial[1],diryinitial[2],color='g',linestyle='--')
ax.quiver(0,0,0,dirzinitial[0],dirzinitial[1],dirzinitial[2],color='b',linestyle='--')
ax.text(dirxinitial[0],dirxinitial[1],dirxinitial[2], '  t = 0ms', size=12, zorder=1, color='r')
ax.text(diryinitial[0],diryinitial[1],diryinitial[2], '  t = 0ms', size=12, zorder=1, color='g')
ax.text(dirzinitial[0],dirzinitial[1],dirzinitial[2], '  t = 0ms', size=12, zorder=1, color='b')

# plot the orientation after each timestep up to the end
for i in range(0,len(At)):
    dirx= At[i] @ [1,0,0]
    diry= At[i] @ [0,1,0]
    dirz= At[i] @ [0,0,1]

    ax.quiver(0,0,0,dirx[0],dirx[1],dirx[2],color='r', alpha=0.3)
    ax.quiver(0,0,0,diry[0],diry[1],diry[2],color='g', alpha=0.3)
    ax.quiver(0,0,0,dirz[0],dirz[1],dirz[2],color='b', alpha=0.3)
    t = (i+1)*t_step
    ax.text(dirx[0],dirx[1],dirx[2], '  t = '+str(t)+'ms', size=12, zorder=1, color='r')
    ax.text(diry[0],diry[1],diry[2], '  t = '+str(t)+'ms', size=12, zorder=1, color='g')
    ax.text(dirz[0],dirz[1],dirz[2], '  t = '+str(t)+'ms', size=12, zorder=1, color='b')

# Axis limits and lables and legend
ax.set_xlim3d(-2,2)
ax.set_ylim3d(-2,2)
ax.set_zlim3d(-2,2)

ax.set_xlabel('Global X')
ax.set_ylabel('Global Y')
ax.set_zlabel('Global Z')

print('Q3 answer', At[-1])
plt.show()
# ----------------------------------------------------------------------------------------------------------------


#### Changes to parameters in Q1 ####


# ----------------------------------------------------------------------------------------------------------------

#### t_end = 0.2  ####

# Initial conditions
t_end = 0.2
t_step = 0.005
A0 = np.array([[1, 0, 0],
               [0, 1, 0],
               [0, 0, 1]])
ω = np.array([0.1, 4, 0.2])

# # Euler integration
At = EulerInt_A(A0,ω,t_step,t_end)

fig = plt.figure(5)
ax = fig.add_subplot(111, projection='3d')
fig.suptitle('4B17 Week 4 tutorial – Kinematic differential equations: Q1 (t_end = 0.2)', fontsize=12)

# Plot the initial coordinate system
dirxinitial=A0 @ [1,0,0]
diryinitial=A0 @ [0,1,0]
dirzinitial=A0 @ [0,0,1]

ax.quiver(0,0,0,dirxinitial[0],dirxinitial[1],dirxinitial[2],color='r',linestyle='--')
ax.quiver(0,0,0,diryinitial[0],diryinitial[1],diryinitial[2],color='g',linestyle='--')
ax.quiver(0,0,0,dirzinitial[0],dirzinitial[1],dirzinitial[2],color='b',linestyle='--')
ax.text(dirxinitial[0],dirxinitial[1],dirxinitial[2], '  t = 0ms', size=12, zorder=1, color='r')
ax.text(diryinitial[0],diryinitial[1],diryinitial[2], '  t = 0ms', size=12, zorder=1, color='g')
ax.text(dirzinitial[0],dirzinitial[1],dirzinitial[2], '  t = 0ms', size=12, zorder=1, color='b')

# plot the orientation after each timestep up to the end
for i in range(0,len(At)):
    dirx= At[i] @ [1,0,0]
    diry= At[i] @ [0,1,0]
    dirz= At[i] @ [0,0,1]

    ax.quiver(0,0,0,dirx[0],dirx[1],dirx[2],color='r', alpha=0.3)
    ax.quiver(0,0,0,diry[0],diry[1],diry[2],color='g', alpha=0.3)
    ax.quiver(0,0,0,dirz[0],dirz[1],dirz[2],color='b', alpha=0.3)
    t = (i+1)*t_step
    ax.text(dirx[0],dirx[1],dirx[2], '  t = '+str(t)+'ms', size=12, zorder=1, color='r')
    ax.text(diry[0],diry[1],diry[2], '  t = '+str(t)+'ms', size=12, zorder=1, color='g')
    ax.text(dirz[0],dirz[1],dirz[2], '  t = '+str(t)+'ms', size=12, zorder=1, color='b')

# Axis limits and lables and legend
ax.set_xlim3d(-2,2)
ax.set_ylim3d(-2,2)
ax.set_zlim3d(-2,2)

ax.set_xlabel('Global X') 
ax.set_ylabel('Global Y')
ax.set_zlabel('Global Z')

plt.show()


# ----------------------------------------------------------------------------------------------------------------

#### t_end = 0.2, t_step = 0.05 ####

# Initial conditions
t_end = 0.2
t_step = 0.05
A0 = np.array([[1, 0, 0],
               [0, 1, 0],
               [0, 0, 1]])
ω = np.array([0.1, 4, 0.2])

# # Euler integration
At = EulerInt_A(A0,ω,t_step,t_end)

fig = plt.figure(6)
ax = fig.add_subplot(111, projection='3d')
fig.suptitle('4B17 Week 4 tutorial – Kinematic differential equations: Q1 (t_end = 0.2, t_step = 0.05)', fontsize=12)

# Plot the initial coordinate system
dirxinitial=A0 @ [1,0,0]
diryinitial=A0 @ [0,1,0]
dirzinitial=A0 @ [0,0,1]

ax.quiver(0,0,0,dirxinitial[0],dirxinitial[1],dirxinitial[2],color='r',linestyle='--')
ax.quiver(0,0,0,diryinitial[0],diryinitial[1],diryinitial[2],color='g',linestyle='--')
ax.quiver(0,0,0,dirzinitial[0],dirzinitial[1],dirzinitial[2],color='b',linestyle='--')
ax.text(dirxinitial[0],dirxinitial[1],dirxinitial[2], '  t = 0ms', size=12, zorder=1, color='r')
ax.text(diryinitial[0],diryinitial[1],diryinitial[2], '  t = 0ms', size=12, zorder=1, color='g')
ax.text(dirzinitial[0],dirzinitial[1],dirzinitial[2], '  t = 0ms', size=12, zorder=1, color='b')

# plot the orientation after each timestep up to the end
for i in range(0,len(At)):
    dirx= At[i] @ [1,0,0]
    diry= At[i] @ [0,1,0]
    dirz= At[i] @ [0,0,1]

    ax.quiver(0,0,0,dirx[0],dirx[1],dirx[2],color='r', alpha=0.3)
    ax.quiver(0,0,0,diry[0],diry[1],diry[2],color='g', alpha=0.3)
    ax.quiver(0,0,0,dirz[0],dirz[1],dirz[2],color='b', alpha=0.3)
    t = (i+1)*t_step
    ax.text(dirx[0],dirx[1],dirx[2], '  t = '+str(t)+'ms', size=12, zorder=1, color='r')
    ax.text(diry[0],diry[1],diry[2], '  t = '+str(t)+'ms', size=12, zorder=1, color='g')
    ax.text(dirz[0],dirz[1],dirz[2], '  t = '+str(t)+'ms', size=12, zorder=1, color='b')

# Axis limits and lables and legend
ax.set_xlim3d(-2,2)
ax.set_ylim3d(-2,2)
ax.set_zlim3d(-2,2)

ax.set_xlabel('Global X')
ax.set_ylabel('Global Y')
ax.set_zlabel('Global Z')

plt.show()

