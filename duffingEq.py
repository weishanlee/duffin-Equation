# -*- coding: utf-8 -*-
"""
Duffing Equation:
ddx + delta * dx + alpha * x + beta * x**3 = gamma * cos(omega * t)

References:
https://en.wikipedia.org/wiki/Duffing_equation
"""
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator

#%% function to plot multiples of pi in x axis
# source:
# https://stackoverflow.com/questions/40642061/how-to-set-axis-ticks-in-multiples-of-pi-python-matplotlib
def multiple_formatter(denominator=2, number=np.pi, latex='\pi'):
    def gcd(a, b):
        while b:
            a, b = b, a%b
        return a
    def _multiple_formatter(x, pos):
        den = denominator
        num = np.int(np.rint(den*x/number))
        com = gcd(num,den)
        (num,den) = (int(num/com),int(den/com))
        if den==1:
            if num==0:
                return r'$0$'
            if num==1:
                return r'$%s$'%latex
            elif num==-1:
                return r'$-%s$'%latex
            else:
                return r'$%s%s$'%(num,latex)
        else:
            if num==1:
                return r'$\frac{%s}{%s}$'%(latex,den)
            elif num==-1:
                return r'$\frac{-%s}{%s}$'%(latex,den)
            else:
                return r'$\frac{%s%s}{%s}$'%(num,latex,den)
    return _multiple_formatter

#%% parameters and rk4 function
   
alpha = -1.0 #1.0 #
beta =  0.25 #5.0 #
gamma = 2.5 #8.0  #
delta = 0.1 #0.02 # 
omega = 2.0  #0.5 # 

#cConst = 0.05
#omega = 0.7
#FConst = 1.0

def f(r,t):
    x = r[0]
    y = r[1]
    z = r[2]
    dx = y
    dy = gamma * np.cos(z) - alpha * x - beta * x**3 - delta * y
    #dy = FConst * np.cos(z) - np.sin(x) - cConst * y
    dz = omega
    return np.array([dx,dy,dz],float)

#%% draw x vs t, v vs t, and the phase diagram.
t0 = 0.0
tn = 100
N = 1000
h = (tn-t0)/N
xIni = 0.5
yIni = 0
zIni = 0

tpoints = np.arange(t0,tn,h)

xpoints = []
ypoints = []
zpoints = []

r = np.array([xIni,yIni,zIni],float)

for t in tpoints:
    xpoints.append(r[0])
    ypoints.append(r[1])
    zpoints.append(r[2])

    k1 = h*f(r,t)
    k2 = h*f(r+0.5*k1,t+0.5*h)
    k3 = h*f(r+0.5*k2,t+0.5*h)
    k4 = h*f(r+k3,t+h)
    r += (k1+2*k2+2*k3+k4)/6

# x(t)
plt.figure("x(t)")
plt.title("x(t)") #(r'Multiples of $\pi$')
ax = plt.gca()
ax.set_xlabel("t",size = 16)
ax.set_ylabel("x(t)",size = 16)
plt.grid()
plt.plot(tpoints,xpoints,'k-')
plt.minorticks_on()
minorLocatorX = AutoMinorLocator(5) #number of minor intervals per major interval
minorLocatorY = AutoMinorLocator(5)
ax.xaxis.set_minor_locator(minorLocatorX) # add minor ticks on x axis
ax.yaxis.set_minor_locator(minorLocatorY) # add minor ticks on y axis
#plt.xlim((0,100))
#plt.ylim((-300,100))
plt.savefig("xt.png")
plt.show()

# dx(t)
plt.figure("$\\dot{x}$")
plt.title(r'Plot of $\dot{x}$ vs t')
ax = plt.gca()
ax.set_xlabel("t",size = 16)
ax.set_ylabel("$\\dot{x}$",size = 16)
plt.grid()    
plt.plot(tpoints,ypoints,'k-')
plt.minorticks_on()
minorLocatorX = AutoMinorLocator(5) #number of minor intervals per major interval
minorLocatorY = AutoMinorLocator(5)
ax.xaxis.set_minor_locator(minorLocatorX) # add minor ticks on x axis
ax.yaxis.set_minor_locator(minorLocatorY) # add minor ticks on y axis
#plt.xlim((0, 500))
plt.savefig("dx.png")
plt.show()

# phase diagrams
plt.figure("phase diagram")
plt.title("phase diagram")
ax = plt.gca()
ax.set_xlabel("x",size = 16)
ax.set_ylabel("$\\dot{x}$",size = 16)
plt.grid()    
plt.plot(xpoints,ypoints,'k-')

ax.xaxis.set_major_locator(plt.MultipleLocator(2*np.pi))
ax.xaxis.set_minor_locator(plt.MultipleLocator(np.pi))
ax.xaxis.set_major_formatter(plt.FuncFormatter(multiple_formatter()))
plt.xlim((-2*np.pi, 2*np.pi))
plt.savefig("phaseDiagram.png")
plt.show()

# 3d plot
from mpl_toolkits.mplot3d import Axes3D
fig = plt.figure("3D plot")
ax = fig.add_subplot(111, projection='3d')
#ax.view_init(90, 270) #elev, azim 
                    # (0,0) -> YZ; (0,-90) -> XZ; (90, 270) -> XY 
ax.plot(xpoints, ypoints, zpoints,color='k')

ax.set_xlabel('X')
ax.set_ylabel("$\\dot{x}$")
ax.set_zlabel('Z=$\\omega$ t')

ax.xaxis.set_major_locator(plt.MultipleLocator(2*np.pi))
ax.xaxis.set_minor_locator(plt.MultipleLocator(np.pi))
ax.xaxis.set_major_formatter(plt.FuncFormatter(multiple_formatter()))
plt.xlim((-3*np.pi, 3*np.pi))
plt.savefig("3DPlot.png")
plt.show()
#%% Poincare Section

# method 1: find out the points in the poincare section
t0 = 0.0
tn = 30000
N = 300000
h = (tn-t0)/N
xIni = 0.5
yIni = 0
zIni = 0

tpoints = np.arange(t0,tn,h)

xpoints = []
ypoints = []
zpoints = []

r = np.array([xIni,yIni,zIni],float)

for t in tpoints:
    xpoints.append(r[0])
    ypoints.append(r[1])
    zpoints.append(r[2])

    k1 = h*f(r,t)
    k2 = h*f(r+0.5*k1,t+0.5*h)
    k3 = h*f(r+0.5*k2,t+0.5*h)
    k4 = h*f(r+k3,t+h)
    r += (k1+2*k2+2*k3+k4)/6

xPS = [] # x values for Poicare Section Plot
yPS = [] # y values for Poicare Section Plot
colorList = []    
k = 0

for xv, yv, zv in zip(xpoints,ypoints,zpoints):
    if zv % (2*np.pi) < h or -zv % (2*np.pi) <h:
        xPS.append(xv)
        yPS.append(yv)
        k += 1
        colorList.append( k*(2*np.pi) )        
"""
# method 2: find out the points in the poincare section    
t0 = 0.0
tn = 30000 #int(1e8)
N = 300000 #int(1e9)
h = (tn-t0)/N

xPS = [] # x values for Poicare Section Plot
yPS = [] # y values for Poicare Section Plot
#colorList = []    
k = 0

ttpoints = np.arange(t0,tn,h)
r = np.array([xIni,yIni,zIni],float)
for t in ttpoints:
    if abs( t-k*(2*np.pi/omega) ) < h:
        xPS.append(r[0])
        yPS.append(r[1])
        k += 1
        #colorList.append( k*(2*np.pi) )
    k1 = h*f(r,t)
    k2 = h*f(r+0.5*k1,t+0.5*h)
    k3 = h*f(r+0.5*k2,t+0.5*h)
    k4 = h*f(r+k3,t+h)
    r += (k1+2*k2+2*k3+k4)/6 
"""
Log = open("Log.txt","w")
Log.write("Number of points in Poicare Section: {}".format(k))
Log.close()
print("Number of points in Poicare Section: {}".format(k))


plt.figure("Poincare Section")
plt.title(r'Poincar$\'e$ Section')
ax = plt.gca()
ax.set_xlabel("x",size = 16)
ax.set_ylabel("$\\dot{x}$",size = 16)
plt.grid()    

plt.plot(xPS,yPS,'k.')
"""
#to plot varicolored poincare section
plt.scatter(xPS, yPS, c=colorList, cmap='gist_rainbow', marker='.')
plt.clim( vmin = min(colorList), vmax = max(colorList) )  
cbar = plt.colorbar(ax=ax)
cbar.set_label(r'z =$\omega$ t',rotation=270,labelpad=30)        
cbar.ax.minorticks_on()
minorLocatorX = AutoMinorLocator(5)
cbar.ax.xaxis.set_minor_locator(minorLocatorX) # add minor ticks on x axis 
"""
ax.xaxis.set_major_locator(plt.MultipleLocator(2*np.pi))
ax.xaxis.set_minor_locator(plt.MultipleLocator(np.pi))
ax.xaxis.set_major_formatter(plt.FuncFormatter(multiple_formatter()))
plt.xlim((-2*np.pi, 2*np.pi))

plt.savefig("PS.png")
plt.show()

#Plot varicolored poincare section
plt.figure("Poincare Section Varicolored")
plt.title(r'Poincar$\'e$ Section')
ax = plt.gca()
ax.set_xlabel("x",size = 16)
ax.set_ylabel("$\\dot{x}$",size = 16)
plt.grid()    


plt.scatter(xPS, yPS, c=colorList, cmap='gist_rainbow', marker='.')
plt.clim( vmin = min(colorList), vmax = max(colorList) )  
cbar = plt.colorbar(ax=ax)
cbar.set_label(r'z =$\omega$ t',rotation=270,labelpad=30)        
cbar.ax.minorticks_on()
minorLocatorX = AutoMinorLocator(5)
cbar.ax.xaxis.set_minor_locator(minorLocatorX) # add minor ticks on x axis 

ax.xaxis.set_major_locator(plt.MultipleLocator(2*np.pi))
ax.xaxis.set_minor_locator(plt.MultipleLocator(np.pi))
ax.xaxis.set_major_formatter(plt.FuncFormatter(multiple_formatter()))
plt.xlim((-2*np.pi, 2*np.pi))

plt.savefig("PS.png")
plt.show()