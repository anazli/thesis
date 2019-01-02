#!/usr/bin/env python3
from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from sys import argv
import numpy as np
import pylab as P 

script, string = argv

L = 7.E-7


def plot_particles():
    ic = np.loadtxt('init_coords.dat')
    lc = np.loadtxt('last_coords.dat')
    d = np.loadtxt('diameters.dat')
    r = 1.1 * d/2. 

    fig1 = plt.figure()
    fig2 = plt.figure()

    ax1 = fig1.add_subplot(111,projection='3d')
    ax2 = fig2.add_subplot(111,projection='3d')

    ax1.scatter(ic[:,0], ic[:,1], ic[:,2],s = np.pi*r**2*1E+18, color = 'blue', alpha = 0.75)
    ax2.scatter(lc[:,0], lc[:,1], lc[:,2],s = np.pi*r**2*1E+18, color = 'red', alpha = 0.75)
    ax1.set_title('Initial configuration of particles.')
    ax1.set_xlabel('x axis')
    ax1.set_ylabel('y axis')
    ax1.set_zlabel('z axis')
    ax1.set_xlim([0., L])
    ax1.set_ylim([0., L])
    ax1.set_zlim([0., L])
    ax2.set_title('Last configuration of particles.')
    ax2.set_xlabel('x axis')
    ax2.set_ylabel('y axis')
    ax2.set_zlabel('z axis')
    ax2.set_xlim([0., L])
    ax2.set_ylim([0., L])
    ax2.set_zlim([0., L])

    plt.show()


def plot_diameters():
    data = np.loadtxt('diameters.dat')
    data = data*1E+9
    plt.figure()
    #P.style.use('dark_background')
    #s = np.random.lognormal(8.9,0.34,10000)
    n,bins,patches = plt.hist(data,15,normed=True,align='mid',histtype='bar',color='blue',label='$N=1000$\n$\sigma=0.34$\n$D_m=8.9 nm$')
    x = np.linspace(min(bins), max(bins),10000)
    y = (np.exp(-(np.log(x/8.9)**2 )/(2*0.34**2))/(x*0.34*np.sqrt(2*np.pi))) 
    plt.plot(x, y,color='r',linewidth=2)
    plt.xlabel('Diameters (m)')
    plt.ylabel('Freq')
    plt.legend(loc=1)
    plt.axis('tight')
    plt.show()


def plot_energy():
    data = np.loadtxt('energy.dat')

    #plt.style.use('dark_background')
    plt.figure()

    plt.plot(data,'b-')
    plt.show()

def plot_magnetization():
    data = np.loadtxt('magnetization.dat')
    plt.figure()
    plt.plot(data,'b-')
    plt.show()

def plot_aggregates():
    data = np.loadtxt('aggregates.dat')
    plt.figure()
    ax = plt.gca()
    plt.plot(data[:,0],'b-',label='Monomers')
    plt.legend(loc=4)
    ax.set_ylim([0,1])
    plt.show()
    plt.figure()
    ax = plt.gca()
    plt.plot(data[:,1],'b-',label='Dimers')
    plt.legend(loc=4)
    ax.set_ylim([0,1])
    plt.show()
    plt.figure()
    ax = plt.gca()
    plt.plot(data[:,2],'b-',label='Agglomerates')
    plt.legend(loc=4)
    ax.set_ylim([0,1])
    plt.show()
   

if string == "particles" or string == "p":

    plot_particles()

elif string == "energy" or string == "e":
    
    plot_energy()

elif string == "diameter" or string == "d":
    
    plot_diameters()

elif string == "magnetization" or string == "m":

    plot_magnetization()

elif string == "aggregates" or string == "a":
    
    plot_aggregates()

else:
 
    print ("????........................")
    print ("What do you want me to plot?")


