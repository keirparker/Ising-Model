"""
Modelling and Visualisation- Checkpoint 1
"""
#import libs
import numpy as np
from numpy.random import rand
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from itertools import chain

#Ask for inputs for initial conditions

N = int(raw_input("please input system size: "))
T = float(raw_input("please input temperature: "))
code = raw_input("please type G or K (for Glauber or Kawasaki): ")
animations = raw_input("Animate? (y/n): ")



sweep = N**2


def initialstate(N):   
    ''' generates a random spin configuration for initial condition'''
    state = 2*np.random.randint(2, size=(N,N))-1
    return state



#glauber dynamics
def mcmove(config,T,N):
    '''Monte Carlo move using Metropolis algorithm '''
    for i in range(N):
        for j in range(N):
                a = np.random.randint(0, N-1)
                b = np.random.randint(0, N-1)
                s =  config[a, b]
                nb = config[(a+1)%(N-1),b] + config[a,(b+1)%(N-1)] + config[(a-1)%(N-1),b] + config[a,(b-1)%(N-1)]
                cost = 2*s*nb
                if cost < 0:
                    s *= -1
                elif rand() < np.exp(-cost*(1/T)):
                    s *= -1
                config[a, b] = s
    return config

#kawasaki dynamics
def mcmoveK(config,T,N):
    '''Monte Carlo move using Metropolis algorithm '''
    for i in range(N):
        for j in range(N):
                a = np.random.randint(0, N-1)
                b = np.random.randint(0, N-1)
                c = np.random.randint(0, N-1)
                d = np.random.randint(0, N-1)
                s1 =  config[a, b]  
                s2 =  config[c, d]  
                
                #if not the same
                if s1!=s2:
                    #random spin config inside matrix
                    #calculate the cost
                    nb1 = config[(a+1)%(N-1),b] + config[a,(b+1)%(N-1)] + config[(a-1)%(N-1),b] + config[a,(b-1)%(N-1)]
                    nb2 = config[(c+1)%(N-1),d] + config[c,(d+1)%(N-1)] + config[(c-1)%(N-1),d] + config[c,(d-1)%(N-1)]
                    cost = (2*s1*nb1 + 2*s2*nb2)
                    
                    if [a, b]== [(c+1)%(N-1),d] or [c,(d+1)%(N-1)] or [(c-1)%(N-1),d] or [c,(d-1)%(N-1)]:
                        cost = (2*s1*nb1 + 2*s2*nb2) +4  
                    
                    #SWAP
                    if cost < 0:
                        m = s1
                        n = s2
                        s2 = m
                        s1 = n
                    elif rand() < np.exp(-cost*(1/T)): 
                        m = s1
                        n = s2
                        s2 = m 
                        s1 = n
                    config[a, b] = s1 
                    config[c, d] = s2 
    return config
#function for animation
def updateim(*args):
    mcmove(config,T,N)
    im = plt.imshow(config,cmap= "CMRmap",animated = True)
    return im,

#function to calculate the average mag
def calcMag(config):
    '''Magnetization of a given configuration'''
    mag = np.sum(config)
    return mag
#function for susceptibility
def sus(mlist,T):
    s=(np.mean(map(lambda x: x**2,maglist))-(np.mean(maglist))**2)/(T*N)
    return s


if code == "K":
    mcmove = mcmoveK

#INITIALIZE
config = initialstate(N)

if animations== "y":

    fig = plt.figure()
    ani = animation.FuncAnimation(fig, updateim, interval=500, blit=True)
    plt.show()



Tlist = np.arange(1,4,0.1)

#SUSCEPTIBILITY
if code == "G":
    avgmaglist = []
    suslist = []
    print Tlist
    p=0
    for j in Tlist:
        maglist = []
        for i in range(10500): # 10000 sweeps, missing the first 100
            mcmove(config,j,N)
            if i >=500:
                if i % 10 == 0:
                    maglist.append(calcMag(config))
                    if i % 100==0:
                        p+=(0.033)
                        print str(p)+"%"
        suslist.append(sus(maglist,j))
        avgmaglist.append(np.mean(maglist))

    suslist = (map(lambda x:x/N,suslist))
    avgmaglist = (map(lambda x:x/N,avgmaglist))
    
    plt.figure()
    plt.plot(Tlist,suslist)
    plt.ylabel("susceptibility N$^-1$")
    plt.xlabel("Temperature")
    plt.savefig("SUS VS TEMP %s" % str(code))
    plt.title("SUS VS TEMP %s" % str(code))
    
    plt.figure()
    plt.plot(Tlist,avgmaglist)
    plt.ylabel("Mag AVG")
    plt.xlabel("Temperature")
    plt.savefig("MAG VS TEMP %s" % str(code))
    plt.title("MAG VS TEMP %s" % str(code))

#SUSCEPTIBILITY K
if code == "K":
    avgmaglist = []
    suslist = []
    print Tlist
    p=0
    for j in Tlist:
        maglist = []
        for i in range(10500): # 10000 sweeps, missing the first 100
            mcmoveK(config,j,N)
            if i >=500:
                if i % 10 == 0:
                    maglist.append(calcMag(config))
                    if i % 100==0:
                        p+=(0.033)
                        print str(p)+"%"
        suslist.append(sus(maglist,j))
        avgmaglist.append(np.mean(maglist))  

    suslist = (map(lambda x:x/N,suslist))
    avgmaglist = (map(lambda x:x/N,avgmaglist))

    plt.figure()
    plt.plot(Tlist,suslist)
    plt.ylabel("susceptibility N$^-1$")
    plt.xlabel("Temperature")
    plt.savefig("SUS VS TEMP %s" % str(code))
    plt.title("SUS VS TEMP %s" % str(code))    
    plt.figure()
    plt.plot(Tlist,avgmaglist)
    plt.ylabel("Mag AVG N$^-1$")
    plt.xlabel("Temperature")
    plt.savefig("MAG VS TEMP %s" % str(code))
    plt.title("MAG VS TEMP %s" % str(code))
    

def energy(config):
    energy = 0
    for i in range(N):
        for j  in range(N):
            S = config[i,j]
            nb = config[(i+1)%N,j] + config[i,(j+1)%N] + config[(i-1)%N,j] + config[i,(j-1)%N]
            energy += -1*nb*S
    return energy/4 #################DIVIDE BY 4#####################

def shc(energylist,T):
    SHC=(np.mean(map(lambda x: x**2,energylist))-(np.mean(energylist))**2)/((T**2)*N)
    return SHC




#ENERGY GLAUBER
if code == "G":
    avgenergylist = []
    shclist = []
    cerror1list = []
    print Tlist
    p=0
    for j in Tlist:
        energylist = []
        cerrorlist = []
        for i in range(10500): # 10000 sweeps, missing the first 100
            mcmove(config,j,N)
            if i >=500:
                if i % 10 == 0:
                    energylist.append(energy(config))
                    
                    if i % 100==0:
                        p+=(0.03)
                        print str(p)
        shclist.append(shc(energylist,j))
        avgenergylist.append(np.mean(energylist))
#errors
        for z in range(len(energylist)):
            y = shc(energylist,j)
            value = energylist.pop(z)
            x = shc(energylist,j)
            energylist.insert(z,value)
            cerrorlist.append((x-y)**2)
            
        cerror1list.append(np.sqrt(np.sum(cerrorlist)))
    
    cerror1list = (map(lambda x:x/N,cerror1list))
    avgenergylist = (map(lambda x:x/N,avgenergylist))
    shclist = (map(lambda x:x/N,shclist))


    outF = open("Kawasaki.txt", "w")

    outF.write("\nTemperature\n")
    for line in Tlist:
        outF.write("%s"%line)
        outF.write("  ")
    outF.write("\nAverage Energy\n")
    for line in avgenergylist :
        outF.write("%s"%line)
        outF.write("  ")
    outF.write("\nErrors\n")
    for line in cerror1list:
        outF.write("%s"%line)
        outF.write("  ")
    outF.write("\nHeatCapacity\n")
    for line in (map(lambda x:x/N,shclist)) :
        outF.write("%s"%line)
        outF.write("  ")  
    outF.write("\nsusceptibility\n")
    for line in (suslist) :
        outF.write("%s"%line)
        outF.write("  ")
    outF.close() 

    
    plt.figure()
    plt.errorbar(Tlist,shclist, yerr = cerror1list)
    plt.ylabel("Specific Heat Capacity N$^-1$")
    plt.xlabel("Temperature")
    plt.savefig("normalised by N SHC VS TEMP %s" % str(code))
    plt.title("normalised by N SHC VS TEMP %s" % str(code))
    plt.figure()
    plt.plot(Tlist,avgenergylist)
    plt.ylabel("Energy AVG")
    plt.xlabel("Temperature")
    plt.savefig("AVG ENERGY VS TEMP %s" % str(code))
    plt.title("AVG ENERGY VS TEMP %s" % str(code))


"""
100 sweeps for equilibrium time
10 sweeps for the auto correlation time

errors:

Magnitization ===> std dev of mean

specific heat and susceptibitly ===>  Resampling:
    
    Resampling looks at the same data many times:
        two methods : 1. Bootstrap(algorithim_stochastic):
                          do n = 1000 measurements (store the results) at a Temperature
                          choose n measurements from 1+n randomly ie. can have multiple of the same value in the set
                          compute avg of value from the set of randomly picked values
                          repeat k times
                          the error is then delta_c = sqrt(avg(c^2) - avg(c)^2)

    
"""


#ENERGY KAWASAKI
if code == "K":
    avgenergylist = []
    shclist = []
    cerror1list = []
    #print Tlist
    p=0
    for j in Tlist:
        energylist = []
        cerrorlist = []
        for i in range(10500): # 10000 sweeps, missing the first 100
            mcmoveK(config,j,N)
            if i >=500:
                if i % 10 == 0:
                    energylist.append(energy(config))

                    if i % 100==0:
                        p+=(0.03)
                        print str(p)
        shclist.append(shc(energylist,j))
        avgenergylist.append(np.mean(energylist))
#errors
        for z in range(len(energylist)):
            y = shc(energylist,j)
            value = energylist.pop(z)
            x = shc(energylist,j)
            energylist.insert(z,value)
            cerrorlist.append((x-y)**2)
            
        cerror1list.append(np.sqrt(np.sum(cerrorlist)))

    cerror1list = (map(lambda x:x/N,cerror1list))
    avgenergylist = (map(lambda x:x/N,avgenergylist))
    shclist = (map(lambda x:x/N,shclist))

#write to file
        
    textList = ["Temperature", "Average Energy", "SpecHeatCapacity Errors", "SpecHeatCapacity"]
    outF = open("GLAUBER.txt", "w")

    outF.write("\nTemperature\n")
    for line in Tlist:
        outF.write("%s"%line)
        outF.write("  ")
    outF.write("\nAverage Energy\n")
    for line in avgenergylist :
        outF.write("%s"%line)
        outF.write("  ")
    outF.write("\nErrors\n")
    for line in cerror1list:
        outF.write("%s"%line)
        outF.write("  ")
    outF.write("\nHeatCapacity\n")
    for line in (map(lambda x:x/N,shclist)) :
        outF.write("%s"%line)
        outF.write("  ")  
    outF.close() 

#plot

    plt.figure()
    print cerror1list
    plt.errorbar(Tlist,shclist,yerr=cerror1list)
    plt.ylabel("Specific Heat Capacity N$^-1$")
    plt.xlabel("Temperature")
    plt.savefig("normalised by N SHC VS TEMP %s" % str(code))
    plt.title("normalised by N SHC VS TEMP %s" % str(code))  
    plt.figure()
    plt.plot(Tlist,avgenergylist)
    plt.ylabel("Energy AVG")
    plt.xlabel("Temperature")
    plt.savefig("AVG ENERGY VS TEMP %s" % str(code))
    plt.title("AVG ENERGY VS TEMP %s" % str(code))





