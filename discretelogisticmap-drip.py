import math
import matplotlib.pyplot as plt
import numpy as num

#parameters
rMin=0.06
rMax=0.2
rResolution=0.0001

#constants
K=1.0
b=0.5
rho=1.0
a=12.0
g=0.33
xc=1.0

#initial conditions
xN=0.0
vN=0.1
mN=0.1

def calculateMass(t):
	return mN+r*t

def calculateOmega(t):
	return math.sqrt(K/calculateMass(t))

def calculateGamma(t):
	return ((b+r)/calculateMass(t))

def calculateB():
	tB = xN-mN*g/K
	return tB

def calculateA():
	tA = ( vN+( (b+r)/(mN) )*calculateB()-(g*r)/(K) )/((math.sqrt(K/mN)))
	return tA

def calculateDOmega(t):
	return (calculateOmega(t)-(K*r*t)/(2.0*calculateMass(t)*calculateMass(t)*calculateOmega(t)))

def findPosition(t):
	omega=calculateOmega(t)
	gamma=calculateGamma(t)
	B=calculateB()
	A=calculateA()
	return (A*math.sin(omega*t)+B*math.cos(omega*t))*math.exp(-gamma*t) + calculateMass(t)*g/K

def findVelocity(t):
	omega=calculateOmega(t)
	gamma=calculateGamma(t)
	B=calculateB()
	A=calculateA()
	return math.exp(-calculateGamma(t)*t)*\
	(\
		calculateA()*(calculateDOmega(t)*math.cos(calculateOmega(t)*t))-\
		calculateB()*(calculateDOmega(t)*math.sin(calculateOmega(t)*t))+\
		((r*t*(b+r))/(calculateMass(t)**2.0)-(b+r)/(calculateMass(t)))*\
		(\
			calculateA()*math.sin(calculateOmega(t)*t)+\
			calculateB()*math.cos(calculateOmega(t)*t)
		)\
	) + g*r/K

def doNRIteration(current,target):
	return current - (findPosition(current)-target)/findVelocity(current)

def findTimeForPosition(target,iterations,start):
	current = start
	for i in range(0,iterations):
		current = doNRIteration(current,target)
	return current

def getNewX(t):
	vc = findVelocity(t)
	return xc - (((3.0*a*vc*calculateMass(t))/(4.0*math.pi*rho))**(1./3.))*(a*vc)

def getNewMass(t):
	return (1.0 - a*findVelocity(t))*calculateMass(t)

def interpolateTimeFromPosition(target,iterations):
	
	time = 0.0
	position = findPosition(0)
	lastPosition = position

	while (time < 1000):
		lastPosition = position
		position = findPosition(time)
		if((lastPosition < target and position > target) or (lastPosition > target and position < target)):
			candidate = findTimeForPosition(target,iterations,time)
			if (candidate > 0):
				return candidate 
		time+=0.1




xs = []
ys = []

r=rMin

#main loop
while (r<=rMax-rResolution):
	r+=rResolution
	for i in range(1,70):
		Tc = interpolateTimeFromPosition(xc,10)
		mC = getNewMass(Tc)
		vC = findVelocity(Tc)
		xC = getNewX(Tc)

		mN=mC
		vN=vC
		xN=xC

		if (i>20):
			xs.append(r)
			ys.append(Tc)


plt.scatter(xs,ys,0.1)
plt.xlabel('Mass Flow Rate')
plt.ylabel('Tc')
plt.savefig('bifurcation_r=[' + str(rMin) + '-' + str(rMax) + '].png')



