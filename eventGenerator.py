import math
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.pylab import *
import time
#define constants
alpha = [2.3e-8,9.3e-9,1e-11,2.29e-3,4.28e-7,0]
beta = [3,5,2.4,3,2.4,2,1]
k = [2.34e-3,1.4,3.2e-5,1.3e-3,3.2e-5,4e-5,1.25e-3]
#values below from channel 1, CCVR, Vignati
s = .015
R0 = 1.15
Rl = 6e10
T0 = 3.35
gamma = 0.5
Vb = 5
Cp = 5e-10
dur=1e-4

def cond(index,start,end):
    return ((k[index-1])*((start)**(beta[index-1]) - (end)**(beta[index-1])))

def Rntd(temp):
    return R0*math.exp((T0/temp)**gamma)

def phonon(a,b,c,d,e,f): #add NTD Power
    return ((Pntd/dur)+cond(1,d,a) + cond(2,b,a) - cond(3,a,s))/(alpha[0]*(a**3))

def electron(a,b,c,d,e,f): #figure out T
    return ((f**2 - ((Vb*Rntd(s))/(Rl+Rntd(s)))**2)/(Rntd(b)) - cond(2,b,a))/(alpha[1]*b)

def heater(a,b,c,d,e,f):
    return (cond(4,d,c) - cond(5,c,s))/(alpha[2]*c)

def crystal(a,b,c,d,e,f):
    return ((Pcrystal/dur)-cond(1,d,a) - cond(4,d,c) - cond(6,d,e))/(alpha[3]*(d**3))

def teflon(a,b,c,d,e,f):
    return (cond(6,d,e) - cond(7,e,s)) / ((alpha[4]*e) + (alpha[5]*(e**3)))

def feedback (a,b,c,d,e,f):
    return (Vb-(f*((Rl+Rntd(b))/(Rntd(b)))))/(Rl*Cp)

# fourth order Runge-Kutta method in 6 dimensions
def rK6(a, b, c, d, e, f, fa, fb, fc, fd, fe, ff, hs):
	a1 = fa(a, b, c, d, e, f)*hs
	b1 = fb(a, b, c, d, e, f)*hs
	c1 = fc(a, b, c, d, e, f)*hs
	d1 = fd(a, b, c, d, e, f)*hs
	e1 = fe(a, b, c, d, e, f)*hs
	f1 = ff(a, b, c, d, e, f)*hs
	ak = a + a1*0.5
	bk = b + b1*0.5
	ck = c + c1*0.5
	dk = d + d1*0.5
	ek = e + e1*0.5
	fk = f + f1*0.5
	a2 = fa(ak, bk, ck, dk, ek, fk)*hs
	b2 = fb(ak, bk, ck, dk, ek, fk)*hs
	c2 = fc(ak, bk, ck, dk, ek, fk)*hs
	d2 = fd(ak, bk, ck, dk, ek, fk)*hs
	e2 = fe(ak, bk, ck, dk, ek, fk)*hs
	f2 = ff(ak, bk, ck, dk, ek, fk)*hs
	ak = a + a2*0.5
	bk = b + b2*0.5
	ck = c + c2*0.5
	dk = d + d2*0.5
	ek = e + e2*0.5
	fk = f + f2*0.5
	a3 = fa(ak, bk, ck, dk, ek, fk)*hs
	b3 = fb(ak, bk, ck, dk, ek, fk)*hs
	c3 = fc(ak, bk, ck, dk, ek, fk)*hs
	d3 = fd(ak, bk, ck, dk, ek, fk)*hs
	e3 = fe(ak, bk, ck, dk, ek, fk)*hs
	f3 = ff(ak, bk, ck, dk, ek, fk)*hs
	ak = a + a3
	bk = b + b3
	ck = c + c3
	dk = d + d3
	ek = e + e3
	fk = f + f3
	a4 = fa(ak, bk, ck, dk, ek, fk)*hs
	b4 = fb(ak, bk, ck, dk, ek, fk)*hs
	c4 = fc(ak, bk, ck, dk, ek, fk)*hs
	d4 = fd(ak, bk, ck, dk, ek, fk)*hs
	e4 = fe(ak, bk, ck, dk, ek, fk)*hs
	f4 = ff(ak, bk, ck, dk, ek, fk)*hs
	a = a + (a1 + 2*(a2 + a3) + a4)/6
	b = b + (b1 + 2*(b2 + b3) + b4)/6
	c = c + (c1 + 2*(c2 + c3) + c4)/6
	d = d + (d1 + 2*(d2 + d3) + d4)/6
	e = e + (e1 + 2*(e2 + e3) + e4)/6
	f = f + (f1 + 2*(f2 + f3) + f4)/6
	return a, b, c, d, e, f

#run algorithm
ntdEventArray = []
crystalEventArray = []
timeSteps = int(9.3e6)
Pcrystal, Pntd = 0,0

def getTemps(eventType):
    global Pcrystal
    global Pntd
    global timeSteps
    TempArray = []
    if eventType == "ntd":
        Pcrystal = 0
        Pntd = 4.80653e-13
    else:
        Pcrystal = 4.80653e-13
        Pntd = 0
    a,b,c,d,e,f,hs = s,s,s,s,s,(Vb*Rntd(s))/(Rl+Rntd(s)),1e-8
    start = time.time()
    for i in range(timeSteps):
        if i>timeSteps/5000:
            Pcrystal=0
            Pntd = 0
        a, b, c, d, e, f = rK6(a, b, c, d, e, f, phonon, electron, heater, crystal, teflon, feedback, hs)
        TempArray.append([a,b,c,d,e,f])
    print eventType
    #print a,b,c,d,e
    end = time.time()
    print "Duration: " + str(end-start)
    TempArray = np.array(TempArray)
    return TempArray

ntdEventArray = getTemps("ntd")
crystalEventArray = getTemps("crystal")

t = linspace(0,5,timeSteps)

plt.plot(t,crystalEventArray[:,1])
plt.xlabel('Time (s)')
plt.ylabel('Voltage (Arb. Units)')
plt.title('Crystal Event Pulse')
crystalFig = plt.gcf()
crystalFig.set_size_inches(6,4)
crystalFig.savefig('crystal_event.png',dpi=150)
amplitude = np.amax(crystalEventArray[:,1])-s
maxTemp = np.argmax(crystalEventArray[:,1])
riseTime =  (5.0/(timeSteps))*len(zip(*np.where(np.logical_and((crystalEventArray[:maxTemp,1]-s)>=0.1*(amplitude),(crystalEventArray[:maxTemp,1]-s)<=0.9*(amplitude)))))
decayTime =  (5.0/(timeSteps))*len(zip(*np.where(np.logical_and((crystalEventArray[maxTemp:,1]-s)>=0.3*(amplitude),(crystalEventArray[maxTemp:,1]-s)<=0.9*(amplitude)))))

print "Amplitude: " + str(amplitude)
print "Rise Time: " + str(riseTime)
print "Decay Time: " + str(decayTime)
print "DecayTime/RiseTime: " + str(decayTime/riseTime)

plt.plot(t,ntdEventArray[:,1])
plt.xlabel('Time (s)')
plt.ylabel('Voltage (Arb. Units)')
plt.title('NTD Event Pulse')
ntdFig = plt.gcf()
ntdFig.set_size_inches(6,4)
ntdFig.savefig('ntd_event.png',dpi=150)
amplitude = np.amax(ntdEventArray[:,1])-s
maxTemp = np.argmax(ntdEventArray[:,1])
riseTime =  (5.0/(timeSteps))*len(zip(*np.where(np.logical_and((ntdEventArray[:maxTemp,1]-s)>=0.1*(amplitude),(ntdEventArray[:maxTemp,1]-s)<=0.9*(amplitude)))))
decayTime =  (5.0/(timeSteps))*len(zip(*np.where(np.logical_and((ntdEventArray[maxTemp:,1]-s)>=0.3*(amplitude),(ntdEventArray[maxTemp:,1]-s)<=0.9*(amplitude)))))
print "Amplitude: " + str(amplitude)
print "Rise Time: " + str(riseTime)
print "Decay Time: " + str(decayTime)
print "DecayTime/RiseTime: " + str(decayTime/riseTime)
