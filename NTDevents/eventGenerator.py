import math
import numpy as np
import csv
import time
import sys
from decimal import Decimal
from prettytable import PrettyTable

#define constants
alpha = [2.7e-8,9.9e-9,1e-11,9.14e3,4.28e-7,0]
beta = [3,4.37,2.4,3,2.4,2,1]
k = [2.34e-3,0.7,3.2e-5,1.3e-3,3.2e-5,4e-5,1.25e-7]

s = .015
R0 = 1.15
Rl = 6e10
T0 = 3.35
gamma = 0.5
Vb = 5
Cp = 5e-10

#thermal model equations
def cond(index,start,end):
    return ((k[index-1]/(beta[index-1]+1))*((start)**((beta[index-1])+1) - (end)**((beta[index-1])+1)))
def Rntd(temp):
    return R0*math.exp((T0/temp)**gamma)
def phonon(a,b,c,d,e,f):
    return ((Entd/dur)+cond(1,d,a) + cond(2,b,a) - cond(3,a,s))/(alpha[0]*(a**3))
def electron(a,b,c,d,e,f):
    return ((f**2 - ((Vb*Rntd(s))/(Rl+Rntd(s)))**2)/(Rntd(b)) - cond(2,b,a))/(alpha[1]*b)
def heater(a,b,c,d,e,f):
    return (cond(4,d,c) - cond(5,c,s))/(alpha[2]*c)
def crystal(a,b,c,d,e,f):
    return ((Ecrystal/dur)-cond(1,d,a) - cond(4,d,c) - cond(6,d,e))/(alpha[3]*((d/232.0)**3))
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

#algorithm parameters
timeSteps = int(5e6)
stepSize = 1e-6
print "Total Time: " + str(stepSize*timeSteps) + " seconds" #Display how many seconds this will simulate

#power input parameters
Ecrystal, Entd = 0,0
eventEnergy = (2615.0)*(1.6022e-16) #2615 KeV
dur = 1e-3

def getConstants():
    t = PrettyTable(['Name','Value'])
    t.add_row([ "Capacitances ",""])
    t.add_row([ "Phonon" , '%.2E' % Decimal(alpha[0]*(s**3))])
    t.add_row([ "Electron" , '%.2E' % Decimal(alpha[1]*(s))])
    t.add_row([ "Heater" , '%.2E' % Decimal(alpha[2]*s)])
    t.add_row([ "Crystal" , '%.2E' % Decimal(alpha[3]*((s/232.0)**3))])
    t.add_row([ "Teflon" , '%.2E' % Decimal((alpha[4]*s) + (alpha[5]*(s**3)))])
    t.add_row([ "",""])
    t.add_row([ "Conductances ",""])
    t.add_row([ "NTD Glue" ,'%.2E' % Decimal(k[0]*(s**beta[0]))])
    t.add_row([ "EP Coupling",'%.2E' % Decimal(k[1]*(s**beta[1]))])
    t.add_row([ "NTD Gold",'%.2E' % Decimal(k[2]*(s**beta[2]))])
    t.add_row([ "Heater Glue",'%.2E' % Decimal(k[3]*(s**beta[3]))])
    t.add_row([ "Heater Gold",'%.2E' % Decimal(k[4]*(s**beta[4]))])
    t.add_row([ "Crystal-Teflon",'%.2E' % Decimal(k[5]*(s**beta[5]))])
    t.add_row([ "Teflon-Sink",'%.2E' % Decimal(k[6]*(s**beta[6]))])
    t.align = "r"
    print t

def getTemps(eventType):
    global Ecrystal, Entd, timeSteps, eventEnergy,stepSize
    TempArray = []
    if eventType == "ntd":
        Ecrystal = 0
        Entd = eventEnergy
    else:
        Ecrystal = eventEnergy
        Entd = 0
    a,b,c,d,e,f,hs = s,s,s,s,s,(Vb*Rntd(s))/(Rl+Rntd(s)),stepSize
    with open('data/output' + eventType + '.csv', 'wb') as g:
        writer = csv.writer(g)
        t0 = time.clock()
        writer.writerow([0,a,b,c,d,e,f])
        for i in range(timeSteps):
            if i>int((1.0*dur)/(stepSize)):
                Ecrystal = 0
                Entd = 0
            a, b, c, d, e, f = rK6(a, b, c, d, e, f, phonon, electron, heater, crystal, teflon, feedback, hs)
            if i%int(1e4)==0:
                writer.writerow([i*(5.0/timeSteps),a,b,c,d,e,f])
                percentage = i*(100.0/timeSteps)
                if i!=0:
                    sys.stdout.write(str(percentage)+"% | ("+ "%02d:%02d" % divmod(time.clock()-t0,60) + "|" + "%02d:%02d" % divmod((time.clock()-t0)*((100-(percentage))/percentage),60) +")\r")
                    sys.stdout.flush()
    #TempArray.append([a,b,c,d,e,f])
    #TempArray = np.asarray(TempArray)
    #return TempArray
