import math
import numpy as np
import time
import sys
from decimal import Decimal
from prettytable import PrettyTable

def getTemps(testNumber, eventType, alpha, beta, k, s, R0, Rl, T0, gamma, Vb, Cp):

    parameters = {
        "alpha": alpha,
        "beta": beta,
        "k": k,
        "s": s,
        "R0": R0,
        "Rl": Rl,
        "T0": T0,
        "gamma": gamma,
        "Vb": Vb,
        "Cp": Cp
    }
    
    t = PrettyTable(['Name','Value'])
    t.add_row([ "Temperature" , '%.2E' % Decimal(s)])
    t.add_row([ "T0" , '%.2E' % Decimal(T0)])
    t.add_row([ "R0" , '%.2E' % Decimal(R0)])
    t.align = "r"
    print t

    #thermal model equations
    def cond(index,start,end):
        return ((k[index-1]/(beta[index-1]+1))*((start)**((beta[index-1])+1) - (end)**((beta[index-1])+1)))
    def Rntd(temp):
        return R0*math.exp((T0/temp)**gamma)
    def phonon(a,b,c,d,e,f):
        return ((Entd/dur)+cond(1,d,a) + cond(2,b,a) - cond(3,a,s))/(alpha[0]*(a**3))
    def electron(a,b,c,d,e,f):
        return ((Eelectron/dur) + (f**2 - ((Vb*Rntd(s))/(Rl+Rntd(s)))**2)/(Rntd(b)) - cond(2,b,a))/(alpha[1]*b)
    def heater(a,b,c,d,e,f):
        return ((Eheater/dur)+cond(4,d,c) - cond(5,c,s))/(alpha[2]*c)
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

    timeSteps = int(5e6)
    stepSize = 1e-6

    #power input parameters
    Ecrystal, Entd, Eheater, Eelectron = 0,0,0,0
    eventEnergy = (2615.0)*(1.6022e-16) #2615 KeV
    dur = 1e-3
    sampling_rate = 125

    info = {
        "test_number":testNumber,
        "event_type": eventType,
        "duration": stepSize*timeSteps,
        "sampling_rate": sampling_rate,
        "event_energy": eventEnergy,
        "power_injection_duration": dur
    }
        
    TempArray = []
    if eventType == "ntd":
        Entd = eventEnergy
        print "NTD Event"
    elif eventType == "heater":
        Eheater = eventEnergy
        print "Heater Event"
    elif eventType == "electron":
        Eelectron = eventEnergy
        print "Electron Event"
    else:
        Ecrystal = eventEnergy
        print "Crystal Event"
        
    a,b,c,d,e,f,hs = s,s,s,s,s,(Vb*Rntd(s))/(Rl+Rntd(s)),stepSize
   
    a_data, b_data, c_data, d_data, e_data, f_data = [],[],[],[],[],[]
    
    t0 = time.clock() 
  
    for i in range(int(timeSteps+1e3)):
        
        #turn off power deposition
        if i>int((1.0*dur)/(stepSize)+(timeSteps/5)):
            Ecrystal = 0
            Entd = 0
            Eheater = 0
            Eelectron = 0
        
        #
        if i%int(1/(sampling_rate*stepSize))==0:
            a_data.append(a)
            b_data.append(b)
            c_data.append(c)
            d_data.append(d)
            e_data.append(e)
            f_data.append(f)
            
        if i>(timeSteps/5):
            a, b, c, d, e, f = rK6(a, b, c, d, e, f, phonon, electron, heater, crystal, teflon, feedback, hs)
            if i%int(dur/stepSize)==0:
                percentage = i*(100.0/timeSteps)
                sys.stdout.write("%.2f" % percentage +"% | ("+ "%02d:%02d" % divmod(time.clock()-t0,60) + "|" + "%02d:%02d" % divmod((time.clock()-t0)*((80-(percentage-20))/(percentage-20)),60) +")\r")   
    data = {
        "a": a_data,
        "b": b_data,
        "c": c_data,
        "d": d_data,
        "e": e_data,
        "f": f_data
    }
    
    alldata = {
        "model_info": info,
        "thermal_parameters": parameters,
        "data": data
    }
    
    return alldata