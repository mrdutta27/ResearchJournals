def cond(index,start,end):
    return ((k[index-1]/(beta[index-1]+1))*(start**(beta[index-1]+1) - end**(beta[index-1]+1)))

def phonon(a,b,c,d,e,f): #add NTD Power
    return (cond(1,d,a) + cond(2,b,a) - cond(3,a,s))/(alpha[0]*(a**3)) + (Pntd/dur)

def electron(a,b,c,d,e,f): #figure out T
    return ((f**2)/(R0*math.exp((T0/np.absolute(b))**gamma)) - cond(2,b,a))/(alpha[1]*b)

def heater(a,b,c,d,e,f):
    return (cond(4,d,c) - cond(5,c,s))/(alpha[2]*c)

def crystal(a,b,c,d,e,f):
    return (-1*cond(1,d,a) - cond(4,d,c) - cond(6,d,e))/(alpha[3]*(d**3))

def teflon(a,b,c,d,e,f):
    return (cond(6,d,e) - cond(7,e,s)) / ((alpha[4]*e) + (alpha[5]*(e**3)))

def feedback (a,b,c,d,e,f):
    return (Vb-(f*((Rl+(R0*math.exp((T0/np.absolute(b))**gamma)))/((R0*math.exp((T0/np.absolute(b))**gamma))))))/(Rl*Cp)
