from eventGenerator import *

alpha = [2.7e-8,9.9e-9,1e-11,9.14e3,4.28e-7,0]
beta = [3,4.37,2.4,3,2.4,2,1]
k = [2.34e-3,1.4,3.2e-5,1.3e-3,3.2e-5,2.5e-5,4e-7]

R0 = 3.2
s = .015
Rl = 6e10
T0 = 5.2
gamma = 0.5
Vb = 5
Cp = 5e-10

getTemps(testNumber, 'crystal', alpha, beta, k, s, R0, Rl, T0, gamma, Vb, Cp)
