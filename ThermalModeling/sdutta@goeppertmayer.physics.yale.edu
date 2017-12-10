from tqdm import tqdm
import numpy as np

from eventGenerator import *

alpha = [2.7e-8,9.9e-9,1e-11,9.14e3,4.28e-7,0]
beta = [3,4.37,2.4,3,2.4,2,1]
k = [2.34e-3,1.4,3.2e-5,1.3e-3,3.2e-5,2.5e-5,4e-7]

R0 = 3.2
s = .015
Rl = 6e10
T0 = 5.2
gamma = 0.5
#Vb = 5
Cp = 5e-10

V_Bias = [2.576, 3.003, 3.432, 3.859, 4.291, 4.718, 5.147, 5.577, 6.006, 6.437, 6.865]
V_Bol = [0.0082, 0.00892, 0.009538, 0.01009, 0.01055, 0.01097, 0.01134, 0.01165, 0.01193, 0.0122, 0.01243]
I_Bol = [4.28e-11, 4.99e-11, 5.7e-11, 6.42e-11, 7.13e-11, 7.84e-11, 8.56e-11, 9.28e-11, 9.99e-11, 1.07e-10, 1.14e-10]

def arrayRandomVars(array,prop):
    return np.random.uniform(np.multiply(array,(1-prop)), np.multiply(array,(1+prop)))

def getchi2(test, data):
    return sum(abs(np.divide(np.power((np.subtract(data,test)),2),test)))

num = 5
prop = 0.2

bestChi2 = 10000
bestAlpha, bestBeta, bestK = [],[],[]

for j in tqdm(range(1,num),ascii=True):

    V_Data, I_Data = [],[]

    alpha = arrayRandomVars(alpha,prop)
    beta = arrayRandomVars(beta,prop)
    k = arrayRandomVars(k,prop)

    for b in tqdm(range(1,len(V_Bias)),ascii=True):
        v,i = getTemps(alpha, beta, k, s, R0, Rl, T0, gamma, V_Bias[b], Cp)
        V_Data.append(v)
        I_Data.append(i)

    v_chi = getchi2(V_Bol,V_Data)

    if v_chi < bestChi2:
        bestChi2 = v_chi
        bestAlpha = alpha
        bestBeta = beta
        bestK = k

print bestAlpha, bestBeta, bestK
