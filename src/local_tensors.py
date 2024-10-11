

import numpy as np
import ncon as nc
from scipy.linalg import sqrtm

def Ising2D(temp, J):

    a = np.zeros((2,2,2,2))
    b = np.zeros((2,2,2,2))
    kron = np.zeros((2,2,2,2))
    kronm = np.zeros((2,2,2,2))

    Q = np.zeros((2,2))

    Q[0,0] = np.exp(J/temp)
    Q[0,1] = np.exp(-J/temp)
    Q[1,0] = np.exp(-J/temp)
    Q[1,1] = np.exp(J/temp)

    for i in range(2):
        kron[i,i,i,i] = 1
        kronm[i,i,i,i] = 1
    kronm[1,1,1,1] = -1

    sqrtQ = sqrtm(Q)

    tensors = [sqrtQ, sqrtQ, sqrtQ, sqrtQ, kron]
    legs = ([-1,1],[-2,2],[-3,3],[-4,4],[1,2,3,4])
    a = nc.ncon(tensors, legs)

    tensors = [sqrtQ, sqrtQ, sqrtQ, sqrtQ, kronm]
    legs = ([-1,1],[-2,2],[-3,3],[-4,4],[1,2,3,4])
    b = nc.ncon(tensors, legs)

    return a, b


