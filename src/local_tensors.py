

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



def qStatePotts(temp, J, q):

    a = np.zeros((q,q,q,q))
    b = np.zeros((q,q,q,q))
    kron = np.zeros((q,q,q,q))
    kronm = np.zeros((q,q,q,q))

    Q = np.ones((q,q))

    for i in range(q):
        Q[i,i] = np.exp(J/temp)
           
    for i in range(q):
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


def AshkinTeller(temp, J, lbda):

    a = np.zeros((4,4,4,4))
    b = np.zeros((4,4,4,4))
    kron = np.zeros((4,4,4,4))
    kronm = np.zeros((4,4,4,4))

    Q = np.ones((4,4))

    x0 = J*(lbda+2)/temp
    x1 = -J*(lbda)/temp
    x2 = J*(lbda-2)/temp
    x3 = -J*(lbda)/temp

    Q[0,0] = np.exp(x0)
    Q[0,1] = np.exp(x1)
    Q[0,2] = np.exp(x2)
    Q[0,3] = np.exp(x3)

    Q[1,0] = np.exp(x2)
    Q[1,1] = np.exp(x0)
    Q[1,2] = np.exp(x3)
    Q[1,3] = np.exp(x1)

    Q[2,0] = np.exp(x1)
    Q[2,1] = np.exp(x3)
    Q[2,2] = np.exp(x0)
    Q[2,3] = np.exp(x2)

    Q[3,0] = np.exp(x3)
    Q[3,1] = np.exp(x2)
    Q[3,2] = np.exp(x1)
    Q[3,3] = np.exp(x0)

        
    for i in range(4):
        kron[i,i,i,i] = 1
        kronm[i,i,i,i] = 1
    kronm[2,2,2,2] = -1

    sqrtQ = sqrtm(Q)

    tensors = [sqrtQ, sqrtQ, sqrtQ, sqrtQ, kron]
    legs = ([-1,1],[-2,2],[-3,3],[-4,4],[1,2,3,4])
    a = nc.ncon(tensors, legs)

    tensors = [sqrtQ, sqrtQ, sqrtQ, sqrtQ, kronm]
    legs = ([-1,1],[-2,2],[-3,3],[-4,4],[1,2,3,4])
    b = nc.ncon(tensors, legs)

    return a, b