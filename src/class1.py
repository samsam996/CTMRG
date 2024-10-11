import numpy as np
import ncon as nc


class CTMRG:
     def __init__(self, chi, d):
        self.chi = chi  
        self.d = d      

class CTMRG_C4v(CTMRG):

    def __init__(self,chi,d):
        super().__init__(chi,d)
        self.C = np.zeros((chi,chi))
        self.T = np.zeros((chi,d,chi))

    def open_boundary_condition(self,a):
        for i in range(self.d):
            self.T[i,i,i] = 1
            self.C[i] = 1

    def fixed_boundary_condition(self,a):
            self.T[0,0,0] = 1
            self.C[0,0] = 1
       
    def evolution(self,a):

        tensors = [self.C, self.T, self.T, a]
        legs = ([1,2],[2,4,-3],[1,3,-1],[3,-2,-4,4])
        Cprime = nc.ncon(tensors, legs)
        Cprime_ = np.reshape(Cprime, (self.chi*self.d, self.chi*self.d)) # dchi dchi

        U,S,V = np.linalg.svd(Cprime_)

        U = U[:,0:self.chi]  # dchi chi_trunc
        U = np.reshape(U, (self.chi, self.d, self.chi) )  # chi d chi_trunc

        tensors = [Cprime, U, U]
        legs = ([1, 2, 4, 3],[1, 2, -1],[4, 3, -2])
        self.C = nc.ncon(tensors, legs)

        tensors = [U, self.T, a, U]
        legs = ([1,2,-1],[1,3,5],[2,-2,4,3],[5,4,-3])
        self.T = nc.ncon(tensors, legs)

        self.C = self.C + np.transpose(self.C)
        self.T = self.T + np.transpose(self.T, (2,1,0))

        self.T = self.T/np.max(self.T)
        self.C = self.C/np.max(self.C)

        return None



class CTMRG_nosymm(CTMRG):
     
    def __init__(self,chi,d):
        super().__init__(chi,d)
        self.T1 = np.zeros(chi,d,chi)
        self.T2 = np.zeros(chi,d,chi)
        self.T3 = np.zeros(chi,d,chi)
        self.T4 = np.zeros(chi,d,chi)
        self.C1 = np.zeros(chi,chi)
        self.C2 = np.zeros(chi,chi)
        self.C3 = np.zeros(chi,chi)
        self.C4 = np.zeros(chi,chi)


    def rightmove(self,a):

        tensors = [self.C1, self.T1] # etc ...
        R = 0 
        Rbar = 0 

        tensors = [self.C1, self.T1]
        legs = 0
        C1prime = nc.ncon(tensors, legs)

        return None
     
     # need to define the left, up and down moves as well. 

     # need to do the evolution made of the four moves. 

