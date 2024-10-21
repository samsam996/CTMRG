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

        self.T = self.T/np.max(np.abs(self.T))
        self.C = self.C/np.max(np.abs(self.C))

        return None
    
    def one_site_obs(self,a,b):

        tensors = [self.C, self.T, self.C, self.T, a, self.T, self.C, self.T, self.C]
        legs = ([1,2],[2,4,5],[5,6],[1,3,8],[3,10,7,4],[6,7,12],[8,9],[9,10,11],[11,12])
        Z = nc.ncon(tensors, legs)

        tensors_up = [self.C, self.T, self.C, self.T, b, self.T, self.C, self.T, self.C]
        Obs = nc.ncon(tensors_up, legs)
        one_site = Obs/Z

        return one_site

    def partition_function(self,a):

        tensors = [self.C, self.T, self.C, self.T, a, self.T, self.C, self.T, self.C]
        legs = ([1,2],[2,4,5],[5,6],[1,3,8],[3,10,7,4],[6,7,12],[8,9],[9,10,11],[11,12])
        Z = nc.ncon(tensors, legs)

        C4 = np.trace(self.C*self.C*self.C*self.C)

        tensors = [self.C, self.C, self.T, self.T, self.C, self.C]
        legs = ([1,2],[2,4],[1,3,5],[4,3,7],[5,6],[6,7])
        CTTC = nc.ncon(tensors, legs)

        return Z*C4/(CTTC*CTTC)

class CTMRG_nosymm(CTMRG):
     
    def __init__(self,chi,d):
        super().__init__(chi,d)
        self.T1 = np.zeros((chi,d,chi))
        self.T2 = np.zeros((chi,d,chi))
        self.T3 = np.zeros((chi,d,chi))
        self.T4 = np.zeros((chi,d,chi))
        self.C1 = np.zeros((chi,chi))
        self.C2 = np.zeros((chi,chi))
        self.C3 = np.zeros((chi,chi))
        self.C4 = np.zeros((chi,chi))

    def invert_sqrt_s(self,s):
        
        invert_sqrt_s = np.zeros((self.chi, self.chi))
        for i in range(len(s)):
            if s[i,i] > 1e-12:
                invert_sqrt_s[i,i] = 1/np.sqrt(s[i,i])
            else:
                invert_sqrt_s[i,i] = 0

        invert_sqrt_s = invert_sqrt_s/np.max(invert_sqrt_s)

        return invert_sqrt_s


    def leftmove(self,a):

        tensors = [self.C1, self.T1, self.T1, self.C2, self.T4, a, a, self.T2]
        legs = ([1,2],[2,4,5],[5,6,8],[8,9],[1,3,-3],[3,-4,7,4],[7,-2,10,6],[9,10,-1])
        R = nc.ncon(tensors, legs) # chi d chi d
        R = np.reshape(R, (self.chi*self.d, self.chi*self.d)) # dchi dchi

        tensors = [self.T4, a, a, self.T2, self.C4, self.T3, self.T3, self.C3]
        legs = ([-1,3,1],[3,4,7,-2],[7,6,8,-4],[-3,8,9],[1,2],[2,4,5],[5,6,10],[10,9]) 
        Rbar = nc.ncon(tensors, legs) 
        Rbar = np.reshape(Rbar, (self.chi*self.d, self.chi*self.d)) # dchi dchi


        # A = np.random.rand(3,3)
        # u,s,v = np.linalg.svd(A)
        # print('AAA', A-u@np.diag(s)@v)

        u,s,v = np.linalg.svd(R@Rbar)
        u = u[:,0:self.chi]
        v = v[0:self.chi,:]
        s = np.diag(s)
        s = s[0:self.chi,0:self.chi]
        
        invert_sqrt_s = self.invert_sqrt_s(s)
        P4bar = Rbar@np.transpose(v)@invert_sqrt_s
        P4 = np.transpose(R)@u@invert_sqrt_s # dchi chit

        # print(P4@np.transpose(P4bar))


        P4 = np.reshape(P4, (self.chi, self.d, self.chi)) # chi d chi_t
        P4bar = np.reshape(P4bar, (self.chi, self.d, self.chi)) # chi d chi_t

        tensors = [self.C1, self.T1, P4bar]
        legs = ([1,2],[2,3,-2],[1,3,-1])
        self.C1 = nc.ncon(tensors, legs)

        tensors = [P4, self.T4, a, P4bar]
        legs = ([1,3,-1],[1,2,4],[2,5,-2,3],[4,5,-3])
        self.T4 = nc.ncon(tensors, legs)

        tensors = [P4, self.C4, self.T3]
        legs = ([1,3,-1],[1,2],[2,3,-2])
        self.C4 = nc.ncon(tensors, legs)

        self.C1 = self.C1/np.max(np.abs(self.C1))
        self.C4 = self.C4/np.max(np.abs(self.C4))
        self.T4 = self.T4/np.max(np.abs(self.T4))

        return None
    
    def rightmove(self,a):
    
        tensors = [self.C1, self.T1, self.T1, self.C2, self.T4, a, a, self.T2]
        legs = ([1,2],[2,4,5],[5,6,8],[8,9],[1,3,-3],[3,-4,7,4],[7,-2,10,6],[9,10,-1])
        Rbar = nc.ncon(tensors, legs) # chi d chi d
        Rbar = np.reshape(Rbar, (self.chi*self.d, self.chi*self.d)) # dchi dchi

        tensors = [self.T4, a, a, self.T2, self.C4, self.T3, self.T3, self.C3]
        legs = ([-1,3,1],[3,4,7,-2],[7,6,8,-4],[-3,8,9],[1,2],[2,4,5],[5,6,10],[10,9]) 
        R = nc.ncon(tensors, legs) 
        R = np.reshape(R, (self.chi*self.d, self.chi*self.d)) # dchi dchi

        u,s,v = np.linalg.svd(R@Rbar)
        
        u = u[:,0:self.chi]
        v = v[0:self.chi,:]
        s = np.diag(s)
        s = s[0:self.chi,0:self.chi]
        invert_sqrt_s = self.invert_sqrt_s(s)

        # tensors = [Rbar,np.transpose(v), invert_sqrt_s]

        P3bar = Rbar@np.transpose(v)@invert_sqrt_s
        P3 = np.transpose(R)@u@invert_sqrt_s

        # print(P3bar@np.transpose(P3))
        # print(invert_sqrt_s)

        P3bar = np.reshape(P3bar, (self.chi, self.d, self.chi))
        P3 = np.reshape(P3, (self.chi, self.d, self.chi))
        
        tensors = [self.T1, self.C2, P3bar]
        legs = ([-1,3,1],[1,2],[2,3,-2])
        self.C2 = nc.ncon(tensors, legs)

        tensors = [P3, a, self.T2, P3bar]
        legs = ([2,1,-1],[-2,4,3,1],[2,3,5],[5,4,-3])
        self.T2 = nc.ncon(tensors, legs)

        tensors = [P3, self.T3, self.C3]
        legs = ([3,1,-2],[-1,1,2],[2,3])
        self.C3 = nc.ncon(tensors, legs)

        self.C2 = self.C2/np.max(np.abs(self.C2))
        self.T2 = self.T2/np.max(np.abs(self.T2))
        self.C3 = self.C3/np.max(np.abs(self.C3))

        return None

    def upmove(self,a):

        tensors = [self.C1, self.T1, self.T4, a, self.T4, a, self.C4, self.T3]
        legs = ([1,2],[2,3,-1],[1,4,5],[4,7,-2,3],[5,6,8],[6,10,-4,7],[8,9],[9,10,-3])
        Rbar = nc.ncon(tensors, legs)
        Rbar = np.reshape(Rbar, (self.d*self.chi, self.chi*self.d))

        tensors = [self.T1, self.C2, a, self.T2, a, self.T2, self.T3, self.C3]
        legs = ([-3,3,1],[1,2],[-4,6,4,3],[2,4,5],[-2,9,7,6],[5,7,8],[-1,9,10],[10,8])
        R = nc.ncon(tensors, legs)
        R = np.reshape(R, (self.d*self.chi, self.d*self.chi))

        u,s,v = np.linalg.svd(R@Rbar)
        u = u[:,0:self.chi]
        v = v[0:self.chi,:]
        s = np.diag(s)
        s = s[0:self.chi,0:self.chi]
        invert_sqrt_s = self.invert_sqrt_s(s)

        P1bar = np.transpose(R)@u@invert_sqrt_s
        P1 = Rbar@np.transpose(v)@invert_sqrt_s

        # P1 = np.transpose(R)@u@invert_sqrt_s
        # P1bar = Rbar@np.transpose(v)@invert_sqrt_s
        # print(P1bar@np.transpose(P1))
        # print(np.shape(P1))
        # print(invert_sqrt_s)

        P1 = np.reshape(P1, (self.chi, self.d, self.chi))
        P1bar = np.reshape(P1bar, (self.chi, self.d, self.chi))

        tensors = [self.C1, self.T4, P1bar]
        legs = ([1,2],[1,3,-1],[2,3,-2])
        self.C1 = nc.ncon(tensors,legs)

        tensors = [P1, self.T1, a, P1bar]
        legs = ([1,2,-1],[1,3,4],[2,-2,5,3],[4,5,-3])
        self.T1 = nc.ncon(tensors,legs)

        tensors = [P1, self.C2, self.T2]
        legs = ([1,3,-1],[1,2],[2,3,-2])
        self.C2 = nc.ncon(tensors, legs)

        self.C1 = self.C1/np.max(np.abs(self.C1))
        self.T1 = self.T1/np.max(np.abs(self.T1))
        self.C2 = self.C2/np.max(np.abs(self.C2))

        return None
    
    def downmove(self,a):

        return None


    def evolution(self,a):

        self.leftmove(a)
        self.rightmove(a)
        self.upmove(a)

        return None
     
    def open_boundary_condition(self,a):
        for i in range(self.d):
            self.T1[i,i,i] = 1
            self.T2[i,i,i] = 1
            self.T3[i,i,i] = 1
            self.T4[i,i,i] = 1
            self.C1[i] = 1
            self.C2[i] = 1
            self.C3[i] = 1
            self.C4[i] = 1

    def fixed_boundary_condition(self,a):
            self.T1[0,0,0] = 1
            self.T2[0,0,0] = 1
            self.T3[0,0,0] = 1
            self.T4[0,0,0] = 1
            self.C1[0,0] = 1
            self.C2[0,0] = 1
            self.C3[0,0] = 1
            self.C4[0,0] = 1

     # need to define the left, up and down moves as well. 

     # need to do the evolution made of the four moves. 

    def one_site_obs(self,a,b):

        tensors = [self.C1, self.T1, self.C2, self.T4, a, self.T2, self.C4, self.T3, self.C3]
        legs = ([1,2],[2,4,5],[5,6],[1,3,8],[3,10,7,4],[6,7,12],[8,9],[9,10,11],[11,12])
        Z = nc.ncon(tensors, legs)

        tensors_up = [self.C1, self.T1, self.C2, self.T4, b, self.T2, self.C4, self.T3, self.C3]
        Obs = nc.ncon(tensors_up, legs)
        one_site = Obs/Z

        return one_site 

    def partition_function(self,a):

        tensors = [self.C1, self.T1, self.C2, self.T4, a, self.T2, self.C4, self.T3, self.C3]
        legs = ([1,2],[2,4,5],[5,6],[1,3,8],[3,10,7,4],[6,7,12],[8,9],[9,10,11],[11,12])
        Z = nc.ncon(tensors, legs)

        tensors = [self.C1, self.C2, self.C3, self.C4]
        legs = ([1,2],[2,3],[4,3],[1,4])
        square = nc.ncon(tensors, legs)

        tensors = [self.C1, self.C2, self.T4, self.T2, self.C4, self.C3]
        legs = ([1,2],[2,4],[1,3,5],[4,3,7],[5,6],[6,7])
        vert = nc.ncon(tensors, legs)

        tensors = [self.C1, self.T1, self.C2, self.C4, self.T3, self.C3]
        legs = ([1,2],[2,4,5],[5,7],[1,3],[3,4,6],[6,7])
        hor = nc.ncon(tensors, legs)

        return Z*square/(hor*vert)