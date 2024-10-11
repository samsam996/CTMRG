
import ncon as nc
import numpy as np

def one_site_obs(a,b,ctm):


    tensors = [ctm.C, ctm.T, ctm.C, ctm.T, a, ctm.T, ctm.C, ctm.T, ctm.C]
    legs = ([1,2],[2,4,5],[5,6],[1,3,8],[3,10,7,4],[6,7,12],[8,9],[9,10,11],[11,12])
    Z = nc.ncon(tensors, legs)

    tensors_up = [ctm.C, ctm.T, ctm.C, ctm.T, b, ctm.T, ctm.C, ctm.T, ctm.C]
    Obs = nc.ncon(tensors_up, legs)


    one_site = Obs/Z

    return one_site

def partition_function(a,ctm):

    tensors = [ctm.C, ctm.T, ctm.C, ctm.T, a, ctm.T, ctm.C, ctm.T, ctm.C]
    legs = ([1,2],[2,4,5],[5,6],[1,3,8],[3,10,7,4],[6,7,12],[8,9],[9,10,11],[11,12])
    Z = nc.ncon(tensors, legs)

    C4 = np.trace(ctm.C*ctm.C*ctm.C*ctm.C)

    tensors = [ctm.C, ctm.C, ctm.T, ctm.T, ctm.C, ctm.C]
    legs = ([1,2],[2,4],[1,3,5],[4,3,7],[5,6],[6,7])
    CTTC = nc.ncon(tensors, legs)

    return Z*C4/(CTTC*CTTC)