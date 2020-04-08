import numpy as np

def pyrBandIndices(pind,band):

    if ((band > np.size(pind,0)) |(band < 1)):
      raise Exception('BAND_NUM must be between 0 and number of pyramid bands')

    if (np.size(pind,1) != 2):
      raise Exception('INDICES must be an Nx2 matrix indicating the size of the pyramid subbands')

    ind = 1

    for l in range(band-1):
        p=1
        for i in range(len(pind[l,:])):
            p=p*pind[l,i]
        ind = ind + p
    p=1
    for i in range(len(pind[band-1,:])):
        p=p*pind[band,i]
    t=[]
    for i in range(ind,ind+p):
       t.append(i)
    #print(np.shape(t))
    return t