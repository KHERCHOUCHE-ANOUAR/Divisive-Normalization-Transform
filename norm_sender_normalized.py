import numpy as np
from scipy.misc import imresize
import math
from skimage.transform import resize
import time
from pyrBand import pyrBand
from pyrBandIndices import pyrBandIndices
def rotate(t,i,j):
  t = np.roll(t, i, axis=0) # row
  t = np.roll(t, j, axis=1) # column
  return t
def norm_sender_normalized(pyro,pind,Nsc,Nor,parent,neighbor,blSzX,blSzY,nbins):
    guardband = 16
    pyro = np.real(pyro)
    Nband = np.size(pind,0)-1
    band = [1 ,3 ,6, 8, 9, 11]
    zRange = 15
    subband=[]
    shape=[]
    p = 0
    for scale in range(Nsc):
        for orien in range(Nor):

            nband = (scale)*Nor+orien+1    #except the ll

            aux = pyrBand(pyro, pind, nband)

            Nsy,Nsx = np.shape(aux)

            prnt = parent & (nband < Nband-Nor)
               # has the subband a parent?
            BL = np.zeros((np.size(aux,0),np.size(aux,1),1 + prnt))
            BL[:,:,0] = aux

            if (prnt):
                auxp = pyrBand(pyro, pind, nband+Nor)
                #print(np.shape(auxp))
                auxp = np.real(resize(auxp, (2*np.size(auxp,0), 2*np.size(auxp,1))))
                #print(np.shape(auxp))


                BL[:,:,1] = auxp[0:Nsy,0:Nsx]

            y=BL
            nv,nh,nb = np.shape(y)
            block = [blSzX,blSzY]

            nblv = nv-block[0]+1	# Discard the outer coefficients
            nblh = nh-block[1]+1   # for the reference (centrral) coefficients (to avoid boundary effects)
            nexp = nblv*nblh			#number of coefficients considered
            N = blSzX*blSzY + prnt # size of the neighborhood

            Ly = (block[0]-1)//2		#block(0) and block(1) must be odd!
            Lx = (block[1]-1)//2

            if (Ly!=math.floor(Ly))|(Lx!=math.floor(Lx)):
                raise Exception('Spatial dimensions of neighborhood must be odd!')

            Y = []		#It will be the observed signal (rearranged in nexp neighborhoods)

            for ny in range(-Ly,Ly+1):	# spatial neighbors
                for nx in range(-Lx,Lx+1):


                    foo = rotate(y[:,:,0],nx,ny)
                    foo = foo[Ly:Ly+nblv,Lx:Lx+nblh]
                    t=np.reshape(foo,(np.size(foo,0)*np.size(foo,1)))

                    Y.append(t)


            if (prnt):

                foo = y[:,:,1]

                foo = foo[Ly:Ly+nblv,Lx:Lx+nblh]

                t=np.reshape(foo,(np.size(foo,0)*np.size(foo,1)))

                Y.append(t)


            #including neighbor
            if (neighbor):

                for neib in range(Nor):
                    if (neib == orien):
                        continue


                    nband1 = (scale)*Nor+neib+1 # except the ll
                    aux1 = pyrBand(pyro, pind, nband1)
                    #print(np.shape(aux1))
                    aux1 = aux1[Ly:Ly+nblv,Lx:Lx+nblh]
                    #print(np.shape(aux1))

                    t=np.reshape(aux1,(np.size(aux1,0)*np.size(aux1,1)))

                    Y.append(t)


            x=0

            for t in range (len(Y)):

                if (np.shape(Y[t])[0]>x):
                    x=np.shape(Y[t])[0]

            E=np.zeros((x,len(Y)))
            for t in range (len(Y)):
                for j in range (len(Y[t])):
                       E[j,t]=Y[t][j]
            Y=E

            C_x = np.inner(np.transpose(Y),np.transpose(Y))/nexp
            # C_x is positive definete covariance matrix
            Q,L = np.linalg.eig(C_x)
            # correct possible negative eigenvalues, without changing the overall variance
            L = np.diag(np.diag(L)*(np.diag(L)>0))*sum(np.diag(L))/(sum(np.diag(L)*(np.diag(L)>0))+(sum(np.diag(L)*(np.diag(L)>0))==0))
            C_x = Q*L*np.transpose(Q)

            o_c = aux[Ly:Ly+nblv,Lx:Lx+nblh]
            o_c = (o_c[:])
            o_c = o_c - np.mean(o_c)

            tempY = np.multiply((np.dot(Y,np.linalg.pinv(C_x))),Y)//N
            z = np.sqrt(np.sum(tempY,1))
            ind = []
            for i in range (len(z)) :
                  if(z[i]!=0):
                   ind.append(i)
            oc=np.reshape(o_c,(1,np.size(o_c,0)*np.size(o_c,1)))
            g_c=[]
            z=np.reshape(z,(1,np.size(z,0)))
            ind1=np.reshape(ind,(1,np.size(ind,0)))
            # print(np.shape(z),np.shape(oc),np.shape(ind),np.size(oc,1))
            for i in range (np.size(ind1,1)-1):
                t=ind1[0,i]
                #print(t)
                if(t>=np.size(oc,1)):
                    break
                #print(t)
                k=oc[0,t]/z[0,t]
                g_c.append(k)



            # consider the guardband

            t=np.sqrt(np.size(g_c))
            g_c=np.reshape(g_c[:(math.floor(t))*(math.floor(t))],(math.floor(t),math.floor(t)))
            gb = guardband/(2^(scale))
            x=np.size(g_c,0)
            y=np.size(g_c,1)

            g_c=g_c[int(gb+1):int(x-gb+1), int(gb+1):int(y-gb+1)]

            g_c= (g_c[:])
            g_c = g_c - np.mean(g_c)


            subband.append(g_c)
            shape.append(np.shape(g_c))





    return subband, np.reshape(shape,(12,2))
#start_time = time.time()
coef=np.load('pyr.npz')
pyr=coef['X']
coef=np.load('pind.npz')
pind=coef['X']

num_or = 6
num_scales = 2
subband, shape=norm_sender_normalized(pyr,pind,num_scales,num_or,1,1,3,3,50)
#print("--- %s seconds ---" % (time.time() - start_time))
print('shape of subbands before DNT transform\n',pind[1:13,:])
print('shape of subbands after DNT transform\n',shape)
