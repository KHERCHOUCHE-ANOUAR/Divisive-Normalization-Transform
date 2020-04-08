import numpy as np
from pyrBandIndices import pyrBandIndices
def pyrBand(pyr, pind, band):
   return np.reshape( pyr[pyrBandIndices(pind,band)],(pind[band,0], pind[band,1]))
