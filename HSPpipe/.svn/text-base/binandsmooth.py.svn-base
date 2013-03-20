##############written by fbb june 2011#####
########utility to bin and tophat smooth arrays of data. 
########june 29 2011 so far it only works with oned array (e.g lighcurves!)
##########################################################################
#function: onedsmooth
#   INPUT: 
#    arr: the imput array numpy array or list, will be converted to 
#         numpy array eitherway
#    bsz: bin size
#    bin: 1 for binning the data, 2 for smoothing the data
##########################################################################

import numpy as np

def onedsmooth(arr,bsz, bin):
    arr=np.array(arr)

    bszm1=bsz-1

    newarr = arr[:-(bszm1)]+arr[(bszm1):]
    for i in range(1,bszm1):
        newarr+=arr[i:-(bszm1-i)]

    if bin>0:
        return newarr[::bsz]

    else:
        return newarr/float(bsz)
        


