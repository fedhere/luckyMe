#Author: Scott Hillberry
#Last updated 8/30/2011
###Tukey Window code credited to Dat Chu of University of Houston. Updated by
###Scott to from a 2d Tukey
###http://leohart.wordpress.com/ << Dat Chu's blog
#############################################################################

import sys
import pyfits as PF
import numpy as np
from scipy.fftpack import ifft2, fftshift

def pad(data, padding):
### padding is the number of pixels to be added to each edge
    pad = np.zeros((data.shape[0]+2*padding,data.shape[1]+2*padding), float)
    #print 'Padding image...'
    #print '%s --padding--> %s' %(str(data.shape),str(pad.shape))
    pad[padding:-padding,padding:-padding] = data
    return pad


def drizzle(data, writefits):
    '''Takes an image, doubles the size, and creates new pixel values using averages of surrounding pixels.
#   +---+---+
#   | 2 | 3 |
#   +---+---+
#   | 1 | 4 |
#   +---+---+
    Illustration of grid, where #1 is the original pixel and 2-4 come from doubling the size.'''

### CREATE LARGER ARRAY ###

    new = np.zeros(tuple(x*2 for x in data.shape), float)
#    padd = pad(data, 1)
#fed: changed padding so that you add a replica of the first and last row and column to the data. this way you do not get crazy values at the edges of the drizzled image

    padd = np.vstack( [ data[0], data , data[-1:] ] )
    padd = np.column_stack( [ padd[:,0], padd, padd[:,-1:] ] )

    new[::2,::2] = data  #setting values for Pxl#1
    new[1::2,::2] = (data + padd[2:,1:-1])/2 #values for Pxl#2
    new[::2,1::2] = (data + padd[1:-1,2:])/2 #values for Pxl#4
    new[1::2,1::2] = (data + padd[2:,1:-1] + padd[1:-1,2:] + padd[2:,2:])/4 #values for Pxl#3

    return new


def hann2d(axis2, axis1):
    #print 'Creating Hann window...'
    #print 'image size: (%d, %d)' %(axis2, axis1)
    n = axis2
    if axis2 > axis1:
        n = axis2
    else:
        n = axis1
## creating square hamming window
    ham2d = np.outer(hamming(n),np.ones(n))
    ham2d = np.sqrt(ham2d * ham2d.T)
## clipping window to image size
    start2 = int(n/2-axis2/2)
    start1 = int(n/2-axis1/2)
    ham2d = ham2d[start2:start2+axis2,start1:start1+axis1]
    #print 'window size: ' + str(ham2d.shape)

    return ham2d


def tukey2d(shp, alpha=0.6):
    '''The Tukey window, also known as the tapered cosine window, can be regarded as a cosine lobe of width alpha * N / 2 that is convolved with a rectangle window of width (1 - alpha / 2). At alpha = 1 it becomes a Hann window, and at alpha = 0 it becomes rectangular.'''

    window_length = shp[0]
    if shp[1] > shp[0]:
        window_length = shp[1] #window_length is longest axis

    # Special cases
    if alpha <= 0:
        window = np.ones(window_length) #rectangular window
    elif alpha >= 1:
        window = np.hanning(window_length)
    else:    # Normal case: 0 < alpha < 1
        x = np.linspace(0, 1, window_length)
        window = np.ones(x.shape)
        # first condition 0 <= x < alpha/2
        first_condition = x<alpha/2
        window[first_condition] = 0.5 * (1 + np.cos(2*np.pi/alpha * (x[first_condition] - alpha/2) ))
        # second condition already taken care of
        # third condition 1 - alpha / 2 <= x <= 1
        third_condition = x>=(1 - alpha/2)
        window[third_condition] = 0.5 * (1 + np.cos(2*np.pi/alpha * (x[third_condition] - 1 + alpha/2)))

###### creating 2d tukey
    window2d = np.outer(window,np.ones(window_length))
    window2d = np.sqrt(window2d * window2d.T)
###### clipping window to image size
    start2 = int(window_length/2-shp[0]/2)
    start1 = int(window_length/2-shp[1]/2)
    window2d = window2d[start2:start2+shp[0],start1:start1+shp[1]]

    return window2d

def correlate(ffts, fast, rad): #NOTE: fits are already fftn'd and tukey'd   
    if fast:
	sys.path.append("~/Downloads/PyFFTW3-0.2.1/build/lib")
           ##path to PyFFTW3##
        import fftw3



    naxis1 = ffts[0].shape[1] #NAXIS1 is x-axis, i.e. 'fast axis' 
    naxis2 = ffts[0].shape[0] #.shape returns (slow, fast)

#### CREATING TUKEY WINDOWS ####
    conj = np.conjugate(ffts[1])

    if fast:
	r = np.zeros(ffts[0].shape, complex)
        rplan = fftw3.Plan(ffts[0]*conj, r, 'backward')
        rplan.execute()
        r = r.real
    else:
        r = ifft2(ffts[0]*conj).real
 

    rshift = fftshift(r)
    reg = rshift[naxis2/2-rad:naxis2/2+rad,naxis1/2-rad:naxis1/2+rad]

    phase = np.where(reg == np.max(reg))

    axis2_shift, axis1_shift = [x - rad for x in phase]
### Checks if image has negative shift ### 
#    if phase[0] > naxis2/2:
#        axis2_shift =  phase[0] - naxis2
#    else:
#        axis2_shift = phase[0]
#
#    if phase[1] > naxis1/2:
#        axis1_shift = phase[1] - naxis1
#    else:
#        axis1_shift = phase[1]

    return [axis2_shift, axis1_shift, r]



def detripling(detrip, luckies,mcoff, hd, inpath):
# gsx, gsy, rad, minsep):
#####     DETRIPLING     #####

    gsx,gsy,rad=hd['gsx'],hd['gsy'],hd['rad']
    coresz = hd['coresz']
    c=hd['coresz']/2
    print 'Separating cores...'
    sys.stdout.flush()
    for i in xrange(mcoff):
        core1=[0,0,0,0,0]
        core2=[0,0,0,0,0]
        name=luckies[i][0]
        tmpfile=PF.getdata(inpath+'/'+name)
        reg=tmpfile[gsy-rad:gsy+rad+1,gsx-rad:gsx+rad+1]
        for y in xrange(1,2*rad-1):
            for x in xrange(1,2*rad-1):
                core=reg[y-c:y+c+1,x-c:x+c+1]
                if core.sum()>=core1[0]:
                    core1=[core.sum(),x,y,core.argmax()%coresz,core.argmax()/coresz]
                    #core1=[core.sum(),x,y,argmax(core)%coresz,argmax(core)/coresz]
                        
                        
                if core.sum()==core1[0]:
                    continue
                elif core.sum()>=core2[0] and np.sqrt((1.*x-core1[1])**2+(1.*y-core1[2])**2)>=hd['minsep']:
                    core2=[core.sum(),x,y,core.argmax()%coresz,core.argmax()/coresz]
                    #core2=[core.sum(),x,y,argmax(core)%coresz,argmax(core)/coresz]

        if detrip=='v':
            d=2
            dtr='VERTICAL'
        elif detrip=='h':
            d=1
            dtr='HORIZONTAL'
        else:
            print 'wrong detripling flag'
            return -1  
        
        if core1[d]>core2[d]:
            luckies[i]=[name,core1[1]-c+core1[3],core1[2]-c+core1[4]]
        else:
            luckies[i]=[name,core2[1]-c+core2[3],core2[2]-c+core2[4]]
        
    return luckies,dtr 
    
def calccents(tmparr):
#caluclates a weighted average centroid within a region of pixels    
    #print tmparr.shape
    if tmparr.shape==(0,0):
        print "empty array passed to calccents, returning (0,0,0)"
        return (0,0,0)
    allx =np.sum(tmparr,axis=0)
    xl = len(allx)
    indx = np.arange(xl)
    allf = np.sum(allx)
    #print 'allx.shape = ' + str(allx.shape)
    #pridt 'indx.shape = ' + str(indx.shape)
    mx = np.sum(allx*indx)/allf

    ally = np.sum(tmparr,axis=1)
    yl = len(ally)
    indy = np.arange(yl)
    my = np.sum(ally*indy)/allf
    
    sky =np.sum(tmparr[0,:])
    sky+=np.sum(tmparr[-1,:])
    sky+=np.sum(tmparr[:,0])
    sky+=np.sum(tmparr[:,-1])
    (lx,ly)=np.shape(tmparr)
    sky*=(lx*ly)/(2.0*(lx+ly))
    allf=allf-sky
    #returning: weighted x centroid, weighted y centroid, sky value
    return (mx,my, allf)

def peakdetect2d(img, lookahead = 30, delta = 0.5):
    """
    Converted from/based on a MATLAB script at http://billauer.co.il/peakdet.html
    Python code pulled from sixtenbe on github.com. Modified for 2d.

    Algorithm for detecting local maximas and minmias in a signal.
    Discovers peaks by searching for values which are surrounded by lower
    or larger values for maximas and minimas respectively
    
    keyword arguments:
    img -- input image
    lookahead -- (optional) distance to look ahead from a peak candidate to
        determine if it is the actual peak (default: 30 pxls) 
        '(sample / period) / f' where '4 >= f >= 1.25' might be a good value
    delta -- (optional) this specifies a minimum difference between a peak and
        the following points, before a peak may be considered a peak. Useful
        to hinder the algorithm from picking up false peaks towards to end of
        the signal. To work well delta should be set to 'delta >= RMSnoise * 5'.
        (default: 0.5)
            Delta function causes a 20% decrease in speed, when omitted
            Correctly used it can double the speed of the algorithm
    
    return -- list [maxtab] containing the positive peaks. Each cell of the listcontains a tupple of:
        (position, peak_value) 
        to get the average peak value do 'np.mean(maxtab, 0)[1]' on the results
    """
    #shpe = img.shape

    maxtab = []
    maxmat = np.zeros(img.shape, bool)
    dump = []   #Used to pop the first hit which always is false
    #length = len(y_axis)
    #if x_axis is None:
    #    x_axis = range(length)
    
    #perform some checks
    #if length != len(x_axis):
    #    raise ValueError, "Input vectors y_axis and x_axis must have same length"
    if lookahead < 1:
        raise ValueError, "Lookahead must be above '1' in value"
    if not (np.isscalar(delta) and delta >= 0):
        raise ValueError, "delta must be a positive number"
    
### PAD TRAILING EDGES ###
    pad = np.zeros((img.shape[0]+lookahead,img.shape[1]+lookahead),float)
    pad.fill(-np.Inf)
    #pad[:-lookahead,:-lookahead].fill(0)
    #print pad.shape
    pad[:-lookahead,:-lookahead] = img
    shpe = pad.shape
    #raw_input()
    #maxima candidates are temporarily stored in mx
    mx = -np.Inf
    
    #Only detect peak if there is 'lookahead' amount of points after it
    for yblock in range(1,shpe[0]-1, lookahead):
        for xblock in range(1, shpe[1]-1, lookahead):
	    #print '--->sublock [%d:%d,%d:%d]' %(yblock, yblock+lookahead, xblock, xblock+lookahead)
	    for y in range(yblock, yblock+lookahead):
	        for x in range(xblock, xblock+lookahead):
		    if not yblock+lookahead  > shpe[0] and \
		       not xblock+lookahead > shpe[1]:
	                f = pad[y,x]
                        #print 'value = %f at (%d,%d)' %(f, y, x)
                        #print 'max = %f' %mx
	                if f > mx:
                            #print 'potential max found...'
		            mx = f
	 	            mxpos = (y,x)

	    y = mxpos[0]
	    x = mxpos[1]
            #print '(y, x) = (%d, %d)' %(y,x)
            forwardx = x + lookahead
	    if forwardx > shpe[1]: forwardx = shpe[1]
	    
	    backx = x - lookahead
	    if backx < 0: backx = 0

	    backy = y - lookahead
	    if backy < 0: backy = 0
            
            forwardy = y + lookahead
	    if forwardy > shpe[0]: forwardy = shpe[0]
	#Maxima peak candidate found
	#look ahead in signal to ensure that this is a peak and not jitter
	    #print pad[y:y+1, x+1:forwardx].max()
            #print pad[y+1:forwardy, x:x+1].max()
            #print pad[y:y+1, backx:x].max()
            #print pad[backy:y, x:x+1].max()

	    if pad[y:y+1, x+1:forwardx].max() < mx - delta and \
	       pad[y+1:forwardy, x:x+1].max() < mx - delta and \
	       pad[y:y+1, backx:x].max() < mx - delta and \
               pad[backy:y, x:x+1].max() < mx - delta and mx != -np.Inf: 
	        #print 'Yes, max found'
	        maxtab.append((mxpos, mx))
	        maxmat[mxpos[0],mxpos[1]] = True
	        dump.append(True)
      
            #print maxtab 
            #raw_input()
	mx = -np.Inf
 
    #Remove the false hit on the first value of the y_axis
    #try:
    #    if dump[0]:
    #        maxtab.pop(0)
    #        print "pop max"
    #    del dump
    #except IndexError:
        #no peaks were found, should the function return empty lists?
    #    pass
   
    maxtab = list(set(maxtab)) 
    maxtab = sorted(maxtab, key=lambda tup: -tup[1])
    return (maxtab, maxmat)
