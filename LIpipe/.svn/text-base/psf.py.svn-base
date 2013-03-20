###############usage: reads in photometry file and for n stars creates x and y average displacements and integrates them over m time steps to create the actual gaussian profile#####################################################
###############



from numpy import *
import sys
import os
#from scipy import *
from scipy import fftpack

class Chdir:         
    def __init__( self, newPath ):  
        self.savedPath = os.getcwd()
        os.chdir(newPath)
        
    def __del__( self ):
        os.chdir( self.savedPath )



def g2d(coords, TENSIG, beta, PSF, nx, ny, c, pixscale):
    data = zeros((nx,ny),float)
    for i in (coords):
        for j in xrange(c[0]-TENSIG,c[0]+TENSIG+1):
            for k in xrange(c[1]-TENSIG,c[1]+TENSIG+1):
                g=10.0/(2.0*3.1415*PSF*PSF/(pixscale*pixscale))*exp(-((j-i[0])*(j-i[0])+(k-i[1])*(k-i[1]))/(2.0*PSF*PSF/(pixscale*pixscale)))
                data[j][k] +=g
   
    return data

def m2d(coords, TENSIG, nx,ny,c, beta, PSF, pixscale):
    data = zeros((nx,ny),float)
    for i in (coords):
        for j in xrange(c[0]-TENSIG,c[0]+TENSIG+1):
            for k in xrange(c[1]-TENSIG,c[1]+TENSIG+1):
                m=10.0/(2.0*3.1415*PSF*PSF/(pixscale*pixscale))*\
                pow((1.0+(((j-i[0])*(j-i[0])+(k-i[1])*(k-i[1]))/(PSF*PSF/(pixscale*pixscale)))),-beta)
                data[j][k] +=m
   
    return data



def mkfits(par, coords1, coords2, newdir, date, PSF, TENSIG, pixscalex, pixscaley,cadence, x, y, nx, ny, c, exposure):
    from pyfits import *
    fitsobj = HDUList()
    # create Primary HDU with minimal header keywords
    hdu = PrimaryHDU()
    # add a 10x5 array of zeros
    h=hdu.header
    h.update('RA', '%s' %par['ra'])
    h.update('Dec', '%s' %par['dec'])
    h.update('DATE', '%s' %date)#par['date'])
    h.update('TELESC', '%s'%par['scope'])
    h.update('CAMERA', '%s'%par['camera'])
    h.update ('IMTYPE', 'LI', 'LI for Lucky Imaging, HSP for high speed photometry.' )
    h.update ('CCDSCLX', '%s' %pixscalex,'arcsec/pixel')
    h.update ('CCDSCLY', '%s' %pixscaley,'arcsec/pixel')

    profile=par['profile'] #'g' for gaussian, 'm' for moffat

    beta = float(par['beta'])
    if profile=='g':
        data1 = g2d(coords1, TENSIG, beta, PSF, nx, ny, c, pixscalex)
        data2 = g2d(coords2, TENSIG, beta, PSF, nx, ny, c, pixscalex)
        hdu.data=concatenate((data1,data2))
        h.update('PROFILE', 'gaussian')
        h.update('PSF', '%s' %PSF, 'arcseconds')
    elif profile == 'm':
        data1 = m2d(coords1, TENSIG, nx,ny,c, beta, PSF, pixscalex)
        data2 = m2d(coords2, TENSIG, nx,ny,c, beta, PSF, pixscalex)
        hdu.data=concatenate((data1,data2))
        h.update('PROFILE', 'moffat')
        h.update('PSFALPHA', '%s' %(par['alpha']))
        h.update('PSFBETA', '%s' %par['beta'])
        
    h.update('DISPLACE',   '%s/coord_list.dat' %newdir, 'photometry file for x and y position')
    h.update('CADENCE',   '%f' %cadence, 'frequency of position update in hz')
    h.update('INTEGRAT', '%s' % par['nsteps'], 'number of integrations')
    exposure =(float(par['nsteps'])*exposure)
    h.update('EXPOSURE', '%f' %exposure, 'exposure in seconds')
    h.update ('NSTARS' , '1', 'number of stars used')
        
        # save to a file, the writeto method will make sure the required
        # keywords are conforming to the data
        
    notes1 = 'if IMTYPE is LI the coordinate refers tot he location of the brightest pixel within a restricted area (typically 25 pix radius) centered on the position of the target at the previous time step. one star is used. coordinate file format is #file   x y  brightest-pixel-counts ----------'
    notes2 = 'if IMTYPE is HSP sextractor and iraf photometry phot package are used to derive x and y position. more then one star can be used.  coordinate file format is #image-index-in-spool  \[x1 y1 flux1 normalized-flux1]*number of stars -----'
    notes =par['notes']

    h.update('REDUCTN', '%s' %(notes1+notes2))
    h.update('NOTES', '%s' %(par['notes']))
    fitsobj.append(hdu)
    

    fname = '%s/psf_%s_%3.1fs.fits'%(newdir,profile,exposure)
    print 'writing fits file to %s'%fname
    if os.path.isfile(fname):
        strg = "rm %s"%fname
        os.system(strg)
    fitsobj.writeto(fname)


###################################################main#######################
def centan(outpath,dispfile, par, nstars, nameroot, newdir):
    from pyfits import open as pfopen
    from pylab import *
    
    if os.path.isfile(dispfile) == 0:
        print "no strehl analysis file ",dispfile,". run analysis first!"
        return -1
    f=open(dispfile,'r')
    allcoordslist=[]
    skip = int(par['nskip'])
    nsteps = int(par['nsteps'])

    #####     HEADER INFO     #####
    firstfits = '%s/unspooled/%s_%05d.fits' %(outpath,nameroot,skip) 
    image=pfopen(firstfits)
    header=image[0].header
    image.close()
    if  'HBIN' in header:
        pixscalex =  float(par['ps'])*float(header['HBIN'])
        pixscaley =  float(par['ps'])*float(header['VBIN'])
    elif 'CCDXBIN' in header:
        pixscalex =  float(par['ps'])*float(header['CCDXBIN'])
        pixscaley =  float(par['ps'])*float(header['CCDYBIN'])


    if 'EXPOSURE' in header:
        exposure =  float(header['EXPOSURE'])
    elif 'EXPTIME' in header:
        exposure =  float(header['EXPTIME'])
    else:
        print "no exposure lenght recognizable key!" 
        return -1

    if 'KCT' in header:
        cadence = float(header['KCT'])
    else: 
        cadence = 1.0/exposure

    if 'FRAME' in header:
        date = header['FRAME']
    elif 'DATE' in header:
        date = header['DATE']

    PSFg=float(par['psf'])
    PSFm=float(par['alpha'])

    nx,ny=100,100
    c=(50,50)
    profile=par['profile'] #'g' for gaussian, 'm' for moffat
    if profile=='g':
        PSF = PSFg
    elif profile == 'm':
        PSF=PSFm        
    else:
        print "unknown profile"
        return -1



    TENSIG=min(int(PSF/pixscalex*5),c[0])
    x,y=arange(nx),arange(ny)
     
    for i in f:
        if i.startswith('#'): continue
        i=i.split()
        allcoordslist.append([i[0], float(i[1]), float(i[2]), 
                              float(i[3]), float(i[4]), float(i[5]), 
                              float(i[6]), float(i[7]), float(i[8])])


        allcoords=sorted(allcoordslist,key=lambda list:list[0])


    if skip>0:
        allcoords=allcoords[skip:]
    
    coordfile = "%s/coord_list.dat"%(newdir)
    f= open(coordfile,'w')
    print >> f ,"#fname 	dx(pix)  dy(pix) dx(arcsec) dy(arcsec) flux(counts, aperture) x(pix,aperture) y(pix, aperture) x(pix, maxflux), y(pix, maxflux) nrightest pixel(counts)"
    x0, y0 = allcoords[0][6],allcoords[0][7]
    for l in allcoords:
        dx=float(l[6])-float(x0)
        dy=float(l[7])-float(y0)
        print >>f, l[0],dx,dy,dx*pixscalex,dy*pixscaley,l[8],\
            l[6],l[7],l[4],l[5],l[3]


#print zip(*allcoordslist)[4]

    mux = []
    muy = []

    for i in xrange(nstars):
        dx=array(zip(*allcoords)[6])
        dy=array(zip(*allcoords)[7])
    
        mux.append(array(dx[:nsteps]-dx[0]+c[1]))
        muy.append(array(dy[:nsteps]-dy[0]+c[0]))
    mx= mean(mux,0)
    my= mean(muy,0)


    xindex = arange(len(dx))

    plt.figure()

#fname = '%s/%s/%s_dx.png'%(LIDIR,par['fits'],par['fits'])

#savefig(fname,dpi=None, facecolor='w', edgecolor='w',
#        orientation='portrait', papertype=None, format=None,
#        transparent=False, bbox_inches=None, pad_inches=0.1)
    subplot(2,1,1)
    plt.xlabel('time (seconds)')
    plt.ylabel('displacement (arcseconds)')
    plt.ylabel('dx (arcseconds)')
    plot (xindex*cadence,(dx-dx[0])*pixscalex, 'o-',label='x')
    subplot(2,1,2)
    plt.ylabel('dy (arcseconds)')
    plot (xindex*cadence,(dy-dy[0])*pixscaley, 'o-',label='y')
    legend(loc=1, ncol=1, shadow=True)
    fname = '%s/%s_dxdy.png'%(newdir,nameroot)
    savefig(fname,dpi=None, facecolor='w', edgecolor='w',
            orientation='portrait', papertype=None, format=None,
            transparent=False, bbox_inches=None, pad_inches=0.1)


    plt.figure()
    plt.xlabel('dx (arcseconds)')
    plt.ylabel('dx (arcseconds)')

#fname = '%s/%s/%s_dx.png'%(LIDIR,par['fits'],par['fits'])

#savefig(fname,dpi=None, facecolor='w', edgecolor='w',
#        orientation='portrait', papertype=None, format=None,
#        transparent=False, bbox_inches=None, pad_inches=0.1)

    plot ((dx-dx[0])*pixscalex,(dy-dy[0])*pixscaley, 'o')
#    legend(loc=1, ncol=1, shadow=True)
    fname = '%s/%s_dxvsdy.png'%(newdir,nameroot)
    savefig(fname,dpi=None, facecolor='w', edgecolor='w',
            orientation='portrait', papertype=None, format=None,
            transparent=False, bbox_inches=None, pad_inches=0.1)



    plt.figure()
    xfft=fft((dx-dx[0])*pixscalex)
    yfft=fft((dy-dy[0])*pixscaley)
    nxfft=len(xfft)
    nyfft=len(yfft)
    powerx = abs(xfft[1:(nxfft/2)])**2
    powery = abs(yfft[1:(nyfft/2)])**2
    nyquist=1./2
    freqx=array(range(nxfft/2))/(nxfft/2.0)*nyquist
    freqy=array(range(nyfft/2))/(nyfft/2.0)*nyquist
    periodx=1./freqx
    periody=1./freqy
    plt.xlabel('period of x and y oscillations [seconds]')
    plt.ylabel('power')
    plot(periodx[1:len(periodx)/2], powerx[0:len(powerx)/2], 'o-',label='x')
    plot(periody[1:len(periody)/2], powery[0:len(powery)/2], 'o-',label='y')
#    plt.xlim(0,max(periodx)/2)
#    xaxis((0,40))
    fname = '%s/%s_fft.png'%(newdir,nameroot)
#    show()
    savefig(fname,dpi=None, facecolor='w', edgecolor='w',
            orientation='portrait', papertype=None, format=None,
            transparent=False, bbox_inches=None, pad_inches=0.1)


    coords1 = array([ zeros(2,float) for i in xrange(nsteps) ]).reshape(nsteps,2)
    coords2 = array([ ones(2,float)*50 for i in xrange(nsteps) ]).reshape(nsteps,2)
    for i in range(nsteps):
        coords1[i][0] = mx[i]
        coords1[i][1] = my[i]
#        coords2[i][0] *=c[0]
#        coords2[i][1] *=c[1]

        
    mkfits(par, coords1, coords2,newdir,date, PSF, TENSIG, pixscalex, pixscaley, cadence,x, y,nx, ny, c, exposure)
    strg = 'cp %s/unspooled/%s_%05d.fits %s'%(outpath, nameroot,skip,newdir)
    os.system(strg)
#    os.chdir(olddir)
#    os.system(strg)
#    strg = 'tar -czvf %s.tgz %s_displacement'%(newdir,nameroot)
#    print strg
#    os.system(strg)
    return 1


if __name__ == '__main__':
    if len(sys.argv) != 2 or sys.argv[1].startswith('-h') or sys.argv[1] == 'h':
        print """Usage. Requires: 
                **name of parameter file conatining :**
                
                Directory containing images
            	#'y' for using displacement, 'n' for just integration
                'disp' : 'y',
                #target coordinates (optional)
                'ra' : '',\
                'dec' : '',
                'profile' : 'm',\
                'alpha' : 1.4,\
                'beta' : 3.0,\
                'psf' : 0.7,\
                #number of steps to use in the psf reconstruction
                'nsteps' : 100,\
                #number of steps images to skip
                'nskip':0,\
                #telescope
                'scope' : 'FTN'
		dark method
		"""
    	sys.exit()


#####     DECLARE VARIABLES     #####
    from mymkdir import mymkdir
    par = readconfig(sys.argv[1])
    print par
    olddir = '%s/%s/' %(LIDIR,par['spool'][0])
    newdir = '%s/%s/%s_displacement' %(LIDIR,par['spool'][0],par['spool'][0])
 
    if mymkdir(newdir)!=0:
        sys.exit(0)
#    strg = 'mkdir %s'%newdir
#    os.system(strg)

    
    dispfile = "%s/%s/strehl_list.dat"%(LIDIR,par['spool'][0])


    centan(doutpath,dispfile, par, 1,nameroot, newdir)
