import pyfits as PF
import sys
import os
from numpy import *
from scipy import *
#sys.path.append("../LIHSPcommon")

#inpath: filepath of input dark spool
#inname: filename of input dark spool
#method: avg med clip[,kappa]
#outpath: filepath of masterdark
#makefits: 1 prints masterdark to .fits


def rflat(filepath,filename):    
    filein = filepath + '/' + filename
    if filename.find('mstrflat') < 0:
        if filename.endswith('.fits'):
            filename=filename[:-5]
        filein=filepath+'/'+filename+'_mstrflat.fits'
    flat=PF.open(filein)
    return flat[0].data

def mkflat(inpath,filenames,masterdark,dexp,outpath):
    from myutils import mymkdir
    from mdark import rdark

    fileout=outpath+'/'+filenames[0].replace('.fits','_mstrflat.fits')
    try:
        masterflat = rflat(outpath,filenames)
        print "Read masterflat :", fileout
    except:
        nimages=len(filenames)
        flat=PF.open(inpath+'/'+filenames[0])
        xarray=flat[0].header['NAXIS1']
        yarray=flat[0].header['NAXIS2']
        cube=zeros((nimages,yarray,xarray),float)
        
        print "flat 1: ", filenames[0]
        header=flat[0].header

        if  'HBIN' in header:
            binning=(float(header['HBIN']),float(header['VBIN']))
        elif 'CCDXBIN' in header:
            binning=(float(header['CCDXBIN']),float(header['CCDYBIN']))
        if binning == (1,1):
                bias=433.5
        elif binning ==(2,2):
            bias=385.5
        else:
            bias=400.0
            print "WARNING: for this binning the bias is just a guess..."

        cube[0]=(flat[0].data-bias)/header['EXPTIME']-(masterdark-bias)/dexp
        i=0
        for f in filenames[1:]:
            print "flat ",str(i),": ", filenames[0]
            i=i+1
            flat=PF.open(inpath+'/'+f)
            header=flat[0].header
            cube[i]=flat[0].data/header['EXPTIME']-masterdark/dexp
        masterflat=zeros((yarray,xarray),int)
    
    
        print 'Taking median with 3-sigma clipping...'
        header.update('FRMTYPE','FLAT')
        header.update('EXPTIME',1.0)
        header.update('METHOD','SIG-CLIP MED','Method of frame combination', after="FRMTYPE")
        sig=std(cube, axis=0)
        mn=cube.mean(axis=0)
        kap=3
        masterflat=zeros((yarray,xarray), float)
        for i in range(yarray):
            for j in range(xarray):
                maskhi=cube[:,i,j]>(mn[i,j]-kap*sig[i,j])
                masklo=cube[:,i,j]<(mn[i,j]+kap*sig[i,j])
                masterflat[i,j]=median(cube[maskhi&masklo,i,j])
        masterflat/=median(masterflat)
        header.update('CLIPPING',kap,'Kappa coefficient of clipping',after='METHOD')
#####     WRITE MASTERFLAT TO .FITS     #####
        
        print 'printing'
        
        PF.writeto(fileout, masterflat, header, clobber=True)               #creates master.fits at fileout
        print 'Master flat written to %s' % fileout
    
    return masterflat
