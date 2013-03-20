################################################################
###############written by dkudrow 08/2010########################
#Reads all of the images in a directory, sorts them by Strehl
#Ratio or consecutive order, and stacks the specified top percent
#either aligned to a guide star or not, weighting if requested.
#This version includes: 
#detripling for compact binaries, drizzle algorithm for 
#subpixel alignment  and dynamic guide region.
#Type $python stacked.py -h for usage
####
#last modified 10/26/2010
#################################################################

print 'Initializing stacking sequence...'
SECINDAY=1.15740741e-5
PLOTIT = 1
import sys,os,time
import pyfits as PF
from numpy import loadtxt, float64, array
from scipy import *
from scipy.interpolate import interp1d
from pylab import figure,subplot,plt,plot,savefig, legend, title
from myutils import readconfig,mymkdir, mjd, mkarray
import types

print 'Loaded all packages'
SATLEVEL= 1400000
ps = 0.13
def calccents(tmparr):
    if shape(tmparr)==(0,0):
        print "empty array passed to calccents, returning (0,0,0)"
        return (0,0,0)
    allx =sum(tmparr,axis=0)
    xl = len(allx)
    indx = arange(xl)
    allf = sum(allx)
    mx = sum(allx*indx)/allf

    ally = sum(tmparr,axis=1)
    yl = len(ally)
    indy = arange(yl)
    my = sum(ally*indy)/allf
    
    sky =sum(tmparr[0,:])
    sky+=sum(tmparr[-1,:])
    sky+=sum(tmparr[:,0])
    sky+=sum(tmparr[:,-1])
    skystd = std (concatenate((tmparr[0,:],tmparr[-1,:],tmparr[:,0],tmparr[:,-1]),axis=0))
    (lx,ly)=shape(tmparr)
    sky /=(2.0*(lx+ly))
    skyall=sky*(lx*ly)
    allf=allf-skyall
    #returning: weighted x centroid, weighted y centroid, sky value
    return (mx,my, allf, sky, skystd)



def myphot(inpath, coordfile, last, rad, follow, nameroot,target):
    files0 = sorted([x for x in os.listdir(inpath) if (x.endswith('.fits') and nameroot+'_' in x)])
    files1 = sorted([x for x in os.listdir(inpath) if (x.endswith('.fits') and nameroot+'._' in x)])
    files=files0+files1
#    target = target-1
    nfiles=len(files)
    if last < nfiles:
        nfiles = last
        files=files[:nfiles]
    image=PF.open(inpath+'/'+files[0])
    header=image[0].header
    image.close()
    xsize=header['NAXIS1']
    ysize=header['NAXIS2']
    if  'HBIN' in header:
        hbin =float(header['HBIN'])
        vbin =float(header['VBIN'])
    elif 'CCDXBIN' in header:
        hbin =float(header['CCDXBIN'])
        vbin =float(header['CCDYBIN'])
    if 'EXPOSURE' in header:
        exp =  (header['EXPOSURE'])
    elif 'EXPTIME' in header:
        exp =  (header['EXPTIME'])
    else:
        print "no exposure lenght recognizable key!" 
        return -1

    if 'MJD' in header:
        tstart = float(header['MJD'])
    elif 'FRAME' in header:
        tmpdate = header['FRAME']
        year = int(tmpdate[0:4])
        month = int(tmpdate[5:7])
        day = int(tmpdate[8:10])    
        hour = float(tmpdate[11:13])/24.0
        minute = float(tmpdate[14:16])/1440.0
        second = float(tmpdate[17:])/86400

        
        tstart=mjd(year,month,day)+hour+minute+second
        

    psx = 2.0/(ps*hbin)
    psy = 2.0/(ps*vbin)
#value in pixel of 2''*2'' to contain centroid update within 2'' , x and y
    psxsq=psx*psx
    psysq=psy*psy

    coords=loadtxt(coordfile, dtype='float')
    coords=coords.transpose()
    if type(coords[0]) is float64:
        bpx = [coords[0]]
        bpy = [coords[1]]
    elif len(coords[0])>1:
        bpx = coords[0]
        bpy = coords[1]
    else: 
        print "missing coordinates in file ",coordfile
    
    if 0 in bpx: bpx[0]=int(xsize/2.0)#where(bpx==0)[0]] = int(xsize/2.0)
    if 0 in bpy: bpy[0]=int(xsize/2.0)#where(bpy==0)[0]] = int(xsize/2.0)



    nstars = len(bpx)
    newbpx = zeros(nstars, float)
    newbpy = zeros(nstars, float)

    print bpx, bpy
    if type(rad) is float or type(rad) is int:
        rad = [rad]
    nrad=len(rad)
    if 'clean' in inpath:
        regdir='%s/../../regions' %(inpath)
    else: regdir='%s/../regions' %(inpath)
    if mymkdir(regdir)!=0:
        return -1

    print "calculating aperture photometry for ",nstars," stars and ",nrad,"apertures"
    phot=[]
    for r in rad:
        phot.append([])
    for name in files:
        tmpfile=PF.open(inpath+'/'+name)
        

        for j,r in enumerate(rad):
            regfile = '%s/%s_%d_%s.reg' %(regdir,name,r,follow[0:3])
#        print "\n\n printing region file for fits ", name,"to ",regfile,"\n\n"
            regf = open(regfile,"w")
            print >> regf, 'global color=green dashlist=8 3 width=1 font="helvetica 10 normal" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1'

            for i in xrange(nstars):
                try:
                    reg=tmpfile[0].data[bpy[i]-r:bpy[i]+r+1,bpx[i]-r:bpx[i]+r+1]
                    if reg==[]:
                        print "something is wrong: did you loose the centroid??"
                        mx,my,apf,sky,skystd=float("nan") ,float("nan") ,float("nan") ,float("nan") ,float("nan")
                    else:
                        maxvalue =reg.max()
                        maxindex =where(reg==reg.max())
                        if maxvalue>=SATLEVEL or len(maxindex[0]) >1:
                            print 'Saturated or repeat brightest pixel found in %s. ' % name
                            sat = 1
                        (mx,my,apf,sky,skystd) =\
                                               calccents(
                            tmpfile[0].data[bpy[i]-r:
                                            bpy[i]+r+1,
                                            bpx[i]-r:
                                            bpx[i]+r+1])
                    
                        newbpx[i]=bpx[i]+mx-r
                        newbpy[i]=bpy[i]+my-r
                except:
                    print "something is wrong: did you loose the centroid??"
                    mx,my,apf,sky,skystd=float("nan") ,float("nan") ,float("nan") ,float("nan") ,float("nan")
            #            print bpx[i],bpy[i],name, reg.shape
                sat = 0
                largeoff=0

                #            print bpx,bpy,mx,my



                print >> regf, "physical"
                print >> regf, "box(",bpx[i]+1,",",bpy[i]+1,",",2*r+1,",",2*r+1,",0)"
                print >> regf,"circle(",int(bpx[i])+1, ",",int(bpy[i])+1,",1)"
                print >> regf,"circle(",int(bpx[i])-r+mx+1,",",int(bpy[i])-r+my+1,",1) # color=red"
            

            
#            print newbpx, newbpy, bpx, bpy, psxsq
            
                if follow=='dynamic':
                    if (newbpx[i]-bpx[i])*(newbpx[i]-bpx[i])<psxsq and (newbpy[i]-bpy[i])*(newbpy[i]-bpy[i])<psysq: 
                
                        if newbpx[i]-2*r > 0 or newbpx[i]+2*r < xsize \
                                or newbpy[i]-2*r > 0 or newbpy[i]+2*r < ysize \
                                or (bpx[i]-newbpx[i])*(bpx[i]-newbpx[i]) > r*r/4.0 \
                                or (bpy[i]-newbpy[i])*(bpy[i]-newbpy[i]) > r*r/4.0  :
                        
                            bpx[i]=newbpx[i]
                            bpy[i]=newbpy[i]
                        else:
                            print 'Cannot update centroid cause it is outside of the allowed region',name
                    else:
                        print 'Cannot update centroid cause it is too far from previous datapoint',name
                        largeoff = 1
        
                phot[j].append({'name':name,
                               'x':bpx[i]-r+mx, 
                               'y':bpy[i]-r+my, 
                               'flux':apf, 
                               'flux_norm':apf, 
                               'sky':sky, 
                               'skystd':skystd,
			       'satflag':sat,
                               'offflag':largeoff})


            regf.close()
            tmpfile.close()

    for j,r in enumerate(rad):
        phot[j]=array(phot[j])
        newshape=(nfiles,nstars)
        phot[j].flatten()

        phot[j]=reshape(phot[j],newshape)


    
#    print len (phot[0]), nfiles*nstars
#    reshape (phot[0],(nfiles,nstars))
#    print phot[0]
#    means = mkarray(0.0, len(rad), nstars)
#    stds = mkarray(0.0, len(rad), nstars)
#    fluxs = mkarray(zeros(nfiles), len(rad), nstars)

    means = zeros(len(rad)*nstars,float)
    means=reshape (means,(len(rad),nstars))

    stds = zeros(len(rad)*nstars,float)
    stds=reshape (stds,(len(rad),nstars))

    ws = zeros(len(rad)*nstars,float)
    ws=reshape (ws,(len(rad),nstars))

    fluxs = zeros((len(rad)*(nstars+2)*nfiles),float)
    fluxs=reshape (fluxs,(len(rad),nstars+2,nfiles))

    xs = zeros((len(rad)*nstars*nfiles),float)
    xs=reshape (xs,(len(rad),nstars,nfiles))

    ys = zeros((len(rad)*nstars*nfiles),float)
    ys=reshape (ys,(len(rad),nstars,nfiles))

    ts = arange(nfiles)
    tsmjd = ts*SECINDAY+tstart

    for j,r in enumerate(rad):
        for k in xrange(nfiles):
            for i in xrange(nstars):
                fluxs[j][i][k] = phot[j][k][i]['flux']
                xs[j][i][k] = phot[j][k][i]['x']
                ys[j][i][k] = phot[j][k][i]['y']
            if target > -1:
                fluxs[j][nstars][k] = phot[j][k][target-1]['flux']

        for i in xrange(nstars):
            means[j][i]= mean(fluxs[j][i])
#            means[j][i]=
            stds[j][i]= std(fluxs[j][i])
            ws[j][i]=1.0/(stds[j][i]*stds[j][i])

            if target>-1:
                if target!=i+1:
                    fluxs[j][nstars+1]+=(fluxs[j][i]/means[j][i])*ws[j][i]
                else:
                    fluxs[j][nstars] /= means[j][target-1] 
        if target>-1:
            fluxs[j][nstars+1]/=sum(ws[j])
            print sum(ws[j])
    print fluxs[0][nstars]
    print fluxs[0][nstars+1]
   
#    plot (ts,fluxs[0][nstars])
#    plot (ts,fluxs[0][nstars+1])
#    show()
    
    
    print "means", means
    print "standard deviations", stds
    print "SNR: ", means/stds

        

    for k,r in enumerate(rad):
        if 'clean' in inpath:
            apfilename = '%s/../../photometry_%d_%s.dat' %(inpath, r, follow[0:3])

        else :        
            apfilename = '%s/../photometry_%d_%s.dat' %(inpath, r, follow[0:3])

        print "here", inpath

        apfile = open(apfilename,'w')
        print >> apfile, "#name    MJD    time(sec after MJD:)",tstart,
        for i in xrange(nstars):
            print >> apfile,"x y flux normalized-flux sky std satration offset",
        print >>apfile, ""

        for j in xrange(nfiles):
            timestamp = j*exp
            print >> apfile, phot[k][j][0]['name'], tstart+timestamp*SECINDAY, timestamp,
            for i in xrange(nstars):
                phot[k][j][i]['flux_norm']/=means[k][i]
                phothere = '%f %f %f %f %f %f %d %d '%(phot[k][j][i]['x'],phot[k][j][i]['y'],phot[k][j][i]['flux'],phot[k][j][i]['flux_norm'],stds[k][i], phot[k][j][i]['sky'],phot[k][j][i]['satflag'],phot[k][j][i]['offflag'])
                print >> apfile, phothere,
            if target>-1:
                phothere = '%f '%(fluxs[k][nstars][j]/fluxs[k][nstars+1][j])
                print >> apfile, phothere
            else:
                print >>apfile, ''
                       


    if PLOTIT==1:
        for k,r in enumerate(rad):
            if 'clean' in inpath:
                plotdir='%s/../../plots_ap%.2f' %(inpath, r)
                if mymkdir(plotdir)!=0:
                    return -1
            else :        
                plotdir='%s//../plots_ap%.2f' %(inpath, r)
                if mymkdir(plotdir)!=0:
                    return -1

            for i in xrange(nstars):
                if follow=='dynamic':
                    fig=figure()
                    ax1 = subplot(111)
                    
                    plt.ylabel('displacement (pixels)')
                
                
                    print "plotting ",r,i
                    lab = "dx %.2f %.2f, ap %.2f"%(xs[k][i][0],ys[k][i][0],r)
                    plot(ts,xs[k][i]-xs[k][i][0], label=lab)
                    lab = "dy %.2f %.2f,ap %.2f"%(xs[k][i][0],ys[k][i][0],r)
                    plot(ts,ys[k][i]-ys[k][i][0], label=lab)
                    ax1.set_xlabel('displacement ( pixels)')
                    ax1.set_xlabel('time (seconds)')
                    
                    legend(loc=1, ncol=1, shadow=True)
                    
                    fname = '%s/displacement_%0.f_%0.2f_%d.png' %(plotdir,xs[k][i][0],ys[k][i][0] , r)
                    savefig(fname,dpi=None, facecolor='w', edgecolor='w',
                            orientation='portrait', papertype=None, format=None,
                            transparent=False, bbox_inches=None, pad_inches=0.1)

            for i in xrange(nstars):
                fig=figure()
                ax1 =subplot(111)
                lab = "flux %.2f %.2f, %.2f\n mean: %.2f, std: %.2f, SNR: %.2f"%(xs[k][i][0],ys[k][i][0],r, means[k][i], stds[k][i], means[k][i]/stds[k][i])
                plot(ts, fluxs[k][i], label=lab)
                #            legend(loc=1, ncol=1, shadow=True)
                title(lab)
                #            ax2 = twiny()
                ax1.set_ylabel('flux (counts)')
                ax1.set_xlabel('time (seconds)')
                fname = '%s/photometry_%0.2f_%0.2f_%d_%s.png' %(plotdir, xs[k][i][0],ys[k][i][0], r, follow[0:3])
                
                savefig(fname,dpi=None, facecolor='w', edgecolor='w',
                        orientation='portrait', papertype=None, format=None,
                        transparent=False, bbox_inches=None, pad_inches=0.1)

                fig=figure()
                ax1 =subplot(111)

                plot(ts, fluxs[k][i]/means[k][i], label=lab)
                title(lab)
                
                #            ax2 = twiny()
                ax1.set_ylabel('normalized flux')
                ax1.set_xlabel('time (seconds)')
                fname = '%s/photometry_%0.2f_%0.2f_%d_%s_norm.png' %(plotdir, xs[k][i][0],ys[k][i][0], r, follow[0:3])
                
                savefig(fname,dpi=None, facecolor='w', edgecolor='w',
                        orientation='portrait', papertype=None, format=None,
                        transparent=False, bbox_inches=None, pad_inches=0.1)
                
            if target > -1:
                fig=figure()
                ax1 =subplot(111)
                
                lab = "flux %6.2f %6.2f, %6.2f, normalized"%(xs[k][i][0],ys[k][i][0],r)
                plot(ts,fluxs[k][nstars],label=lab)
                lab = "flux %6.2f %6.2f, %6.2f, corrected"%(xs[k][i][0],ys[k][i][0],r)
                
                plot(ts, fluxs[k][nstars]/fluxs[k][nstars+1], label=lab)
                
                legend(loc=1, ncol=1, shadow=True)
                
                #            axrfluxs2 = twiny()
                
                ax1.set_ylabel('flux (counts)')
                ax1.set_xlabel('time (seconds)')
                fname = '%s/differential_photometry_%0.2f_%0.2f_%d_%s.png' %(plotdir, xs[k][i][0],ys[k][i][0], r, follow[0:3])
                
                savefig(fname,dpi=None, facecolor='w', edgecolor='w',
                        orientation='portrait', papertype=None, format=None,
                        transparent=False, bbox_inches=None, pad_inches=0.1)
                
                
                
#    plt.xlabel('time (seconds)')
#    plt.ylabel('')
#    plot (xindex*par['cadence'],(dx-dx[0])*0.16, label='x')
#fname = '%s/%s/%s_dx.png'%(LIDIR,par['fits'],par['fits'])

#savefig(fname,dpi=None, facecolor='w', edgecolor='w',
#        orientation='portrait', papertype=None, format=None,
#        transparent=False, bbox_inches=None, pad_inches=0.1)

    

    return 1

##########################################################################################################################################################


if __name__ == '__main__':
    if len(sys.argv) != 2 or sys.argv[1].startswith('-h') or sys.argv[1] == 'h':
        print """Usage. Requires: 
                **name of parameter file conatining :**
                
                #original spool name 
                'spool':'whatever.fits'
                #Directory containing images
                #target coordinates (optional)
                'coords':'coords.tmp',\
                #number of steps images to skip
                'nskip':0,\
		"""
    	sys.exit()


#####     DECLARE VARIABLES     #####
    from mymkdir import mymkdir,mygetenv
    par = readconfig(sys.argv[1])
    SPEEDYOUT=mygetenv('SPEEDYOUT')
#    strg = 'mkdir %s'%newdir
#    os.system(strg)

    if (cosmic):
	if (isinstance(par['spool'],types.StringTypes)):
            nameroot=par['spool']
            if nameroot.endswith('.fits'):
        	nameroot = nameroot[:-5]
            inpath = '%s/%s//clean/'%(SPEEDYOUT,nameroot)

	else:
            nameroot=par['spool'][0]
            if nameroot.endswith('.fits'):
        	nameroot = nameroot[:-5]
	
            inpath = '%s/%s_all/clean/'%(SPEEDYOUT,nameroot)
	
    else:
	if (isinstance(par['spool'],types.StringTypes)):
            nameroot=par['spool']
            if nameroot.endswith('.fits'):
        	nameroot = nameroot[:-5]
            inpath = '%s/%s/unspooled/'%(SPEEDYOUT,nameroot)

	else:
            nameroot=par['spool'][0]
            if nameroot.endswith('.fits'):
        	nameroot = nameroot[:-5]
            inpath = '%s/%s_all/unspooled/'%(SPEEDYOUT,nameroot)
	
    myphot(inpath, par['coords'], int(par['last']), par['ap'], follow, nameroot, par['target'])


