
#################################################################
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

import sys,os,time
import pyfits as PF
from numpy import *
from scipy import *
from scipy.interpolate import interp1d
from scipy.fftpack import fftn, ifftn
from pylab import plt, axvline, savefig, subplot,  figure, plot, legend, title, hist
from myutils import readconfig
import mysciutils
sys.path.append("../LIHSPcommon")

print 'Loaded all packages'
satlevel = 1400000

def calccents(tmparr):
    if shape(tmparr)==(0,0):
        print "empty array passed to calccents, returning (0,0,0)"
        return (0,0,0)
    allx =sum(tmparr,axis=0)
    xl = len(allx)
    indx = arange(xl)
    allf = sum(allx)
    #print 'allx.shape = ' + str(allx.shape)
    #print 'indx.shape = ' + str(indx.shape)
    mx = sum(allx*indx)/allf

    ally = sum(tmparr,axis=1)
    yl = len(ally)
    indy = arange(yl)
    my = sum(ally*indy)/allf
    
    sky =sum(tmparr[0,:])
    sky+=sum(tmparr[-1,:])
    sky+=sum(tmparr[:,0])
    sky+=sum(tmparr[:,-1])
    (lx,ly)=shape(tmparr)
    sky*=(lx*ly)/(2.0*(lx+ly))
    allf=allf-sky
    #returning: weighted x centroid, weighted y centroid, sky value
    return (mx,my, allf)

def createstack(inpath,gsx,gsy,rad,select,pc,shift,detrip,minsep,outpath,coresz,follow,ps):

    print "searching for files in %s"%inpath
    files = sorted([x for x in os.listdir(inpath) if x.endswith('.fits')])
    nfiles=len(files)
    if type(pc) is float or type(pc) is int:
        pc = [pc]
        
    cutoff = zeros(len(pc), float)
    for i,c in enumerate(pc): 
        cutoff[i] = int(c/100.0*nfiles)
        if cutoff[i]==0:
            cutoff[i]=1
#####     HEADER INFO     #####
            
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


    psx = 2.0/(ps*hbin)
    psy = 2.0/(ps*vbin)
#value in pixel of 2''*2'' to contain centroid update within 2'' , x and y
    psxsq=psx*psx
    psysq=psy*psy

    align='UNALIGN'
    if shift == 'align':
        align='ALIGN'

    if select=='lucky':
        imtype='LUCKY'
        if shift == 'none':
            print """Will not perform lucky imaging without aligning frames
                 Changing parameter 'shift' to '1'"""
            shift = 'align'
            align='ALIGN'
    elif select == 'coadded':
        imtype='CONSEC'
    elif select == 'weighted':
        imtype='WEIGHTED_TIPTILT'
        if shift == 'none':
            print """Will not perform weighted tip/tilt imaging without 
                 aligningframes
                 Changing parametes 'shift' to '1'"""
            shift = 'align'
            align='ALIGN'
    elif select == 'corr' or select == 'correlate':
	imtype='LUCKY_CORRELATED'
	#NOTE: do I need to change shift and align? 

#	weights = []
#	for coff in cutoff:
#        	weights.append(zeros(coff,float))
    else:
        print "invalid 'select' paramter (%s): use 'coadded' for standard stacking, 'lucky' for lucky, 'weighted' for weighted average exposures, 'corr' for lucky with phase correlation"%select
        return -1  

    if gsx == 0:
        gsx = int(xsize/2)
    if gsy == 0:
        gsy = int(ysize/2)

    for pcn,coff in enumerate(cutoff):
        coff=int(coff)
	print "Generating ", imtype, " image using ",coff," images, using", align,"method"
        if align == 'ALIGN': print "centroid: ",gsx," ",gsy
        if detrip!='none': print "detripling core ",coresz, "min distance",minsep

        sys.stdout.flush()


#####     SELECTION     #####

################creating name files that will be useful later
        outname=files[0][:-9].replace('.fits','')
        if detrip == 'v':
            outfile=outpath+'/'+outname+'_%3.2f_%s_%s_%s_%s.fits' % (float(pc[pcn]),imtype,align,follow, 'dv')
        elif detrip == 'h':
            outfile=outpath+'/'+outname+'_%3.2f_%s_%s_%s_%s.fits' % (float(pc[pcn]),imtype,align,follow, 'dh')
        else:
            outfile=outpath+'/'+outname+'_%3.2f_%s_%s_%s.fits' % (float(pc[pcn]),imtype,align,follow)
        histfile = outpath+'/'+outname+'_'+follow+'.hist.png' 
        listfilename = '%s/strehl_list_%s.dat' %(outpath,follow)


        best=[]
        if imtype == 'CONSEC' and align == 'UNALIGN':
            for name in files:
                best.append([name,
                             0,
                             0,
                             0,
                             gsx,
                             gsy,
                             gsx,gsy, 0])
                
        else:
            print '#Evaluating %s images...' % nfiles
            sys.stdout.flush()
            if os.path.isfile(listfilename) == True:
		print listfilename + ' exists'
                print 'reading list of strehl sorted files'
                listfile = open(listfilename)
                sys.stdout.flush()

                lsplit=listfile.readline().split()
                if len(lsplit)<9:
                    print 'old version of strehl ratio file. remove ', listfilename, 'first and rerun reduction.'
                    return -1
                best.append([lsplit[0], int(lsplit[1]), int(lsplit[2]), 
                             float(lsplit[3]), 
                             int(lsplit[4]), int(lsplit[5]),
                             lsplit[6], lsplit[7],lsplit[8]])
                for l in listfile:
                    lsplit=l.split()
                    best.append([lsplit[0], int(lsplit[1]), int(lsplit[2]), 
                                 float(lsplit[3]), 
                                 int(lsplit[4]), int(lsplit[5]),
                                 float(lsplit[6]), float(lsplit[7]),
                                 float(lsplit[8])])
                    
            else: #strehl list does not exist
                bpx,bpy=gsx,gsy
                nrejected = 0
                name = files[0]
                regfile = '%s/../%s_%s.reg' %(outpath,name,follow)
                print "\n\n printing region file for fits ", name,"to ",regfile,"\n\n"
                regf = open(regfile,"w")
                tmpfile=PF.open(inpath+'/'+name)
                reg=tmpfile[0].data[bpy-rad:bpy+rad+1,bpx-rad:bpx+rad+1]
		print reg
		print '[%d:%d,%d:%d]' %(bpy-rad,bpy+rad+1,bpx-rad,bpx+rad+1)

                maxvalue =reg.max()
                if maxvalue<satlevel:
                    maxindex =where(reg==reg.max())
                    if len(maxindex[0]) ==1:

                        (mx,my,apf) =\
                        calccents(tmpfile[0].data[bpy-rad+maxindex[0]-rad:\
						  bpy-rad+maxindex[0]+rad+1,\
						  bpx-rad+maxindex[1]-rad:\
						  bpx-rad+maxindex[1]+rad+1])
                        print >> regf, 'global color=green dashlist=8 3 width=1 font="helvetica 10 normal" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1'
                        print >>regf, "physical"
                        print >> regf, "box(",(maxindex[1])[0]+bpx-rad+1,",",(maxindex[0])[0]+bpy-rad+1,",",2*rad+1,",",2*rad+1,",0)"
                        print >> regf,"circle(",maxindex[1][0]+int(bpx)-rad+1, ",",maxindex[0][0]+int(bpy)-rad+1,",1)"
                        print >> regf,"circle(",bpx+maxindex[1][0]-2*rad+mx+1,",",bpy+maxindex[0][0]-2*rad+my+1,",1) # color=red"
                        regf.close()

                        best.append([name,
                                     maxindex[1],
                                     maxindex[0],
                                     maxvalue,
                                     int(bpx),
                                     int(bpy),
                                     bpx+maxindex[1][0]-2*rad+mx, 
                                     bpy+maxindex[0][0]-2*rad+my, apf])
                        
                        if follow=='dynamic':
                            print "dynamic"
                            #only updating the centroid if the new brightest pixel is within 2'' of the old one
                            newbpx=bpx+maxindex[1][0]-rad
                            newbpy=bpy+maxindex[0][0]-rad
                            if (newbpx-bpx)*(newbpx-bpx)<2.0/psx and (newbpy-bpy)*(newbpy-bpy)<2.0/psy: 

                                if newbpx-2*rad > 0 or newbpx+2*rad < xsize \
                                    or newbpy-2*rad > 0 or newbpy+2*rad < ysize \
                                    or (bpx-newbpx)*(bpx-newbpx) > rad*rad/4.0 \
                                    or (bpy-newbpy)*(bpy-newbpy) > rad*rad/4.0  :
                                
                                    bpx=newbpx
                                    bpy=newbpy
                            else:
                                print 'Cannot update centroid cause it is outside of the allowed region'
                    else:
                        print 'Repeat brightest pixel found in %s. Frame skipped.' % name
                else:
                    print 'Saturated  pixel found in %s. Frame skipped.' % name
                    tmpfile.close()

                for name in files[1:]:
		    #print 'Running image %s ...' %name
                    tmpfile=PF.open(inpath+'/'+name)
                        
                    reg=tmpfile[0].data[bpy-rad:bpy+rad+1,bpx-rad:bpx+rad+1]
                    maxvalue =reg.max()
                    if maxvalue<satlevel:
                        maxindex =where(reg==reg.max())
                        if len(maxindex[0]) ==1: 
                            (mx,my,apf) =\
                            mysciutils.calccents(
                                tmpfile[0].data[bpy-rad+maxindex[0][0]-rad:
                                                    bpy-rad+maxindex[0][0]+rad+1,
                                                bpx-rad+maxindex[1][0]-rad:
                                                    bpx-rad+maxindex[1][0]+rad+1])


                            best.append([name,
                                         maxindex[1],
                                         maxindex[0],
                                         maxvalue,
                                         int(bpx),
                                         int(bpy),
                                         bpx+maxindex[1][0]-2*rad+mx, 
					 bpy+maxindex[0][0]-2*rad+my, apf])
			    l = best[-1]
            		    #print '%s %d %d %f %d %d %f %f %f\n'\
                		#%(l[0],l[1],l[2],l[3],int(l[4]+l[1])-rad,
                    		#int(l[5]+l[2])-rad,l[6],l[7],l[8])
			    #raw_input()
                            if follow=='dynamic':
                                newbpx=bpx+maxindex[1][0]-rad
                                newbpy=bpy+maxindex[0][0]-rad
                       #         print newbpx, newbpy, bpx, bpy
                                if (newbpx-bpx)*(newbpx-bpx)<psxsq and (newbpy-bpy)*(newbpy-bpy)<psysq: 

                                    if newbpx-2*rad > 0 or newbpx+2*rad < xsize \
                                            or newbpy-2*rad > 0 or newbpy+2*rad < ysize \
                                            or (bpx-newbpx)*(bpx-newbpx) > rad*rad/4.0 \
                                            or (bpy-newbpy)*(bpy-newbpy) > rad*rad/4.0  :
                                    
                                        bpx=newbpx
                                        bpy=newbpy

                                else:
                                    print 'Cannot update centroid cause it is outside of the allowed region, if this is common choose another star or shrink the search region radius ', name
                        else:
                            print 'Repeat brightest pixel found in %s. Frame skipped.' % name
                    else:
                        print 'Saturated  pixel found in %s. Frame skipped.' % name
                    tmpfile.close()
                    
                    
            if imtype!='CONSEC':
                best=sorted(best,key=lambda list:list[3] ,reverse=True)
            else :
                best=sorted(best,key=lambda list:list[0] ,reverse=False)                
#        print best[:10]
        nrejected = nfiles-len(best)

        if coff > nfiles-nrejected: 
            coff=nfiles-nrejected
        luckies=best[:coff]
        exposure =  0.0
##########PRINTING STREHL LIST######################


        print '%d images selected.' % coff
        sys.stdout.flush()
        if imtype !='CONSEC':
            if os.path.isfile(listfilename) != True:
                print '#printing list of selected images to '+listfilename
                sys.stdout.flush()
                listout=open(listfilename,'w')
            

                for l in best:
                    strhere = '%s %d %d %f %d %d %f %f %f\n'%(
                        l[0],l[1],
                        l[2],l[3],
                        int(l[4]+l[1])-rad,
                        int(l[5]+l[2])-rad,
                        l[6],l[7],l[8])
                    
                    listout.write(strhere)
                listout.close()
            if select == 'weighted':
            #########creating weights as AS^2.4 according to Robert Tubbs###
                weights = array(zip(*luckies)[3])
                strehlnorm = 0.3/weights[0]
                weights = pow(weights*strehlnorm,2.4)
                a =  coff/sum(weights)
                weights = weights*a
###########CREATING HISTOGRAM##################

            if os.path.isfile(histfile) != 1:
                
                bins = linspace(best[-1][3],best[0][3],50)
                
                strehl=[i[3] for i in best]
                events, edges, patches = hist(strehl,bins, normed=False)
            
                lower = resize(edges, len(edges)-1)
                binmid =  lower + 0.5*diff(edges)

                lenbins=len(bins)-2
                allevents = sum(events)
                nevents = 0
                xsaved = 0
                for i in range(lenbins):
                    nevents +=events[lenbins-i]
                    if nevents>allevents/100.0 and xsaved == 0:
                        xsaved = binmid[lenbins-i]
                
                plt.xlabel('counts')
                plt.ylabel('frequency')
                axvline (xsaved)
                savefig(histfile,dpi=None, facecolor='w', edgecolor='w',
                        orientation='portrait', papertype=None, format=None,
                        transparent=False, bbox_inches=None, pad_inches=0.1)
#    show()
#    sys.exit()

########################################

#####     DETRIPLING     #####
        dtr='NONE'
        if detrip!='none':
            c=coresz/2
            print 'Separating cores...'
            sys.stdout.flush()
            for i in xrange(coff):
                core1=[0,0,0,0,0]
                core2=[0,0,0,0,0]
                name=luckies[i][0]
                tmpfile=PF.getdata(inpath+'/'+name)
                reg=tmpfile[gsy-rad:gsy+rad+1,gsx-rad:gsx+rad+1]
                for y in xrange(1,2*rad-1):
                    for x in xrange(1,2*rad-1):
                        core=reg[y-c:y+c+1,x-c:x+c+1]
                        if core.sum()>=core1[0]:
                            core1=[core.sum(),x,y,argmax(core)%coresz,argmax(core)/coresz]


                        if core.sum()==core1[0]:
                            continue
                        elif core.sum()>=core2[0] and sqrt((1.*x-core1[1])**2+(1.*y-core1[2])**2)>=minsep:
                            core2=[core.sum(),x,y,argmax(core)%coresz,argmax(core)/coresz]

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

#####     ALIGNMENT AND STACKING     #####

        print 'Compiling Lucky Image...'
        sys.stdout.flush()

        if select=='corr' or select=='correlate':
	    print 'Using phase correlation...'
            ntotal = coff
            print '\nlength of fits list: %d' %ntotal
            mastername = luckies[0][0]
            print 'masterfits: %s' %mastername
            try: masterdata = PF.getdata(inpath+'/'+mastername)
	    except:
		print """Cannot find image %s\nAborting...""" %(inpath+'/'+mastername)
		return -1
            
            masterdata = mysciutils.drizzle(masterdata, 0)
	    _PAD_ = 50
            coaddedtmp = mysciutils.pad(masterdata, _PAD_)
############ TUKEY WINDOW MASTER #####
	    lrg = int(coaddedtmp.shape[0] < coaddedtmp.shape[1]) #longest edge
	    alpha = double(masterdata.shape[lrg])/coaddedtmp.shape[lrg]
            print "alpha = %f" %alpha
            tukey = mysciutils.tukey2d(tuple(int(x*alpha-50) for x in coaddedtmp.shape), alpha)
            padding =  int(round((float(coaddedtmp.shape[0])-tukey.shape[0])/2))
            tukeypad = mysciutils.pad(tukey, padding)
            masterfft = fftn(coaddedtmp*tukeypad)

            for k in xrange(1,coff):
		print '\nRunning %d/%d' %(k, coff-1)
                name = luckies[k][0]
		print '...' + name
		fits = PF.getdata(inpath+'/'+name)
		
		fits = mysciutils.drizzle(fits, 0)
		fitspad = mysciutils.pad(fits, _PAD_)

		lrg = int(fits.shape[0] < fits.shape[1]) #longest edge
	        alpha = double(fits.shape[lrg])/fitspad.shape[lrg]
        	tukey = mysciutils.tukey2d(tuple(int(x*alpha-50) for x in fitspad.shape), alpha)
        	padding =  int(round((float(fitspad.shape[0])-tukey.shape[0])/2))
        	tukeypad = mysciutils.pad(tukey, padding)
        	fitsfft = fftn(fitspad*tukeypad)

############## FINDING PHASE ############
		print 'Finding phase...'
		axis2_shift, axis1_shift = mysciutils.correlate([masterfft, fitsfft])
		
		if axis2_shift >= 0:
	            if axis2_shift == 0: axis2_shift = -fitspad.shape[0]
        	    if axis1_shift >= 0:
                	if axis1_shift == 0: axis1_shift = -fitspad.shape[1]
                	coaddedtmp[axis2_shift:,axis1_shift:] += fitspad[:-axis2_shift,:-axis1_shift]
            	    else: #axis1_shift < 0
                	coaddedtmp[axis2_shift:,:-abs(axis1_shift)] += fitspad[:-axis2_shift,abs(axis1_shift):]
        	else: #axis2_shift < 0
            	    if axis1_shift >= 0:
                	if axis1_shift == 0: axis1_shift = -fitspad.shape[1]
                	coaddedtmp[:-abs(axis2_shift),axis1_shift:] += fitspad[abs(axis2_shift):,:-axis1_shift]
            	    else: #axis1_shift < 0
                	coaddedtmp[:-abs(axis2_shift),:-abs(axis1_shift)] += fitspad[abs(axis2_shift):,abs(axis1_shift):]

	    stack = coaddedtmp

        else: 
            dx,dy=rad,rad
            stack=zeros((2*(ysize-2)-7*rad,2*(xsize-2)-7*rad))
            for k in xrange(coff):
		print 'Running %d/%d...' %(k, coff-1)
                name=luckies[k][0]
                tmp=PF.getdata(inpath+'/'+name)
            
 	        frame = mysciutils.drizzle(tmp, 0)
            
                if shift=='align':
                    subpix=[[]]*4
                    if follow=='dynamic':
#                        print luckies[k], gsx
                        dx=2*(luckies[k][4]-gsx+luckies[k][1])
                        dy=2*(luckies[k][5]-gsy+luckies[k][2])
                        x=luckies[k][4]-rad+luckies[k][1]
                        y=luckies[k][5]-rad+luckies[k][2]
                    else:
                        dx=2*luckies[k][1]
                        dy=2*luckies[k][2]
                        x=gsx-rad+luckies[k][1]
                        y=gsy-rad+luckies[k][2]

                    subpix[0]=[0,0,3*tmp[y,x]+tmp[y-1,x-1]+
                           tmp[y-1,x]+tmp[y,x-1]]   # top left
                    subpix[1]=[0,1,3*tmp[y,x]+tmp[y-1,x]+
                           tmp[y-1,x+1]+tmp[y,x+1]] # top right
                    subpix[2]=[1,0,3*tmp[y,x]+tmp[y,x-1]+
                           tmp[y+1,x-1]+tmp[y+1,x]] # bot left
                    subpix[3]=[1,1,3*tmp[y,x]+tmp[y,x+1]+
                           tmp[y+1,x]+tmp[y+1,x+1]] # bot right
                    offset=sorted(subpix, key=lambda list:list[2], reverse=True)
                    dx+=offset[0][1]
                    dy+=offset[0][0]
#                    print dy,2*(ysize-2)-8*rad+dy,dx,2*(xsize-2)-8*rad+dx
            
                crop=frame[dy:2*(ysize-2)-7*rad+dy,dx:2*(xsize-2)-7*rad+dx]
                tmpcounter=0
        
                if shape(crop)!=shape(stack):
                    print "skipped image %s(if too many images are skipped try running with static centroid) "%name
                    continue
                if select == 'weighted':
                    crop = crop*weights[k]
                    exposure +=  exp*weights[k]
                else:
                    exposure +=exp
            
                stack+=crop


#####     UPDATE HEADER     #####
        if 'EXPOSURE' in header: header['EXPOSURE']=exposure
        elif 'EXPTIME' in header: header['EXPTIME']=exposure
        header.update('IMAGTYPE', '%s' % imtype, 'Consecutive, Lucky, Weighted, or Lucky_Correlated')
        header.update('SELECTED', '%s' % pc[pcn], 'Percentage of images selected')
        header.update('IMGALIGN', '%s' % align, 'Image realignment around guide star')
        header.update('DETRIPLE', '%s' % dtr, 'Axis of detripling')
        if detrip!='none':
            header.update('DETRCORE', '%s' % coresz, 'Size of detripling core')
            header.update('DETRSEP', '%s' % minsep, 'Minimum core separation')
        else:
            header.update('DETRCORE', '0', 'Size of detripling core')
            header.update('DETRSEP', '0', 'Minimum core separation')
        header.update('ALIGNCX', '%d' %gsx, 'x-pixel initial centroid')
        header.update('ALIGNCY', '%d' %gsy, 'y-pixel initial centroid')
        header.update('ALIGNR', '%d' %rad, 'search region radius')
        header.update('CENTROID', '%s' %follow, 'Dynamic or static centroid')
    
    

    
#####     WRITE OUT    #####
    
        if os.path.isfile(outfile):
            os.remove(outfile)

        PF.writeto(outfile, stack, header)
    
        print 'Composite image written to %s' %  outfile
    	strg ='chmod 777 '+outpath+'/'+outname+'/*'
	os.system(strg)
    return 1
########################################################
########################################################
########################################################
########################################################




if __name__ == '__main__':
    if len(sys.argv) != 1 or sys.argv[1].startswith('-h') or sys.argv[1] == 'h':
        print """Usage. Requires: 
                **name of parameter file conatining :**
		
                Directory containing images
		dark file
		dark method
		Guide star x coordinate
		Guide star y coordinate
		region dimensions x
		percentage of images to be selected (0-100)
		lucky: \'lucky\',\'weighted\', \'coadded\', \'corr\' 
		shift: \'align\' or \'none\'
                detripling: \'v\' or \'h\' or \'none\'
                minimum separation or cores
                core size for detripling (odd integer)
                dynamic guide region; 1 or 0
            """
        sys.exit()


#####     DECLARE VARIABLES     #####
    pars = readconfig(sys.argv[1])
    inpath=pars['impath']+'/'+pars['imgdir']
    gsx=pars['x']
    gsy=pars['y']
    rad=pars['r']
    select=pars['sel']
    pc = pars['percent']
    ps = pars['ps']
    shift=pars['align']
    detrip=pars['detrip']
    minsep=float(pars['separation'])
    outpath=pars['outpath']
    coresz=pars['core']
    follow=pars['centroid']
    print inpath
    for sel in select:

        createstack(inpath,gsx,gsy,rad,sel,pc,shift,detrip,minsep,outpath,coresz,follow,ps)
