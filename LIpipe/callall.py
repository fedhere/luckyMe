import sys,os,time
from numpy import ones

sys.path.append("../LIHSPcommon")

from myutils import mygetenv,readconfig





#unspool
def unspool(par,outpath):
    from myutils import mymkdir,mygetenv
    from unspool import unspoolit

    print 'creating ',outpath
    
    if mymkdir(outpath)!= 0:
        sys.exit()
    if mymkdir(outpath+'/unspooled')!= 0:
	sys.exit()
#    strg = 'mkdir %s'%outpath
#    os.system(strg)
    heredir=mygetenv('SPEEDYOUT')+'/darks/'
    if mymkdir(heredir)!= 0:
        sys.exit()
#    strg = 'mkdir mygetenv('SPEEDYOUT')//darks'
#    os.system(strg)
    dodarks = [par['dodark']]
    flatten = [par['flatten']]

    for i in range(1,len(par['spool'])):
        if dodarks[0] == 1 or dodarks[0] == 0:
		dodarks.append(dodarks[0])
	else:
		dodarks.append(3)

                
    inpath=par['impath'].strip()+'/'+par['imdir'].strip()+'/'
    print inpath
    print par['spool'], dodarks
    for  i,img in enumerate(par['spool']):
        fname = '%s/%s/%s'%(par['impath'].strip(),par['imdir'],  img)
        if os.path.isfile(fname):
	    ret = unspoolit(inpath, img, inpath,par['dark'],'avg', dodarks[i], heredir,outpath,0,flatten,par['flats'])
	    #ret = unspoolit(inpath, img, inpath,par['dark'],'clip', dodarks[i], heredir,outpath,0)
 	    if ret !=1:
		print "\n\n\n!!!!!!!!!!!!!!!PANIC: unspooling failed. exiting!!!!!!!!!!!\n\n\n"
		sys.exit(0)		
        else: 
            print 'no spool %s to be found'%fname
   	    sys.exit(0)
###################             
#############creating strehl ratio###########

def runstack(pars, outpath):
#    from weightedstack_corr import createstack
    outpath=outpath.replace('\n','').replace('\\n','').strip()
    from makemelucky import createstack
    inpath='%s/unspooled/'%(outpath)
    gsx=pars['x']
    gsy=pars['y']
    rad=pars['r']
    ps =pars['ps']
    select=pars['sel']
    pc = pars['percent']
    shift=pars['align']
    detrip=pars['detrip']
    minsep=float(pars['separation'])
    coresz=pars['core']
    follow=pars['centroid']
    saturation=float(pars['saturation'])
    fast = pars['fftw3'].startswith('y')


    if len(select[0]) == 1:
        print "\n\none selection method:"
        print "\nprocessing %s\n\n\n"%select
#	ret = createstack(inpath,gsx,gsy,rad,select,pc,shift,detrip,minsep,outpath,coresz,follow,ps)
        ret = createstack(inpath,gsx,gsy,rad,select,pc,shift,detrip,minsep,outpath,coresz,follow,ps,saturation, fast)
    else:
        for sel in select: 
            print "\n\n\nprocessing %s\n\n\n"%sel
#            ret = createstack(inpath,gsx,gsy,rad,sel,pc,shift,detrip,minsep,outpath,coresz,follow,ps)
            ret = createstack(inpath,gsx,gsy,rad,sel,pc,shift,detrip,minsep,outpath,coresz,follow,ps, saturation, fast)
            if ret !=1:
		print "\n\n\n!!!!!!!!!!!!!!!PANIC: LI reducing   failed. exiting!!!!!!!!!!!\n\n\n"
		sys.exit(0)		

def runcentan(pars, outpath, nameroot):
    from psf import centan
    newdir = '%s/%s_displacement' %(outpath,nameroot)
    from mymkdir import mymkdir
    if mymkdir(newdir)!= 0:
        sys.exit()
#    strg = 'mkdir %s'%newdir
#    os.system(strg)
    dispfile = "%s/strehl_list_%s.dat"%(outpath,pars['centroid'])
    profile=par['profile'] #'g' for gaussian, 'm' for moffat
    PSFg=float(par['psf'])
    PSFm=float(par['alpha'])
    displacement=par['disp']
    beta=float(par['beta'])
    nx,ny=100,100
    c=(50,50)
    if profile=='g':
        PSF = PSFg
    elif profile == 'm':
        PSF=PSFm
        
    else:
        print "unknown profile"
        sys.exit()
        
    ret =centan(outpath,dispfile, par, 1,nameroot, newdir)
    if ret !=1:
		print "\n\n\n!!!!!!!!!!!!!!!PANIC: centroid analysis failed. exiting!!!!!!!!!!!\n\n\n"
		sys.exit(0)	
######################################################################################

print sys.argv[1]    
par = readconfig(sys.argv[1])
print par, par['impath']
print "\n\n\n ABOUT TO RUN THE LCOGT LI PIPELINE TO\n"
if par['unspool'].startswith('y'):
    print "***UNSPOOL\n"
if par['reduce'].startswith('y'):
    print "***REDUCE: finding strehl ratio and generating composite imaged\n"    
if par['centan'].startswith('y'):
    print "***ANALYZE: the spool (x-y motion, strehl distribition)\n\n\n"





nameroot=par['spool'][0]
if nameroot.endswith('.fits'):
    nameroot = nameroot[:-5]

if len(par['spool'])>1:	
	outpath = '%s/%s_all'%(mygetenv('SPEEDYOUT'),nameroot)
else:
        outpath = '%s/%s'%(mygetenv('SPEEDYOUT'),nameroot)

inpath = par['impath']+'/'+par['imdir']

if par['unspool'].startswith('y'):
    print "UNSPOOLING"
    print par['impath']
    unspool(par, outpath)
    
cosmic = 0
#if par['cosmic'].startswith('y'):  
#    cosmic = 1


if par['reduce'].startswith('y'):
	runstack(par, outpath)

if par['centan'].startswith('y'):
	runcentan(par, outpath, nameroot)


print "all done"
