import sys,os,time
import pyfits as PF
sys.path.append("../LIHSPcommon")
from myutils import mygetenv,readconfig, mymkdir, mjd

speedyout=mygetenv('SPEEDYOUT')
def readconfig(configfile):
    f = open(configfile, 'r')
    config_string = f.read()
    parameters = eval(config_string)
    return parameters



#unspool
def unspool(par,outpath):
    from myutils import mymkdir,mygetenv
    from unspool import unspoolit
    if mymkdir(outpath)!= 0:
        sys.exit()
    if mymkdir(outpath+'/unspooled')!= 0:
	sys.exit()
#    strg = 'mkdir %s'%outpath
#    os.system(strg)
    heredir=mygetenv('SPEEDYOUT')+'/darks/'
    if mymkdir(heredir)!= 0:
        sys.exit()


    dodarks = [par['dodark']]



    if (isinstance(par['spool'],types.StringTypes)):
        print "only 1 file to unspool"
    else:    
        
        for i in range(1,len(par['spool'])):
            if dodarks[0] == 3 or dodarks[0] == 0:
		dodarks.append(dodarks[0])
            else:
		dodarks.append(2)

                
    inpath=par['impath']+'/'+par['imdir']+'/'
    print inpath
    print par['spool'], dodarks

    heredir=speedyout+'/darks/'
    if mymkdir(heredir)!= 0:
        sys.exit()

    if (isinstance(par['spool'],types.StringTypes)):
    	nameroot=par['spool']
	if nameroot.endswith('.fits'):
       		nameroot = nameroot[:-5]
        fname = '%s/%s/%s.fits'%(par['impath'],par['imdir'], nameroot)
        if os.path.isfile(fname):
	    ret = unspoolit(inpath, nameroot+'.fits', inpath,par['dark'],'avg', dodarks[0], heredir,outpath,0)
 	    if ret !=1:
		print "\n\n\n!!!!!!!!!!!!!!!PANIC: unspooling failed. exiting!!!!!!!!!!!\n\n\n"
		sys.exit(0)		
        else: 
            print 'no spool %s to be found'%fname
   	    sys.exit(0) 
    else :
        for  i,img in enumerate(par['spool']):
	    if img.endswith('.fits'):
       		img = img[:-5]
            fname = '%s/%s/%s.fits'%(par['impath'],par['imdir'],  img)
            if os.path.isfile(fname):
                ret = unspoolit(inpath, img+'.fits', inpath,par['dark'],'avg', dodarks[i], heredir,outpath)
                if ret !=1:
                    print "\n\n\n!!!!!!!!!!!!!!!PANIC: unspooling failed. exiting!!!!!!!!!!!\n\n\n"
                    sys.exit(0)		
            else: 
                print 'no spool %s to be found'%fname
                sys.exit(0) 
################### aperture photometry
def myapphot(par, cosmic):
    from myutils import mymkdir,mygetenv
    from myapphot import *

    SPEEDYOUT=mygetenv('SPEEDYOUT')
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
             
        

    
    
    ret =myphot(inpath, par['coords'], int(par['last']), par['ap'], par['centroid'], nameroot, par['target'])
    if ret !=1:
        print "\n\n\n!!!!!!!!!!!!!!!PANIC: photometry failed. exiting!!!!!!!!!!!\n\n\n"
        sys.exit(0)	
######################################################################################

################### cosmics and sextractor
def runsex(par,cosmic):
    if (cosmic == 1):        
        strg = 'python runcosmic.py /science/fbianco/HSPdata/unspooled/%s %s_ 0 %d'%(par['spool'],par['spool'],par['last'])
        os.system(strg)
        
        strg='for  i in /science/fbianco/HSPdata/%s/clean/unspooled/*fits ; do python runsextractor.py $i; done'%par['spool']
        os.system(strg)
        
    if (cosmic == 0):
        strg='for  i in /science/fbianco/HSPdata/%s/unspooled/*fits ; do python runsextractor.py $i; done'%par['spool']
        print strg
	os.system(strg)

#photometry
def photometry(par, cosmic):
    if (cosmic):
        strg = 'python phot.py /science/fbianco/HSPdata/%s/unspooled/clean/ %s %d %s'%(par['spool'],par['coordat'], par['last'], par['spool'])
        os.system(strg)

    else:
        strg = 'python phot.py /science/fbianco/HSPdata/%s/unspooled/ %s %d %s'%(par['spool'],par['coordat'], par['last'], par['spool'])
        
        os.system(strg)

                                                
######################################################################################
from myutils import mymkdir,mygetenv
import types
par = readconfig(sys.argv[1])
print par
SPEEDYOUT=mygetenv('SPEEDYOUT')

if (isinstance(par['spool'],types.StringTypes)):
    nameroot=par['spool']
    if nameroot.endswith('.fits'):
        nameroot = nameroot[:-5]

    outpath = '%s/%s/'%(SPEEDYOUT,nameroot)
else :
    nameroot=par['spool'][0]
    if nameroot.endswith('.fits'):
        nameroot = nameroot[:-5]

    outpath = '%s/%s_all/'%(SPEEDYOUT,nameroot)


tmp = '%s/%s/%s.fits'%(par['impath'],par['imdir'],  nameroot)
image=PF.open(tmp)
header=image[0].header


print "last image in spool: ",par['last']

if par['last'] == 0 or par['last'] >image[0].header['NAXIS3']:
        par['last']=int(image[0].header['NAXIS3'])
        print "last image in spool: ",par['last']


if par['unspool'].startswith('y'):
    print "\n\n\nUNSPOOLING\n\n\n"
    unspool(par, outpath)

cosmic = 0
if par['cosmic'].startswith('y'):  
    print "\n\n\nremoving cosmics...\n\n\n"
    cosmic = 1

if par['sextract'].startswith('y'):
    print "\n\n\nextracting (sex)...\n\n\n"
    runsex(par, cosmic)

if par['phot'].startswith('y'):
    print "\n\n\nrunning iraf photometry...\n\n\n"    
    photometry(par, cosmic)

if par['createlc'].startswith('y'):
    print "\n\n\ncreating lcvs...\n\n\n"
    if(cosmic):
        strg = 'python extractlc.py /science/fbianco/HSPdata/%s/unspooled/clean/ %s %s %d  %d'%(nameroot,par['coordat'],nameroot,par['last'],par['ap'])
        os.system(strg)
    else:	
        strg = 'python extractlc.py /science/fbianco/HSPdata/%s/unspooled %s %s %d  %d'%(nameroot,par['coordat'],nameroot,par['last'],par['ap'])
	os.system(strg)



if par['myapphot'].startswith('y'):
    print "\n\n\nrunning my aperture photometry...\n\n\n"    
    myapphot(par, cosmic)


