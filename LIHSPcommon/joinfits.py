import sys,os
import pyfits as PF
from myutils import mygetenv, mymkdir

inpath=mygetenv('SPEEDYIN')

def modheader(fl):
    image=PF.open(inpath+'/'+fl)
    header=image[0].header
    try:
        dataimage=PF.open(inpath+'/tmp-'+fl)
    except:
        print "no tmp file"
        return()
    dataheader=dataimage[0].header

    for i in header.ascard:
        if  'NAXIS' in i.key or i.key in ['BITPIX','END']:
            continue
        elif 'COMMENT' in i.key:
            dataheader.add_comment(i.comment)
#        else:
        else: 
            try:
                dataheader.update(i.key, i.value, i.comment)
            except:
                try:
                    dataheader.update('hierarch '+i.key, i.value, i.comment)
                except:
                    continue
                

    os.system('mv '+inpath+'/'+fl+' '+inpath+'/'+fl+'.save')
    #os.remove(inpath+'/'+fl)
    dataimage.writeto(inpath+'/'+fl)
    return()    
    
if __name__ == '__main__':
    try :
        sys.argv[1]
    except:
        print "needs names of header files"
        sys.exit()

    for i in sys.argv[1:]:
        modheader(i)
    
