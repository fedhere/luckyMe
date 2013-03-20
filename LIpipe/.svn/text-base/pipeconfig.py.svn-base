import sys
from enthought.traits.api import *
from enthought.traits.ui.api import *
from enthought.traits.ui.menu import OKButton, CancelButton
import os
from numpy import where
sys.path.append("../LIHSPcommon")
from myutils import readconfig

class Configfile(HasTraits):
    """ configuration file object """

    impath = String('/science/fbianco/bpldata',
        desc="the path to the directory containing your data spools",
        label="path to data spools ", )

    dodark = Enum('no dark','compute it','create and save it','use existing',
        desc="the dark method ",
        label="dark subtraction ", )

    darkpath = String('dark_2x2_10hz.fits',
        desc="the name of the spool used to generate your masterdark fits ",
        label="dark spool ", )

    spool0 = String('spool_22.fits',
        desc="the data spool fits name\n(all data spool inserted will be processed as 1 continuous serie of data) ",
        label="data spool ", )
    spool1 = String('',
        desc="the data spool fits name\n(all data spool inserted will be processed as 1 continuous serie of data) ",
        label="data spool ", )

    spool2 = String('',
        desc="the data spool fits name\n(all data spool inserted will be processed as 1 continuous serie of data) ",
        label="data spool ", )

    spool3 = String('',
        desc="the data spool fits name\n(all data spool inserted will be processed as 1 continuous serie of data) ",
        label="data spool ", )

    spool4 = String('',
        desc="the data spool fits name\n(all data spool inserted will be processed as 1 continuous serie of data) ",
        label="data spool ", )

    spool5 = String('',
        desc="the data spool fits name\n(all data spool inserted will be processed as 1 continuous serie of data) ",
        label="data spool ", )

    xc = CInt(0, label="x centroid (pixel) ", desc="the x pixel center of the brightest pixel search region")
    yc = CInt(0, label="y centroid (pixel) ", desc="the y pixel center of the brightest pixel search region")
    r = CInt(25, label="search radius (pixels) ", desc="the radius in pixel of the brightest pixel search region")

    sel = Enum('lucky','weighted','coadded',
        desc="the processing method: lucky imaging, weighted (images in the average are coadded with a weight proportional to their strehl ratio), or coadded (vanilla image stack equivalent to a longer exposure). Specifying the align option with the weighted and coadded method is equivalent to a software tip-tilt.",
        label="selection method ")

    pc1 = Range(low=0.0, high = 100, value = 1, desc="the percentage of images to use to create the combined image", label="selection percentage")
    pc2 = Range(low=0.0, high = 100, value = 10, desc="the percentage of images to use to create the combined image", label="selection percentage")
    pc3 = Range(low=0.0, high = 100, value = 30, desc="the percentage of images to use to create the combined image", label="selection percentage")
    

    align = Enum('yes','no',
        desc="the option to align images while stacking",
        label="aligning images ", )

    detripling = Enum('none','horizontal','vertical',
        desc="the option to detriple the stack, for close binary targets",
        label="detripling ", )

    centroid = Enum('dynamic','static',
        desc="the centroid method: static preserves original coordinated for the brightest pixel search region, dynamic updates them",
        label="centroid method ", )

    ra = String('',
        label="ra ", desc="the ra of the center of the field (optional)")
    dec  = String('',
        label="dec ", desc="the dec of the center of the field (optional)")

    scope  = Enum('FTN','FTS','1m','0.4m',
        label="telescope ", desc="the telescope name or class (optional)")
    camera  = Enum('LucaR','iXon888',
        label="camera ", desc="the camera model (optional)")
    psf = CFloat(1.0, label="expected psf at site ('')",desc="the expected 'instantaneous' psf")
    ps = CFloat(0.26, label="plate scale (''/pixel) ", desc="the pixel scale (''/pixel)")

    nskip = CInt(0, label="number of files to skip ", desc='the number of files to skip at the beginning of the first spool for the centroid analysis')

    unspool = Enum('y','n',
                   label ="perform unspooling ",)
    red = Enum('y','n',
                   label ="perform image reduction ",)
    centan = Enum('y','n',
                   label ="perform centroid analysis ",)

    view = View(Item('impath'),
                Item('spool0'),
                Item('spool1'),
                Item('spool2'),
                Item('spool3'),
                Item('spool4'),
                Item('spool5'),
                Item('dodark', style='custom'),
                Item('darkpath'),#, orientation='horizontal'),
                Item('xc'),Item('yc'),Item('r'),#, orientation='horizontal'),
                Item('sel', style='custom'),Item('pc1'),
                Item('pc2'),
                Item('pc3'),
                Item('detripling', style='custom'),
                Item('align', style='custom'),
                Item('centroid', style='custom'),#,orientation='horizontal'),
                Item('ra'),Item('dec'),Item('scope', style='custom'),Item('camera', style='custom'),#,orientation='horizontal'),
                Item('ps'),Item('psf'),Item('nskip'),#,orientation='horizontal'),
                Item('unspool', style='custom'),Item('red', style='custom'),Item('centan', style='custom'),
                resizable = True,
                title ="wanna get lucky?",
                scrollable=True,
                buttons = [OKButton, CancelButton])

#                buttons = [ 'OK', 'Cancel' ] )
    

    def creatit(self):
        """ creates a configuration file and runs pipeline with it """
        spool = []
        for i in [self.spool0,self.spool1,self.spool2,self.spool3,self.spool4,self.spool5]:
            if len(i)>0:
                spool.append(i)
        pc = []
        for i in [self.pc1,self.pc2,self.pc3]:
            if i==0:
                continue
            pc.append(i)

        
#        spoolstring = ''
#        for i in spool:
#            if len(i)>0:
#                if spoolstring != '':
#                    spoolstring = spoolstring+','
#                spoolstring = spoolstring + "'"+i + "'"
#        spoolstring +="]

        if self.dodark=='no dark':
            dodark = 0
        elif self.dodark=='compute it':
            dodark = 1
        elif self.dodark=='create and save it':
            dodark = 2
        elif self.dodark=='use existing':
            dodark = 3            

        f=open('config.txt','w')
        print >>f, """{'imdir' : './',\\
'separation' : 8,\\
'core' : 3,\\
'profile' : 'm',\\
'alpha' : 1.4,\\
'beta' : 3.0,\\
'nsteps' : 100,\\
'disp' : 'y',\\
'notes':'',\\"""
        print >>f,"'impath' : '"+self.impath+"',\\"
        print >>f,"'dodark' : %i,\\"%dodark
        print >>f,"'dark' : '"+self.darkpath+"',\\"
        print >>f,"'spool' : ",spool, ",\\"
        print >>f,"'x' : %i,\\"%self.xc
        print >>f,"'y' : %i,\\"%self.yc
        print >>f,"'r' : %i,\\"%self.r
        print >>f,"'sel' : '"+self.sel+"',\\"
        print >>f,"'percent' : ",pc,",\\"
        print >>f,"'psf' : %3.1f,\\"%self.psf
        if self.ps == 0:
            if self.camera == 'LucaR':
                self.ps = 0.08
            elif self.camera == 'iXon888':
                self.ps = 0.13
            else:
                print "please select a camera or a pixel scale"
                sys.exit()

        print >>f,"'ps' : %3.1f,\\"%self.ps
        print >>f,"'centroid' : '"+self.centroid+"',\\"
        print >>f,"'ra' : '"+self.ra+"',\\"
        print >>f,"'dec' : '"+self.dec+"',\\"
        print >>f,"'scope' : '"+self.scope+"',\\"
        print >>f,"'camera' : '"+self.camera+"',\\"
        print >>f,"'nskip' : %i,\\"%self.nskip
        print >>f,"'unspool' : '"+self.unspool+"',\\"
        print >>f,"'reduce' : '"+self.red+"',\\"
        print >>f,"'centan' : '"+self.centan+"',\\"
        if self.align == 'yes':
          print >>f,"'align' : 'align',\\"  
        else:
          print >>f,"'align' : 'none',\\"  
        if self.detripling == 'none':
            print >>f,"'detrip' : 'none'}"  
        elif self.detripling.startswith('h'):
            print >>f,"'detrip' : 'h'}"
        else:
            print >>f,"'detrip' : 'v'}"

        f.close()
        strg = "head config.txt"
        os.system(strg)

        strg = "python callall.py config.txt "
        os.system(strg)


if __name__ == "__main__":
    if os.path.isfile('config.txt'):
       par = readconfig('config.txt')
    
       configfile = Configfile()
       spool = par['spool']
       configfile.impath =par['impath']
       if len(spool[0])==1:
           configfile.spool0 =par['spool']
       else:
           configfile.spool0 =par['spool'][0]
           if len(spool)>1:
               configfile.spool1 =par['spool'][1]
               if len(spool)>2:
                   configfile.spool2 =par['spool'][2]
                   if len(spool)>3:
                       configfile.spool3 =par['spool'][3]
                       if len(spool)>4:
                           configfile.spool4 =par['spool'][4]
                           if len(spool)>5:
                               configfile.spool5 =par['spool'][5]

           
               
       if par['dodark']==0:
           configfile.dodark = 'no dark'
       elif par['dodark']==1:
           configfile.dodark ='compute it'
       elif par['dodark']==2:
           configfile.dodark ='create and save it'
       elif par['dodark']==3:
           configfile.dodark = 'use existing'
       
#       configfile.dark =par['dark']
#       configfile.spool0 =par['spool']
       configfile.darkpath=par['dark']
       configfile.xc =int(par['x'])
       configfile.yc =int(par['y'])
       configfile.r =int(par['r'])
       configfile.sel =par['sel']  
       print par['percent'], len(par['percent'])
       configfile.pc1 =float(par['percent'][0])
       configfile.pc2 =0.0
       configfile.pc3 =0.0
       if len(par['percent'])==2:
           configfile.pc2 =float(par['percent'][1])
       if len(par['percent'])==3:
           configfile.pc2 =float(par['percent'][1])
           configfile.pc3 =float(par['percent'][2])
       configfile.psf =par['psf']
       configfile.ps =par['ps']
       configfile.centroid =par['centroid']
       if par['align'].startswith('a'):
           configfile.align = 'yes'
       else :
           configfile.align = 'no'
       configfile.ra =par['ra']
       configfile.dec =par['dec']
       configfile.scope =par['scope']
       configfile.camera =par['camera']
       configfile.nskip =int(par['nskip'])
       configfile.unspool =par['unspool']
       configfile.reduce =par['reduce']
       configfile.centan =par['centan']
    
    configfile.configure_traits(view = 'view')
    configfile.creatit()
