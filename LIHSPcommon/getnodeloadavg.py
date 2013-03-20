##
## GETNODELOADAVG
##
## Function to get the load average of a beowulf cluster node
##
## Version    Date    	  Author  Comments
## 0.1	      2009-07-24  RAS 	  Initial version
##


## Import required modules:
from commands import getstatusoutput

## Function: getnodeloadavg
def getnodeloadavg(node):

    ## Read the system's /proc/loadavg file
    tstr = "%s %s %s" % ( 'bpsh', node, 'cat /proc/loadavg')
    (iexec,coutput) = getstatusoutput(tstr)

    ## Parse the resulting string of numbers, with catch for those nodes which
    ## return 'Node N is down'.  These get their loadavg set to 100% to prevent them 
    ## being used.  
    if coutput.__contains__('is down') == True:
      	loadavg = 100.0
    else:
      	loadavg = float(coutput.split()[0])

    ## Return the load average:
    return loadavg
