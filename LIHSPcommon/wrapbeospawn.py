##
## WRAPBEOSPAWN
##
## Function to wrap a spawning of an external process depending on the
## machine available.  If the machine is stand-alone, then the call should
## execute normally.  If the machine is a cluster, then the function should
## find a node with a low load average and send the job to that node using 
## bpsh.  
##
## Version    Date    	  Author  Comments
## 0.1	      2009-07-24  RAS 	  Initial version
## 0.2	      2010-03-17  RAS 	  Added flag for debug output
## 0.3	      2011-03-10  RAS 	  Updated to call the subprocess module if the code 
##    	      	      	      	  is not running on a Beowulf
## 0.4	      2011-06-22  RAS 	  Updated to use the subprocess Popen object for
##    	      	      	      	  bpsh calls also in the process of enabling this
##    	      	      	      	  function to handle redirected input and output
##    	      	      	      	  for cluster calls. 


############################
# IMPORT OF REQUIRED MODULES
from commands import getstatusoutput
from getnodeloadavg import *
from os import spawnv, P_NOWAIT
from sys import exit
from subprocess import Popen

############################
# FUNCTION: WRAPCALL

def wrapbeospawn(comstr,AllowedNodesList=['9999']):

    '''Function to wrap a spawning of an external process depending on
    the machine available.  If the machine is stand-alone, then the call should
    execute normally.  If the machine is a Beowulf cluster, then the function should
    find a node with a low load average and send the job to that node using  bpsh.
    Inputs: 
    comstr: Single concatenated string containing shell commandline and commandline 
    arguments
    AllowedNodesList: List of node numbers as strings which the function is allowed
    to use.  If no list is given, or if list = [ '9999' ], then all available nodes will be used.
    Outputs: Tuple of the node number used for the process and the 
    process ID number (PID).
    '''  
    
    # Debug flag, set to 1 for debug output:
    dbg = 0

    # First step is to find out whether the current machine is part of a
    # Beowulf cluster or not.  
    UseNode = 0
    try:
      	(istat,NodesOutput) = getstatusoutput('beostat | grep Node')
	Nnodes = len(NodesOutput.split('\n')) - 1
      	if istat != 0:
     	    Nnodes = -1
     	    UseNode = -1
    except:
      	Nnodes = -1
      	UseNode = -1
    if dbg==1:
        print 'WRAPBEOSPAWN Nnodes='+repr(Nnodes)+' UseNode='+repr(UseNode)+'\n'
      	print AllowedNodesList
	
    # If this is on a Beowulf cluster, and AllowedNodesList = [9999], then
    # we're allowed to run on all nodes, so find out which ones those are:
    if Nnodes > 0 and AllowedNodesList[0] == '9999':
      	AllowedNodesList = []
	for node in range(0,Nnodes,1):
	    AllowedNodesList.append(repr(node))
      	if dbg==1:  
	    print 'WRAPBEOSPAWN Allowed nodes are\n'
      	    print AllowedNodesList
      	    print '\n'

    # If this is a beowulf cluster (Nnodes > 0), then find a node that's not
    # busy from the list of usable nodes.  If all nodes are occupied, this
    # finds the one with the lowest load average.  
    if Nnodes > 0:
      	UseNode = -1
      	minload = 10000.0
      	for node in AllowedNodesList:
      	    load = (getnodeloadavg(node))
      	    if load < minload:
	      	UseNode = node
	      	minload = load
      	if dbg ==1:
	    print 'WRAPBEOSPAWN Using node '+repr(UseNode)+' with current load '+repr(minload)+'\n'
    
    # Now fire the command given at the shell, using bpsh if we're running on 
    # a Beowulf cluster. 
    if Nnodes > 0 and UseNode > -1:
	scriptname = '/usr/bin/bpsh'
	inputfile = 'NONE'
	outputfile = 'NONE'
      	if dbg==1:
	    print 'WRAPBEOSPAWN Submitting bpsh call\n'
	    print 'Original command: ',comstr
	
	# If the commandstring indicates an input commandfile
	# should be used, parse the original command string:
	if comstr.__contains__('<') == True:
	    codename = comstr.split()[0]
	    idx = comstr.split().index('<')
	    inputfile = comstr.split()[idx+1]
	    if path.isfile(inputfile) == True:
	      	infileobj = open(inputfile,'r')
	      	args = ['/usr/bin/bpsh', UseNode, codename]
      	    else:
     	      	Nnodes = -1
     	      	UseNode = -1
	
	# If the commandstring indicates an output commandfile
	# should be used, parse the original command string:
	if comstr.__contains__('>') == True:
	    codename = comstr.split()[0]
	    idx = comstr.split().index('>')
	    outputfile = comstr.split()[idx+1]
	    (tfile,tpath) = path.split(outputfile)
	    if path.isdir(tpath) == True:
	      	outfileobj = open(outputfile,'w')
	      	args = ['/usr/bin/bpsh', UseNode, codename]
	    else:
	      	Nnodes = -1
     	      	UseNode = -1
	
	# Otherwise, just split the commandstring into a list:
	if comstr.__contains__('<') == False and \
	  comstr.__contains__('>') == False:
	    comlist = comstr.split()
	    args = ['/usr/bin/bpsh', UseNode] + comlist
	    
      	if dbg==1:
	    print 'Scriptname: ',scriptname
	    print 'Args: ',args
	    print 'Inputfile: ',inputfile
	    print 'Outputfile: ',outputfile
      	
	# Contine to issue command only if we've been able to 
	# parse the commandline properly:
	if Nnodes > 0 and UseNode > 0:
	    if inputfile == 'NONE' and outputfile=='NONE':
            	PID = Popen(args).pid
	    elif inputfile != 'NONE' and outputfile=='NONE':
	    	PID = Popen(args,stdin=infileobj).pid
	    elif inputfile == 'NONE' and outputfile!='NONE':
	    	PID = Popen(args,stdout=outfileobj).pid
	    elif inputfile != 'NONE' and outputfile!='NONE':
	    	PID = Popen(args,stdin=infileobj,stdout=outfileobj).pid
      	    if dbg==1:
	      print 'PID: ',PID
	
	    if inputfile != 'NONE':
	    	infileobj.close()
	    if outputfile != 'NONE':
	    	outfileobj.close()
	
    # Otherwise fire the command at the shell normally. 
    else:
      	if dbg==1:
	    print 'WRAPBEOSPAWN Submitting normal call\n'
	comlist = comstr.split()
        PID = Popen(comlist).pid

    # Return the execution status code and output from the command
    return UseNode, PID
    
    
    
############################
# COMMANDLINE OPERATION

if __name__ == '__main__':

    # Import required functions:
    from sys import argv
    from sys import exit
    
    # Customisable list of allowed nodes:
    List = [ '9999' ]
    
    # Read commandline arguments and concatenate as a comstr
    comstr = ' '.join(argv[1:])
    
    # Read list of allowed nodes if given:
    if len(argv) == 10:
      	List = []
	List.append(argv[3])

    print comstr
    # Call function
    (inode,pid) = wrapbeospawn(comstr,List)
    
    # Print output
    #print 'Execution code: '+repr(iexec)+'\n'
    #print 'Output: '+coutput+'\n'
    print inode,pid
