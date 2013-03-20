from os import path,system, getenv


def mymkdir(mypath):
    print 
    if path.isfile(mypath):
        print "\n\ncannot create a directory: the path points to a file\n\n"
        return -1
    elif path.isdir(mypath):
        print "\n\ndirectory " +mypath+ " exists\n\n"
        return 0
    else: strg = 'mkdir %s'%mypath 
    print "here",strg
    system(strg)
    return 0


def mygetenv(myenv):
    import sys
    env = getenv(myenv)   
    if env is not None:
        print env
    else :
        print "no environmental variable ",myenv," set. set your environmental variable to the appropriate path and then try again!"
        sys.exit()
    return(env)

def readconfig(configfile):
    f = open(configfile, 'r')
    config_string = f.read()
    parameters = eval(config_string)
    return parameters


def mjd(year, month, day):
    a=(14-month)//12
    y=year+4800-a
    m=month+(12*a)-3
    p=day+(((153*m)+2)//5)+(365*y)
    q=(y//4)-(y//100)+(y//400)-32045
    return p+q-2400000.5

def mkarray(value, nrows, ncols):
  return [[value]*ncols for _ in range(nrows)]
