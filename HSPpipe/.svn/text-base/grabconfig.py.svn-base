#!/usr/bin/env python

def grabconfig(configfile):
    import ConfigParser
    #Read config file and set variables
    config = ConfigParser.SafeConfigParser()
    config.read(configfile)
    out = {}
    out['imdir'] = config.get('environment','imdir')
    out['archive'] = config.get('environment','archive')
    out['rundir'] = config.get('environment','rundir')
    out['apertures'] = config.get('parameters','apertures')
    out['ensemble'] = config.get('parameters','ensemble')
    out['stacklist'] = config.get('parameters','stacklist')
    out['coords'] = 'stack.coo'
    out['psfd'] = config.get('parameters','psfd')
    out['target'] = int(config.get('parameters','target'))
    return out

if __name__ == '__main__':
    import sys
    out = grabconfig(sys.argv[1])
    print out
