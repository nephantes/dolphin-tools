#!/bin/env python
import os,sys
from optparse import OptionParser
from ConfigParser import SafeConfigParser

def getConfig(filename='../default_params/config.ini'):
   parser = SafeConfigParser()
   filename = sys.path[0]+"/../default_params/config.ini"

   parser.read(filename)

   params_section="Biocore"		
   if (os.environ.has_key('DOLPHIN_PARAMS_SECTION')):		
      params_section=os.environ['DOLPHIN_PARAMS_SECTION']

   config={}
   for name, value in parser.items(params_section):
      config[name]=value
   return config


