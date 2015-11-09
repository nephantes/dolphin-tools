#!/usr/bin/python
import logging
from optparse import OptionParser
import urllib,urllib2
import json
import os
import re
import cgi
import warnings
import sys
import time
sys.path.insert(0, sys.path[0])
from config import *
from funcs import *

class jobStatus:
    url=""
    f=""
    def __init__(self, url, f ):
        self.url = url
        self.f = f
    
    def checkAllJobsFinished(self, username, wkey, servicename, logging):
        data = urllib.urlencode({'func':'checkAllJobsFinished', 'username':username, 
                                 'wkey':wkey, 'servicename':servicename})
        self.f.queryAPI(self.url, data, servicename, logging) 
    
    def insertJob(self,  username, wkey, com, jobname, servicename, jobnum, result, logging):
        data = urllib.urlencode({'func':'insertJob', 'username':username, 
                                 'wkey':wkey, 'servicename':servicename, 
                                 'com':com , 'jobname':jobname, 'jobnum':str(jobnum), 
                                 'result':str(result)})
        self.f.queryAPI(self.url, data, "INSERT"+jobname, logging) 
   
    def updateJob(self,  username, wkey, jobname, servicename, field, jobnum, result, logging):
        data = urllib.urlencode({'func':'updateJob', 'username':username, 
                                 'wkey':wkey, 'servicename':servicename, 
                                 'field':field , 'jobname':jobname, 'jobnum':str(jobnum), 
                                 'result':str(result)})
        self.f.queryAPI(self.url, data, "UPDATE"+jobname, logging) 
    
    def insertJobOut(self, username, wkey, jobnum, outdir, edir, logging):
        file=str(outdir)+"/tmp/lsf/"+str(jobnum)+".std"
        print "insertJobOut file:"+file
        if os.path.isfile(file) and os.access(file, os.R_OK):
            command="python " + edir  + "/readLSFout.py -f "+file
            child = os.popen(command)
            jobout = child.read().rstrip()
            err = child.close()
         
            if err:
                 logging.info('ERROR: %s failed w/ exit code %d' % (command, err))
                 raise RuntimeError, 'ERROR: %s failed w/ exit code %d' % (command, err) 
            jobout=self.f.plaintext2html(re.sub('\'', "\"", jobout))
           
            data = urllib.urlencode({'func':'insertJobOut', 'username':username, 
                                 'wkey':wkey, 'jobnum':str(jobnum), 'jobout':jobout})
            self.f.queryAPI(self.url, data, "jobnum:"+str(jobnum), logging) 

       
def main():
   try:
        #python finishJob.py -u kucukura -k nVy1THnthvrRWfXcj187KeDDQrNAkY -s splitFastQ
        parser = OptionParser()
        parser.add_option('-u', '--username', help='defined user in the cluster', dest='username')
        parser.add_option('-k', '--key', help='defined key for the workflow', dest='wkey')
        parser.add_option('-c', '--com', help='bash script of the command', dest='com')
        parser.add_option('-j', '--jobname', help='name of of the job', dest='jobname')
        parser.add_option('-s', '--servicename', help='service name', dest='servicename')
        parser.add_option('-t', '--type', help='type of the operation', dest='type')
        parser.add_option('-n', '--jobnum', help='submitted job number', dest='jobnum')
        parser.add_option('-m', '--message', help='resulting message of the job', dest='message')
        parser.add_option('-o', '--outdir', help='output directory', dest='outdir')
        parser.add_option('-f', '--config', help='config parameters section', dest='config')

        (options, args) = parser.parse_args()
   except:
        #parser.print_help()
        print "for help use --help"
        sys.exit(2)


   USERNAME                = options.username   
   WKEY                    = options.wkey 
   COM                     = options.com
   JOBNAME                 = options.jobname
   SERVICENAME             = options.servicename
   TYPE                    = options.type
   JOBNUM                  = options.jobnum
   MESSAGE                 = options.message
   OUTDIR                  = options.outdir
   CONFIG                  = options.config
        

   f = funcs()
   config = getConfig(CONFIG)
   jobStat = jobStatus(config['url'], f)
   
   edir=os.path.dirname(sys.argv[0])
   print OUTDIR
   print TYPE
   logfile="%s/tmp/lsf/%s.jobStatus.log"%(OUTDIR,str(JOBNUM))
   logging.basicConfig(filename=logfile, filemode='a',format='%(asctime)s %(message)s', datefmt='%m/%d/%Y %I:%M:%S %p', level=logging.DEBUG)
   logging.info("File Path:%s"%os.getcwd())
   logging.info("WKEY:"+str(WKEY))
   logging.info("TYPE:"+str(TYPE))

   if (TYPE == "dbSubmitJob"):
        jobStat.insertJob(USERNAME, WKEY, COM, JOBNAME, SERVICENAME, JOBNUM, MESSAGE, logging)
   elif (TYPE == "dbSetStartTime"):
        jobStat.updateJob(USERNAME, WKEY, JOBNAME, SERVICENAME, "start_time", JOBNUM, MESSAGE, logging)  	
   elif (TYPE == "dbSetEndTime"):
        jobStat.updateJob(USERNAME, WKEY, JOBNAME, SERVICENAME, "end_time", JOBNUM, MESSAGE, logging) 	
        jobStat.insertJobOut(USERNAME, WKEY, JOBNUM, OUTDIR, edir, logging)

   jobStat.checkAllJobsFinished(USERNAME, WKEY, SERVICENAME, logging)

if __name__ == "__main__":
    main()


