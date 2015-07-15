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

class jobStatus:
    config=""
    url=""
    re_string = re.compile(r'(?P<htmlchars>[<&>])|(?P<space>^[ \t]+)|(?P<lineend>\r\n|\r|\n)|(?P<protocal>(^|\s)((http|ftp)://.*?))(\s|$)', re.S|re.M|re.I)
    
    def plaintext2html(self, text, tabstop=4):
        def do_sub(m):
            c = m.groupdict()
            if c['htmlchars']:
                return cgi.escape(c['htmlchars'])
            if c['lineend']:
                return '<br>'
            elif c['space']:
                t = m.group().replace('\t', '&nbsp;'*tabstop)
                t = t.replace(' ', '&nbsp;')
                return t
            elif c['space'] == '\t':
                return ' '*tabstop;
            else:
                url = m.group('protocal')
                if url.startswith(' '):
                    prefix = ' '
                    url = url[1:]
                else:
                    prefix = ''
                last = m.groups()[-1]
                if last in ['\n', '\r', '\r\n']:
                    last = '<br>'
                return '%s<a href="%s">%s</a>%s' % (prefix, url, url, last)
        return re.sub(self.re_string, do_sub, text)
    
    def queryAPI(self, data, name, logging):
        opener = urllib2.build_opener(urllib2.HTTPHandler())
        trials=0
        while trials<5:
           try:
              mesg = opener.open(self.url, data=data).read()
              trials=10
           except:
              print "Couldn't connect to dolphin server (%s)"%trials
              logging.info("Couldn't connect to dolphin server (%s)"%trials)
              time.sleep(15)
           trials=trials+1
        ret=str(json.loads(mesg))
        logging.info("%s:%s"%(name,ret))
    
        if (ret.startswith("ERROR")):
              logging.info("%s:%s"%(name,ret))
              print name + ":" + ret + "\n"
              sys.exit(2);
        
    
    def checkAllJobsFinished(self, username, wkey, servicename, logging):
        data = urllib.urlencode({'func':'checkAllJobsFinished', 'username':username, 
                                 'wkey':wkey, 'servicename':servicename})
        
        self.queryAPI(data, servicename, logging) 
    
    def insertJob(self,  username, wkey, com, jobname, servicename, jobnum, result, logging):
        data = urllib.urlencode({'func':'insertJob', 'username':username, 
                                 'wkey':wkey, 'servicename':servicename, 
                                 'com':com , 'jobname':jobname, 'jobnum':str(jobnum), 
                                 'result':str(result)})
        print data
        self.queryAPI(data, "INSERT"+jobname, logging) 
    
    def updateJob(self,  username, wkey, jobname, servicename, field, jobnum, result, logging):
        data = urllib.urlencode({'func':'updateJob', 'username':username, 
                                 'wkey':wkey, 'servicename':servicename, 
                                 'field':field , 'jobname':jobname, 'jobnum':str(jobnum), 
                                 'result':str(result)})
        self.queryAPI(data, "UPDATE"+jobname, logging) 
    
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
            jobout=self.plaintext2html(re.sub('\'', "\"", jobout))
           
            data = urllib.urlencode({'func':'insertJobOut', 'username':username, 
                                 'wkey':wkey, 'jobnum':str(jobnum), 'jobout':jobout})
            self.queryAPI(data, "jobnum:"+str(jobnum), logging) 

       
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
        
   jobStat = jobStatus()
   jobStat.config = getConfig(CONFIG)
   jobStat.url = jobStat.config['url']
   
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


