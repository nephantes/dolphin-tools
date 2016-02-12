#!/usr/bin/env python
 
import os, re, string, sys
import warnings
import json
import time
import urllib,urllib2
from sys import argv, exit, stderr
from optparse import OptionParser

sys.path.insert(0, sys.path[0]+"/../../src")
from config import *
from funcs import *

class dbcomm:
  url=""
  f=""
  def __init__(self, url, f ):
      self.url = url
      self.f = f
        
  def getJobNums(self, wkey):
    data = urllib.urlencode({'func':'getJobNums', 'wkey':wkey})
    ret = json.loads(self.f.queryAPI(self.url, data, "getJobNums:"+wkey))
    return ret
    
  def insertJobStats(self, username, wkey, jobnum, outdir):
    file=str(outdir)+"/tmp/logs/"+str(jobnum)+".out"
    if os.path.isfile(file) and os.access(file, os.R_OK):
        lines = [line.rstrip('\n') for line in open(file)]
        stats = {
                 'CPU time' : 0,
                 'Max Memory' : 0,
                 'Average Memory': 0,
                 'Total Requested Memory' : 0,
                 'Delta Memory' : 0,
                 'Max Processes' : 0,
                 'Max Threads' : 0
                 }
     
        for line in lines:
           if re.match("\s*(.*)\s:\s*(.*)\s(sec.|MB|)", line):
              m = re.match("\s*(.*)\s:\s*([^\s]*)\s?(sec.|MB|)", line)
              stats[m.groups()[0]] = m.groups()[1]
        
        data = urllib.urlencode({'func':'insertJobStats', 'username':username, 'wkey':wkey, 'jobnum':jobnum, 'stats':str(json.dumps(stats)) })
        print data
        self.f.queryAPI(self.url, data, "insertJobStats:"+wkey ) 

  def insertReportTable(self, reportfile):
 
  
    if os.path.isfile(reportfile): 
      with open(reportfile,'r') as source:
        for line in source:
          wkey, version, type, file=re.split(r'\t+', line.rstrip())
          data = urllib.urlencode({'func':'insertReportTable', 'wkey':wkey, 'version':version, 'type':type, 'file':file})
          self.f.queryAPI(self.url, data, "insertReportTable:"+wkey )

def main():
    try:
        parser = OptionParser()
        parser.add_option('-r', '--runparamsid', help='runparamsid', dest='runparamsid')
        parser.add_option('-w', '--wkey', help='wkey', dest='wkey')
        parser.add_option('-i', '--insertreport', help='insert report table', dest='insertreport')
        parser.add_option('-f', '--func', help='function', dest='func')
        parser.add_option('-u', '--username', help='username', dest='username')
        parser.add_option('-o', '--outdir', help='output directory', dest='outdir')
        parser.add_option('-c', '--config', help='config', dest='config')

        (options, args) = parser.parse_args()
    except:
        parser.print_help()
        print "for help use --help"
        sys.exit(2)

    RUNPARAMSID             = options.runparamsid
    INSERTREPORT            = options.insertreport
    FUNC                    = options.func
    WKEY                    = options.wkey
    OUTDIR                  = options.outdir
    USERNAME                = options.username
    CONFIG                  = options.config
   
    f = funcs()
    config = getConfig(CONFIG)
    print config['url']
    dbcon = dbcomm(config['url'], f)

    if (FUNC == "insertreport"):
       dbcon.insertReportTable(INSERTREPORT)
    elif (FUNC == "insertJobStats"):
         jobnums = dbcon.getJobNums(WKEY)
         for job in jobnums:
            dbcon.insertJobStats(USERNAME, WKEY, job['job_num'], OUTDIR)
    sys.exit(0)

if __name__ == "__main__":
    main()
