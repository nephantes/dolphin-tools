#!/usr/bin/env python
 
import os, re, string, sys, commands
import warnings
import MySQLdb
from parameters import *
from ZSI.client import NamedParamBinding as NPBinding, Binding
import json
import time

from sys import argv, exit, stderr
from optparse import OptionParser
 
warnings.filterwarnings('ignore', '.*the sets module is deprecated.*',
                        DeprecationWarning, 'MySQLdb')

url="http://localhost/dolphin_webservice/service.php"

def runSQL(sql):

    db = MySQLdb.connect(
      host = DBHOST,
      user = DBUSER,
      passwd = DBPASS,
      db = DB,
      port = DBPORT)
    try:
        cursor = db.cursor()
        cursor.execute(sql)
        print sql
        results = cursor.fetchall()
        cursor.close()
        del cursor

    except Exception, err:
        print "ERROR DB:for help use --help"
        db.rollback()
        sys.stderr.write('ERROR: %s\n' % str(err))
        sys.exit(2)
    finally:
        db.commit()
        db.close()
    return results

def getJobNums(wkey):
    kw = {'url':url}
    b = NPBinding(**kw)
     
    trials=0
    while trials<5:
        try:
            mesg=b.getJobNums(b=wkey)
            res=json.loads(mesg['return'])
            return res
            trials=10
        except:
            print "Couldn't connect to dolphin server (%s)"%trials
            time.sleep(15)         
        trials=trials+1
        
def insertJobStats(username, wkey, jobnum, outdir):
    file=str(outdir)+"/tmp/lsf/"+str(jobnum)+".out"
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
        
        kw = {'url':url}
        b = NPBinding(**kw)
     
        trials=0
        while trials<5:
            try:
               mesg=b.insertJobStats(a=username, c=wkey, b=jobnum, d=str(json.dumps(stats)))
               trials=10
            except:
               print "Couldn't connect to dolphin server (%s)"%trials
               logging.info("Couldn't connect to dolphin server (%s)"%trials)
               time.sleep(15)         
            trials=trials+1

def updateRunParams(runparamsid, wkey):
   
    sql = "UPDATE biocore.ngs_runparams set run_status=2, wkey='"+str(wkey)+"' where id="+str(runparamsid)
 
    return runSQL(sql)

def insertReportTable(reportfile):
  
  if os.path.isfile(reportfile): 
    with open(reportfile,'r') as source:
      for line in source:
         wkey, version, type, file=re.split(r'\t+', line.rstrip())
         sql = "INSERT INTO report_list(wkey, version, type, file) VALUES ('%s', '%s','%s','%s')"%(wkey, version, type, file)
         runSQL(sql)

def main():
    try:
        parser = OptionParser()
        parser.add_option('-r', '--runparamsid', help='runparamsid', dest='runparamsid')
        parser.add_option('-w', '--wkey', help='wkey', dest='wkey')
        parser.add_option('-i', '--insertreport', help='insert report table', dest='insertreport')
        parser.add_option('-f', '--func', help='function', dest='func')
        parser.add_option('-u', '--username', help='username', dest='username')
        parser.add_option('-o', '--outdir', help='output directory', dest='outdir')

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
    
    if (FUNC == "running"):
       updateRunParams(RUNPARAMSID, WKEY)
    elif (FUNC == "insertreport"):
       insertReportTable(INSERTREPORT)
    elif (FUNC == "insertJobStats"):
         jobnums = getJobNums(WKEY)
         for job in jobnums:
            print job
            insertJobStats(USERNAME, WKEY, job['job_num'], OUTDIR)
       
       
    sys.exit(0)

if __name__ == "__main__":
    main()
