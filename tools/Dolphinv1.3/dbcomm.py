#!/usr/bin/env python
 
import os, re, string, sys, commands
import warnings
import MySQLdb
from parameters import *

from sys import argv, exit, stderr
from optparse import OptionParser
 
warnings.filterwarnings('ignore', '.*the sets module is deprecated.*',
                        DeprecationWarning, 'MySQLdb')

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
        #print sql
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

def updateRunParams(runparamsid, wkey):
   
    sql = "UPDATE biocore.ngs_runparams set run_status=2, wkey='"+str(wkey)+"' where id="+str(runparamsid)
 
    return runSQL(sql)

def insertReportTable(reportfile):
   
  with open(reportfile,'r') as source:
     for line in source:
         wkey, version, type, file=re.split(r'\t+', line.rstrip())
         sql = "INSERT INTO report_list(wkey, version, type, file) VALUES ('%s', '%s','%s','%s')"%(wkey, version, type, file)
         print sql
         runSQL(sql)



def main():
    try:
        parser = OptionParser()
        parser.add_option('-r', '--runparamsid', help='runparamsid', dest='runparamsid')
        parser.add_option('-w', '--wkey', help='wkey', dest='wkey')
        parser.add_option('-i', '--insertreport', help='insert report table', dest='insertreport')
        parser.add_option('-f', '--func', help='function', dest='func')
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
    
    if (FUNC == "running"):
       updateRunParams(RUNPARAMSID, WKEY)
    elif (FUNC == "insertreport"):
       insertReportTable(INSERTREPORT)
    
    sys.exit(0)

if __name__ == "__main__":
    main()
