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
   
    sql = "UPDATE biocore.ngs_runparams set wkey='"+str(wkey)+"' where id="+str(runparamsid)
 
    return runSQL(sql)


def main():
    try:
        parser = OptionParser()
        parser.add_option('-r', '--runparamsid', help='runparamsid', dest='runparamsid')
        parser.add_option('-w', '--wkey', help='wkey', dest='wkey')
        parser.add_option('-o', '--outdir', help='output directory', dest='outdir')

        (options, args) = parser.parse_args()
    except:
        parser.print_help()
        print "for help use --help"
        sys.exit(2)

    RUNPARAMSID             = options.runparamsid
    WKEY                    = options.wkey
    OUTDIR                  = options.outdir

    updateRunParams(RUNPARAMSID, WKEY)
    
    sys.exit(0)

if __name__ == "__main__":
    main()
