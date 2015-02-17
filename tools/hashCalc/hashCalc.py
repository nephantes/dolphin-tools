#!/share/bin/python

from optparse import OptionParser
import warnings

warnings.filterwarnings('ignore', '.*the sets module is deprecated.*',
                        DeprecationWarning, 'MySQLdb')
import os
import re
import MySQLdb
import sys
import re
import cgi
   
def runSQL(dbhostname, sql): 
    port=3306
    db = MySQLdb.connect(
      host = dbhostname, 
      user = 'biocore', 
      passwd = 'biocore2013', 
      db = 'biocore', 
      port = port)
    #db = MySQLdb.connect(dbhostname,"biocore","biocore2013","biocore" )
    try:
        cursor = db.cursor()

        #print sql
        cursor.execute(sql)
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

def calcHash(inputfile):
   command="md5sum "+inputfile
   child = os.popen(command)
   jobout = child.read().rstrip()
   err = child.close()
   hashstr = jobout.split(" ")[0]
   print hashstr 
   return hashstr

def runCalcHash(dbhostname,  input, dirname, outdir):
   files=input.split(',')
   if (len(files)>1):
     hashstr=calcHash(dirname+"/"+files[0].strip())
     hashstr=hashstr+","+calcHash(dirname+"/"+files[1].strip())
   else:
     hashstr=calcHash(dirname+"/"+input.strip())
    
   sql  = "UPDATE  biocore.ngs_fastq_files nff, "
   sql += "(SELECT nff.id FROM biocore.ngs_fastq_files nff, biocore.ngs_dirs nd "
   sql += "where nff.dir_id = nd.id AND "
   sql += "nff.file_name='"+str(input)+"' and "
   sql += "nd.fastq_dir='"+str(dirname)+"') a "
   sql += "set nff.checksum='"+str(hashstr)+"' "
   sql += "where a.id=nff.id"

   results = runSQL(dbhostname, sql)

def main():
   try:
        parser = OptionParser()
        parser.add_option('-b', '--dbhostname', help='defined dbhostname for the db', dest='dbhostname')
        parser.add_option('-d', '--dirname', help='dirname', dest='dirname')
        parser.add_option('-i', '--inputfiles', help='inputfiles', dest='inputfiles')
	parser.add_option('-o', '--outdir', help='output directory', dest='outdir')

        (options, args) = parser.parse_args()
   except:
	parser.print_help()
        print "for help use --help"
        sys.exit(2)

   DBHOSTNAME              = options.dbhostname
   DIRNAME                 = options.dirname
   INPUT                   = options.inputfiles
   OUTDIR                  = options.outdir
  
      

   if (DBHOSTNAME == None):
        DBHOSTNAME="galaxy.umassmed.edu"
   if (INPUT == None or DIRNAME == None):
     print "for help use --help"
     sys.exit(2)

   print "INPUT:"+INPUT
   print "DIRNAME:"+DIRNAME        
   runCalcHash(DBHOSTNAME, INPUT, DIRNAME, OUTDIR)

if __name__ == "__main__":
    main()


