#!/usr/bin/env python
 
import os, re, string, sys, commands
import warnings
import MySQLdb

from sys import argv, exit, stderr
from optparse import OptionParser
 
from boto.s3.connection import S3Connection
from boto.s3.key import Key

warnings.filterwarnings('ignore', '.*the sets module is deprecated.*',
                        DeprecationWarning, 'MySQLdb')

def runSQL(sql):
    port=3306
    db = MySQLdb.connect(
      host = 'galaxy.umassmed.edu',
      user = 'biocore',
      passwd = 'biocore2013',
      db = 'biocore',
      port = port)
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

def getFileList():
    try:
        sql = "SELECT d.fastq_dir, f.file_name, d.amazon_bucket FROM biocore.ngs_fastq_files f, biocore.ngs_dirs d where d.id=f.dir_id;"
        results = runSQL(sql)
    except Exception, err:
        sys.stderr.write('ERROR: %s\n' % str(err))
        sys.exit(2)
    return results

def uploadFile(pb, amazon_bucket, fastq_dir, filename ):
    k = Key(pb)
    inputfile = "%s/%s"%(fastq_dir.strip(), filename.strip())
    s3file = "%s/%s"%(amazon_bucket.strip(), filename.strip())
   
    print inputfile
    print s3file
    m=pb.get_key(s3file)
    
    if m and m.exists():
        print "Already uploaded %s" % s3file
    else:
        print 'Upload started=>'
        k.name = s3file
        k.set_contents_from_filename(inputfile)

def getSampleList(run_group_id, barcode):
    if (barcode != "NONE"):
        sql = "SELECT DISTINCT ns.id, ns.name, d.fastq_dir, d.backup_dir, d.amazon_bucket FROM biocore.ngs_runlist nr, ngs_samples ns, ngs_temp_lane_files sf, ngs_dirs d where sf.lane_id=ns.lane_id and d.id=sf.dir_id and ns.id=nr.sample_id and nr.run_Group_id='"+str(run_group_id)+"';"
    else:
        sql = "SELECT DISTINCT ns.id, ns.name, d.fastq_dir, d.backup_dir, d.amazon_bucket FROM biocore.ngs_runlist nr, ngs_samples ns, ngs_temp_sample_files sf, ngs_dirs d where sf.sample_id=ns.id and d.id=sf.dir_id and ns.id=nr.sample_id and nr.run_Group_id='"+str(run_group_id)+"';"

    return runSQL(sql)

def getAmazonCredentials(clusteruser):
    sql = 'SELECT DISTINCT ac.* FROM biocore.amazon_credentials ac, biocore.group_amazon ga, biocore.users u where ac.id=ga.amazon_id and ga.group_id=u.group_id and u.clusteruser="'+clusteruser+'";'
    results = runSQL(sql)

    return results

def main():
    try:
        parser = OptionParser()
        parser.add_option('-b', '--barcode', help='barcode', dest='barcode')
        parser.add_option('-j', '--jobsubmit', help='jobsubmit', dest='jobsubmit')
        parser.add_option('-g', '--groupid', help='group id', dest='groupid')
        parser.add_option('-u', '--username', help='username', dest='username')
        parser.add_option('-p', '--pairedend', help='pairedend', dest='paired')
        parser.add_option('-o', '--outdir', help='output directory', dest='outdir')

        (options, args) = parser.parse_args()
    except:
        parser.print_help()
        print "for help use --help"
        sys.exit(2)

    BARCODE                 = options.barcode
    PAIRED                  = options.paired
    USERNAME                = options.username
    RUNGROUPID              = options.groupid
    JOBSUBMIT               = options.jobsubmit
    OUTDIR                  = options.outdir



    if (OUTDIR == None or JOBSUBMIT == None):
        print "for help use --help"
        sys.exit(2)

    print BARCODE
    print PAIRED
    print JOBSUBMIT
    print OUTDIR
    print USERNAME
    
    
    amazon = getAmazonCredentials(USERNAME)

    if (amazon!=()):   
       conn = S3Connection(amazon[0][1], amazon[0][2])
       pb = conn.get_bucket(amazon[0][3])

    samplelist=getSampleList(RUNGROUPID, BARCODE)
    print samplelist


    sys.exit(0)
    
    for sample in samplelist:
        
        filename=sample[1]
        backup_dir=sample[3]
        amazon_bucket=sample[4]
        calcHashFilesInCluster(OUTDIR, BARCODE, filename)
        #amazon_bucket=result[2]
        #if (filename.find(',')!=-1):
        #    files=filename.split(',')
        #    uploadFile(pb, amazon_bucket, fastq_dir, files[0] )
        #    uploadFile(pb, amazon_bucket, fastq_dir, files[1] )
        #else:
        #    uploadFile(pb, amazon_bucket, fastq_dir, filename )
        sys.exit(0)
    sys.exit(0)

if __name__ == "__main__":
    main()
