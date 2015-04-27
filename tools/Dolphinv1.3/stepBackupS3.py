#!/usr/bin/env python
 
import os, re, string, sys, commands
import warnings
import MySQLdb
import itertools

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

def getSQLval(sql):
    results = runSQL(sql)
    for row in results:
        return row[0]

  

def getSampleList(runparams_id, barcode):
    if (barcode != "NONE"):
        sql = "SELECT DISTINCT ns.id, sf.id, d.id, ns.lane_id, ns.name, sf.file_name, d.backup_dir, d.amazon_bucket, ns.owner_id, ns.group_id, ns.perms FROM biocore.ngs_runlist nr, ngs_samples ns, ngs_temp_lane_files sf, ngs_dirs d where sf.lane_id=ns.lane_id and d.id=sf.dir_id and ns.id=nr.sample_id and nr.run_id='"+str(runparams_id)+"';"
    else:
        sql = "SELECT DISTINCT ns.id, sf.id, d.id, ns.lane_id, ns.name, sf.file_name, d.backup_dir, d.amazon_bucket, ns.owner_id, ns.group_id, ns.perms FROM biocore.ngs_runlist nr, ngs_samples ns, ngs_temp_sample_files sf, ngs_dirs d where sf.sample_id=ns.id and d.id=sf.dir_id and ns.id=nr.sample_id and nr.run_id='"+str(runparams_id)+"';"

    return runSQL(sql)

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


def getAmazonCredentials(clusteruser):
    sql = 'SELECT DISTINCT ac.* FROM biocore.amazon_credentials ac, biocore.group_amazon ga, biocore.users u where ac.id=ga.amazon_id and ga.group_id=u.group_id and u.clusteruser="'+clusteruser+'";'
    results = runSQL(sql)

    return results

def updateInitialFileCounts(file_id, tablename, inputdir, filename, paired):
    
    count=-1
    if (paired != "None"):
    	files=filename.split(',')
    	firstCount=getValfromFile(inputdir+'/tmp/'+files[0]+'.count')
    	secondCount=getValfromFile(inputdir+'/tmp/'+files[1]+'.count')

        #print "f:[%s]:s:[%s]"%(firstCount, secondCount)
    	if (firstCount == secondCount):
    	    count=firstCount
    	else:
    	    print "ERROR 85: The # of reads in paired end libary not equal"
            print "%s"%(inputdir+'/tmp/'+files[0]+'.count')
            print "%s"%(inputdir+'/tmp/'+files[1]+'.count')
    	    sys.exit(85)
    else:
    	count=getValfromFile(inputdir+'/tmp/'+filename+'.count')
    sql = "UPDATE "+tablename+" set total_reads='"+count+"' where id="+str(file_id)
    runSQL(sql)


def getValfromFile(filename):
    val=''
    if (os.path.isfile(filename)):
    	infile = open(filename)
    	val = infile.readlines()[0].rstrip()
    else:
    	sys.exit(84)
    return val

def getFastqFileId(sample):
    sample_id=sample[0]
    sql="select id from ngs_fastq_files where sample_id="+str(sample_id)
    return getSQLval(sql)

def upadateFastqFile(fastq_id, md5sum, count, filename, sample):
    sample_id=sample[0]
    owner_id=sample[8]
    sql="update ngs_fastq_files set checksum='"+md5sum+"', total_reads='"+count+"', date_modified=now(), last_modified_user="+str(owner_id)+" where sample_id="+str(sample_id)+" and id="+str(fastq_id)
    runSQL(sql)
    
def insertFastqFile(checksum, total_reads, filename, sample):
    (sample_id, sf_id, dir_id, lane_id, libname, org_file_name, backup_dir, amazon_bucket, owner_id, group_id, perms) = sample
    sql="INSERT INTO `biocore`.`ngs_fastq_files` ( `file_name`, `total_reads`, `checksum`, `sample_id`, `lane_id`,`dir_id`,`owner_id`,`group_id`,`perms`,`date_created`,`date_modified`,`last_modified_user`) VALUES('"+filename+"','"+total_reads+"','"+checksum+"','"+str(sample_id)+"','"+str(lane_id)+"','"+str(dir_id)+"','"+str(owner_id)+"', '"+str(group_id)+"', '"+str(perms)+"', now(), now(), '"+str(owner_id)+"');"
    runSQL(sql)
    
def checkReadCounts(sample_id, tablename):
    sql='SELECT sum(total_reads) FROM biocore.ngs_temp_sample_files where sample_id='+str(sample_id)
    print sql
    temp_count=getSQLval(sql)
    print temp_count
    sql='SELECT total_reads FROM biocore.ngs_fastq_files where sample_id='+str(sample_id)
    merged_count=getSQLval(sql)
    print merged_count
    if (temp_count==merged_count):
    	return 1
    else:
        return 0

    
def processFastqFiles(sample, paired):
    md5sum=''
    count=''
    libname=sample[4]
    backup_dir=sample[6]
    filename=''
 
    if (paired != 'None'):
        firstCount=getValfromFile(backup_dir+'/'+libname+'.1.fastq.gz.count')
        secondCount=getValfromFile(backup_dir+'/'+libname+'.2.fastq.gz.count')
        if (firstCount == secondCount):
            count=firstCount
            filename= libname+'.1.fastq.gz,'+libname+'.2.fastq.gz'
            firstmd5sum=getValfromFile(backup_dir+'/'+libname+'.1.fastq.gz.md5sum').split(' ')[0]
            secondmd5sum=getValfromFile(backup_dir+'/'+libname+'.2.fastq.gz.md5sum').split(' ')[0]
            md5sum=firstmd5sum+','+secondmd5sum
        else:
            print "ERROR 85: The # of reads in paired end libary not equal"
            print "%s"%(backup_dir+'/'+libname+'.1.fastq.gz.count')
            sys.exit(85)
    else:
        filename= libname+'.fastq.gz'
        count=getValfromFile(backup_dir+'/'+libname+'.fastq.gz.count')
        md5sum=getValfromFile(backup_dir+'/'+libname+'.fastq.gz.md5sum')
    fastq_id=getFastqFileId(sample)
    if (fastq_id>0):
    	upadateFastqFile(fastq_id, md5sum, count, filename, sample)
    else:
    	insertFastqFile(md5sum, count, filename, sample)

def main():
    try:
        parser = OptionParser()
        parser.add_option('-b', '--barcode', help='barcode', dest='barcode')
        parser.add_option('-j', '--jobsubmit', help='jobsubmit', dest='jobsubmit')
        parser.add_option('-r', '--runparamsid', help='group id', dest='runparamsid')
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
    RUNPARAMSID             = options.runparamsid
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

    samplelist=getSampleList(RUNPARAMSID, BARCODE)

    inputdir=OUTDIR+"/input"
    tablename="ngs_temp_sample_files"
    if (BARCODE != "NONE"):
        inputdir=OUTDIR+"/seqmapping/barcode"
        tablename="ngs_temp_lane_files"
    
    processedLibs=[]
    for sample in samplelist:
        sample_id=sample[0]
        file_id=sample[1]
        libname=sample[4]  
        filename=sample[5]
        backup_dir=sample[6]
        amazon_bucket=sample[7]
        amazon_bucket = re.sub('s3://'+amazon[0][3]+'/', '', amazon_bucket)

        updateInitialFileCounts(file_id, tablename, inputdir, filename, PAIRED)
        if (not [libname, sample_id] in processedLibs):
            processFastqFiles(sample, PAIRED)
            processedLibs.append([libname, sample_id])

    for libname, sample_id in processedLibs:
        if (checkReadCounts(sample_id, tablename)):

            if (filename.find(',')!=-1):
                files=filename.split(',')
                uploadFile(pb, amazon_bucket, backup_dir, libname+'.1.fastq.gz' )
                uploadFile(pb, amazon_bucket, backup_dir, libname+'.2.fastq.gz' )
            else:
                uploadFile(pb, amazon_bucket, backup_dir, libname+'.fastq.gz' )
        else:
            print "ERROR 86: The # of read counts doesn't match: %s",libname
            sys.exit(86)
    sys.exit(0)

if __name__ == "__main__":
    main()
