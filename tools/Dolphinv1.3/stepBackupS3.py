#!/usr/bin/env python
 
import os, re, string, sys, commands
import warnings
import itertools
import json

from boto.s3.connection import S3Connection
from boto.s3.key import Key
from sys import argv, exit, stderr
from optparse import OptionParser

sys.path.insert(0, sys.path[0]+"/../../src")
from config import *
from funcs import *

class stepBackup:
  url=""
  f=""
  def __init__(self, url, f ):
      self.url = url
      self.f = f
      
  def getSampleList(self, runparamsid, barcode):  
    data = urllib.urlencode({'func':'getSampleList', 'runparamsid':str(runparamsid), 'barcode':str(barcode)})
    ret = eval(self.f.queryAPI(self.url, data, "getSampleList:"+runparamsid))
    return ret
  
  def uploadFile(self, pb, amazon_bucket, fastq_dir, filename ):
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


  def getAmazonCredentials(self, username):
    data = urllib.urlencode({'func':'getAmazonCredentials', 'username':username})
    ret = self.f.queryAPI(self.url, data, "getAmazonCredentials:"+str(username))
    if (ret):
      ret=json.loads(ret)
    return ret

  def updateInitialFileCounts(self, file_id, tablename, inputdir, filename, paired):
   try: 
      count=-1
      if (paired):
          files=filename.split(',')
          firstCount=self.getValfromFile(inputdir+'/tmp/'+files[0]+'.count')
          secondCount=self.getValfromFile(inputdir+'/tmp/'+files[1]+'.count')
  
          if (firstCount == secondCount):
              count=firstCount
          else:
              print "ERROR 85: The # of reads in paired end libary not equal"
              print "%s"%(inputdir+'/tmp/'+files[0]+'.count')
              print "%s"%(inputdir+'/tmp/'+files[1]+'.count')
              sys.exit(85)
      else:
          count=self.getValfromFile(inputdir+'/tmp/'+filename+'.count')
          data = urllib.urlencode({'func':'updateInitialFileCounts', 'tablename':str(tablename), 'total_reads':str(count), 'file_id':str(file_id)})
          ret = eval(self.f.queryAPI(self.url, data, "updateInitialFileCounts:"+str(count)))
   except Exception, ex:
      self.stop_err('Error (line:%s)running stepBackupS3.py updateInitialFileCounts\n%s'%(format(sys.exc_info()[-1].tb_lineno), str(ex)))



  def getValfromFile(self,filename):
    val=''
    if (os.path.isfile(filename)):
    	infile = open(filename)
    	val = infile.readlines()[0].rstrip()
    else:
    	sys.exit(84)
    return val

  def getFastqFileId(self, sample):
    data = urllib.urlencode({'func':'getFastqFileId', 'file_id':str(sample['file_id'])})
    ret = json.loads(self.f.queryAPI(self.url, data, "getFastqFileId:"+str(sample['file_id'])))
    return ret['id']

  def upadateFastqFile(self, fastq_id, md5sum, count, filename, sample):
    sample_id=sample['sample_id']
    owner_id=sample['owner_id']
    data = urllib.urlencode({'func':'upadateFastqFile', 'sample_id':str(sample_id), 'owner_id':str(owner_id), 'md5sum':str(md5sum), 'total_reads':str(count), 'fastq_id':str(fastq_id)})
    ret = eval(self.f.queryAPI(self.url, data, "upadateFastqFile:"+str(sample_id)))

    
  def insertFastqFile(self, checksum, total_reads, filename, sample):
    data = urllib.urlencode({'func':'insertFastqFile', 'filename':str(filename),'total_reads':str(total_reads),'checksum':str(checksum), 'sample_id':str(sample['sample_id']),'lane_id':str(sample['lane_id']),'dir_id': str(sample['dir_id']),'owner_id':str(sample['owner_id']), 'group_id':str(sample['group_id']), 'perms':str(sample['perms'])})
    ret = eval(self.f.queryAPI(self.url, data, "insertFastqFile:"+str(sample['sample_id'])))

    
  def checkReadCounts(self, sample_id, tablename):
    data = urllib.urlencode({'func':'checkReadCounts', 'sample_id':str(sample_id), 'tablename':str(tablename)})
    print data
    ret = self.f.queryAPI(self.url, data, "checkReadCounts:"+str(sample_id))
    ret = json.loads(ret)
    print ret
    return ret['Result']
    
    
  def processFastqFiles(self, sample, paired):
    try:
      md5sum=''
      count=''
      libname=sample['samplename']
      backup_dir=sample['backup_dir']
      filename=''

      if (paired):
          firstCount=self.getValfromFile(backup_dir+'/'+libname+'.1.fastq.gz.count')
          secondCount=self.getValfromFile(backup_dir+'/'+libname+'.2.fastq.gz.count')
          if (firstCount == secondCount):
              count=firstCount
              filename= libname+'.1.fastq.gz,'+libname+'.2.fastq.gz'
              firstmd5sum=self.getValfromFile(backup_dir+'/'+libname+'.1.fastq.gz.md5sum').split(' ')[0]
              secondmd5sum=self.getValfromFile(backup_dir+'/'+libname+'.2.fastq.gz.md5sum').split(' ')[0]
              md5sum=firstmd5sum+','+secondmd5sum
          else:
              print "ERROR 85: The # of reads in paired end libary not equal"
              print "%s"%(backup_dir+'/'+libname+'.1.fastq.gz.count')
              sys.exit(85)
      else:
          filename= libname+'.fastq.gz'
          count=self.getValfromFile(backup_dir+'/'+libname+'.fastq.gz.count')
          md5sum=self.getValfromFile(backup_dir+'/'+libname+'.fastq.gz.md5sum')
      fastq_id=self.getFastqFileId(sample)
      
      if (fastq_id>0):
          self.upadateFastqFile(fastq_id, md5sum, count, filename, sample)
      else:
          self.insertFastqFile(md5sum, count, filename, sample)

    except Exception, ex:
      self.stop_err('Error (line:%s)running stepBackupS3.py processFastqFiles\n%s'%(format(sys.exc_info()[-1].tb_lineno), str(ex)))


  # error
  def stop_err(self, msg ):
        sys.stderr.write( "%s\n" % msg )
        sys.exit(2)

def main():
  try:
        parser = OptionParser()
        parser.add_option('-b', '--barcode', help='barcode', dest='barcode')
        parser.add_option('-j', '--jobsubmit', help='jobsubmit', dest='jobsubmit')
        parser.add_option('-r', '--runparamsid', help='group id', dest='runparamsid')
        parser.add_option('-u', '--username', help='username', dest='username')
        parser.add_option('-p', '--pairedend', help='pairedend', dest='paired')
        parser.add_option('-a', '--amazonupload', help='amazonupload', dest='amazonupload')
        parser.add_option('-o', '--outdir', help='output directory', dest='outdir')
        parser.add_option('-c', '--config', help='config parameters section', dest='config')

        (options, args) = parser.parse_args()
  except:
        parser.print_help()
        print "for help use --help"
        sys.exit(2)


  try:
    BARCODE                 = options.barcode
    PAIRED                  = options.paired
    USERNAME                = options.username
    RUNPARAMSID             = options.runparamsid
    JOBSUBMIT               = options.jobsubmit
    OUTDIR                  = options.outdir
    AMAZONUPLOAD            = options.amazonupload
    CONFIG                  = options.config

    f = funcs()
    config = getConfig(CONFIG)
    backup = stepBackup(config['url'], f)
   
    if (OUTDIR == None or JOBSUBMIT == None):
        print "for help use --help"
        sys.exit(2)

    print BARCODE
    print PAIRED
    print JOBSUBMIT
    print OUTDIR
    print USERNAME
    
    amazon = backup.getAmazonCredentials(USERNAME)
    
    if (amazon):
       conn = S3Connection(amazon['aws_access_key_id'], amazon['aws_secret_access_key'])
       pb = conn.get_bucket(amazon['bucket'])
     
    samplelist=backup.getSampleList(RUNPARAMSID, BARCODE)
    
    inputdir=OUTDIR+"/input"
    tablename="ngs_temp_sample_files"
    if (BARCODE != "NONE"):
        inputdir=OUTDIR+"/seqmapping/barcode"
        tablename="ngs_temp_lane_files"
    
    processedLibs=[]
    amazon_bucket=""
    for sample in samplelist:
        if (type(sample) is str):
            sample = samplelist
       
        sample_id=sample['sample_id']
        file_id=sample['file_id']
        libname=sample['samplename']  
        filename=sample['file_name']
        backup_dir=sample['backup_dir']
        amazon_bucket=sample['amazon_bucket']
        PAIRED=None
        if (filename.find(',')!=-1):
           PAIRED="Yes"
    
        backup.updateInitialFileCounts(file_id, tablename, inputdir, filename, PAIRED)
        if (not [libname, sample_id] in processedLibs):
            backup.processFastqFiles(sample, PAIRED)
            processedLibs.append([libname, sample_id])

    if (AMAZONUPLOAD.lower() != "no" and amazon!=() and amazon_bucket!=""):
      amazon_bucket = re.sub('s3://'+amazon['bucket']+'/', '', amazon_bucket)
      print amazon_bucket
      for libname, sample_id in processedLibs:
        print libname + ":" + str(sample_id)
        if (backup.checkReadCounts(sample_id, tablename)):
            if (filename.find(',')!=-1):
                files=filename.split(',')
                backup.uploadFile(pb, amazon_bucket, backup_dir, libname+'.1.fastq.gz' )
                backup.uploadFile(pb, amazon_bucket, backup_dir, libname+'.2.fastq.gz' )
            else:
                backup.uploadFile(pb, amazon_bucket, backup_dir, libname+'.fastq.gz' )
        else:
            print "ERROR 86: The # of read counts doesn't match: %s",libname
            sys.exit(86)
  except Exception, ex:
        backup.stop_err('Error (line:%s)running stepBackupS3.py\n%s'%(format(sys.exc_info()[-1].tb_lineno), str(ex)))

  sys.exit(0)

if __name__ == "__main__":
    main()
