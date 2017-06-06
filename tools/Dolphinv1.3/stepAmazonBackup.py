#!/usr/bin/env python

import ConfigParser, os, re, string, sys, commands
import warnings
import json

import boto3
from boto3.s3.transfer import S3Transfer
from botocore.client import Config
from sys import argv, exit, stderr
from optparse import OptionParser
from binascii import hexlify, unhexlify
from simplecrypt import encrypt, decrypt

sys.path.insert(0, sys.path[0]+"/../../src")
from config import *
from funcs import *

class stepBackup:
  url=""
  f=""
  def __init__(self, url, f ):
      self.url = url
      self.f = f
      
  def getAmazonBackupList(self):  
    data = urllib.urlencode({'func':'getAmazonBackupList'})
    ret = self.f.queryAPI(self.url, data, "getAmazonBackupList")
    if (ret):
       ret=json.loads(ret)
    return ret
  
  def uploadFile(self, amazon, amazon_bucket, filename ):
    try:
       passgrab = ConfigParser.ConfigParser()
       passgrab.readfp(open(sys.path[0]+'/../default_params/.salt'))
       s3 = boto3.resource('s3', 'us-east-1',
       aws_access_key_id=decrypt(passgrab.get('Dolphin', 'AMAZON'), unhexlify(amazon['aws_access_key_id'])),
       aws_secret_access_key=decrypt(passgrab.get('Dolphin', 'AMAZON'), unhexlify(amazon['aws_secret_access_key'])),
       config=Config(signature_version='s3v4'))

       p = amazon_bucket.split("/")
       s3_file_name = "%s/%s"%(str(re.sub(p[0]+'/', '', amazon_bucket)), os.path.basename(filename))
       amazon_bucket = p[0]

       print 'Upload started[%s]=>[%s]'%(filename, s3_file_name)
       s3_bucket = s3.Bucket(str(amazon_bucket))

       s3_object = s3_bucket.Object(s3_file_name)
       mpu = s3_object.initiate_multipart_upload()

       # chunk , min size 5mb
       CHUNK = 5 * 1024 * 1024

       part_info = {
           'Parts': []
       }

       with open(filename, 'rb') as f:

          while True:

            # read chunk from file
            chunk = f.read(CHUNK)

            if not chunk:
                # complete
                mpu.complete(MultipartUpload=part_info)
                break

            # part number
            part_nr = len(part_info['Parts'])+1
            # create part
            part = mpu.Part(part_nr)

            # upload
            response = part.upload(Body=chunk)

            # add part
            part_info['Parts'].append({
                'PartNumber': part_nr,
                'ETag': response['ETag'] 
            })

    except Exception, ex:
        self.stop_err('Error (line:%s)running stepAmazonBackup.py\n%s'%(format(sys.exc_info()[-1].tb_lineno), str(ex)))

  def getAmazonCredentials(self, username):
    data = urllib.urlencode({'func':'getAmazonCredentials', 'username':username})
    ret = self.f.queryAPI(self.url, data, "getAmazonCredentials:"+str(username))
    if (len(ret)>2):
      ret=json.loads(ret)[0]
    else:
      ret=''
    return ret

  def updateMD5sum(self, table, field, file_id, filename, md5sum):
   try: 
      data = urllib.urlencode({'func':'updateMD5sum', 'table':str(table), 'field': str(field), 'id':file_id, 'md5sum':str(md5sum)})
      ret = eval(self.f.queryAPI(self.url, data, "updateMD5sum:"+str(md5sum)))
   except Exception, ex:
      self.stop_err('Error (line:%s)running stepAmazonBackup.py updateInitialFileCounts\n%s'%(format(sys.exc_info()[-1].tb_lineno), str(ex)))

  def getValfromFile(self,filename):
    val=''
    if (os.path.isfile(filename)):
    	infile = open(filename)
    	val = infile.readlines()[0].rstrip()
    else:
        print "Filename:"+filename
        sys.exit(84)
    return val

  def processUpload(self, amazon, amazon_bucket, file_id, filename):
    try:
      md5sum = self.getValfromFile(filename+'.md5sum').split(" ", 1)[0]
      self.uploadFile(amazon, amazon_bucket, filename )     
      self.updateMD5sum("amazon_backup", "checksum", file_id, filename, md5sum)
      
    except Exception, ex:
      self.stop_err('Error (line:%s)running stepBackupS3.py processFastqFiles\n%s'%(format(sys.exc_info()[-1].tb_lineno), str(ex)))

  # error
  def stop_err(self, msg ):
        sys.stderr.write( "%s\n" % msg )
        sys.exit(2)

def main():
  try:
        parser = OptionParser()
        parser.add_option('-u', '--username', help='username', dest='username')
        parser.add_option('-c', '--config', help='config parameters section', dest='config')

        (options, args) = parser.parse_args()
  except:
        parser.print_help()
        print "for help use --help"
        sys.exit(2)


  USERNAME                = options.username
  CONFIG                  = options.config

  f = funcs()
  config = getConfig(CONFIG)
  backup = stepBackup(config['url'], f)
  
  try:
    if (CONFIG == None):
        print "for help use --help"
        sys.exit(2)

    print USERNAME
    
    amazon = backup.getAmazonCredentials(USERNAME)
     
    filelist=backup.getAmazonBackupList()
    
    for file_info in filelist:
        file_id=file_info['id']
        filename=file_info['file_name']
        print filename
        amazon_bucket=file_info['s3bucket']
        amazon_bucket = str(re.sub('s3://', '', amazon_bucket))
        backup.processUpload(amazon, amazon_bucket, file_id, filename)

  except Exception, ex:
        backup.stop_err('Error (line:%s)running stepBackupS3.py\n%s'%(format(sys.exc_info()[-1].tb_lineno), str(ex)))

  sys.exit(0)

if __name__ == "__main__":
    main()
