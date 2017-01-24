import os, re, string, sys, commands
import logging
import urllib,urllib2
import warnings
import json
import boto3
import MySQLdb
import ConfigParser
from binascii import hexlify, unhexlify
from simplecrypt import encrypt, decrypt
from boto3.s3.transfer import S3Transfer
from botocore.client import Config
from sys import argv, exit, stderr
from optparse import OptionParser
sys.path.insert(0, sys.path[0]+"/../../src")
from config import *
from funcs import *

class botoDownload:
    f=""
    config = ConfigParser.ConfigParser()
    url = ""
    
    def __init__(self, url, f):
        self.f = f
        self.url = url 
        
    def getAmazonCredentials(self, username):
        data = urllib.urlencode({'func':'getAmazonCredentials', 'username':username})
        ret = self.f.queryAPI(self.url, data, "getAmazonCredentials:"+str(username))
        if (len(ret)>2):
            ret=json.loads(ret)[0]
        else:
            ret=''
        return ret
    
    def downloadFile(self, amazon, s3_filename, outfile):
        try:
            saltconfig = ConfigParser.ConfigParser()
            saltconfig.readfp(open(sys.path[0]+'/../default_params/.salt'))
            saltpass = saltconfig.get('Dolphin','AMAZON')
            s3 = boto3.client('s3', 'us-east-1',
            aws_access_key_id=decrypt(saltpass, unhexlify(amazon['aws_access_key_id'])),
            aws_secret_access_key=decrypt(saltpass, unhexlify(amazon['aws_secret_access_key'])))
            
            print 'Download started[%s]=>[%s]'%(s3_filename, outfile)
            inputarr = s3_filename.split("/")
            amazon_bucket = inputarr[2]
            inputfile = "/".join(inputarr[3:])
            print 'Bucket[%s]=>[%s]'%(amazon_bucket, inputfile)
            s3.download_file(amazon_bucket, inputfile, outfile)	
            
        except Exception, ex:
            print ex
            
def main():
    f = funcs()
    
    #define options
    parser = OptionParser()
    parser.add_option("-i", "--inputfile", dest="filename")
    parser.add_option("-o", "--outputfile", dest="output")
    parser.add_option("-u", "--username", dest="username")
    parser.add_option("-c", "--config", dest="config")
    # parse
    options, args = parser.parse_args()
    # retrieve options
    FILE  = options.filename
    OUTPUT    = options.output
    USERNAME  = options.username
    CONFIG = options.config

    config = getConfig(CONFIG)
    boto = botoDownload(config['url'], f)
    amazon = boto.getAmazonCredentials(USERNAME)

    boto.downloadFile(amazon, FILE, OUTPUT)
    
    sys.exit(0)

if __name__ == "__main__":
    main()
    
