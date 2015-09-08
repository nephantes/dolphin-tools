#!/usr/bin/env python
 
import os, re, string, sys, commands
import warnings
import MySQLdb
from sys import argv, exit, stderr
from optparse import OptionParser
 
warnings.filterwarnings('ignore', '.*the sets module is deprecated.*',
                        DeprecationWarning, 'MySQLdb')

sys.path.insert(0, sys.path[0]+"/../../src")
from config import *

class stepGetTotalReads:
  config=''
  url=''
  def runSQL(self, sql):

    db = MySQLdb.connect(
      host = self.config['dbhost'],
      user = self.config['dbuser'],
      passwd = self.config['dbpass'],
      db = self.config['db'],
      port = int(self.config['dbport']))
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

  def getFileList(self, runparamsid, barcode):
    if (barcode != "NONE"):
        sql = "SELECT DISTINCT ns.id, ns.samplename, sf.file_name, d.fastq_dir, d.backup_dir FROM ngs_runlist nr, ngs_samples ns, ngs_temp_lane_files sf, ngs_dirs d where sf.lane_id=ns.lane_id and d.id=sf.dir_id and ns.id=nr.sample_id and nr.run_id='"+str(runparamsid)+"';"
    else:
        sql = "SELECT DISTINCT ns.id, ns.samplename, sf.file_name, d.fastq_dir, d.backup_dir FROM ngs_runlist nr, ngs_samples ns, ngs_temp_sample_files sf, ngs_dirs d where sf.sample_id=ns.id and d.id=sf.dir_id and ns.id=nr.sample_id and nr.run_id='"+str(runparamsid)+"';"
    print sql
    return self.runSQL(sql)

  def submitJob(self, JOBSUBMIT, name, command):
    child = os.popen(JOBSUBMIT+' -n stepGetTotalReads_'+name+' -c "'+command+'"')
    print JOBSUBMIT+' -n stepGetTotalReads_'+name+' -c "'+command+'"'+"\n\n"
    jobout = child.read().rstrip()
    err = child.close()
    return jobout
        
  def gzipFileAndGetCount(self, JOBSUBMIT, inputdir, filename, backupdir):
    command = "mkdir -p "+backupdir+" && gzip -c "+inputdir+"/"+filename+" > "+backupdir+"/"+filename+".gz && s=\$(zcat "+backupdir+"/"+filename+".gz|wc -l) && echo \$((\$s/4)) > "+backupdir+"/"+filename+".gz.count && md5sum "+backupdir+"/"+filename+".gz> "+backupdir+"/"+filename+".gz.md5sum"
    self.submitJob(JOBSUBMIT, "gzip_"+filename, command)

  def getCount(self, JOBSUBMIT, outputdir, inputdir, filename):
    cat = "cat"
    if ('.gz' in filename):
        cat = "zcat"
    command = "mkdir -p "+outputdir+" && s=\$("+cat+" "+inputdir+"/"+filename+"|wc -l) && echo \$((\$s/4))  > "+outputdir+"/"+filename+".count"
    self.submitJob(JOBSUBMIT, "count_"+filename, command)
    
def main():
    try:
        parser = OptionParser()
        parser.add_option('-b', '--barcode', help='barcode', dest='barcode')
        parser.add_option('-j', '--jobsubmit', help='jobsubmit', dest='jobsubmit')
        parser.add_option('-r', '--runparamsid', help='run params id', dest='runparamsid')
        parser.add_option('-u', '--username', help='username', dest='username')
        parser.add_option('-p', '--pairedend', help='pairedend', dest='paired')
        parser.add_option('-o', '--outdir', help='output directory', dest='outdir')
        parser.add_option('-c', '--config', help='config parameters section', dest='config')

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
    CONFIG                  = options.config

    totalReads = stepGetTotalReads()
    totalReads.config = getConfig(CONFIG)
    totalReads.url = totalReads.config['url']


    if (OUTDIR == None or JOBSUBMIT == None):
        print "for help use --help"
        sys.exit(2)

    print BARCODE
    print PAIRED
    print JOBSUBMIT
    print OUTDIR
    print USERNAME
    print RUNPARAMSID
    
    filelist=totalReads.getFileList(RUNPARAMSID, BARCODE)
    
    inputdir=OUTDIR+"/input"
    if (BARCODE != "NONE"):
        inputdir=OUTDIR+"/seqmapping/barcode"
    
    processedLibs=[]
    for sample in filelist:
        
        libid=sample[0]
        libname=sample[1]
        filename=sample[2]
        fastq_dir=sample[3]
        backup_dir=sample[4]
        
        print libid
        print libname
        print filename
        print fastq_dir
        print backup_dir
        
        if (filename.find(',')!=-1):
            files=filename.split(',')
            totalReads.getCount(JOBSUBMIT, inputdir + "/tmp", fastq_dir,  files[0])
            totalReads.getCount(JOBSUBMIT, inputdir + "/tmp", fastq_dir,  files[1])
            if (not libname in processedLibs):
                totalReads.gzipFileAndGetCount(JOBSUBMIT, inputdir, libname+".1.fastq", backup_dir)
                totalReads.gzipFileAndGetCount(JOBSUBMIT, inputdir, libname+".2.fastq", backup_dir)
        else:
            totalReads.getCount(JOBSUBMIT, inputdir + "/tmp", fastq_dir,  filename)
            if (not libname in processedLibs):
                totalReads.gzipFileAndGetCount(JOBSUBMIT, inputdir, libname+".fastq", backup_dir)
        processedLibs.append(libname)
        
    sys.exit(0)

if __name__ == "__main__":
    main()
