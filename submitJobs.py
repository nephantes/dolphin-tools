# !/bin/env python
        
from optparse import OptionParser
import sys
import os
import re
import string
import subprocess
#import time

#def getKey(num):
#       return ''.join(random.choice(string.ascii_letters) for x in range(num))
   
def main():

   try:
        parser = OptionParser()
        parser.add_option('-u', '--username', help='defined user in the cluster', dest='username')
        parser.add_option('-d', '--dbhostname', help='defined hostname for the db', dest='dbhostname')
        parser.add_option('-k', '--key', help='defined key for the workflow', dest='wkey')
        parser.add_option('-s', '--servicename', help='service name', dest='servicename')
        parser.add_option('-c', '--command', help='command that is goinf to be run', dest='com')
        parser.add_option('-n', '--name', help='name of the run', dest='name')
        parser.add_option('-p', '--cpu', help='the # of cpu', dest='cpu')
        parser.add_option('-o', '--outdir', help='output directory', dest='outdir')
 	(options, args) = parser.parse_args()
   except:
        print "OptionParser Error:for help use --help"
        sys.exit(2)
   USERNAME    = options.username
   DBHOSTNAME  = options.dbhostname
   WKEY        = options.wkey 
   OUTDIR      = options.outdir 
   SERVICENAME = options.servicename
   COM         = options.com
   NAME        = options.name
   CPU         = options.cpu

   python      = "module load openssl/1.0.1e;module load python/2.7.5; python ";
   
   print WKEY
   #if (WKEY==None):
   #     WKEY=getKey(30);
   if (USERNAME==None):
        USERNAME=subprocess.check_output("whoami", shell=True).rstrip()
   
   print "USER:"+str(USERNAME)+"\n";

   if (NAME == None):
        NAME="job";
   if (OUTDIR == None):
        OUTDIR="~/out";
  
   if (DBHOSTNAME == None):
        DBHOSTNAME="galaxy.umassmed.edu"
        
   #print "COMMAND: [" + com + "]\n"
   #print "NAME: [" + name + "]\n"
   #print "cpu: [" + cpu + "]\n"


   exec_dir=os.path.dirname(os.path.abspath(__file__))
   #print "EXECDIR" + exec_dir
   sdir="/project/umw_biocore/bin/workflow/scripts"
   src=OUTDIR + "/tmp/src"
   lsf=OUTDIR + "/tmp/lsf"
  
   os.system("mkdir -p "+src)
   os.system("mkdir -p "+lsf)

   f=open(src+"/"+NAME+".submit.bash", 'w')
   f.write("#!/bin/bash\n")
   f.write("#BEGINING-OF-FILE\n")
   f.write("cd " + exec_dir + "\n")
   f.write("echo \""+COM+"\"\n")
   f.write(python+" " + sdir + "/jobStatus.py -d " + str(DBHOSTNAME) + " -u " + str(USERNAME) + " -k " + str(WKEY) + " -s " + str(SERVICENAME) + " -t dbSetStartTime -n $LSB_JOBID -j "+ str(NAME)+ " -m 2\n")
   f.write("\n\n"+ COM +"\n\n")
   f.write("retval=$?\necho \"[\"$retval\"]\"\nif [ $retval -eq 0 ]; then\n")
   f.write(python+" " + sdir + "/jobStatus.py -d " + str(DBHOSTNAME) + " -u " + str(USERNAME) + " -k " + str(WKEY) + " -s " + str(SERVICENAME) + " -t dbSetEndTime -n $LSB_JOBID -j "+ str(NAME)+ " -m 3\n")
   f.write("  echo success\nelse\n  echo failed\n  exit 127\nfi\ndate\n")

   f.write("#END-OF-FILE\n")
   f.close();
   os.system("chmod 755 "+src+"/"+NAME+".submit.bash")

   #f=open(pbs+"/"+NAME+".submit.tmp", 'w')
   #f.write("#!/bin/bash\n")

     
   #f.write("NODE:"+ str(nodename) + " part[0]:" + part[0] + "\n")
   #CHANGE this submition script according to the system.
   command="bsub -c 100 -R \"rusage[mem=1024]\" -J "+NAME+" -o "+lsf+" < "+src+"/"+NAME+".submit.bash"
   print command
   #f.write("SUBMIT SCRIPT[" + command +"]\n\n")
   #while True:
   output = os.popen(command).readlines()
   print str(output)+"\n"
        #f.write("SUBMIT OUT:[" + str(output) + "]\n")
   words = re.split('[\<\>]+', str(output))
   num = words[1]
   #     if num>0:
   #         time.sleep( 5 )
   #         break
    
   command = python+" " + sdir + "/jobStatus.py -d " + str(DBHOSTNAME) + " -u " + str(USERNAME) + " -k " + str(WKEY) + " -s " + str(SERVICENAME) + " -t dbSubmitJob -n "+ str(num) + " -j "+ str(NAME) + " -m 1 -c \"" + src+"/"+NAME+".submit.bash\"" 
   print command
   #f.write("RUN COMMAND:\n" + str(command) + "\n")
   #f.write("NUM:[" + str(num) + "]\n")
   if num>0:
      os.system(command)   
   #f.close()
 
if __name__ == "__main__":
    main()
