#!/bin/env python
        
from optparse import OptionParser
import sys
import os
import re
import string
import subprocess
#import time
import math
import json

sys.path.insert(0, sys.path[0])
from config import *
from funcs import *

#def getKey(num):
#       return ''.join(random.choice(string.ascii_letters) for x in range(num))
sys.path.insert(0, sys.path[0])
from config import *
from funcs import *


class submitJobs:
    url = ""
    f=""
    def __init__(self, url, f ):
        self.url = url
        self.f = f
    
    def checkJob(self, jobname, wkey, logging):
        data = urllib.urlencode({'func':'checkJob', 'jobname':jobname, 'wkey':wkey})
        return json.loads(self.f.queryAPI(self.url, data, jobname, logging))['Result']
                
    def runcmd(self, command, logging): 
       print command
       logging.info("\n\n\nCOM:["+command+"]\n\n\n")
       child = os.popen(command)
       data = child.read()
       logging.info("\n\n\ndata:["+data+"]\n\n\n")
       print data
       err = child.close()
       if err:
           logging.info('ERROR: %s failed w/ exit code %d' % (command, err))
           return 'ERROR: %s failed w/ exit code %d' % (command, err)
           
       return data

def main():

   try:
        parser = OptionParser()
        parser.add_option('-u', '--username', help='defined user in the cluster', dest='username')
        parser.add_option('-k', '--key', help='defined key for the workflow', dest='wkey')
        parser.add_option('-s', '--servicename', help='service name', dest='servicename')
        parser.add_option('-c', '--command', help='command that is goinf to be run', dest='com')
        parser.add_option('-n', '--name', help='name of the run', dest='name')
        parser.add_option('-o', '--outdir', help='output directory', dest='outdir')
        parser.add_option('-f', '--config', help='configuration parameter section', dest='config')
        (options, args) = parser.parse_args()
   except:
        print "OptionParser Error:for help use --help"
        sys.exit(2)
   USERNAME    = options.username
   WKEY        = options.wkey 
   OUTDIR      = options.outdir 
   SERVICENAME = options.servicename
   COM         = options.com
   NAME        = options.name
   CONFIG      = options.config
   python      = "python"
   
   config=getConfig(CONFIG)
   f = funcs()
   submitjobs = submitJobs(config['url'], f)
   

   if (USERNAME==None):
        USERNAME=subprocess.check_output("whoami", shell=True).rstrip()
   
   print "USER:"+str(USERNAME)+"\n";

   sdir=os.path.dirname(sys.argv[0])
   exec_dir=os.path.dirname(os.path.abspath(__file__))
   #print "EXECDIR" + exec_dir
   track=OUTDIR + "/tmp/track"
   src=OUTDIR + "/tmp/src"
   logs=OUTDIR + "/tmp/logs"
  
   os.system("mkdir -p "+track)
   os.system("mkdir -p "+src)
   os.system("mkdir -p "+logs)
   
   logfile="%s/JOB.%s.log"%(logs, NAME)
   logging.basicConfig(filename=logfile, filemode='a',format='%(asctime)s %(message)s', datefmt='%m/%d/%Y %I:%M:%S %p', level=logging.DEBUG)
   logging.info("File Path:%s"%os.getcwd())
   
   success_file = track+"/"+str(NAME)+".success";
   jobstatus_cmd = "python %(sdir)s/jobStatus.py -f %(CONFIG)s -u %(USERNAME)s -k %(WKEY)s -s %(SERVICENAME)s -t %(TYPE)s -o %(OUTDIR)s -j %(NAME)s -m %(MESSAGE)s"
   if not os.path.exists(success_file):
     f=open(src+"/"+NAME+".tmp.bash", 'w')
     f.write("#!/bin/bash\n")
     f.write("#BEGINING-OF-FILE\n")

     f.write("JOB_NUM=$1\n")
     f.write("sleep 1\n")
     f.write("cd " + exec_dir + "\n")
     COMSTR=re.sub(r"'",r"''", COM)
     f.write("echo '"+str(COMSTR)+"'\n")
     MESSAGE="2"
     TYPE="dbSetStartTime"
     f.write(jobstatus_cmd % locals() + " -n $JOB_NUM")
     f.write("\n   retval=$?\n   if [ $retval -ne 0 ]; then\n     exit 66\n   fi\n")
     f.write("\n\n"+ str(COM) +"\n\n")
     f.write("retval=$?\necho \"[\"$retval\"]\"\nif [ $retval -eq 0 ]; then\n")
     if (str(NAME) != str(SERVICENAME)):
       f.write("touch "+success_file+"\n")
     
     MESSAGE="3"
     TYPE="dbSetEndTime"
     f.write(jobstatus_cmd%locals() + " -n $JOB_NUM")
     f.write("\n    retval=$?\n   if [ $retval -ne 0 ]; then\n     exit 66\n   fi\n")
     f.write("  echo success\nelse\n  echo failed\n")
     MESSAGE="0"
     f.write(jobstatus_cmd%locals() + " -n $JOB_NUM")
     f.write("\n    retval=$?\n   if [ $retval -ne 0 ]; then\n     exit 66\n   fi\n")
     f.write("  exit 127\nfi\ndate\n")
     f.write("#END-OF-FILE\n")
     f.close();
     os.system("chmod 755 "+src+"/"+NAME+".tmp.bash")
     f=open(src+"/"+NAME+".submit.bash", 'w')
     f.write(src + "/" +NAME + ".tmp.bash $1> " + logs + "/$1.std 2>&1")
     f.close();

     os.system("chmod 755 "+src+"/"+NAME+".submit.bash")

     command=src+"/"+NAME+".submit.bash $$"
     print "\n\n\nCOM:"+command+"\n\n\n"
     pid = subprocess.Popen(command, shell=True).pid
     print "\n\n\nPID:"+str(pid)+"\n\n\n"

     MESSAGE="1"
     TYPE="dbSubmitJob"
     jobstatus_cmd = jobstatus_cmd + " -n %(pid)s -c '%(src)s/%(NAME)s.submit.bash'"
     command = jobstatus_cmd % locals()
     #PUT TRY CATCH HERE 
     if pid>0:
        return submitjobs.runcmd(command, logging)
 
if __name__ == "__main__":
    main()
