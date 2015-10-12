#!/bin/env python
        
from optparse import OptionParser
import sys
import os
import re
import string
import subprocess
#import time

#def getKey(num):
#       return ''.join(random.choice(string.ascii_letters) for x in range(num))
   
def runcmd(command): 
    print command
    child = os.popen(command)
    data = child.read()
    print data
    err = child.close()
    if err:
        return 'ERROR: %s failed w/ exit code %d' % (command, err)
    return data

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
        parser.add_option('-t', '--time', help='time', dest='time')
        parser.add_option('-m', '--memory', help='memory', dest='memory')
        parser.add_option('-o', '--outdir', help='output directory', dest='outdir')
        parser.add_option('-q', '--queue', help='queue', dest='queue')
        parser.add_option('-f', '--config', help='configuration parameter section', dest='config')
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
   TIME        = options.time
   MEMORY      = options.memory
   QUEUE       = options.queue
   CONFIG      = options.config
   python      = "python"

   if (USERNAME==None):
        USERNAME=subprocess.check_output("whoami", shell=True).rstrip()
   
   print "USER:"+str(USERNAME)+"\n";

   if (NAME == None):
        NAME="job";
   if (OUTDIR == None):
        OUTDIR="~/out";
   if (CPU == None):
        CPU="1";
   if (TIME == None):
        TIME="600";
   if (QUEUE == None):
        queue="-q short"
   else: 
        queue="-q "+str(QUEUE)
   if (MEMORY == None):
        MEMORY="4096";
   if (DBHOSTNAME == None):
        DBHOSTNAME="localhost"
   sdir=os.path.dirname(sys.argv[0])
   exec_dir=os.path.dirname(os.path.abspath(__file__))
   #print "EXECDIR" + exec_dir
   track=OUTDIR + "/tmp/track"
   src=OUTDIR + "/tmp/src"
   lsf=OUTDIR + "/tmp/lsf"
  
   os.system("mkdir -p "+track)
   os.system("mkdir -p "+src)
   os.system("mkdir -p "+lsf)
   success_file = track+"/"+str(NAME)+".success";
   jobstatus_cmd = "python %(sdir)s/jobStatus.py -f %(CONFIG)s -u %(USERNAME)s -k %(WKEY)s -s %(SERVICENAME)s -t %(TYPE)s -o %(OUTDIR)s -j %(NAME)s -m %(MESSAGE)s"
   if not os.path.exists(success_file):
     f=open(src+"/"+NAME+".tmp.bash", 'w')
     f.write("#!/bin/bash\n")
     f.write("#BEGINING-OF-FILE\n")

     f.write("JOB_NUM=$1\n")
     f.write("sleep 1\n")
     f.write("cd " + exec_dir + "\n")
     f.write("echo '"+str(COM)+"'\n")
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
     f.write(src + "/" +NAME + ".tmp.bash $1> " + lsf + "/$1.std 2>&1")
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
        return runcmd(command)
 
if __name__ == "__main__":
    main()
