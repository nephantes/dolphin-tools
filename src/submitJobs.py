#!/bin/env python
        
from optparse import OptionParser
import sys
import os
import re
import string
import subprocess

sys.path.insert(0, sys.path[0])
from config import *

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
 	(options, args) = parser.parse_args()
   except:
        print "OptionParser Error:for help use --help"
        sys.exit(2)
   USERNAME    = options.username
   DBHOST      = options.dbhostname
   WKEY        = options.wkey 
   OUTDIR      = options.outdir 
   SERVICENAME = options.servicename
   COM         = options.com
   NAME        = options.name
   CPU         = options.cpu
   TIME        = options.time
   MEMORY      = options.memory
   QUEUE       = options.queue
   python      = "python"

   config=getConfig()

   com="module list 2>&1 |grep python/2.7.5"
   pythonload=str(os.popen(com).readline().rstrip())
   if (len(pythonload)<5):
       com     = "module load python/2.7.5;";
       pythonload=str(os.popen(com).readline().rstrip())
       python      = "module load python/2.7.5;" + python;

   com="module list 2>&1 |grep openssl/1.0.1g"
   sslload=str(os.popen(com).readline().rstrip())
   if (len(sslload)<5):
       com = "module load openssl/1.0.1g;"
       sslload=str(os.popen(com).readline().rstrip())
   
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
        TIME="240";

   if (QUEUE == None):
        queue="-q short"
   else: 
        queue="-q "+str(QUEUE)

   if (MEMORY == None):
        MEMORY="4096";
  
   COM.replace('\"{','\'{')
   COM.replace('}\"','}\'')
   print "COMMAND: [" + COM + "]\n"
   print "NAME: [" + NAME + "]\n"
   print "cpu: [" + CPU + "]\n"


   exec_dir=os.path.dirname(os.path.abspath(__file__))
   print "EXECDIR" + exec_dir
   sdir=config['tooldir']+"/src"
   track=OUTDIR + "/tmp/track"
   src=OUTDIR + "/tmp/src"
   lsf=OUTDIR + "/tmp/lsf"
  
   os.system("mkdir -p "+track)
   os.system("mkdir -p "+src)
   os.system("mkdir -p "+lsf)
   success_file = track+"/"+str(NAME)+".success";
   if not os.path.exists(success_file):
     f=open(src+"/"+NAME+".tmp.bash", 'w')
     f.write("#!/bin/bash\n")
     f.write("#BEGINING-OF-FILE\n")
     f.write("a=$( module list 2>&1)\n")
     f.write("echo $a\n")

     f.write("if [[ $a != *python/2.7.5* ]];\n")
     f.write("then\n")
     f.write("  module load python/2.7.5\n")
     f.write("fi\n")

     f.write("if [[ $a != *openssl/1.0.1g* ]];\n")
     f.write("then\n")
     f.write("  module load openssl/1.0.1g\n")
     f.write("fi\n")

     f.write("cd " + exec_dir + "\n")
     #f.write("echo '"+str(COM)+"'\n")
     f.write("python " + sdir + "/jobStatus.py -u " + str(USERNAME) + " -k " + str(WKEY) + " -s " + str(SERVICENAME) + " -t dbSetStartTime -o "+str(OUTDIR)+" -n $LSB_JOBID -j "+ str(NAME)+ " -m 2\n")
     f.write("echo \"python " + sdir + "/jobStatus.py -u " + str(USERNAME) + " -k " + str(WKEY) + " -s " + str(SERVICENAME) + " -t dbSetStartTime -o "+str(OUTDIR)+" -n $LSB_JOBID -j "+ str(NAME)+ " -m 2\"\n")
     f.write("   retval=$?\n   if [ $retval -ne 0 ]; then\n     exit 66\n   fi\n")
     f.write("\n\n"+ str(COM) +"\n\n")
     f.write("retval=$?\necho \"[\"$retval\"]\"\nif [ $retval -eq 0 ]; then\n")
     if (str(NAME) != str(SERVICENAME)):
       f.write("touch "+success_file+"\n")
     f.write("echo \"python " + sdir + "/jobStatus.py -u " + str(USERNAME) + " -k " + str(WKEY) + " -s " + str(SERVICENAME) + " -t dbSetEndTime -o "+str(OUTDIR)+" -n $LSB_JOBID -j "+ str(NAME)+ " -m 3\"\n")
     f.write("  python " + sdir + "/jobStatus.py -u " + str(USERNAME) + " -k " + str(WKEY) + " -s " + str(SERVICENAME) + " -t dbSetEndTime -o "+str(OUTDIR)+" -n $LSB_JOBID -j "+ str(NAME)+ " -m 3\n")
     f.write("    retval=$?\n   if [ $retval -ne 0 ]; then\n     exit 66\n   fi\n")
     f.write("  echo success\nelse\n  echo failed\n")
     f.write("  python " + sdir + "/jobStatus.py -u " + str(USERNAME) + " -k " + str(WKEY) + " -s " + str(SERVICENAME) + " -t dbSetEndTime -o "+str(OUTDIR)+" -n $LSB_JOBID -j "+ str(NAME)+ " -m 0\n")
     f.write("    retval=$?\n   if [ $retval -ne 0 ]; then\n     exit 66\n   fi\n")
     f.write("  exit 127\nfi\ndate\n")

     f.write("#END-OF-FILE\n")
     f.close();
     os.system("chmod 755 "+src+"/"+NAME+".tmp.bash")
     f=open(src+"/"+NAME+".submit.bash", 'w')
     f.write(src + "/" +NAME + ".tmp.bash > " + lsf + "/$LSB_JOBID.std 2>&1")
     #f.write(src + "/" +NAME + ".tmp.bash 1> " + lsf + "/dene.std 2> >(tee -a " + lsf + "/dene.std >&2)")
     f.close();

     os.system("chmod 755 "+src+"/"+NAME+".submit.bash")

     f=open(src+"/"+NAME+".submit.log", 'w')
     #f.write("#!/bin/bash\n")
     
     #f.write("NODE:"+ str(nodename) + " part[0]:" + part[0] + "\n")
     #CHANGE this submition script according to the system.
     #PUT TRY CATCH HERE 
     command="bsub "+queue+" -m blades -P dolphin -R \"span[hosts=1]\" -n "+str(CPU)+" -W "+str(TIME)+" -R \"rusage[mem="+str(MEMORY)+"]\" -J "+NAME+" -o "+lsf+" < "+src+"/"+NAME+".submit.bash"
     print command
     f.write("SUBMIT SCRIPT[" + command +"]\n\n")
     #while True:
     output = runcmd(command)
     #print str(output)+"\n"
     f.write("SUBMIT OUT:[" + str(output) + "]\n")
     words = re.split('[\<\>]+', str(output))
     num = words[1]
     #     if num>0:
     #         time.sleep( 5 )
     #         break
     
     
     command = python+" " + sdir + "/jobStatus.py -u " + str(USERNAME) + " -k " + str(WKEY) + " -s " + str(SERVICENAME) + " -t dbSubmitJob -o "+str(OUTDIR)+" -n "+ str(num) + " -j "+ str(NAME) + " -m 1 -c \"" + src+"/"+NAME+".submit.bash\"" 
     #print command
     f.write("RUN COMMAND:\n" + str(command) + "\n")
     #f.write("NUM:[" + str(num) + "]\n")
     #PUT TRY CATCH HERE 
     if num>0:
        return runcmd(command)
     f.close()
 
if __name__ == "__main__":
    main()
