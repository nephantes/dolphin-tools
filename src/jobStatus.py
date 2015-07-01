#!/usr/bin/python
import logging
from optparse import OptionParser
from ZSI.client import NamedParamBinding as NPBinding, Binding
import json
import os
import re
import cgi
import warnings
import sys
import time
   
url="http://biocore.umassmed.edu/pipeline/service.php"

re_string = re.compile(r'(?P<htmlchars>[<&>])|(?P<space>^[ \t]+)|(?P<lineend>\r\n|\r|\n)|(?P<protocal>(^|\s)((http|ftp)://.*?))(\s|$)', re.S|re.M|re.I)

def plaintext2html(text, tabstop=4):
    def do_sub(m):
        c = m.groupdict()
        if c['htmlchars']:
            return cgi.escape(c['htmlchars'])
        if c['lineend']:
            return '<br>'
        elif c['space']:
            t = m.group().replace('\t', '&nbsp;'*tabstop)
            t = t.replace(' ', '&nbsp;')
            return t
        elif c['space'] == '\t':
            return ' '*tabstop;
        else:
            url = m.group('protocal')
            if url.startswith(' '):
                prefix = ' '
                url = url[1:]
            else:
                prefix = ''
            last = m.groups()[-1]
            if last in ['\n', '\r', '\r\n']:
                last = '<br>'
            return '%s<a href="%s">%s</a>%s' % (prefix, url, url, last)
    return re.sub(re_string, do_sub, text)



def checkAllJobsFinished(username, wkey, service_name, logging):
    kw = {'url':url}
    b = NPBinding(**kw)
    trials=0
    while trials<5:
       try:
          mesg=b.checkAllJobsFinished(a=username , c=wkey , b=service_name)
          trials=10
       except:
          print "Couldn't connect to dolphin server (%s)"%trials
          logging.info("Couldn't connect to dolphin server (%s)"%trials)
          time.sleep(15)
       trials=trials+1
    
    data=json.dumps(mesg)
    wkey=json.loads(data)

    ret=str(wkey['return'])
    #print "JOBs checked:"+ret+"\n"
    if (ret.startswith("ERROR")):
          logging.info("%s:%s"%(jobname,ret))
          logging.info("Check the job#:%s"%jobnum)
          print service_name + ":" + ret + "\n"
          print "Check the service:"+service_name+"\n"
          sys.exit(2);

def insertJob( username, wkey, com, jobname, service_name, jobnum, result, logging):
    
    #print "username="+username+", wkey="+wkey+", com="+com+", jobname="+jobname+", ser="+service_name+", jobnum="+jobnum+", res="+result 
    #kw = {'url':url, 'tracefile':sys.stdout}
    kw = {'url':url}
    b = NPBinding(**kw)
    trials=0
    while trials<5:
       try:
          mesg=b.insertJob(a=username , c=wkey , b=com , e=jobname, d=service_name , f=jobnum, h=result)
          trials=10
       except:
          print "Couldn't connect to dolphin server (%s)"%trials
          logging.info("Couldn't connect to dolphin server (%s)"%trials)
          time.sleep(15)
       trials=trials+1
    
    data=json.dumps(mesg)
    wkey=json.loads(data)

    ret=str(wkey['return'])
    #print "JOB inserted:"+ret+"\n"
    if (ret.startswith("ERROR")):
          logging.info("%s:%s"%(jobname,ret))
          logging.info("Check the job#:%s"%jobnum)
          print jobname + ":" + ret + "\n"
          print "Check the job#:"+jobnum+"\n"
          sys.exit(2);

def updateJob( username, wkey, jobname, service_name, field, jobnum, result, logging):
    #print "username="+str(username)+", wkey="+str(wkey)+", jobname="+str(jobname)+", ser="+str(service_name)+", field="+str(field)+", jobnum="+str(jobnum)+", res="+str(result) 
    kw = {'url':url}
    b = NPBinding(**kw)
    trials=0
    while trials<5:
       try:
          mesg=b.updateJob(a=username , c=wkey , b=jobname, e=service_name , d=field, g=jobnum, f=result)
          trials=10
       except:
          print "Couldn't connect to dolphin server (%s)"%trials
          logging.info("Couldn't connect to dolphin server (%s)"%trials)
          time.sleep(15)         
       trials=trials+1
    logging.info(mesg)
    data=json.dumps(mesg)
    wkey=json.loads(data)

    ret=str(wkey['return'])
    #print "JOB updated:"+ret+"\n"
    if (ret.startswith("ERROR")):
          logging.info("%s:%s"%(jobname,ret))
          logging.info("Check the job#:%s"%jobnum)
          print jobname + ":" + ret + "\n"
          print "Check the job#:"+jobnum+"\n"
          sys.exit(2);

def insertJobOut(username, wkey, jobnum, outdir, edir, logging):
    file=str(outdir)+"/tmp/lsf/"+str(jobnum)
    if os.path.isfile(file) and os.access(file, os.R_OK):
        command="python " + edir  + "/readLSFout.py -f "+file
        child = os.popen(command)
        jobout = child.read().rstrip()
        err = child.close()
     
        if err:
             logging.info('ERROR: %s failed w/ exit code %d' % (command, err))
             raise RuntimeError, 'ERROR: %s failed w/ exit code %d' % (command, err) 
        jobout=plaintext2html(re.sub('\'', "\"", jobout))
        kw = {'url':url}
        b = NPBinding(**kw)
     
        trials=0
        while trials<5:
            try:
               mesg=b.insertJobOut(a=username , c=wkey , b=jobnum, f=jobout)
               trials=10
            except:
               print "Couldn't connect to dolphin server (%s)"%trials
               logging.info("Couldn't connect to dolphin server (%s)"%trials)
               time.sleep(15)         
            trials=trials+1
       
def main():
   try:
        #python finishJob.py -u kucukura -k nVy1THnthvrRWfXcj187KeDDQrNAkY -s splitFastQ
        parser = OptionParser()
        parser.add_option('-u', '--username', help='defined user in the cluster', dest='username')
        parser.add_option('-k', '--key', help='defined key for the workflow', dest='wkey')
        parser.add_option('-c', '--com', help='bash script of the command', dest='com')
        parser.add_option('-j', '--jobname', help='name of of the job', dest='jobname')
        parser.add_option('-s', '--servicename', help='service name', dest='servicename')
        parser.add_option('-t', '--type', help='type of the operation', dest='type')
        parser.add_option('-n', '--jobnum', help='submitted job number', dest='jobnum')
        parser.add_option('-m', '--message', help='resulting message of the job', dest='message')
        parser.add_option('-o', '--outdir', help='output directory', dest='outdir')



        (options, args) = parser.parse_args()
   except:
        #parser.print_help()
        print "for help use --help"
        sys.exit(2)


   USERNAME                = options.username   
   WKEY                    = options.wkey 
   COM                     = options.com
   JOBNAME                 = options.jobname
   SERVICENAME             = options.servicename
   TYPE			           = options.type
   JOBNUM                  = options.jobnum
   MESSAGE                 = options.message
   OUTDIR                  = options.outdir
        
   edir=os.path.dirname(sys.argv[0])
   print OUTDIR
   print TYPE
   logfile="%s/tmp/lsf/%s.jobStatus.log"%(OUTDIR,str(JOBNUM))
   logging.basicConfig(filename=logfile, filemode='a',format='%(asctime)s %(message)s', datefmt='%m/%d/%Y %I:%M:%S %p', level=logging.DEBUG)
   logging.info("File Path:%s"%os.getcwd())
   logging.info("WKEY:"+str(WKEY))
   logging.info("TYPE:"+str(TYPE))

   if (TYPE == "dbSubmitJob"):
        insertJob(USERNAME, WKEY, COM, JOBNAME, SERVICENAME, JOBNUM, MESSAGE, logging) # MESSAGE has jobnum here
   elif (TYPE == "dbSetStartTime"):
        updateJob(USERNAME, WKEY, JOBNAME, SERVICENAME, "start_time", JOBNUM, MESSAGE, logging)  	
   elif (TYPE == "dbSetEndTime"):
        updateJob(USERNAME, WKEY, JOBNAME, SERVICENAME, "end_time", JOBNUM, MESSAGE, logging) 	
        insertJobOut(USERNAME, WKEY, JOBNUM, OUTDIR, edir, logging)

   checkAllJobsFinished(USERNAME, WKEY, SERVICENAME, logging)

if __name__ == "__main__":
    main()


