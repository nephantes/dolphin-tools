#!/share/bin/python

from optparse import OptionParser
import warnings

warnings.filterwarnings('ignore', '.*the sets module is deprecated.*',
                        DeprecationWarning, 'MySQLdb')
import MySQLdb
import sys
   
def runSQL(dbhostname, sql): 
    db = MySQLdb.connect(
      host = 'ghpcc06', 
      user = 'biocore', 
      passwd = 'biocore2013', 
      db = 'biocore', 
      port = 33066)
    #db = MySQLdb.connect(dbhostname,"biocore","biocore2013","biocore" )
    try:
        cursor = db.cursor()

        #print sql
        cursor.execute(sql)
        results = cursor.fetchall()
        cursor.close()
        del cursor

    except Exception, err:
        print "DB Error:for help use --help"
        db.rollback()
        sys.stderr.write('ERROR: %s\n' % str(err))
        sys.exit(2)
    finally: 
        db.commit()
    	db.close()
    return results

def getWorkflowId(dbhostname, username, wkey):
    try: 
        sql = "select `workflow_id` from `workflow_run` where `username`='" + str(username) + "' and  `wkey`='" + str(wkey) + "'"
        results = runSQL(dbhostname, sql)
        row     = results[0]
        id      = row[0]
    except:
        id      = 1
    return id

def getId(dbhostname, name, val, username):
    try:
        sql = "select `" + str(name) + "_id` from `"+ str(name) + "s` where `" + str(name) + "name`='" + str(val) + "' and `username`='" + str(username) + "'"
        results = runSQL(dbhostname, sql)
        row     = results[0]
        id      = row[0]
    except:
        id      = 1
    return id

def checkAllJobsFinished(dbhostname, username, wkey, workflow_id, service_id):
    select  = "select count(job_id) c from jobs "
    where1 = " where `username`= '" + str(username) + "' and `wkey`='" + str(wkey) + "' and `workflow_id`=" + str(workflow_id) + " and `service_id`=" + str(service_id)
    where2 = " and `result`=3 "
    sql = "select s1.c, s2.c from (" + select + where1 + ") s1,  (" + select + where1 + where2 + ") s2"
    #f = open("~/scratch/sonuc.txt", "w")
    #f = open("/home/svcgalaxy/scratch/sonuc.txt", "w")
    #f.write(sql+"\n")
    results = runSQL(dbhostname, sql)
    row     = results[0]
    s1     = row[0]
    s2     = row[1]

    if(s1==s2):
    #     f.write("EQUAL\n")
        updateService(dbhostname, wkey, service_id, 1)
    #    f.write(str(sql))
 
 
def updateService(dbhostname,  wkey, service_id, result):
   sql = "update service_run set `end_time`=now(), `result`=" + str(result) + " where `wkey`='" + str(wkey) + "' and `service_id`=" + str(service_id) 
   results = runSQL(dbhostname, sql)

def insertJob(dbhostname, username, wkey, com, jobname, workflow_id, service_id, jobnum, result):
   sql = "insert into jobs(`username`, `wkey`, `run_script`, `jobname`, `workflow_id`, `service_id`, `result`, `submit_time`, `job_num`) values ('" + str(username) + "','" + str(wkey) + "','" + str(com) + "','" + str(jobname) + "'," +str(workflow_id) + "," + str(service_id) + ", '"+ result+"', now(), "+ jobnum +")"
   results = runSQL(dbhostname, sql)

def updateJob(dbhostname, username, wkey, jobname, workflow_id, service_id, field, jobnum, result):
   sql = "update jobs set `" + str(field) + "`=now(), `result`=" + str(result) + " where `username`= '" + str(username) + "' and `wkey`='" + str(wkey) + "' and `jobname`='" + str(jobname) + "' and `workflow_id`=" + str(workflow_id) + " and `service_id`=" + str(service_id) + " and `job_num`=" + str(jobnum)
   results = runSQL(dbhostname, sql)

       
def main():
   try:
	#python finishJob.py -u kucukura -k nVy1THnthvrRWfXcj187KeDDQrNAkY -s splitFastQ
        parser = OptionParser()
        parser.add_option('-u', '--username', help='defined user in the cluster', dest='username')
        parser.add_option('-d', '--dbhostname', help='defined dbhostname for the db', dest='dbhostname')
        parser.add_option('-k', '--key', help='defined key for the workflow', dest='wkey')
        parser.add_option('-c', '--com', help='bash script of the command', dest='com')
        parser.add_option('-j', '--jobname', help='name of of the job', dest='jobname')
        parser.add_option('-s', '--servicename', help='service name', dest='servicename')
	parser.add_option('-t', '--type', help='type of the operation', dest='type')
	parser.add_option('-n', '--jobnum', help='submitted job number', dest='jobnum')
	parser.add_option('-m', '--message', help='resulting message of the job', dest='message')

        (options, args) = parser.parse_args()
   except:
	#parser.print_help()
        print "for help use --help"
        sys.exit(2)


   USERNAME                = options.username   
   DBHOSTNAME              = options.dbhostname
   WKEY                    = options.wkey 
   COM                     = options.com
   JOBNAME                 = options.jobname
   SERVICENAME             = options.servicename
   TYPE			   = options.type
   JOBNUM                  = options.jobnum
   MESSAGE                 = options.message

   if (DBHOSTNAME == None):
        DBHOSTNAME="galaxy.umassmed.edu"
        
   workflow_id = getWorkflowId(DBHOSTNAME, USERNAME, WKEY)
   service_id = getId(DBHOSTNAME, "service", SERVICENAME, USERNAME)
   
   #print str(workflow_id) + ":" + str(service_id)
   if (TYPE == "dbSubmitJob"):
	insertJob(DBHOSTNAME, USERNAME, WKEY, COM, JOBNAME, workflow_id, service_id, JOBNUM, MESSAGE) # MESSAGE has jobnum here
   elif (TYPE == "dbSetStartTime"):
        updateJob(DBHOSTNAME, USERNAME, WKEY, JOBNAME, workflow_id, service_id, "start_time", JOBNUM, MESSAGE)  	
   elif (TYPE == "dbSetEndTime"):
        updateJob(DBHOSTNAME, USERNAME, WKEY, JOBNAME, workflow_id, service_id, "end_time", JOBNUM, MESSAGE) 	

   checkAllJobsFinished(DBHOSTNAME, USERNAME, WKEY, workflow_id, service_id)

if __name__ == "__main__":
    main()


