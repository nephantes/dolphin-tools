*********
Quick Start Guide
*********

dolphin-tools

This guide will walk you through creating new pipelines in dolphin step by step.

Getting Started
===============

dolphin-tools need to be located in the processing units. It can be either a HPC system or a standalone machine. All required tools need to be installed in these systems.
First we are going to create a simple test pipeline. To be able to create a test pipeline and run it in the cluster we need a workflow file.

Installation
============
Please check quick start guide for Docker or for HPC please use installation guide for the appropriate cluster section. 

Workflow file
=============

Workflow file is a text file and includes the directives about the workflow you want to create. This workflow file can reside either in your client machine or in your cluster. 
runWorkflow.py script uses this workflow file to run the steps in your target processing units. runWorkflow script is basically a client script to connect a web-api to run the scripts in a distributed environment.
In addition to this, this script runs in the client system until all workflow ends succcessfully. This workflow file consists of three columns. First column is a name for the step. Second column the command you want to run in your remote machine  and the last column is the time interval that the client checks the run status.
If the run finished sucessfully or killed in the cluster or your host machine, the client senses this and continue or exit accordingly.
There are examples in tests/pipeline folder.

* Ex:
.. code-block:: bash
   w1.txt

   step1	ls -l	1
   step2	echo 1	1
   step3	ls	1

Here we are going to run three steps in the cluster and runWorkflow script will check the steps in every second, if they finished sucessfully or not.

Running a workflow
==================

To run a workflow runWorkflow.py file need to be called with the parameters below.

* Usage: runWorkflow.py [options]

Options:
  -h, --help            show this help message and exit
  -i INPUTPARAM, --inputparam=INPUTPARAM
                        input parameters for the workflow
  -p DEFAULTPARAM, --defaultparam=DEFAULTPARAM
                        defined parameter file that will be run on cluster
  -u USERNAME, --username=USERNAME
                        defined user in the cluster
  -k WKEY, --wkey=WKEY  defined key for the workflow
  -w WORKFLOWFILE, --workflowfile=WORKFLOWFILE
                        workflow filename
  -d DBHOST, --dbhost=DBHOST
                        dbhost name
  -o OUTDIR, --outdir=OUTDIR
                        output directory in the cluster
  -f CONFIG, --config=CONFIG
                        configuration parameter section
                        
please chose your -f option according to your installation. If you are running this on Docker and made your definitions right on your Docker section right. The command should be something like below;

.. code-block:: bash
   python path/to/installation_dir/src/runWorkflow.py -f Docker -w path/to/w1.txt -o /export/TEST2

Standart output of a run in Docker
==================================

If everything is successfull you need to see an output something like below;

.. code-block:: bash
   Docker
   http://localhost/dolphin/api/service.php
   localhost
   /export/tmp/logs
   WORKFLOW STARTED:hiB4l9c3DZNw7YS4ZzjAdmxyyXPzHS
   
   step1:hiB4l9c3DZNw7YS4ZzjAdmxyyXPzHS:ls -l
   "RUNNING(2):python \/usr\/local\/share\/dolphin_tools\/src\/runService.py -f Localhost -d localhost   -o \/export\/TEST -u kucukura -k hiB4l9c3DZNw7YS4ZzjAdmxyyXPzHS -c \"ls -l\" -n step1 -s step1 2>&1"
   RUNNING(2):python /usr/local/share/dolphin_tools/src/runService.py -f Localhost -d localhost   -o /export/TEST -u kucukura -k hiB4l9c3DZNw7YS4ZzjAdmxyyXPzHS -c "ls -l" -n step1 -s step1 2>&1
   
   step1:hiB4l9c3DZNw7YS4ZzjAdmxyyXPzHS:ls -l
   "Service ended successfully!!!"
   Service ended successfully!!!
   
   step2:hiB4l9c3DZNw7YS4ZzjAdmxyyXPzHS:echo 1
   "RUNNING(2):python \/usr\/local\/share\/dolphin_tools\/src\/runService.py -f Localhost -d localhost   -o \/export\/TEST -u kucukura -k hiB4l9c3DZNw7YS4ZzjAdmxyyXPzHS -c \"echo 1\" -n step2 -s step2 2>&1"
   RUNNING(2):python /usr/local/share/dolphin_tools/src/runService.py -f Localhost -d localhost   -o /export/TEST -u kucukura -k hiB4l9c3DZNw7YS4ZzjAdmxyyXPzHS -c "echo 1" -n step2 -s step2 2>&1
   
   step2:hiB4l9c3DZNw7YS4ZzjAdmxyyXPzHS:echo 1
   "Service ended successfully!!!"
   Service ended successfully!!!
   
   step3:hiB4l9c3DZNw7YS4ZzjAdmxyyXPzHS:ls
   "RUNNING(2):python \/usr\/local\/share\/dolphin_tools\/src\/runService.py -f Localhost -d localhost   -o \/export\/TEST -u kucukura -k hiB4l9c3DZNw7YS4ZzjAdmxyyXPzHS -c \"ls\" -n step3 -s step3 2>&1"
   RUNNING(2):python /usr/local/share/dolphin_tools/src/runService.py -f Localhost -d localhost   -o /export/TEST -u kucukura -k hiB4l9c3DZNw7YS4ZzjAdmxyyXPzHS -c "ls" -n step3 -s step3 2>&1
   
   step3:hiB4l9c3DZNw7YS4ZzjAdmxyyXPzHS:ls
   "Service ended successfully!!!"
   Service ended successfully!!!
   
All    the services Ended

The directory structure:
==========================================

* For each step you want to run will be a script under OUTDIR/scripts directory.
* The standard output will be logged under tmp/lsf folder with its PID.std.
There are other log files are about communication with mySQL and LSF logs if you are running them in LSF cluster
* Intermediate submission scripts are in tmp/src folder
* If there are other jobs submitted in the steps, they are going to be tracked under track folder to be able to resumed the jobs. But in this test, there is no such jobs. 

.. code-block:: bash
   
   /export/TEST
   |-- scripts
   |   |-- step1.bash
   |   |-- step2.bash
   |   `-- step3.bash
   `-- tmp
       |-- lsf
       |   |-- 862.jobStatus.log
       |   |-- 862.std
       |   |-- 895.jobStatus.log
       |   |-- 895.std
       |   |-- 927.jobStatus.log
       |   `-- 927.std
       |-- src
       |   |-- step1.submit.bash
       |   |-- step1.tmp.bash
       |   |-- step2.submit.bash
       |   |-- step2.tmp.bash
       |   |-- step3.submit.bash
       |   `-- step3.tmp.bash
       `-- track







