#!/bin/csh -f
#  TestLB.cmd
#
#  UGE job for TestLB built Sat Nov 12 13:41:20 PST 2016
#
#  The following items pertain to this script
#  Use current working directory
#$ -cwd
#  input           = /dev/null
#  output          = /u/home/s/sanjaydm/voom2/src/Applications/Landau-Brazovskii/gauss/eps_mp08/TestLB.joblog.$JOB_ID
#$ -o /u/home/s/sanjaydm/voom2/src/Applications/Landau-Brazovskii/gauss/eps_mp08/TestLB.joblog.$JOB_ID
#  error           = Merged with joblog
#$ -j y
#  The following items pertain to the user program
#  user program    = /u/home/s/sanjaydm/voom2/src/Applications/Landau-Brazovskii/gauss/eps_mp08/TestLB
#  arguments       = InputGauss.LB
#  program input   = Specified by user program
#  program output  = Specified by user program
#  Resources requested
#
#$ -l h_data=1024M,h_rt=8:00:00
#
#  Name of application for log
#$ -v QQAPP=job
#  Email address to notify
#$ -M sanjaydm@mail
#  Notify at beginning and end of job
#$ -m bea
#  Job is not rerunable
#$ -r n
#
# Initialization for serial execution
#
  unalias *
  set qqversion = 
  set qqapp     = "job serial"
  set qqidir    = /u/home/s/sanjaydm/voom2/src/Applications/Landau-Brazovskii/gauss/eps_mp08
  set qqjob     = TestLB
  set qqodir    = /u/home/s/sanjaydm/voom2/src/Applications/Landau-Brazovskii/gauss/eps_mp08
  cd     /u/home/s/sanjaydm/voom2/src/Applications/Landau-Brazovskii/gauss/eps_mp08
  source /u/local/bin/qq.sge/qr.runtime
  if ($status != 0) exit (1)
#
  echo "UGE job for TestLB built Sat Nov 12 13:41:20 PST 2016"
  echo ""
  echo "  TestLB directory:"
  echo "    "/u/home/s/sanjaydm/voom2/src/Applications/Landau-Brazovskii/gauss/eps_mp08
  echo "  Submitted to UGE:"
  echo "    "$qqsubmit
  echo "  SCRATCH directory:"
  echo "    "$qqscratch
#
  echo ""
  echo "TestLB started on:   "` hostname -s `
  echo "TestLB started at:   "` date `
  echo ""
#
# Run the user program
#
  source /u/local/Modules/default/init/modules.csh
  module load intel/13.cs
#
  echo TestLB "InputGauss.LB" \>\& TestLB.output.$JOB_ID
  echo ""
  time /u/home/s/sanjaydm/voom2/src/Applications/Landau-Brazovskii/gauss/eps_mp08/TestLB InputGauss.LB >& /u/home/s/sanjaydm/voom2/src/Applications/Landau-Brazovskii/gauss/eps_mp08/TestLB.output.$JOB_ID
#
  echo ""
  echo "TestLB finished at:  "` date `
#
# Cleanup after serial execution
#
  source /u/local/bin/qq.sge/qr.runtime
#
  echo "-------- /u/home/s/sanjaydm/voom2/src/Applications/Landau-Brazovskii/gauss/eps_mp08/TestLB.joblog.$JOB_ID --------" >> /u/local/apps/queue.logs/job.log.serial
  if (`wc -l /u/home/s/sanjaydm/voom2/src/Applications/Landau-Brazovskii/gauss/eps_mp08/TestLB.joblog.$JOB_ID  | awk '{print $1}'` >= 1000) then
	head -50 /u/home/s/sanjaydm/voom2/src/Applications/Landau-Brazovskii/gauss/eps_mp08/TestLB.joblog.$JOB_ID >> /u/local/apps/queue.logs/job.log.serial
	echo " "  >> /u/local/apps/queue.logs/job.log.serial
	tail -10 /u/home/s/sanjaydm/voom2/src/Applications/Landau-Brazovskii/gauss/eps_mp08/TestLB.joblog.$JOB_ID >> /u/local/apps/queue.logs/job.log.serial
  else
	cat /u/home/s/sanjaydm/voom2/src/Applications/Landau-Brazovskii/gauss/eps_mp08/TestLB.joblog.$JOB_ID >> /u/local/apps/queue.logs/job.log.serial
  endif
  exit (0)
