#!/usr/bin/env python
# Submit jobs on CRAB.
# This file is part of https://github.com/hh-italian-group/h-tautau.

import argparse
import sys
import re
from sets import Set

parser = argparse.ArgumentParser(description='Submit jobs on CRAB.',
                  formatter_class = lambda prog: argparse.HelpFormatter(prog,width=90))
parser.add_argument('--workArea', required=False, default="work_area", help="Work area")
parser.add_argument('--cfg', required=True, help="CMSSW configuration file")
parser.add_argument('--scriptExe', required=False, default="", help="Custom script to execute")
parser.add_argument('--site', required=True, help="Site for stage out.")
parser.add_argument('--dryrun', action="store_true", help="Submission dryrun.")
parser.add_argument('--output', required=True, help="output path after /store/user/USERNAME")
parser.add_argument('--blacklist', required=False, default="", help="list of sites where the jobs shouldn't run")
parser.add_argument('--whitelist', required=False, default="", help="list of sites where the jobs can run")
parser.add_argument('--jobNames', required=False, default="",
					help="list of job names to submit (if not specified - submit all)")
parser.add_argument('--lumiMask', required=False, default="",
					help="json file with a lumi mask (default: apply lumi mask from the config file)")
parser.add_argument('--jobNameSuffix', required=False, default="", help="suffix that will be added to each job name")
parser.add_argument('--unitsPerJob', required=False, type=int, default=-1,
					help="number of units per job (default: use values from the config file)")
parser.add_argument('--maxMemory', required=False, type=int, default=2500,
					help="maximum amount of memory (in MB) a job is allowed to use (default: 2500 MB )")
parser.add_argument('--numCores', required=False, type=int, default=1, help="number of cores per job (default: 1)")
parser.add_argument('--allowNonValid', action="store_true", help="Allow nonvalid dataset as an input.")
parser.add_argument('job_file', nargs='+', help="text file with jobs descriptions")
args = parser.parse_args()

from CRABClient.UserUtilities import config, ClientException, getUsernameFromSiteDB
from CRABAPI.RawCommand import crabCommand
from httplib import HTTPException

config = config()

config.General.workArea = args.workArea

config.JobType.pluginName = 'Analysis'
config.JobType.psetName = args.cfg
config.JobType.maxMemoryMB = args.maxMemory
config.JobType.numCores = args.numCores
if len(args.scriptExe) > 0:
    config.JobType.scriptExe = args.scriptExe

config.Data.inputDBS = 'global'
config.Data.allowNonValidInputDataset = args.allowNonValid
config.General.transferOutputs = True
config.General.transferLogs = True
config.Data.publication = False


config.Site.storageSite = args.site
config.Data.outLFNDirBase = "/store/user/{}/{}".format(getUsernameFromSiteDB(), args.output)

if len(args.blacklist) != 0:
	config.Site.blacklist = re.split(',', args.blacklist)
if len(args.whitelist) != 0:
	config.Site.whitelist = re.split(',', args.whitelist)

job_names = Set(filter(lambda s: len(s) != 0, re.split(",", args.jobNames)))

from crab_tools import JobCollection
try:
    for job_file in args.job_file:
        job_collection = JobCollection(job_file, job_names, args.lumiMask, args.jobNameSuffix, args.unitsPerJob)
        print job_file
        print job_collection
        job_collection.submit(config,args.dryrun)
except RuntimeError as err:
    print >> sys.stderr, "ERROR:", str(err)
    sys.exit(1)
