#
# Copyright 2017 National Renewable Energy Laboratory
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
#

"""
    This library provides tools for executing cases with drivers contained in the
    OpenFAST framework. Any new drivers should have a corresponding public driver
    function called `def run[NewDriver]Case` in this library.
"""

import argparse
import os
import sys
import shutil
import subprocess
import rtestlib as rtl

def _runCase(executable, inputFile, logFile, stdout, restart=False, ExtraFlags=""):
    if logFile is None:
        command = f"{executable} {inputFile} {ExtraFlags}"
    elif restart:
        command = f"{executable} -restart {os.path.splitext(inputFile)[0]} > {logFile}"
    else:
        command = f"{executable} {inputFile} {ExtraFlags} > {logFile}"
    print(command)
    return subprocess.call(command, stdout=stdout, shell=True)
    
def _runGenericCase(inputFile, executable, verbose=False, restart=False, ExtraFlags=""):
    stdout = sys.stdout if verbose else open(os.devnull, 'w')
    
    rtl.validateFileOrExit(inputFile)
    rtl.validateExeOrExit(executable)
    
    if verbose:
        logFile = None
    else:
        caseparent = os.path.sep.join(inputFile.split(os.path.sep)[:-1])
        casebase = caseparent.split(os.path.sep)[-1]  # assumes that the directory structure name is the name of the .log file. (for consistent driver + glue-code names)
        if restart:
            logFile = caseparent + os.path.sep + casebase + '_2.log'
        else:
            logFile = caseparent + os.path.sep + casebase + '.log'
    
    returnCode = _runCase(executable, inputFile, logFile, stdout, restart, ExtraFlags)
    print("COMPLETE with code {}".format(returnCode), flush=True)    
    
    return returnCode

def _runUACase(inputFile, executable, verbose=False):
    stdout = sys.stdout if verbose else open(os.devnull, 'w')
    
    rtl.validateFileOrExit(inputFile)
    rtl.validateExeOrExit(executable)

    if verbose:
        logFile = None
    else:
        logFile = os.path.splitext(inputFile)[0]+'.log'
    
    returnCode = _runCase(executable, inputFile, logFile, stdout)
    print("COMPLETE with code {}".format(returnCode), flush=True)    
    
    return returnCode


def runOpenfastCase(inputFile, executable, verbose=False, restart=False):
    return _runGenericCase(inputFile, executable, verbose, restart)

def runAeromapCase(inputFile, executable, verbose=False):
    return _runGenericCase(inputFile, executable, verbose, restart=False, ExtraFlags="-steadystate")

def runAerodynDriverCase(inputFile, executable, verbose=False):
    caseDirectory = os.path.sep.join(inputFile.split(os.path.sep)[:-1])
    os.chdir(caseDirectory)
    return _runGenericCase(inputFile, executable, verbose)

def runUnsteadyAeroDriverCase(inputFile, executable, verbose=False):
    caseDirectory = os.path.dirname(inputFile)
    os.chdir(caseDirectory)
    return _runUACase(inputFile, executable, verbose)

def runBeamdynDriverCase(inputFile, executable, verbose=False):
    caseDirectory = os.path.sep.join(inputFile.split(os.path.sep)[:-1])
    os.chdir(caseDirectory)
    return _runGenericCase(inputFile, executable, verbose)

def runHydrodynDriverCase(inputFile, executable, verbose=False):
    caseDirectory = os.path.sep.join(inputFile.split(os.path.sep)[:-1])
    os.chdir(caseDirectory)
    return _runGenericCase(inputFile, executable, verbose)

def runMoordynDriverCase(inputFile, executable, verbose=False):
    print("MoorDyn run: {}".format(inputFile))
    return _runGenericCase(inputFile, executable, verbose)

def runSubdynDriverCase(inputFile, executable, verbose=False):
    caseDirectory = os.path.sep.join(inputFile.split(os.path.sep)[:-1])
    os.chdir(caseDirectory)
    return _runGenericCase(inputFile, executable, verbose)

def runInflowwindDriverCase(inputFile, executable, verbose=False):
    caseDirectory = os.path.sep.join(inputFile.split(os.path.sep)[:-1])
    os.chdir(caseDirectory)
    return _runGenericCase(inputFile, executable, verbose)

def runMoordynDriverCase(inputFile, executable, verbose=False):
    caseDirectory = os.path.sep.join(inputFile.split(os.path.sep)[:-1])
    os.chdir(caseDirectory)
    return _runGenericCase(inputFile, executable, verbose)

def runSeaStateDriverCase(inputFile, executable, verbose=False):
    caseDirectory = os.path.sep.join(inputFile.split(os.path.sep)[:-1])
    os.chdir(caseDirectory)
    return _runGenericCase(inputFile, executable, verbose)

def runAerodiskDriverCase(inputFile, executable, verbose=False):
    caseDirectory = os.path.sep.join(inputFile.split(os.path.sep)[:-1])
    os.chdir(caseDirectory)
    return _runGenericCase(inputFile, executable, verbose)

def runSimpleElastodynDriverCase(inputFile, executable, verbose=False):
    caseDirectory = os.path.sep.join(inputFile.split(os.path.sep)[:-1])
    os.chdir(caseDirectory)
    return _runGenericCase(inputFile, executable, verbose)
