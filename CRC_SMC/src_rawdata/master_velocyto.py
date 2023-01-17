#!/usr/bin/python2

# 2022.12.04; version 1

import argparse
import os
import commands


def run_cellranger(processedArg, annoArg, startArg, endArg, coreArg, memArg):
    if not processedArg[-1] == '/': processedArg += '/'
    
    rawDir = os.listdir(processedArg)
    rawDir.sort()
    if startArg:
        if endArg: rawDir = rawDir[startArg:endArg]
        else: rawDir = rawDir[startArg:]
    else:
        if endArg: rawDir = rawDir[:endArg]
        else: pass

    veloBase = '/home/sangho/miniconda3/envs/velocyto/bin/velocyto '
    veloBase += 'run10x '
    if coreArg:
        veloBase += '-@ ' + coreArg + ' '
    else:
        veloBase += '-@ 8 '
    
    for rawdataname in rawDir:
        print ('\n### running ' + rawdataname)
        
        cmdVelo = veloBase + processedArg + rawdataname + '/ '
        cmdVelo += annoArg
        		
        print (cmdVelo + '\n')
        commands.getoutput(cmdVelo)


def main(args):
	run_cellranger(args.Processed, args.Anno, args.start, args.end, args.core, args.mem)

if __name__ == '__main__':
    parser = argparse.ArgumentParser('')
    parser.add_argument('-p', '--Processed', help = 'run10x; cellranger output directory (absolute path)', required = True)
    parser.add_argument('-a', '--Anno', help = 'Annotation meta files (absolute path)', required = True)
    parser.add_argument('--start', help = 'start dir index', required = False, type = int)
    parser.add_argument('--end', help = 'end dir index', required = False, type = int)
    parser.add_argument('--core', help = '# of cores; (*default is 8)', required = False)
    parser.add_argument('--mem', help = 'use memory', required = False)
    args = parser.parse_args()
    main(args)
