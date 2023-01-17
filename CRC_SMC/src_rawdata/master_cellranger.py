#!/usr/bin/python2

# 2020.07.09; version 2
# local core and memory

import argparse
import os
import commands

def run_cellranger(rawArg, outArg, annoArg, startArg, endArg, coreArg, memArg):
    if not rawArg[-1] == '/': rawArg += '/'
    if not outArg[-1] == '/': outArg += '/'
    
    rawDir = os.listdir(rawArg)
    rawDir.sort()
    if startArg:
        if endArg: rawDir = rawDir[startArg:endArg]
        else: rawDir = rawDir[startArg:]
    else:
        if endArg: rawDir = rawDir[:endArg]
        else: pass

    pwd = commands.getoutput('pwd')
    crgbase = commands.getoutput('which cellranger')
    
    for rawdataname in rawDir:
        print ('### running ' + rawdataname + ' ###')
        os.chdir(outArg)
		
        cmd_crg = crgbase + ' count'
        if coreArg:
            cmd_crg += ' --localcores=' + coreArg
        if memArg:
            cmd_crg += ' --localmem=' + memArg
        cmd_crg += ' --id=' + rawdataname.replace('-', '_')

        samplePrefix = list(set(map(lambda x: '_'.join(x.split('_')[:2]), os.listdir(rawArg + rawdataname))))
        samplePrefix.sort()
        if len(samplePrefix) > 1:
            cmd_crg += ' --sample=' + ','.join(samplePrefix)
        
        cmd_crg += ' --fastqs=' + rawArg + rawdataname + '/'
        cmd_crg += ' --transcriptome=' + annoArg
        print (cmd_crg + '\n')
        commands.getoutput(cmd_crg)
        os.chdir(pwd)


def main(args):
	run_cellranger(args.Raw, args.Out, args.Anno, args.start, args.end, args.core, args.mem)

if __name__ == '__main__':
    parser = argparse.ArgumentParser('')
    parser.add_argument('-r', '--Raw', help = 'Raw data directory (absolute path)', required = True)
    parser.add_argument('-o', '--Out', help = 'Output directory (absolute path)', required = True)
    parser.add_argument('-a', '--Anno', help = 'Annotation meta files (absolute path)', required = True)
    parser.add_argument('--start', help = 'start dir index', required = False, type = int)
    parser.add_argument('--end', help = 'end dir index', required = False, type = int)
    parser.add_argument('--core', help = '# of cores', required = False)
    parser.add_argument('--mem', help = 'use memory', required = False)
    args = parser.parse_args()
    main(args)
