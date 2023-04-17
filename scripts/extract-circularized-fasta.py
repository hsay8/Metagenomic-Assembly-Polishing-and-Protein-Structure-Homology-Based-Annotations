#!/usr/bin/python3.6
import sys

if len(sys.argv) == 6:
    circ_contigs=[]
    contoanalyze = open(sys.argv[3], "w")
    infofh = open(sys.argv[1])
    for line in infofh:
        if "#" not in line:
            info = line.split("\t")
            if (info[3] == "Y" or info[3] == "+") and int(info[1]) > int(sys.argv[5]):
                contoanalyze.write(info[0] + " \n")
                circ_contigs.append(info[0])
    infofh.close()
    contoanalyze.close()
    # Write fasta entry to new file if circularized
    writeflag=False
    with open(sys.argv[2]) as allfh:
        for line in allfh:
            if ">" in line:
                writeflag=False
            if line[1:].strip() in circ_contigs:
                writeflag=True
                filename = line[1:].strip()
            if writeflag:
                temp = sys.argv[4]
                with open("%s/%s" %(temp, filename), "a") as circfh:
                    circfh.write(line)

elif len(sys.argv) == 7:
    circ_contigs=[]
    contoanalyze = open(sys.argv[3], "w")
    infofh = open(sys.argv[1])
    for line in infofh:
        if "#" not in line:
            info = line.split("\t")
            if info[3] == "+" and int(info[1]) <= int(sys.argv[5]):
                contoanalyze.write(info[0] + " \n")
                circ_contigs.append(info[0])
    infofh.close()
    contoanalyze.close()
    # Write fasta entry to new file if circularized
    writeflag=False
    with open(sys.argv[2]) as allfh:
        for line in allfh:
            if ">" in line:
                writeflag=False
            if line[1:].strip() in circ_contigs:
                writeflag=True
                filename = line[1:].strip()
            if writeflag:
                temp = sys.argv[4]
                with open("%s/%s" %(temp, filename), "a") as circfh:
                    circfh.write(line)

else:
    with open(sys.argv[1]) as allfh:
        writeflag=True
        for line in allfh:
            if ">" in line:
                filename = line[1:].strip()
            if writeflag:
                directory = sys.argv[2]
                with open("%s/%s/%s-final-assembly.fasta" %(directory, filename, filename), "a") as circfh:
                    circfh.write(line)
