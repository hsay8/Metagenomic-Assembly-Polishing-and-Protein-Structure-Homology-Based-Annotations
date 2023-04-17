#!/usr/bin/python3.8
import re
import sys
import getopt

def flye_convert_to_miniasm():

    print('Converting {} to pseudo - miniasm format'.format(snakemake.input[0]))
    with open(snakemake.input[0],'r') as file:
        filedata = file.read()
        formatted = re.sub(r'(edge_\d+)('	')', r'\1l\2',filedata)
        formatted = re.sub(r'(	\d+	)', r'	',formatted)
        if formatted:
            print(len(formatted),'changes to make')

    with open(snakemake.output[0], 'w') as file:
        file.write(formatted)
        print('Done! Saved {} as in GFA1 format'.format(snakemake.output[0]))

flye_convert_to_miniasm()
