#!/usr/bin/python3.8
import re
import sys
import getopt

def flye_convert_to_miniasm():
    file_input = None
    file_output = None

    argv = sys.argv[1:]

    try:
        opts, args = getopt.getopt(argv, "i:o:", ["input =", "output ="])
    except:
        print("Error")

    for opt, arg in opts:
        if opt in ['-i', '--input']:
            file_input = arg
        elif opt in ['-o', '--output']:
            file_output = arg

    print('Converting {0} to pseudo - miniasm format'.format(file_input))
    with open(file_input,'r') as file:
        filedata = file.read()
        formatted = re.sub(r'(edge_\d+)('	')', r'\1l\2',filedata)
        formatted = re.sub(r'(	\d+	)', r'	',formatted)
        if formatted:
            print(len(formatted),'changes to make')

    with open(file_output, 'w') as file:
        file.write(formatted)
        print('Done! Saved as {} in GFA1 format'.format(file_output))


flye_convert_to_miniasm()
