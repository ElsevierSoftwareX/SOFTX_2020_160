
# script used to convert posterior samples from json to dat format
# python json2dat posteriors.json
# output -> xyz.dat
# before to run it -> source /home/vedovato/virtualenv/pycwb/bin/activate

import sys
try:
    from pesummary.gw.file.read import read
except:
    print("")
    print("error loading pesummary")
    print("try")
    print("source /home/vedovato/virtualenv/pycwb/bin/activate")
    print("")
    sys.exit(1)

arguments = len(sys.argv) - 1
if arguments == 0:
    # output error, and return with an error code
    print("")
    print("missing json file")
    print("")
    sys.exit(1)

print("")
print(sys.argv[1])
print("")

f = read(sys.argv[1])
f.to_dat(outdir="./")

