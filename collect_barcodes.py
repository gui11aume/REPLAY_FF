import sys

seen = set()
for fname in sys.argv[1:]:
   with open(fname) as f:
      for line in f:
         barcode = line.split()[0]
         seen.update([barcode])

for brcd in seen:
   print brcd
