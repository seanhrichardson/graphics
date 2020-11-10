import sys
import os

print("prefix low hi delay(4-10ish)")

name= sys.argv[1]
lo = int(sys.argv[2])
hi = int(sys.argv[3])
delay = sys.argv[4]

os.system("mkdir pngs")
print("Making pngs...")
for i in range(lo,hi+1):
#convert [old_name].xwd [new_name].png
    os.system("convert " + name + str(i).zfill(4) + ".xwd pngs/pic" + str(i).zfill(4) + ".png")
print("converting to gif...")
os.system("convert -delay " + delay + " " + "pngs/*.png -loop 0 movie.gif")
