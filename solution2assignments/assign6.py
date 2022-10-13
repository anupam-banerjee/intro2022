#!/usr/bin/env python3

import sys
import prody as pd
import numpy as np

pd.LOGGER._setverbosity('none')

pdb = sys.argv[1]
prot = pd.parsePDB(pdb)

stage = 100
if len(sys.argv) > 2:
    stage = int(sys.argv[2])
    
gnm, a = pd.calcGNM(prot)
anm, a = pd.calcANM(prot)

gflucts = pd.calcSqFlucts(gnm)
aflucts = pd.calcSqFlucts(anm)

print("Max GNM SqFluct: %.3f" % gflucts.max())

if stage >= 65:
    print("Max ANM SqFluct: %.3f" % aflucts.max())

if stage < 70:
    sys.exit(0)
    
    
print("Max Abs Difference: %.3f" % np.abs(gflucts-aflucts).max())

if stage == 70:
    sys.exit(0)
    
if stage > 70:
	gflucts /= np.max(gflucts)
	aflucts /= np.max(aflucts)

diff = np.abs(gflucts-aflucts)
print("Max Abs Norm Difference: %.3f" % diff.max())

if stage < 80:
    sys.exit(0)

res = list(enumerate(diff))
res.sort(key=lambda iv: (-iv[1],iv[0]))


if stage > 90:
	for r in res[:10]:
		i = prot.ca.getIndices()[r[0]]
		a = prot[i]
		print('%.3f %s%s' % (r[1],a.getResname(),a.getResnum()))
elif stage == 90:
	for r in res[:10]:
		print('%0.3f %d' % (r[1],r[0]))


sys.exit(0)

prot.setBetas(0)
prot.ca.setBetas(diff)
pd.writePDB('%s_diffs.pdb'%pdb,prot)


import matplotlib.pylab as plt
plt.plot(prot.ca.getResnums(),gflucts,label="GNM")
plt.plot(prot.ca.getResnums(),aflucts,label="ANM")
plt.legend(loc='best')
plt.xlabel("Resnum")
plt.ylabel("Fluctuation")
plt.savefig('modes.png',bbox_inches=None)
