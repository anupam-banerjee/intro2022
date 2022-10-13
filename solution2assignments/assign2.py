#!/usr/bin/python3
from pysb import * 
from pysb.macros import catalyze
from pysb.macros import catalyze_one_step
from pysb.macros import bind
from pysb.integrate import odesolve
import numpy as np
import sys

#set inh_kf to the first commandline argument interpretted as a float
inh_kf = float(sys.argv[1])

#these are our default reaction rates
kf_bind = 1e-5
kr_bind = 1e-1
kcat = 1e-1

stage = 5
if len(sys.argv) > 2:
  stage = int(sys.argv[2])

Model()


Monomer('PG3',['b'])
Monomer('PG2',['b'])
Monomer('PGM',['b'])
Monomer('INH',['b'])
Monomer('ENO',['b'])
Monomer('PEP',['b'])
Monomer('PKM',['b','ba'])
Monomer('PYR',['b'])
Monomer('ADP',['ba'])
Monomer('ATP')
Monomer('LDH',['b','bn'])
Monomer('LAC',['b'])
Monomer('NADH',['bn'])
Monomer('NAD',['bn'])


Parameter('kc', kcat)
Parameter('INH_0', 5000)
Parameter('PG3_0', 20000)
Parameter('ADP_0', 10000)
Parameter('NADH_0', 2000)
Parameter('PGM_0',5000)
Parameter('ENO_0',5000)
Parameter('PKM_0',5000)
Parameter('LDH_0',1000)

Initial(INH(b=None), INH_0)
Initial(PG3(b=None), PG3_0)
Initial(ADP(ba=None), ADP_0)
Initial(NADH(bn=None), NADH_0)
Initial(PGM(b=None), PGM_0)
Initial(ENO(b=None), ENO_0)
Initial(PKM(b=None,ba=None), PKM_0)
Initial(LDH(b=None,bn=None), LDH_0)

Observable('oNAD',NAD(bn=None))
Observable('oNADH',NADH(bn=None))
Observable('oADP',ADP(ba=None))
Observable('oPYR',PYR(b=None))
Observable('oLAC',LAC(b=None))
Observable('oPEP',PEP(b=None))
Observable('oPG2',PG2(b=None))
Observable('oPG3',PG3(b=None))

brates = [kf_bind, kr_bind]
rates = [kf_bind, kr_bind, kcat]

#need at least one rule
bind(PGM,'b',INH,'b',[inh_kf,kr_bind])

if stage >= 1:
  catalyze(PGM(), 'b', PG3(), 'b', PG2(b=None), rates)

if stage >= 2:
  catalyze(ENO(), 'b', PG2(), 'b', PEP(b=None), rates)

if stage >= 3:
  bind(PKM,'b',PEP,'b',brates)
  bind(PKM,'ba',ADP,'ba',brates)
  Rule('pkm_catalyze',PKM(b=1,ba=2) % ADP(ba=2) % PEP(b=1) >> PKM(b=None,ba=None) + ATP() + PYR(b=None), kc)

if stage >= 4:
  bind(LDH,'b',PYR,'b',brates)
  bind(LDH,'bn',NADH,'bn',brates)
  Rule('ldh_catalyze_forward',LDH(b=1,bn=2) % NADH(bn=2) % PYR(b=1) >> LDH(b=None,bn=None) + NAD(bn=None) + LAC(b=None), kc)
  
if stage >= 5:
  bind(LDH,'b',LAC,'b',brates)
  bind(LDH,'bn',NAD,'bn',brates)
  Rule('ldh_catalyze_backward',LDH(b=1,bn=2) % NAD(bn=2) % LAC(b=1) >> LDH(b=None,bn=None) + NADH(bn=None) + PYR(b=None), kc)


#solve for several datapoints from 0 to 600 seconds
t = np.linspace(0,600)
out = odesolve(model, t)

#print out observables
for i in range(0,len(out),7):
  print('NAD+: %.0f  LAC: %.0f  PYR: %.0f  PEP: %.0f  PG2: %.0f  PG3: %.0f' % \
  (out['oNAD'][i],out['oLAC'][i],out['oPYR'][i],out['oPEP'][i],out['oPG2'][i],out['oPG3'][i]))


