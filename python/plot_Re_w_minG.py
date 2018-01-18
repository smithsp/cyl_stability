#! /usr/pppl/python/2.3.5x/bin/python
import get_EV, Gnuplot, Numeric, sys, sets
minn = 1; maxn = 3100
if 'minn' in sys.argv:
  ind = sys.argv.index('minn')
  minn = int(sys.argv[ind+1])
if 'maxn' in sys.argv:
  ind = sys.argv.index('maxn')
  maxn = int(sys.argv[ind+1])
d = get_EV.get_all_list(minn=minn,maxn=maxn)
g = Gnuplot.Gnuplot(persist=1)
g.itemlist = []
r0 = .84
if 'r0' in sys.argv:
  ind = sys.argv.index('r0')
  r0 = float(sys.argv[ind+1])
ref = int(sys.argv[1])
var_key = sys.argv[2]
var_key2 = sys.argv[3]
matching = get_EV.get_matching(ref,d=d,nomatch_keys=get_EV.get_numeric_nomatch_keys()+[var_key,var_key2])
var_list = []
var2_list = []
min_eval = []
nums = []
for i in range(len(matching)):
  var_list.append(matching[i].data[var_key]) 
  var2_list.append(matching[i].data[var_key2])
  nums.append(matching[i].num)
  evals = matching[i].evals
  if 'neg' in sys.argv:
    inds = Numeric.nonzero(evals.real<1e-5)
    evals = Numeric.take(evals,inds)
  ind = Numeric.argmin(evals.imag)
  min_eval.append(evals[ind].real)
print nums   
for u2 in Numeric.sort(list(sets.Set(var2_list))):
  inds = Numeric.nonzero(Numeric.array(var2_list)==u2)
  v1 = Numeric.take(var_list,inds)
  ev1 = Numeric.take(min_eval,inds)
  d1 = matching[inds[0]]
  g.itemlist.append(Gnuplot.Data(v1,ev1,title='%s=%g'%(var_key2,d1.data[var_key2])))
  g.itemlist.append(Gnuplot.Data(v1,d1.data['k']*Numeric.array([matching[j].Vz([r0])[0] for j in inds])+ev1[0]-d1.data['k']*d1.Vz([r0]),title='V_z(r=%g)'%r0))
#print g.itemlist
g('set key left top Left reverse')
g.xlabel('Vz0')   
g.ylabel('{/Symbol w}_{min{/Symbol G}}')
g.title('V_z(r=a)=%g'%matching[0].Vz([matching[0].data['a']]))
fnb = 'plots/%d_var_%s_%s_Re_minG'%(min(nums),var_key,var_key2)
print fnb
g('set key out')
g.hardcopy(fnb+'.eps',eps=True)
g.hardcopy(fnb+'.ps')
g.refresh()
#plot_most_unstable.py 2656 Vr ngrid r0 .84 ly minn 2400 maxn 2700 yr .0001 .01 pk b tw Vz0 epsVz
