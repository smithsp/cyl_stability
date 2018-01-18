#! /usr/pppl/python/2.3.5x/bin/python
import Gnuplot,sys,Numeric,sets,cmath
from get_EV import *
ref = int(sys.argv[1])
var_key3 = 'ngrid'#sys.argv[2]
var_key2 = 'Vz0'
var_key = 'b'#sys.argv[3]
sfi = 1e5
if 'sfi' in sys.argv:
  ind = sys.argv.index('sfi')
  sfi = float(sys.argv[ind+1])
d = get_all_list(maxn=2500)
nomatch_keys = get_numeric_nomatch_keys()
nomatch_keys+=[var_key,var_key2,var_key3]
matching_files = get_matching(ref,d,nomatch_keys)
d2 = matching_files[0]
var_list1 = []
var_list2 = []
var_list3 = []
match_num = []
for d1 in matching_files:
  var_list1.append(d1.data[var_key])
  var_list2.append(d1.data[var_key2])
  var_list3.append(d1.data[var_key3])
  match_num.append(d1.data['num'])
(b,Vz0,ev)=get_Vmin_gamma(matching_files)
#for i in range(len(b)):
  #print b[i],Vz0[i]
uv1 = list(sets.Set(var_list1).intersection(sets.Set(b)))
uv2 = list(sets.Set(var_list2))
hexp = 1
if 'hexp' in sys.argv:
  ind = sys.argv.index('hexp')
  hexp = float(sys.argv[ind+1])
if 'var2r' in sys.argv:
  ind = sys.argv.index('var2r')
  uv2=[]
  for var2 in sys.argv[ind+1:]:
    uv2.append(float(var2))
uv2.sort()
g = Gnuplot.Gnuplot()
g._clear_queue()
g.xlabel(var_key)
if var_key=='Vz0':
  g.xlabel('V_{z0}')
if 'norm' in sys.argv:
  g.ylabel('{/Symbol G} / {/Symbol G}_{no-flow}')
else: 
  g.ylabel('{/Symbol G}')
tw = 1.0
if 'WN' in sys.argv:#Switch for normalizing growth rates to the wall time
  tw = matching_files[0].data['tw']
  g.ylabel('{/Symbol G t}_w')
for v2 in uv2:
  inds = Numeric.nonzero(Numeric.array(var_list2)==v2)
  evi = []
  evn = []
  v1list = []
  for v1 in uv1:
    try:
      ind = list(Numeric.take(var_list1,inds)).index(v1)
    except:
      continue
    evi.append(get_h0_extrap(Numeric.take(match_num,inds)[ind],matching_files,hexp=hexp))
    ind = b.index(v1)
    evn.append(evi[-1]/ev[ind])
    v1list.append(v1)
  inds = Numeric.argsort(v1list)
  v1list = Numeric.take(v1list,inds)
  evi = Numeric.take(evi,inds)*tw
  evn = Numeric.take(evn,inds)*tw
  if 'norm' in sys.argv:
    evi = evn
  if 'ASH' in sys.argv: #This is the arcsinh scaling
    for j in range(len(evi)):
      evi[j] = cmath.asinh(sfi*evi[j]).real
  if len(v1list)>0:
    ttl1 = '%s=%g'%(var_key2,v2)
    if var_key2=='tw':
      ttl = '{/Symbol t}_w'
      if v2==-1:
        ttl1 = "%s={/Symbol \245}"%(ttl)
      elif v2==0:
        ttl1 = "%s=0"%(ttl)
      else:
        ttl1 = '%s=10^%d'%(ttl,int((cmath.log10(v2)).real))
    elif var_key2=='Vz0':
        if v2==0:
          ttl1 = 'V_{z0}=0'
        else:
          wA = abs(d2.get_wA())
          ttl1 = 'V_{z0}=%2.2g v_{A}'%(v2/abs(wA/d2.data['k']))
        g('set key box spacing 1.5')
    g._add_to_queue([Gnuplot.Data(v1list,evi,with='linespoints',title=ttl1)])
y0 = 0*v1list
#print min(v1list),max(v1list),y0
g._add_to_queue([Gnuplot.Data(v1list,y0,with='lines 3')])
#print g.itemlist
if 'ASH' in sys.argv:
  set_ytics = 'set ytics ('
  for i in range(2):
    for j in range(-int(cmath.log10(sfi).real)+0,13,1):
      scale = (-1)**i*10**j
      if j%2 != 0:
        set_ytics  = set_ytics + '"" %g 1, '%(cmath.asinh(scale*sfi).real)
      else:
        set_ytics  = set_ytics + '"%g" %g, '%(scale,cmath.asinh(scale*sfi).real)

  set_ytics = set_ytics + '"0" 0 )'
  g(set_ytics)
else:
  g('set log y')
if 'xr' in sys.argv:
  ind = sys.argv.index('xr')
  g('set xrange [%s:%s]'%(sys.argv[ind+1],sys.argv[ind+2]))
  g('set xtics %s,.01,%s'%(sys.argv[ind+1],sys.argv[ind+2]))
  if len(g.itemlist)>4:
    x0 = [float(sys.argv[ind+1]),1.64]
    y0 = [0,0]
    g.itemlist[-1]=Gnuplot.Data(x0,y0,with='lines 3')
g('set key center right Left reverse')
g('set key left top')
#g.refresh()
fn = 'plots/%s_var_%s_%s_extrap.eps'%(min(match_num),var_key,var_key2)
print 'Output to %s'%fn
g.hardcopy(fn[:-3]+'ps',eps=False)
g.hardcopy(fn,eps=True,fontsize=20,solid=True, size=(5,3))
sys.exit()
