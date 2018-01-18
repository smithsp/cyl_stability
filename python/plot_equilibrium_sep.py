#! /usr/pppl/python/2.3.5x/bin/python
import Gnuplot, Numeric, sys, os, get_EV
num = int(sys.argv[1])
try:
  equilib_path = 'equilibria_vcyl'
  filename = 'plots/%d_equilibrium.eps'%(num)
  infile = '%s/%d.txt'%(equilib_path,num)
  f= file(infile,'r')
  s = f.readline()
  f.close()
except IOError:
  equilib_path = 'equilibria_cyl'
  filename = 'plots/%d_equilibrium.eps'%(num)
  infile = '%s/%d.txt'%(equilib_path,num)
  f= file(infile,'r')
  s = f.readline()
  f.close()
  #print infile + ' does not exist'
  #sys.exit()
d = get_EV.get_EV(num)
g = Gnuplot.Gnuplot()
#g.xlabel('r')
g('set xtics ("0" 0, "a" 1)')
g._clear_queue()
g._add_to_queue([Gnuplot.File(infile,using='1:2',title=None,with='lines')])
g('set ytics ("p_0" %g)'%(d.data['P0']))
g('set y2tics ("0" 0)')
#g.title('pressure')
g.hardcopy('plots/%s_pressure.eps'%num,eps=True,color=False,fontsize=24,linewidth=3)
g._clear_queue()
g._add_to_queue([Gnuplot.File(infile,using='1:8',title=None,with='lines')])
q0 = (d.data['Bz0']**2+2*(d.data['Bt0']**2-d.data['P0']))**.5*d.data['k']/d.data['Bt0']
qa = d.data['Bz0']*d.data['k']/d.data['Bt0']
g('set ytics ("q_0" %g)'%(q0))
g('set y2tics ("q_a" %g)'%(qa))
#g.title('q')
g('set yrange [%g:%g]'%(.95*q0,1.05*qa))
g.hardcopy('plots/%s_q.eps'%num,eps=True,color=False,fontsize=24,linewidth=3)


g._clear_queue()
g._add_to_queue([Gnuplot.File(infile,using='1:($3/%s)'%d.data['Bz0'],title=None,with='lines')])
Bzmax=((d.data['Bz0']**2+2*(d.data['Bt0']**2-d.data['P0']))**.5/d.data['Bz0'])
g('set label "[B_{za}^2+2B_{{/Symbol q}a}^2-2p_0]^{1/2}" at 0.02,%g'%(Bzmax*1.005))
g('set ytics (" " %g)'%(Bzmax))
g('set y2tics ("B_{za}" 1)')
g('set yrange [%g:%g]'%(.99,Bzmax*1.01))
g.title('B_z')
g.hardcopy('plots/%s_Bz.eps'%num,eps=True,color=False,fontsize=24,linewidth=3)

g('set yrange [*:*]')
g('unset label')
g._clear_queue()
g._add_to_queue([Gnuplot.File(infile,using='1:4',with='lines',title=None)])
g('set ytics ("0" 0)')
g('set y2tics ("B_{{/Symbol q}a}" %g)'%(d.data['Bt0']))
g.title('B_{/Symbol q}(r)')
g.hardcopy('plots/%s_Bt.eps'%num,eps=True,color=False,fontsize=24,linewidth=3)

g._clear_queue()
g._add_to_queue([Gnuplot.File(infile,using='1:5',title=None,with='lines')])
g._add_to_queue([Gnuplot.File(infile,using='1:6',title=None,with='lines')])
g.title('{/Symbol r}, V_z')
g('set ytics ("V_{z0}" 0, "{/Symbol r}_0" %d)'%d.data['rho0'])
g('set yrange [-.1:1.1]')
g.hardcopy('plots/%s_rhoVz.eps'%num,eps=True,color=False,fontsize=24,linewidth=3)
sys.exit()

g._clear_queue()
g._add_to_queue([Gnuplot.File(infile,using='1:7')])
g.title('V_z')
g.hardcopy('plots/%s_Vz.eps'%num,eps=True)

sys.exit()

#g('set title "B_z(r)"')
#g('plot "%s/%d.txt" using 1:3 notitle with lines'%(equilib_path,num))

#g('set title "B_{/Symbol q}(r)"')
#g('plot "%s/%d.txt" using 1:4 notitle with lines'%(equilib_path,num))

#g('set title "q(r)"')
#g('plot "%s/%d.txt" using 1:8 notitle with lines'%(equilib_path,num))

#g('set title "{/Symbol r}(r)"')
#g('plot "%s/%d.txt" using 1:5 notitle with lines'%(equilib_path,num))

#g('set title "V_{/Symbol q}(r)"')
#g('plot "%s/%d.txt" using 1:6 notitle with lines'%(equilib_path,num))

#g('set title "V_z(r)"')
#g('plot "%s/%d.txt" using 1:7 notitle with lines'%(equilib_path,num))

#os.system('gv %s &'%(filename))
