from pathlib import Path
exec(open(Path(__file__).with_name('wg_green_audit.py')).read().split('for k in [')[0])
import csv
path=str(Path(__file__).resolve().parents[3] / 'src/beyn_theorems_validations/theorem_2_1_publication_ready/paper_a_sweep_M32.csv')
out=[]
for row in csv.DictReader(open(path)):
 if row['status']!='one': continue
 a,e,k=map(float,[row['a'],row['epsilon'],row['kb']])
 for H in [20,40,60]:
  S=l.lattice_sums(4.,k,0.,200,H)
  # Distinct boundary locations, including the near-wall arc.
  errors=[]
  for t,s in [(.4,1.3),(1.1,1.8),(.2,2.4)]:
   X=e*(np.cos(t)-np.cos(s)); y=1+a+e*np.sin(t); eta=1+a+e*np.sin(s)
   v=l.greens_dirichlet(X,y,0.,eta,S,k,2.); ref=modal(X,y,eta,complex(k),N=2000)
   errors.append(abs(v-ref)/max(abs(ref),1e-30))
  out.append(dict(a=a,e=e,H=H,max_relative_error=max(errors)))
for H in [20,40,60]:
 rows=[r for r in out if r['H']==H]; print('H',H,'worst',max(rows,key=lambda r:r['max_relative_error']))
 print('above1e-4',sum(r['max_relative_error']>1e-4 for r in rows),'of',len(rows))
json.dump(out,open('/tmp/wg_sweep_green_audit.json','w'),indent=2)
print('elapsed',time.time()-start)
