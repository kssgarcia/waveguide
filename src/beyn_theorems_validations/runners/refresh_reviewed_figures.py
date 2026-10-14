"""Rebuild reviewed figures from signed historical tables, never call solvers."""
import sys,json,csv,typing,types,shutil,hashlib,ast,subprocess
from dataclasses import fields,replace
from pathlib import Path
from datetime import datetime
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
import first_theorem_beyn as f
import second_theorem_X as x
import second_theorem_Y as y
import execution_harness as h

BASE=Path(__file__).resolve().parents[1]

def records(path,cls):
    hints=typing.get_type_hints(cls);out=[]
    for row in csv.DictReader(Path(path).open()):
        values={}
        for field in fields(cls):
            if field.name not in row:continue
            v=row[field.name];t=hints[field.name]
            if typing.get_origin(t) in (typing.Union,types.UnionType):
                t=next(a for a in typing.get_args(t) if a is not type(None))
                if v=='':values[field.name]=None;continue
            values[field.name]=(v=='True') if t is bool else t(v)
        out.append(cls(**values))
    return out

def config(module,run):
    data=json.loads((run/'run_manifest.json').read_text())['config']
    return module.Config(**{k:tuple(v) if isinstance(v,list) else v for k,v in data.items()})

def latest(name):return Path(json.loads((BASE/'results'/name/'runs/latest_run.json').read_text())['directory'])

def nonplot_ast(source):
    tree=ast.parse(source)
    tree.body=[n for n in tree.body if not isinstance(n,(ast.FunctionDef,ast.AsyncFunctionDef)) or not n.name.startswith('plot_')]
    return ast.dump(tree,include_attributes=False)

def copy_pngs(source,dest):
    dest.mkdir(parents=True,exist_ok=True)
    for p in source.glob('*.png'):shutil.copy2(p,dest/p.name)

def main():
    output=BASE/'results'/'reviewed_figures'/datetime.now().strftime('%Y-%m-%d_%H-%M-%S')
    output.mkdir(parents=True,exist_ok=False)
    provenance={'purpose':'presentation-only regeneration; no solves','sources':{},'files':{},'numerical_code_unchanged':{}}
    for module in (f,x,y):
        path=Path(module.__file__);relative=path.relative_to(BASE.parents[1])
        # Compare to the staged scientific version (including pre-existing work).
        original=subprocess.check_output(['git','show',f':{relative}'],text=True)
        unchanged=nonplot_ast(original)==nonplot_ast(path.read_text())
        provenance['numerical_code_unchanged'][path.name]=unchanged
        if not unchanged:raise RuntimeError(f'Non-plot code changed: {path.name}')
    for name,module in [('theorem_2_3_iii',x),('theorem_2_3_iv',y)]:
        run=latest(name);dest=output/name;copy_pngs(run,dest);cfg=config(module,run)
        summaries=records(run/'summary.csv',module.ValidationResult)
        convergence=records(run/'internal_convergence.csv',module.InternalConvergenceRow)
        module.plot_internal_convergence(convergence,dest,cfg)
        module.plot_summary(summaries,dest,cfg)
        if module is x:module.plot_mu_convergence(records(run/'mu_mfs_convergence.csv',module.MuMFSRow),dest,cfg)
        else:
            module.plot_mfs_convergence(records(run/'mfs_convergence.csv',module.MFSRow),dest,cfg)
            module.plot_placement_a1_fit(summaries,summaries[0].a1_formula_2_7,dest,cfg)
            for result in summaries:
                cp=run/'epsilon_cases'/f'epsilon_{result.epsilon:.8f}'/'completed.json'
                payload=h.decode(json.loads(json.loads(cp.read_text())['payload']),module)
                module.plot_uniqueness_scan(result.epsilon,payload[3],payload[4],result,dest,cfg)
            # Effective physical diagnostic with PASS/ambiguous separated.
            fig,axes=y.plt.subplots(1,2,figsize=(12,5))
            for parameter,marker in [('finite_difference_step','o'),('lattice_terms','s'),('harmonic_order','^')]:
                sub=[r for r in convergence if r.parameter==parameter]
                for r in sub:
                    color='green' if r.physical_bic_verified else 'red'
                    axes[0].scatter(r.normalized_residual,r.open_channel_relative_amplitude,marker=marker,color=color)
                    axes[1].scatter(r.value,r.a_shift_from_reference,marker=marker,color=color)
            axes[0].set_xscale('log');axes[0].set_yscale('log');axes[0].axvline(cfg.singular_vector_residual_tolerance,color='black',ls='--');axes[0].axhline(cfg.physical_open_channel_relative_tolerance,color='black',ls='--')
            axes[0].set_xlabel('normalized residual');axes[0].set_ylabel('open-channel ratio')
            axes[1].set_xscale('log');axes[1].set_yscale('symlog',linthresh=1e-12);axes[1].set_ylim(bottom=0);axes[1].set_xlabel('parameter value (o=h, square=lattice, triangle=H)');axes[1].set_ylabel('|a - a_reference|')
            fig.suptitle('Retuned sensitivity: green=physical PASS; red=ambiguous; absent metrics are not PASS')
            fig.tight_layout();fig.savefig(dest/'sensitivity_physical_certificates.png',dpi=220);y.plt.close(fig)
            detail=json.loads((run/'internal_physical_certificates.json').read_text())
            fig,ax=y.plt.subplots(figsize=(12,5))
            labels=[]
            for j,record in enumerate(detail):
                labels.append(f"{record['parameter']}={record['value']:g}")
                physical=record.get('physical')
                if physical is None:
                    ax.scatter(j,1.,marker='x',color='red');ax.annotate('not evaluated',(j,1.),rotation=90,fontsize=8)
                    continue
                metrics=[physical['wall_relative_residual']/cfg.physical_wall_relative_tolerance,
                         physical['boundary_integral_relative_residual']/cfg.physical_boundary_residual_tolerance,
                         physical['first_open_channel_relative_amplitude']/cfg.physical_open_channel_relative_tolerance,
                         physical['modal_projection_discrepancy']/cfg.physical_modal_projection_tolerance,
                         max(physical['decay_relative_error_left'],physical['decay_relative_error_right'])/cfg.physical_decay_relative_tolerance]
                ax.scatter(j,max(metrics),marker='o' if record['verified'] else 'x',color='green' if record['verified'] else 'red',s=55)
            ax.axhline(1.,color='black',ls='--',label='physical tolerance limit')
            ax.set_yscale('log');ax.set_xticks(range(len(labels)),labels,rotation=65,ha='right',fontsize=8)
            ax.set_ylabel('Largest physical metric / tolerance')
            ax.set_title('Complete physical sensitivity: green=certified; red=ambiguous')
            ax.legend();fig.tight_layout();fig.savefig(dest/'sensitivity_physical_limits.png',dpi=220);y.plt.close(fig)
        provenance['sources'][name]=str(run)
        for p in run.glob('*.csv'):provenance['files'][str(p)]=hashlib.sha256(p.read_bytes()).hexdigest()
    name='theorem_2_1';run=latest(name);baseline=run if (run/'summary.csv').exists() else BASE/'results'/name/'runs/2026-10-07_11-06-56_full_832cf4c8c71d';dest=output/name;copy_pngs(baseline,dest)
    cfg=config(f,baseline)
    old=records(baseline/'paper_a_sweep_M24.csv',f.PaperSweepRow);new=records(run/'paper_a_sweep_M24.csv',f.PaperSweepRow)
    merged={(r.a,r.epsilon):r for r in old};merged.update({(r.a,r.epsilon):r for r in new});rows=list(merged.values())
    f.write_dataclass_csv(dest/'consolidated_paper_a_sweep_M24.csv',rows)
    f.plot_paper_kb_vs_epsilon_individual(rows,cfg,dest);f.plot_paper_kb_vs_epsilon_all_a(rows,cfg,dest)
    f.plot_critical_height(records(run/'critical_height.csv',f.CriticalHeightRow),dest)
    f.plot_internal_convergence(records(baseline/'internal_convergence.csv',f.InternalConvergenceRow),dest)
    # Fixed-M candidate panels are historical localization, not mesh convergence.
    for p in dest.glob('candidate_*.png'):p.unlink()
    # Replace single-M candidate plots with an explicitly historical certificate panel.
    summary=list(csv.DictReader((baseline/'summary.csv').open()));fig,axes=f.plt.subplots(1,2,figsize=(11,4.5))
    axes[0].plot([float(r['epsilon']) for r in summary],[float(r['kb_numerical']) for r in summary],'o-');axes[0].set_ylabel('kb');axes[0].set_xlabel('epsilon')
    axes[1].semilogy([float(r['epsilon']) for r in summary],[float(r['normalized_residual']) for r in summary],'o-');axes[1].set_ylabel('normalized residual');axes[1].set_xlabel('epsilon')
    version=json.loads((baseline/'run_manifest.json').read_text())['assembly_version']
    fig.suptitle(f'2.1 {version}: fixed-M=24 localization; not spatial convergence')
    fig.tight_layout();fig.savefig(dest/'baseline_localization_certificates.png',dpi=220);f.plt.close(fig)
    provenance['sources'][name]={'baseline':str(baseline),'supplement':str(run),'limitation':'Baseline from final run when available; otherwise historical baseline with explicit supplemental table overrides.'}
    for source in (baseline,run):
        for p in source.glob('*.csv'):provenance['files'][str(p)]=hashlib.sha256(p.read_bytes()).hexdigest()
    for folder in output.iterdir():
        if folder.is_dir():
            (folder/'README.md').write_text('Figures regenerated from saved tables; no spectral computation. See ../manifest.json for source versions.\nGlobal uniqueness and Big-O are not established by these figures.\n')
    provenance['generated_pngs']=len(list(output.rglob('*.png')))
    h.atomic_json(output/'manifest.json',provenance)
    (output/'README.md').write_text('# Figuras revisadas\n\nSin repetir barridos. Originales preservados. Las tablas fuente y su versión están en manifest.json.\n\nX: existencia y unicidad impar; unicidad global no establecida. Y: existencia física; screen finito no certifica unicidad global. 2.1: baseline histórico y suplemento final explícitamente separados. Las curvas asintóticas son términos principales, no pruebas de Big-O.\n')
    h.atomic_json(output.parent/'latest_figures.json',{'directory':str(output.resolve())})
    print(output.resolve(),flush=True)

if __name__=='__main__':main()
