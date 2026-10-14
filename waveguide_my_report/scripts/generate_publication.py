#!/usr/bin/env python3
"""Frozen-data publication generation; no spectral searches or eigenvalue tuning.

Run from any directory: python waveguide_my_report/scripts/generate_publication.py
Existing publication artifacts are replaced deterministically; solver outputs stay intact.
--reconstruct-fields computes two saved-root examples and archives their complex arrays.
Subsequent default invocations reuse the checked latest reconstruction.
"""
from __future__ import annotations
import argparse, dataclasses, hashlib, importlib, json, shutil, sys, platform
from pathlib import Path
from datetime import datetime, timezone
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.path import Path as MplPath

ROOT=Path(__file__).resolve().parents[2]
PAPER=ROOT/'waveguide_my_report'
LOCK=json.loads((PAPER/'scripts/source_lock.json').read_text())
RUNS={k:ROOT/v for k,v in LOCK['runs'].items()}
COLORS=['#0072B2','#D55E00','#009E73']
LABELS={'theorem_2_1':'Theorem 2.1','theorem_X':'Theorem 2.3(iii)','theorem_Y':'Theorem 2.3(iv)'}
MARKERS=['o','s','^']
ARTIFACTS=[]
plt.rcParams.update({'font.size':9,'axes.labelsize':9,'legend.fontsize':7.5,
    'xtick.labelsize':8,'ytick.labelsize':8,'pdf.fonttype':42,'ps.fonttype':42,
    'figure.dpi':140,'savefig.facecolor':'white','axes.spines.top':False,
    'axes.spines.right':False,'lines.linewidth':1.25,'lines.markersize':4})


def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def csv(k,n):return pd.read_csv(RUNS[k]/n)
def check_sources():
    for rel,expected in LOCK['sources'].items():
        if sha(ROOT/rel)!=expected:raise RuntimeError(f'Frozen scientific source changed: {rel}')
    for k in RUNS:
        if len(csv(k,'summary.csv'))!=10:raise RuntimeError(f'Incomplete baseline: {k}')
        if not csv(k,'summary.csv').mode_exists_supported.all():raise RuntimeError(f'Uncertified baseline: {k}')

def export(fig,name,sources,transformations):
    fig.savefig(PAPER/'figures'/f'{name}.pdf')
    fig.savefig(PAPER/'figures'/f'{name}.png',dpi=300)
    plt.close(fig)
    ARTIFACTS.append({'name':name,'files':[f'figures/{name}.pdf',f'figures/{name}.png'],
        'sources':sources,'transformations':transformations,'uncertainty':'none unless deterministic brackets explicitly displayed'})

def panel(ax,letter):ax.text(-.14,1.04,f'({letter})',transform=ax.transAxes,fontweight='bold')
def vals(s,col):return s[col].to_numpy(float)
def coeff(s):return float(s.asymptotic_coefficient.iloc[0])
def texnum(x,digits=4):
    x=float(x)
    if not np.isfinite(x):return r'\textemdash'
    if x==0:return '$0$'
    if abs(x)<.001 or abs(x)>=10000:
        m,e=f'{x:.{digits-1}e}'.split('e');return f'${m}\\times10^{{{int(e)}}}$'
    return f'${x:.{digits}g}$'

def table(name,head,rows):
    n=len(head);lines=[r'\begin{tabular}{'+'r'*n+'}',r'\toprule',' & '.join(head)+r' \\',r'\midrule']
    lines.extend(' & '.join(row)+r' \\' for row in rows)
    lines.extend([r'\bottomrule',r'\end{tabular}'])
    (PAPER/'tables'/f'{name}.tex').write_text('\n'.join(lines)+'\n')

def make_tables(S):
    # Retain all columns and rows in companion CSVs; compact typeset subsets.
    for k,s in S.items():
        p=csv(k,'physical_diagnostics.csv')
        open_col='first_open_channel_relative_amplitude'
        rows=[]
        for _,r in s.iterrows():
            if k=='theorem_Y':
                rows.append([texnum(r.epsilon,3),texnum(r.a_leading,5),texnum(r.a_numerical,5),texnum(r.scaled_a_correction,4),texnum(r.sigma_numerical,5),texnum(r.normalized_residual,3)])
            else:
                rows.append([texnum(r.epsilon,3),str(int(r.beyn_estimated_rank)),texnum(r.kb_numerical,8),texnum(r.sigma_numerical,5),texnum(r.sigma_over_epsilon_squared,5),texnum(r.normalized_residual,3)])
        heads=([r'$\varepsilon$',r'$\varepsilon a_1$',r'$a_{\rm BIC}/b$',r'$c$',r'$\sigma$',r'$r_A$'] if k=='theorem_Y' else
               [r'$\varepsilon$',r'$r_\Gamma$',r'$kb$',r'$\sigma$',r'$\sigma/\varepsilon^2$',r'$r_A$'])
        table(k,heads,rows)
    diagnostics=[]
    for k,s in S.items():
        p=csv(k,'physical_diagnostics.csv')
        for eps in [.02,.10,.20]:
            r=s[np.isclose(s.epsilon,eps)].iloc[0];q=p[np.isclose(p.epsilon,eps)].iloc[0]
            channel=q.get('first_open_channel_relative_amplitude',np.nan)
            dec=max(q.decay_relative_error_left,q.decay_relative_error_right)
            diagnostics.append([{'theorem_2_1':'2.1','theorem_X':'2.3(iii)','theorem_Y':'2.3(iv)'}[k],texnum(eps,3),texnum(r.normalized_residual,3),texnum(q.boundary_integral_relative_residual,3),texnum(q.wall_relative_residual,3),texnum(channel,3),texnum(dec,3)])
    table('diagnostics',['Case',r'$\varepsilon$',r'$r_A$',r'$r_{\rm BIE}$',r'$r_{\rm wall}$',r'$r_1$',r'$r_{\rm decay}$'],diagnostics)
    configs={k:json.loads((RUNS[k]/'run_manifest.json').read_text())['config'] for k in RUNS}
    rows=[]
    for k,c in configs.items():
        rows.append([{'theorem_2_1':'2.1','theorem_X':'2.3(iii)','theorem_Y':'2.3(iv)'}[k],'24',str(c['harmonic_order']),str(c['lattice_terms']),texnum(c['finite_difference_step'],2),texnum(c['physical_boundary_residual_tolerance'],2)])
    table('configuration',['Case',r'$M$',r'$H$',r'$N_{\rm lat}$',r'$h/b$',r'BIE tolerance'],rows)
    # In TeX use decimal macros directly inside math; no thousands separators.
    macros={}
    for k,s in S.items():
        prefix={'theorem_2_1':'TOne','theorem_X':'TX','theorem_Y':'TY'}[k]
        macros[prefix+'Leading']=coeff(s)
        macros[prefix+'MaxResidual']=s.normalized_residual.max()
        phys=csv(k,'physical_diagnostics.csv');macros[prefix+'MaxBIE']=phys.boundary_integral_relative_residual.max()
        macros[prefix+'MaxWall']=phys.wall_relative_residual.max()
        macros[prefix+'MaxDecay']=float(max(phys.decay_relative_error_left.max(),phys.decay_relative_error_right.max()))
        if 'first_open_channel_relative_amplitude' in phys:macros[prefix+'MaxChannel']=phys.first_open_channel_relative_amplitude.max()
        if 'mu' in s:macros[prefix+'Mu']=s.mu.iloc[0];macros[prefix+'Nu']=s.nu.iloc[0]
        for eps,tag in [(.02,'Small'),(.10,'Middle'),(.20,'Large')]:
            r=s[np.isclose(s.epsilon,eps)].iloc[0]
            macros[prefix+'Error'+tag]=100*r.relative_error_sigma
    sy=S['theorem_Y'];macros['TYAOne']=sy.a1_formula_2_7.iloc[0]
    macros['TOneCriticalLeading']=csv('theorem_2_1','critical_height.csv').a0_star.iloc[0]
    macros['TOneBracketWidth']=csv('theorem_2_1','critical_height.csv').bracket_width.max()
    fit=csv('theorem_Y','placement_a1_fit.csv');
    for _,r in fit.iterrows():
        if r.point_count==5:
            macros['TYAOneFit']=r.a1_fit_linear;macros['TYAOneFitPercent']=100*r.linear_relative_difference;break
    lines=[]
    for k,v in macros.items():
        if abs(v)<.001:
            m,e=f'{v:.5e}'.split('e');t=f'{m}\\times10^{{{int(e)}}}'
        else:t=f'{v:.7g}'
        lines.append(f'\\newcommand{{\\{k}}}{{{t}}}')
    (PAPER/'tables/numeric_macros.tex').write_text('\n'.join(lines)+'\n')
    (PAPER/'data/numeric_macros.json').write_text(json.dumps(macros,indent=2)+'\n')


def geometry(S):
    theta=np.linspace(0,2*np.pi,600)
    fig,axs=plt.subplots(1,3,figsize=(7.1,2.6),layout='constrained')
    for i,(k,ax) in enumerate(zip(S,axs)):
        eps=.10;beta=.5
        if i==0:x=eps*np.cos(theta);y=.6+eps*np.sin(theta)
        elif i==1:x=eps*(np.cos(theta)-beta/2*np.cos(2*theta));y=eps*(np.sin(theta)-beta/2*np.sin(2*theta))
        else:
            a=float(S[k][np.isclose(S[k].epsilon,eps)].a_numerical.iloc[0]);x=eps*(np.sin(theta)-beta/2*np.sin(2*theta));y=a+eps*(-np.cos(theta)+beta/2*np.cos(2*theta))
        ax.fill(x,y,facecolor='#999999',edgecolor='black')
        ax.axhline(1,color='black');ax.axhline(-1,color='black');ax.axhline(0,color='#777777',lw=.65,ls=':')
        ax.set(xlim=(-1,1),ylim=(-1.1,1.1),xlabel=r'$x/b$',ylabel=r'$y/b$',title=LABELS[k]);ax.set_aspect('equal');panel(ax,'abc'[i])
    export(fig,'geometry',[f'{v}/run_manifest.json' for v in LOCK['runs'].values()]+[f'{LOCK["runs"]["theorem_Y"]}/summary.csv'],['Current parametrizations at epsilon=.10; obstacle scale not magnified; b=1'])


def spectral(S):
    fig,axs=plt.subplots(2,3,figsize=(7.1,4.5),layout='constrained')
    for i,(k,s) in enumerate(S.items()):
        e=vals(s,'epsilon');c=coeff(s);ratio=vals(s,'sigma_numerical')/e**2/c
        axs[0,i].plot(e,ratio,color=COLORS[i],marker=MARKERS[i],label='BEM')
        axs[0,i].axhline(1,color='black',ls='--',lw=1,label='Leading term')
        axs[0,i].set(title=LABELS[k],ylabel=r'$\sigma/(C\varepsilon^2)$',xlabel=r'$\varepsilon$',ylim=(.45,1.04));panel(axs[0,i],'abc'[i])
        axs[1,i].plot(e,(vals(s,'sigma_numerical')-c*e**2)/(e**3*np.abs(np.log(e))),color=COLORS[i],marker=MARKERS[i])
        axs[1,i].set(xlabel=r'$\varepsilon$',ylabel=r'$(\sigma-C\varepsilon^2)/(\varepsilon^3|\log\varepsilon|)$');panel(axs[1,i],'def'[i])
    axs[0,0].legend()
    export(fig,'spectral_scaling',[f'{v}/summary.csv' for v in LOCK['runs'].values()],['Dimensionless leading coefficient normalization; signed scaled remainder, all ten points; no fitted error exponent'])


def beyn():
    fig,axs=plt.subplots(3,2,figsize=(7.1,6.0),layout='constrained')
    for j,k in enumerate(['theorem_2_1','theorem_X']):
        d=csv(k,'beyn_diagnostics.csv');c=csv(k,'beyn_convergence.csv')
        for i,eps in enumerate([.02,.10,.20]):
            q=d[np.isclose(d.epsilon,eps)].sort_values('quadrature_points');z=c[np.isclose(c.epsilon,eps)].sort_values('quadrature_points')
            lab=fr'$\varepsilon={eps:.2f}$'
            axs[0,j].plot(q.quadrature_points,q.threshold_rank,color=COLORS[i],marker=MARKERS[i],label=lab+' threshold')
            axs[0,j].plot(q.quadrature_points,q.gap_rank,color=COLORS[i],marker='x',ls='--',label=lab+' gap')
            # First comparison is undefined; preserve missing values, never fill zero.
            axs[1,j].plot(q.quadrature_points,q.s0_relative_change,color=COLORS[i],marker=MARKERS[i],label=lab+r' $S_0$')
            axs[1,j].plot(q.quadrature_points,q.s1_relative_change,color=COLORS[i],marker='x',ls='--',label=lab+r' $S_1$')
            axs[2,j].plot(z.quadrature_points,z.max_candidate_shift_from_previous,color=COLORS[i],marker=MARKERS[i],label=lab)
        config=json.loads((RUNS[k]/'run_manifest.json').read_text())['config'];tol=config['beyn_moment_convergence_tolerance']
        axs[1,j].axhline(tol,color='black',lw=.8,ls=':',label='Moment tolerance')
        axs[0,j].set(ylim=(-.1,3.2),yticks=[0,1,2,3],ylabel='Estimated rank',title=LABELS[k])
        axs[1,j].set(yscale='log',ylabel='Relative moment change')
        axs[2,j].set_yscale('symlog',linthresh=1e-12);axs[2,j].set(ylabel='Candidate shift in $kb$',xlabel=r'Contour points $N_q$')
        for n in range(3):panel(axs[n,j],'abcdef'[2*n+j])
    axs[0,0].legend(ncol=2,fontsize=6);axs[1,0].legend(ncol=2,fontsize=6);axs[2,0].legend(fontsize=7)
    export(fig,'beyn_convergence',[f'{LOCK["runs"][k]}/{f}' for k in ['theorem_2_1','theorem_X'] for f in ['beyn_diagnostics.csv','beyn_convergence.csv','run_manifest.json']],['Representative epsilons .02/.10/.20; raw threshold/gap ranks and both moment changes; first undefined changes remain missing; candidate shift symlog retains exact zeros'])


def critical():
    d=csv('theorem_2_1','paper_a_sweep_M24.csv');q=csv('theorem_2_1','critical_height.csv')
    fig,axs=plt.subplots(1,2,figsize=(7.1,3.2),layout='constrained')
    styles={'one':('o','#0072B2','One certified mode'),'zero-supported':('x','#D55E00','Exclusion supported'),'ambiguous':('s','#777777','Unresolved')}
    for st,(m,c,l) in styles.items():
        z=d[d.status==st]
        if len(z):axs[0].scatter(z.epsilon,z.a,marker=m,color=c,s=18,label=l)
    axs[0].axhline(q.a0_star.iloc[0],color='black',ls='--',lw=1,label=r'$a_0^*$')
    axs[0].set(xlabel=r'$\varepsilon$',ylabel=r'$a/b$');axs[0].legend(fontsize=6.5,ncol=2,loc='lower center',bbox_to_anchor=(.5,1.07))
    lo=q.a_lower_zero_mode.to_numpy();hi=q.a_upper_one_mode.to_numpy();mid=(lo+hi)/2
    axs[1].errorbar(q.epsilon,mid,yerr=np.vstack([mid-lo,hi-mid]),fmt='o',capsize=4,color=COLORS[0],label='Numerical bracket')
    axs[1].axhline(q.a0_star.iloc[0],color='black',ls='--',lw=1,label=r'$a_0^*$');axs[1].set(xlabel=r'$\varepsilon$',ylabel=r'$a_{\rm critical}/b$');axs[1].legend(fontsize=7)
    for i,ax in enumerate(axs):panel(ax,'ab'[i])
    export(fig,'critical_height',[f'{LOCK["runs"]["theorem_2_1"]}/{f}' for f in ['paper_a_sweep_M24.csv','critical_height.csv']],['All99 geometry classifications retained; four exact lower/upper deterministic brackets, no interpolation or confidence intervals'])


def kb_by_height():
    d=csv('theorem_2_1','paper_a_sweep_M24.csv').sort_values(['a','epsilon'])
    cutoff=float(csv('theorem_2_1','summary.csv').kb_cutoff.iloc[0])
    colors=['#0072B2','#D55E00','#009E73','#CC79A7','#E69F00','#000000']
    markers=['o','s','^','D','v','P']
    fig,ax=plt.subplots(figsize=(7.1,3.8),layout='constrained')
    grid,axs=plt.subplots(2,3,figsize=(7.1,5.2),layout='constrained')
    records=[];theory_records=[]
    for i,(height,q) in enumerate(d.groupby('a',sort=True)):
        valid=(q.status=='one') & q.mode_exists_supported & q.physical_checks_passed & q.diagnostic_complete & np.isfinite(q.kb)
        y=q.kb.where(valid,np.nan).to_numpy(float)
        if valid.any():
            ax.plot(q.epsilon,y,color=colors[i],marker=markers[i],label=fr'$a/b={height:g}$')
        target=axs.flat[i]
        target.plot(q.epsilon,y,color=colors[i],marker=markers[i],label='BEM')
        coefficient=float(q.asymptotic_coefficient.iloc[0])
        if coefficient>0:
            e=np.linspace(float(q.epsilon.min()),float(q.epsilon.max()),250)
            sigma=coefficient*e**2
            radicand=cutoff**2-sigma**2
            predicted=np.sqrt(np.where(radicand>0,radicand,np.nan))
            ax.plot(e,predicted,color=colors[i],ls='--',lw=1.15)
            target.plot(e,predicted,color=colors[i],ls='--',lw=1.15,label='Leading asymptotic')
            theory_records.extend({'a_over_b':height,'epsilon':float(x),'sigma_asymptotic':float(z),'kb_asymptotic':float(k)} for x,z,k in zip(e,sigma,predicted))
            target.legend(fontsize=6.2)
        target.axhline(cutoff,color='#777777',ls=':',lw=.9)
        target.set(title=fr'$a/b={height:g}$',xlabel=r'$\varepsilon$',ylabel=r'$kb$')
        panel(target,'abcdef'[i])
        if not valid.any():
            target.text(.5,.5,'Exclusion supported\nNo certified eigenvalue',ha='center',va='center',transform=target.transAxes,fontsize=8)
            target.set_ylim(cutoff-.01,cutoff+.001)
            target.set_xlim(0,float(q.epsilon.max())*1.04)
        for (_,row),v in zip(q.iterrows(),valid):
            records.append({'a_over_b':height,'epsilon':row.epsilon,'kb':row.kb if v else np.nan,'status':row.status,'plotted':bool(v)})
    ax.plot([],[],color='#555555',ls='--',label='Leading asymptotic (same color)')
    ax.axhline(cutoff,color='#777777',ls=':',label=r'First threshold $\pi/2$')
    ax.set(xlabel=r'$\varepsilon$',ylabel=r'$kb$');ax.legend(fontsize=7,ncol=2)
    source=[f'{LOCK["runs"]["theorem_2_1"]}/paper_a_sweep_M24.csv',f'{LOCK["runs"]["theorem_2_1"]}/summary.csv']
    transform=['All99 rows retained in exported figure data; only finite physically certified complete one-mode rows plotted; absent or ambiguous roots NaN, never zero; lines connect observed values only, no fit; b1 M24; dashed leading theory sigma=C(a)*epsilon², kb=sqrt((pi/2)²-sigma²), drawn only for positive C and real kb; dotted line is threshold']
    export(fig,'kb_vs_epsilon_by_a',source,transform)
    export(grid,'kb_vs_epsilon_each_a',source,transform+['Separate y scales by height to resolve near-threshold branches; excluded a.35 has no eigenvalue curve'])
    pd.DataFrame(records).to_csv(PAPER/'data/kb_vs_epsilon_by_a.csv',index=False)
    pd.DataFrame(theory_records).to_csv(PAPER/'data/kb_vs_epsilon_theory.csv',index=False)


def placement(S):
    s=S['theorem_Y'];e=vals(s,'epsilon');a=vals(s,'a_numerical');a1=float(s.a1_formula_2_7.iloc[0]);f=csv('theorem_Y','placement_a1_fit.csv');f=f[f.point_count==5].iloc[0]
    fig,axs=plt.subplots(1,2,figsize=(7.1,2.8),layout='constrained')
    axs[0].plot(e,a/e,color=COLORS[2],marker='^',label='Tuned BIC')
    axs[0].axhline(a1,color='black',ls='--',lw=1,label=r'$a_1$ from (2.7)')
    z=np.linspace(0,.10,80);axs[0].plot(z,f.a1_fit_linear+f.a2_fit_linear*z,color='#777777',ls=':',label=r'Linear fit, $\varepsilon\leq .10$')
    axs[0].set(xlabel=r'$\varepsilon$',ylabel=r'$a_{\rm BIC}/(b\varepsilon)$');axs[0].legend(fontsize=7)
    axs[1].plot(e,(a-e*a1)/e**2,color=COLORS[2],marker='^');axs[1].set(xlabel=r'$\varepsilon$',ylabel=r'$(a_{\rm BIC}/b-\varepsilon a_1)/\varepsilon^2$')
    for i,ax in enumerate(axs):panel(ax,'ab'[i])
    export(fig,'placement',[f'{LOCK["runs"]["theorem_Y"]}/{f}' for f in ['summary.csv','placement_a1_fit.csv']],['All ten tuned placements; direct a1 (not -beta/12); fitted five-point epsilon<=.10 intercept displayed separately, not a theorem'])


def physical(S):
    fig,axs=plt.subplots(2,2,figsize=(7.1,4.8),layout='constrained')
    for i,(k,s) in enumerate(S.items()):
        p=csv(k,'physical_diagnostics.csv');e=vals(s,'epsilon')
        axs[0,0].semilogy(e,s.normalized_residual,color=COLORS[i],marker=MARKERS[i],label=LABELS[k])
        axs[0,1].semilogy(p.epsilon,p.boundary_integral_relative_residual,color=COLORS[i],marker=MARKERS[i])
        axs[1,1].semilogy(p.epsilon,np.maximum(p.decay_relative_error_left,p.decay_relative_error_right),color=COLORS[i],marker=MARKERS[i])
        if i>0:axs[1,0].semilogy(p.epsilon,p.first_open_channel_relative_amplitude,color=COLORS[i],marker=MARKERS[i],label=LABELS[k])
    axs[0,0].axhline(1e-8,color='black',ls=':',label='Spectral tolerance')
    axs[0,0].set(ylabel=r'$\|Av\|/(\|A\|_2\|v\|)$',xlabel=r'$\varepsilon$');axs[0,0].legend(fontsize=6.8)
    axs[0,1].axhline(.005,color='black',ls=':',label='BIE tolerance');axs[0,1].set(ylabel='Off-grid BIE relative residual',xlabel=r'$\varepsilon$');axs[0,1].legend(fontsize=7)
    axs[1,0].axhline(1e-5,color='black',ls=':',label='Open-channel tolerance');axs[1,0].set(ylabel='First-channel relative amplitude',xlabel=r'$\varepsilon$');axs[1,0].legend(fontsize=6.8)
    axs[1,1].axhline(.2,color='black',ls=':',label='Decay tolerance');axs[1,1].set(ylabel='Max. relative decay-rate error',xlabel=r'$\varepsilon$');axs[1,1].legend(fontsize=7)
    for i,ax in enumerate(axs.flat):panel(ax,'abcd'[i])
    export(fig,'physical_diagnostics',[f'{v}/{f}' for v in LOCK['runs'].values() for f in ['summary.csv','physical_diagnostics.csv']],['All baseline points, separate physical dimensions/acceptance tolerances; first-channel amplitude applies to embedded cases only; max of left/right decay error'])


def sensitivity():
    fig,axs=plt.subplots(1,3,figsize=(7.1,2.8),layout='constrained')
    pars=['harmonic_order','lattice_terms','finite_difference_step'];names=[r'Harmonic order $H$',r'Lattice terms $N_{\rm lat}$',r'Difference step $h/b$']
    for i,k in enumerate(RUNS):
        s=csv(k,'internal_convergence.csv')
        for j,p in enumerate(pars):
            q=s[s.parameter==p].sort_values('value');col='relative_sigma_shift_from_baseline' if 'relative_sigma_shift_from_baseline' in q else 'relative_sigma_shift_from_reference'
            axs[j].plot(q.value,q[col],color=COLORS[i],marker=MARKERS[i],label=LABELS[k])
            if k=='theorem_Y':
                bad=q[~q.physical_bic_verified.astype(bool)];axs[j].scatter(bad.value,bad[col],marker='x',color='black',s=55,label='Y certificate rejected' if j==0 else None,zorder=4)
            axs[j].set_yscale('symlog',linthresh=1e-12);axs[j].set(xlabel=names[j])
            if j==2:axs[j].set_xscale('log')
    axs[0].set_ylabel(r'Relative shift in $\sigma$');axs[0].legend(fontsize=6)
    for j,ax in enumerate(axs):panel(ax,'abc'[j])
    export(fig,'sensitivity',[f'{v}/internal_convergence.csv' for v in LOCK['runs'].values()],['All15 variants per theorem, exact-zero shifts preserved with symlog; rejected Y H=5 and h=1e-5/1e-4 shown with crosses; fixed M24 internal sensitivity, not mesh convergence'])


def mfs_coefficients():
    fig,axs=plt.subplots(2,2,figsize=(7.1,4.7),layout='constrained')
    for i,k in enumerate(['theorem_X','theorem_Y']):
        d=csv(k,'mu_mfs_convergence.csv' if i==0 else 'mfs_convergence.csv')
        for j,(scale,q) in enumerate(d.groupby('source_scale')):
            q=q.sort_values('order');reference=float(q.mu.iloc[-1]);axs[0,i].plot(q.order,(q.mu-reference)/abs(reference),marker=['o','s','^'][j%3],color=COLORS[j%3],label=f'Source scale {scale:g}')
        axs[0,i].set_yscale('symlog',linthresh=1e-15);axs[0,i].set(xlabel='MFS collocation points',ylabel=r'$(\mu_N-\mu_{160})/|\mu_{160}|$',title=LABELS[k]);axs[0,i].legend(fontsize=7)
    d=csv('theorem_Y','a1_robustness.csv')
    for j,(scale,q) in enumerate(d.groupby('source_scale')):
        # Several quadrature variants may share order; retain individual points.
        axs[1,0].plot(q.quadrature_points,q.relative_difference_from_primary,marker=['o','s','^'][j%3],color=COLORS[j%3],label=f'Source scale {scale:g}')
    axs[1,0].set(xlabel='Boundary quadrature points',ylabel=r'Relative difference in $a_1$');axs[1,0].ticklabel_format(axis='y',style='sci',scilimits=(0,0));axs[1,0].legend(fontsize=7)
    q=csv('theorem_Y','a1_small_beta_diagnostic.csv').sort_values('beta')
    axs[1,1].plot(q.beta,q.a1_formula_2_7/q.beta,marker='o',color=COLORS[2],label='Direct (2.7)')
    axs[1,1].axhline(-.5,color='black',ls='--',label=r'$-1/2$ from substitution')
    axs[1,1].axhline(-1/12,color='#777777',ls=':',label=r'$-1/12$ printed Example 5.4')
    axs[1,1].set(xlabel=r'$\beta$',ylabel=r'$a_1/\beta$');axs[1,1].legend(fontsize=6.5)
    for i,ax in enumerate(axs.flat):panel(ax,'abcd'[i])
    export(fig,'mfs_coefficients',[f'{LOCK["runs"]["theorem_X"]}/mu_mfs_convergence.csv']+[f'{LOCK["runs"]["theorem_Y"]}/{n}' for n in ['mfs_convergence.csv','a1_robustness.csv','a1_small_beta_diagnostic.csv']],['MFS mu relative to order160 reference, exact zeros retained with symlog; a1 robustness at fixed order160 versus boundary quadrature, all source scales; small-beta discrepancy shown explicitly without replacing computed a1 by printed example'])


def read_config(mod,k):
    saved=json.loads((RUNS[k]/'run_manifest.json').read_text())['config'];allowed={f.name for f in dataclasses.fields(mod.Config)}
    kw={n:v for n,v in saved.items() if n in allowed}
    for n,v in list(kw.items()):
        if isinstance(v,list):kw[n]=tuple(v)
        elif n.endswith('dir'):kw[n]=Path(v)
    return mod.Config(**kw)


def modal(points,vector,geom,kb,b,terms=128):
    """Same separated-section double-layer modal representation as Y _modal_field.
    Source derivatives are analytic; near field retains current finite differences.
    """
    xi,eta,xp,yp=geom[:4];points=np.asarray(points,float)
    p=np.arange(1,terms+1)*np.pi/(2*b);k=kb/b;q=np.sqrt(p*p-k*k+0j);q=np.where(p<k,-1j*np.sqrt(k*k-p*p+0j),q)
    dx=points[:,0,None]-xi[None,:]
    if np.any(abs(dx)<b):raise ValueError('Modal section insufficiently separated')
    ex=np.exp(-abs(dx)[:,:,None]*q[None,None,:])
    source=(xp[:,None]*p[None,:]*np.cos((eta[:,None]+b)*p[None,:])/q[None,:]-np.sign(dx)[:,:,None]*yp[None,:,None]*np.sin((eta[:,None]+b)*p[None,:]))
    out=-(2*np.pi/len(vector))/(2*b)*np.sum(ex*source*vector[None,:,None],axis=1)*np.sin((points[:,1,None]+b)*p[None,:])
    tail=np.sum(abs(out[:,terms//2:]),axis=1);scale=np.maximum(np.sum(abs(out),axis=1),np.finfo(float).tiny)
    if np.any(tail>1e-12*scale):raise ValueError('Modal tail failed')
    return out.sum(axis=1)


def reconstruct(S):
    sys.path.insert(0,str(ROOT/'src/beyn_theorems_validations'))
    X=importlib.import_module('second_theorem_X');Y=importlib.import_module('second_theorem_Y')
    stamp=datetime.now(timezone.utc).strftime('%Y-%m-%d_%H-%M-%S')
    folder=PAPER/'data/field_reconstructions'/stamp;folder.mkdir(parents=True,exist_ok=False)
    meta={'stage':'in-progress','created_utc':stamp,'source_lock_sha256':sha(PAPER/'scripts/source_lock.json'),'cases':{}}
    (folder/'manifest.json').write_text(json.dumps(meta,indent=2)+'\n')
    for label,mod in [('theorem_X',X),('theorem_Y',Y)]:
        cfg=read_config(mod,label);r=S[label][np.isclose(S[label].epsilon,.10)].iloc[0];eps=float(r.epsilon);kb=float(r.kb_numerical);a=float(r.get('a_numerical',0));theta=mod.boundary_nodes(24)
        A=mod.assemble_full_matrix(kb,eps,24,cfg) if label=='theorem_X' else mod.assemble_full_matrix(kb,eps,a,24,cfg)
        _,ss,vh=np.linalg.svd(A);v=vh[-1].conj();v/=np.linalg.norm(v);v*=np.exp(-1j*np.angle(v[np.argmax(abs(v))]))
        res=float(np.linalg.norm(A@v)/(np.linalg.norm(A,2)*np.linalg.norm(v)))
        if res>1e-8:raise ValueError('Unrestricted saved-root singular vector failed')
        geom=mod.obstacle_geometry(theta,eps,cfg) if label=='theorem_X' else mod.obstacle_geometry(theta,eps,a,cfg)
        dense=np.linspace(0,2*np.pi,800,endpoint=False);boundary=mod.obstacle_geometry(dense,eps,cfg) if label=='theorem_X' else mod.obstacle_geometry(dense,eps,a,cfg)
        gx=np.linspace(-3,3,91);gy=np.linspace(-cfg.b,cfg.b,51);xx,yy=np.meshgrid(gx,gy);pts=np.column_stack([xx.ravel(),yy.ravel()])
        dist=np.min(np.hypot(pts[:,0,None]-boundary[0][None,:],pts[:,1,None]-boundary[1][None,:]),axis=1)
        mask=MplPath(np.column_stack(boundary[:2])).contains_points(pts)|(dist<.025*cfg.b)
        field=np.full(len(pts),np.nan+1j*np.nan,dtype=complex)
        valid=np.flatnonzero(~mask);far=valid[np.min(abs(pts[valid,0,None]-geom[0][None,:]),axis=1)>=cfg.b];near=np.setdiff1d(valid,far)
        for chunk in np.array_split(far,max(1,int(np.ceil(len(far)/250)))):
            if len(chunk):field[chunk]=modal(pts[chunk],v,geom,kb,cfg.b)
        # Existing Y batched off-grid kernel depends only on source coordinates,
        # normals and Green parameters; it is also applicable to X geometry.
        for chunk in np.array_split(near,max(1,int(np.ceil(len(near)/80)))):
            if len(chunk):
                kernels=Y._offgrid_kernel_array(pts[chunk,0],pts[chunk,1],*geom[:4],kb,cfg)
                field[chunk]=(2*np.pi/24)*(kernels@v)
        # Verify batching against the current independent scalar reconstruction.
        picks=near[np.linspace(0,len(near)-1,7).astype(int)];G=mod.make_extended_field_green(kb,cfg)
        direct=(mod.reconstruct_field_at_points(pts[picks].tolist(),v,theta,eps,cfg,G) if label=='theorem_X' else mod.reconstruct_field_at_points(pts[picks].tolist(),v,theta,eps,a,cfg,G))
        check=float(np.max(abs(direct-field[picks]))/max(np.max(abs(direct)),1e-30))
        if check>1e-6:raise ValueError(f'Batched field disagrees with scalar: {check}')
        yn,yw=np.polynomial.legendre.leggauss(cfg.physical_quadrature_points);yn*=cfg.b;yw*=cfg.b
        # Match production sections; extra distances extend visualization over
        # meaningful decay lengths without changing the eigenvalue or vector.
        production=np.asarray(cfg.physical_x_over_b)*cfg.b
        distances=np.unique(np.r_[production,np.linspace(production[-1],production[-1]+4*cfg.b/float(r.sigma_numerical),18)])
        sections=np.zeros((2,len(distances),len(yn)),complex)
        for j,sgn in enumerate([-1,1]):
            for n,d in enumerate(distances):sections[j,n]=modal(np.column_stack([np.full(len(yn),sgn*d),yn]),v,geom,kb,cfg.b)
        norms=np.sqrt(np.sum(abs(sections)**2*yw,axis=2));phi=np.cos(np.pi*yn/(2*cfg.b));channel=np.sum(sections*phi*yw,axis=2);rel=abs(channel)/(np.sqrt(np.sum(phi**2*yw))*norms)
        prod_ix=np.flatnonzero(np.isin(distances,production));rates=[-np.polyfit(production,np.log(norms[i,prod_ix]),1)[0] for i in range(2)]
        stored=csv(label,'physical_diagnostics.csv');stored=stored[np.isclose(stored.epsilon,.1)].iloc[0]
        parity=float(np.linalg.norm(v+v[::-1])/np.linalg.norm(v)) if label=='theorem_X' else float(np.linalg.norm(v-v[::-1])/np.linalg.norm(v))
        details={'epsilon':eps,'kb':kb,'a':a,'sigma':float(r.sigma_numerical),'M':24,'normalized_full_residual':res,'scalar_batch_relative_discrepancy':check,'decay_rates':rates,'stored_decay_rates':[float(stored.decay_rate_left),float(stored.decay_rate_right)],'max_production_channel_ratio':float(np.max(rel[:,prod_ix])),'boundary_parity_relative_residual':parity,'normalization':'unit boundary-vector Euclidean norm; max component phase positive real','near_boundary_exclusion_width_over_b':.025,'modal_terms':128,'sources':[str(RUNS[label].relative_to(ROOT)/'summary.csv')],'method':'saved-root unrestricted full-matrix SVD; no eigenvalue refinement; existing centered-difference lattice kernel near source, analytic modal separated sections'}
        if max(abs(np.asarray(rates)-float(r.sigma_numerical)))/float(r.sigma_numerical)>.001:raise ValueError('Reconstructed decay failed')
        np.savez_compressed(folder/f'{label}.npz',x=gx,y=gy,field=field.reshape(xx.shape),mask=mask.reshape(xx.shape),theta=theta,boundary_vector=v,boundary_x=boundary[0],boundary_y=boundary[1],section_y=yn,section_weights=yw,distances=distances,section_fields=sections,section_norms=norms,first_channel_coefficients=channel,first_channel_relative_amplitudes=rel)
        details['npz_sha256']=sha(folder/f'{label}.npz');meta['cases'][label]=details
        (folder/'manifest.json').write_text(json.dumps(meta,indent=2)+'\n');print('Field reconstruction:',label,details,flush=True)
    meta['stage']='completed';(folder/'manifest.json').write_text(json.dumps(meta,indent=2)+'\n')
    (PAPER/'data/field_reconstructions/latest.json').write_text(json.dumps({'directory':str(folder.relative_to(PAPER)),'manifest_sha256':sha(folder/'manifest.json')},indent=2)+'\n')
    return folder


def fields(folder):
    meta=json.loads((folder/'manifest.json').read_text());fig,axs=plt.subplots(2,2,figsize=(7.1,4.9),layout='constrained');section_rows=[]
    for i,k in enumerate(['theorem_X','theorem_Y']):
        d=np.load(folder/f'{k}.npz');detail=meta['cases'][k]
        if sha(folder/f'{k}.npz')!=detail['npz_sha256']:raise ValueError('Field data changed')
        z=d['field'];scale=float(np.nanmax(abs(z)));re=z.real/scale
        im=axs[0,i].pcolormesh(d['x'],d['y'],np.ma.array(re,mask=d['mask']),cmap='RdBu_r',vmin=-1,vmax=1,shading='nearest',rasterized=True)
        axs[0,i].fill(d['boundary_x'],d['boundary_y'],color='#777777');axs[0,i].set(xlabel=r'$x/b$',ylabel=r'$y/b$',title=LABELS[k]);axs[0,i].set_aspect('equal');fig.colorbar(im,ax=axs[0,i],label=r'$\Re(u)/\max|u|$',shrink=.6)
        dist=d['distances'];norm=d['section_norms'];s=detail['sigma'];scaled=(dist-dist[0])*s
        for sign in range(2):
            for n,distance in enumerate(dist):
                channel=d['first_channel_coefficients'][sign,n];section_rows.append({'theorem':k,'side':['left','right'][sign],'distance_over_b':float(distance),'scaled_distance':float(scaled[n]),'section_L2_norm':float(norm[sign,n]),'normalized_section_L2_norm':float(norm[sign,n]/norm[sign,0]),'first_channel_real':float(channel.real),'first_channel_imag':float(channel.imag),'first_channel_relative_amplitude':float(d['first_channel_relative_amplitudes'][sign,n])})
        axs[1,i].semilogy(scaled,norm[0]/norm[0,0],color=COLORS[0],marker='o',ms=3,label='Left sections')
        axs[1,i].semilogy(scaled,norm[1]/norm[1,0],color=COLORS[1],marker='s',ms=3,ls='--',label='Right sections')
        axs[1,i].semilogy(scaled,np.exp(-scaled),color='black',ls=':',label=r'$e^{-\sigma (|x|-x_0)}$')
        axs[1,i].set(xlabel=r'$\sigma(|x|-x_0)$',ylabel=r'$\|u(x,\cdot)\|_2/\|u(x_0,\cdot)\|_2$');axs[1,i].legend(fontsize=6.5)
        panel(axs[0,i],'ab'[i]);panel(axs[1,i],'cd'[i])
    pd.DataFrame(section_rows).to_csv(PAPER/'data/field_sections.csv',index=False)
    export(fig,'fields_decay',[str((folder/'manifest.json').relative_to(ROOT))]+[str((folder/f'{k}.npz').relative_to(ROOT)) for k in ['theorem_X','theorem_Y']],['Actual saved-root unrestricted full BEM field; independently phased each case and normalized by its own peak; masked interior and .025b near-boundary band; measured L2 cross-section norms normalized separately at first section, plotted in scaled decay distance'])


def main():
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--reconstruct-fields',action='store_true');args=ap.parse_args()
    check_sources()
    for d in ['figures','tables','data/source']: (PAPER/d).mkdir(parents=True,exist_ok=True)
    for k,p in RUNS.items():
        dest=PAPER/'data/source'/k;dest.mkdir(exist_ok=True)
        for f in p.glob('*.csv'):shutil.copyfile(f,dest/f.name)
        shutil.copyfile(p/'run_manifest.json',dest/'run_manifest.json')
    S={k:csv(k,'summary.csv') for k in RUNS}
    make_tables(S);geometry(S);spectral(S);beyn();critical();kb_by_height();placement(S);physical(S);sensitivity();mfs_coefficients()
    field_folder=None
    if args.reconstruct_fields:field_folder=reconstruct(S)
    else:
        pointer=PAPER/'data/field_reconstructions/latest.json'
        if pointer.exists():
            entry=json.loads(pointer.read_text());field_folder=PAPER/entry['directory']
            if sha(field_folder/'manifest.json')!=entry['manifest_sha256']:raise ValueError('Field manifest changed')
            m=json.loads((field_folder/'manifest.json').read_text())
            if m['source_lock_sha256']!=sha(PAPER/'scripts/source_lock.json') or m['stage']!='completed':raise ValueError('Field reconstruction stale or incomplete')
    if field_folder:fields(field_folder)
    manifest={'schema':1,'environment':{'python':platform.python_version(),'numpy':np.__version__,'pandas':pd.__version__,'matplotlib':matplotlib.__version__,'platform':platform.platform()},'generated_utc':datetime.now(timezone.utc).isoformat(),'generator_sha256':sha(Path(__file__)),'source_lock':LOCK,'artifacts':ARTIFACTS,'tables':[str(p.relative_to(PAPER)) for p in (PAPER/'tables').glob('*.tex')],'field_reconstruction':str(field_folder.relative_to(PAPER)) if field_folder else None,'scope':'fixed M24; finite contours; X count only odd-y sector; no Y global uniqueness; physical tolerance diagnostics not exact PDE certification','publisher_requirements':'provisional general journal figures; exact target journal unspecified'}
    outputs=list((PAPER/'figures').glob('*.pdf'))+list((PAPER/'figures').glob('*.png'))+list((PAPER/'tables').glob('*.tex'))
    manifest['output_sha256']={str(p.relative_to(PAPER)):sha(p) for p in outputs}
    (PAPER/'data/publication_manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
    print('Generated',len(ARTIFACTS),'figures and',len(manifest['tables']),'tables/macros',flush=True)

if __name__=='__main__':main()
