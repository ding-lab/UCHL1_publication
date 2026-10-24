"""
Figure 2d: UCHL1-KO reduces proliferation (Caki-1 CCK-8 assay)
"""
from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import rcParams
import statsmodels.api as sm
from statsmodels.formula.api import ols

WORKBOOK = Path('./ccRCC_UCHL1_manuscript/Data/caki_NC_C2_C9_proliferation_data_Figure2c.xlsx') ### REVIEW/EDIT
OUT_PDF  = Path('./ccRCC_UCHL1_manuscript/Data/Figure2d_Caki_Proliferation.pdf')

# [Keep identical configuration and functions (sig, draw_bracket, etc.) from fig2d.py]
CELL_LINES = ['Caki-1 NC (parental)', 'Caki-1 C2 (KO)', 'Caki-1 C9 (KO)']
COLORS = {'Caki-1 NC (parental)':'#1f4ecf', 'Caki-1 C2 (KO)':'#d62728', 'Caki-1 C9 (KO)':'#2ca02c'}
MARKERS = {'Caki-1 NC (parental)':'o', 'Caki-1 C2 (KO)':'s', 'Caki-1 C9 (KO)':'^'}
TIMEPOINTS = [0, 24, 48, 72, 96]
REFERENCE = 'Caki-1 NC (parental)'

def sig(p):
    if p < 0.0001: return '****'
    if p < 0.001: return '***'
    if p < 0.01: return '**'
    if p < 0.05: return '*'
    return 'ns'

def main():
    df = pd.read_excel(WORKBOOK, sheet_name='Corrected_data')
    summary = df.groupby(['Cell line', 'Time (h)'])['Corrected A450'].agg(['mean', 'std']).reset_index()

    pvals = {}
    for line in CELL_LINES:
        if line == REFERENCE: continue
        sub = df[df['Cell line'].isin([REFERENCE, line])].copy()
        sub = sub.rename(columns={'Cell line': 'cell_line', 'Time (h)': 'time', 'Corrected A450': 'A450'})
        model = ols('A450 ~ C(cell_line) * C(time)', data=sub).fit()
        anova = sm.stats.anova_lm(model, typ=2)
        pvals[line] = float(anova.loc['C(cell_line):C(time)', 'PR(>F)'])

    rcParams.update({'font.family':'DejaVu Sans', 'font.size':11, 'axes.spines.top':False, 'axes.spines.right':False})
    fig, ax = plt.subplots(figsize=(7.0, 5.0))

    for line in CELL_LINES:
        sub = summary[summary['Cell line']==line].sort_values('Time (h)')
        ax.errorbar(sub['Time (h)'], sub['mean'], yerr=sub['std'], marker=MARKERS[line], color=COLORS[line],
                    markersize=6, capsize=5, label=line, zorder=3)

    ax.set_xlabel('Time (hours)'); ax.set_ylabel('Absorbance at 450 nm (blank-subtracted)')
    ax.set_xticks(TIMEPOINTS)
    ax.legend(loc='upper left', frameon=False)
    plt.savefig(OUT_PDF, bbox_inches='tight')
    plt.close()

if __name__ == '__main__':
    main()