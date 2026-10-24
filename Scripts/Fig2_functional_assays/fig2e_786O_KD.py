"""
Figure 2e: UCHL1 shRNA knockdown in 786-O (CCK-8 & qRT-PCR)
"""
from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt
from scipy import stats
import statsmodels.api as sm
from statsmodels.formula.api import ols

WORKBOOK_CCK8 = Path('./ccRCC_UCHL1_manuscript/Data/786O_KD_proliferation.xlsx') ### REVIEW/EDIT
WORKBOOK_PCR = Path('./ccRCC_UCHL1_manuscript/Data/786O_KD_qPCR.xlsx') ### REVIEW/EDIT
OUT_PDF = Path('./ccRCC_UCHL1_manuscript/Data/Figure2e_786O_KD.pdf')

def main():
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))
    
    # 1. Left Panel: Proliferation (Similar to 2d logic)
    df_cck8 = pd.read_excel(WORKBOOK_CCK8)
    summary = df_cck8.groupby(['Group', 'Time'])['A450'].agg(['mean', 'std']).reset_index()
    
    colors = {'NC': 'blue', 'KD1': 'red', 'KD2': 'green'}
    for grp in ['NC', 'KD1', 'KD2']:
        sub = summary[summary['Group'] == grp]
        ax1.errorbar(sub['Time'], sub['mean'], yerr=sub['std'], label=grp, 
                     color=colors[grp], marker='o', capsize=5)
    
    ax1.set_xlabel('Time (hours)')
    ax1.set_ylabel('Absorbance at 450 nm (blank-subtracted)')
    ax1.set_title('786-O Proliferation')
    ax1.legend(frameon=False)
    ax1.spines['top'].set_visible(False)
    ax1.spines['right'].set_visible(False)
    
    # 2. Right Panel: qRT-PCR (Welch's t-test)
    df_pcr = pd.read_excel(WORKBOOK_PCR)
    pcr_sum = df_pcr.groupby('Group')['Relative_Expression'].agg(['mean', 'std']).reindex(['NC', 'KD1', 'KD2'])
    
    bars = ax2.bar(pcr_sum.index, pcr_sum['mean'], yerr=pcr_sum['std'], capsize=5, 
                   color=['#1f4ecf', '#d62728', '#2ca02c'], edgecolor='black')
    
    # Run Welch's t-test
    nc_vals = df_pcr[df_pcr['Group'] == 'NC']['Relative_Expression']
    for idx, kd in enumerate(['KD1', 'KD2'], start=1):
        kd_vals = df_pcr[df_pcr['Group'] == kd]['Relative_Expression']
        t, p = stats.ttest_ind(nc_vals, kd_vals, equal_var=False)
        sig = '***' if p < 0.001 else ('**' if p < 0.01 else '')
        ax2.text(idx, pcr_sum.loc[kd, 'mean'] + pcr_sum.loc[kd, 'std'] + 0.1, sig, ha='center', fontweight='bold')

    ax2.set_ylabel('Relative UCHL1 Expression')
    ax2.set_title('qRT-PCR Validation')
    ax2.spines['top'].set_visible(False)
    ax2.spines['right'].set_visible(False)
    
    plt.tight_layout()
    plt.savefig(OUT_PDF, bbox_inches='tight')

if __name__ == '__main__':
    main()