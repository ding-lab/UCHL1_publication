"""
Figure 2f: Inhibitor Dose-Response (24h)
"""
from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt

WORKBOOK = Path('./ccRCC_UCHL1_manuscript/Data/inhibitor_dose_response.xlsx') ### REVIEW/EDIT
OUT_PDF = Path('./ccRCC_UCHL1_manuscript/Data/Figure2f_Inhibitors.pdf')

def main():
    df = pd.read_excel(WORKBOOK)
    summary = df.groupby(['CellLine', 'Inhibitor', 'Dose_uM'])['Viability'].agg(['mean', 'std']).reset_index()
    
    fig, axes = plt.subplots(1, 2, figsize=(10, 4), sharey=True)
    doses = [0.1, 1, 5, 10, 50, 100]
    
    for ax, cell in zip(axes, ['Caki-1', '786-O']):
        sub_cell = summary[summary['CellLine'] == cell]
        for inh, col in zip(['6RK73', 'IMP-1710', 'LDN-57444'], ['#a6cee3', '#e31a1c', '#b2df8a']):
            sub_inh = sub_cell[sub_cell['Inhibitor'] == inh]
            ax.errorbar(sub_inh['Dose_uM'], sub_inh['mean'], yerr=sub_inh['std'], 
                        label=inh, color=col, marker='o', capsize=4)
            
        ax.set_xscale('log')
        ax.set_xticks(doses)
        ax.set_xticklabels(doses)
        ax.set_xlabel('Concentration (µM)')
        ax.set_title(cell)
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        
    axes[0].set_ylabel('Cell Viability (%)')
    axes[-1].legend(frameon=False)
    plt.tight_layout()
    plt.savefig(OUT_PDF, bbox_inches='tight')

if __name__ == '__main__':
    main()