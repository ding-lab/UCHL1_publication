"""
Figure 2g: IMP-1710 Timecourse (96h)
"""
from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt

WORKBOOK = Path('./ccRCC_UCHL1_manuscript/Data/IMP1710_timecourse.xlsx') ### REVIEW/EDIT
OUT_PDF = Path('./ccRCC_UCHL1_manuscript/Data/Figure2g_IMP1710_Timecourse.pdf')

def main():
    df = pd.read_excel(WORKBOOK)
    summary = df.groupby(['CellLine', 'Treatment', 'Time'])['Viability'].agg(['mean', 'std']).reset_index()
    
    fig, axes = plt.subplots(1, 2, figsize=(10, 4), sharey=True)
    
    for ax, cell in zip(axes, ['Caki-1', '786-O']):
        sub_cell = summary[summary['CellLine'] == cell]
        for tx, col in zip(['DMSO', 'IMP-1710'], ['black', 'red']):
            sub_tx = sub_cell[sub_cell['Treatment'] == tx]
            ax.errorbar(sub_tx['Time'], sub_tx['mean'], yerr=sub_tx['std'], 
                        label=tx, color=col, marker='s', capsize=4)
            
        ax.set_xlabel('Time (hours)')
        ax.set_xticks([0, 24, 48, 72, 96])
        ax.set_title(cell)
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        
    axes[0].set_ylabel('Cell Viability (%)')
    axes[-1].legend(frameon=False)
    plt.tight_layout()
    plt.savefig(OUT_PDF, bbox_inches='tight')

if __name__ == '__main__':
    main()