from matplotlib.pyplot import rc_context
import matplotlib.pyplot as plt
import seaborn as sns
import scanpy as sc
import numpy as np
import pandas as pd


cols_Cell_Line = {"WTC":"#2088C9", "H9":"#FFC502", "WIBJ2":"#332288",
                  "H1":"#44AA99", "ESC":"lightgray","iPSC":"lightgray", "HES3":"#BC0F4D",
                 '?': 'lightgray'} # WTC-HES3-H9


cols_SingleMorph = {"Ctrl":"#000000", 
                    'Control':"#000000", 
                     "SHH":   "#F8A220", 
                     "RA": "#D9E70D", 
                     "FGF8": "#883EB6", 
                     "CHIR": "#098684",
                 
                     "BMP4": "#E6376C", # BMP4
                     "BMP7": "#A34030", # BMP7
                     
                     
                     "CycA": "#1B3BB1", # CycA
                     "XAV939":"#57FFF9",
                       "?":"lightgray" }

Morphogen_full = {"Ctrl":"#000000",
                  'SHH_A':"#FFDDAD",
                     'SHH_E':"#C66611",
                     'RA_A':"#E9F17E",
                     'RA_E':"#686F0C",
                     'CHIR_A':"#D8E9E9",
                     'CHIR_E':"#1B7472",
                     'FGF8_late_A':"#C7AED6",
                     'FGF8_late_E':"#8452A3",
                     'CHIR_tA':"#aadce0", 
                     'CHIR_tC':"#528fad",
                     'CHIR_tE':"#1e466e"}

Morphogen_full2 = {
    "Ctrl": "#000000",
    "SHH_A": "#FFC370",
    "SHH_E": "#FF990A",
    "RA_A": "#E0EB47",
    "RA_E": "#9BA512",
    "CHIR_A": "#B0D4D4",
    "CHIR_E": "#569F9E",
    "FGF8_late_A": "#C7AED6",
    "FGF8_late_E": "#8452A3",
    "CHIR_tA": "#FFB3C9",
    "CHIR_tC": "#EC5B86",
    "CHIR_tE": "#C7184C"
}


Experiment = {'time':"#574571",
                'conc':"#dec5da",
                'comb': "#7fa074"}
                  

          
cols_Batch = {'B1':"#88CCEE",
                'B2':'#DDCC77'}


Medium = {'NIM':'#117733',
          'NPM':'#882255'}

cell_type_colors = {
    'Telencephalic Progenitors': '#e0762f',
    'Spinal Cord Progenitors': '#a66d9b',
    'Hindbrain Progenitors': '#e78e97',
    'Neuroectoderm/Neuroepithelium': '#f4c617',
    'Hypothalamic Progenitors': '#e83d63',
    'CNS Neurons': '#3a5a89',
    'Cortical Hem/Diencephalon': '#7A0000',
    'PNS Neurons': '#2983aa',
    'Floor Plate': '#8d4c88',
    'Midbrain Progenitors': '#30d8c8',
    'Retinal Progenitors': '#ac383f',
    'Mesenchyme': '#92c051',
    'Neural Crest': '#474a82',
    'Extraembryonic Tissue': '#A5BC7B',
    'Endoderm/mesoderm-derived tissues': '#67BD1B',
    'Non-neurectodermal tissues': '#92c051',
    'Hindbrain_NPM?': '#f3c2c7',
    "Hindbrain Progenitors" :"#f3c2c7"

    
}




neuron_colors = {
    'MGE-derived GABA Neurons': '#e69138',
    'PNS GLY Neurons': '#d0c7ff',
    'Hypothalamic Neurons SIM1-': '#ff8cbf',
    'Hypothalamic Neurons SIM1+': '#f30d73',
    'Cortical Neurons': '#ff6d01',
    'PNS CHOL Neurons': '#a66d9b',
    'PNS Neurons': '#2983aa',
    'Retinal Neurons': '#ac383f',
    'MGE-derived GABA Immature Neurons': '#f6b26b'
}

palette_OSMGT = {'Pluripotent Stem Cells':"#D9CB96", # iPSCs
                         'Neuroectoderm': "#FBE183", # Neuroectoderm
                         'Neuroepithelium': "#F4C617", # Neuroepithelium
                         'Telencephalic Progenitors': "#E0762F", # Telencephalic Progenitors
                         'Diencephalic Progenitors': "#7A0000", # Diencephalic
                         'Retinal Progenitors': "#AC383F", # Retinal
                           'Hypothalamic Progenitors':"#E83D63", # Hypothalamic
                         'Hindbrain Progenitors': "#E78E97", # Hindbrain
                          'Spinal Cord Progenitors':"#A66D9B", # SC
                          'Posterior Floor Plate':"#8D4C88", # FOXA2
                          'Neural Crest': "#543F7A", # NC
                         'PNS Progenitors': "#20B6D4", # PNS progenitors
                         'PNS Neurons': "#2983AA", # PNS neurons
                         'CNS Neurons': "#3A5A89", # CNS Neurons
                         'Non-Neural Ectoderm': "#008C8D", 
                'Mesenchymal Cells': "#C4E255", 
                 'Endoderm-derived Tissue':"#67BD1B", 
                 'Extraembryonic Tissue':"#A5BC7B"}
 





def plot_separately(adata, col_name, save=False, s=1, palette=None, basis='X_umap'):

    
    # Save is by default False, otherwise - pathname to save fig
    # s - dot size
    # Get unique categories
    unique_categories = np.array(adata.obs[col_name].unique().tolist()).astype(str)
    unique_categories = [x for x in unique_categories if x != 'nan']
    unique_categories = [x for x in unique_categories if x != 'No']

    nrows = int(np.ceil(len(unique_categories) / 3))
    
    # Set up the subplots
    fig, axes = plt.subplots(nrows, 3, figsize=(15, 5 * nrows))
    
    # Loop through each category and create a UMAP plot
    for i, category in enumerate(unique_categories):
        if nrows == 1:
            ax = axes[i]
        else:
            colpos = i % 3
            rowpos = i // 3
            ax = axes[rowpos, colpos]
        
        # Filter data for the current category
        category_data = adata[adata.obs[col_name] == category, :]
        
        # Plot all other cells in grey
        other_data = adata[adata.obs[col_name] != category, :]
        ax.scatter(other_data.obsm[basis][:, 0], other_data.obsm[basis][:, 1], 
                   color='lightgrey', s=s, alpha=0.5)
        
        # Plot the UMAP for the selected category in color
        if palette is None:
            sc.pl.embedding(category_data, basis = basis,color=col_name, ax=ax, show=False, 
                       frameon=False, size=s + 1, legend_loc=None)
        else:
            
            #category_color = palette.get(category, 'black')
            #print(category_color)
            sc.pl.embedding(category_data, basis = basis, ax=ax, show=False, 
                       frameon=False, size=s + 1, legend_loc=None, palette=palette, color=col_name)
        
        ax.set_title(f'{category}')
        
        # Turn off grid and axis for empty subplots
        if not category_data.shape[0]:
            ax.set_axis_off()

    # Adjust layout
    plt.tight_layout()
    
    # Save the plot if a save path is provided
    if save:
        plt.rcParams['pdf.fonttype'] = 42
        plt.savefig(save, dpi=300)

    plt.show()








def plot_stacked_barplot(celltable, conditions,condition_column,column_celltype,
                         cell_types_order,palette, 
                         cond_order,fig_name=None):
    fig, axs = plt.subplots(2, 2, figsize=(40, 20))  # Adjust figsize to your needs
    for i, ax in enumerate(axs.flatten()):
       
        sns.set_style("white")
        sns.set_style('ticks')
        celltable_subset = celltable[celltable[condition_column]==conditions[i]]
       
        dfc = pd.crosstab(celltable_subset.morphogen_full, celltable_subset[column_celltype], normalize='index').mul(100).round(100)
        
        try:
            dfc = dfc.loc[cond_order]
        except:
            pass
    
        # Unique cell types for ordering
        celltype_set = celltable_subset[column_celltype].unique().tolist()
        cell_types_order_pic = [x for x in cell_types_order if x in celltype_set]
        dfc = dfc.loc[:, cell_types_order_pic]
    
        legend_order = [i for i in range(len(cell_types_order_pic))]
        legend_order.reverse()
        
        
        dfc.plot(
            kind='bar', 
            stacked=True, 
            rot=0, 
            ax=ax, 
            color=palette, 
            width=0.95, 
            ylabel='Percent(%)'
        )
        
        handles, labels = ax.get_legend_handles_labels()
        ax.legend(
            [handles[i] for i in legend_order], [labels[i] for i in legend_order],
            loc='center left', 
            ncol=1, 
            fancybox=False, 
            shadow=False, 
            fontsize=16, 
            bbox_to_anchor=(1.05, 0.5)
        )
        ax.set_xticks(ax.get_xticks())
        #ax.set_xticklabels(dfc_subset.index, weight='bold')
        ax.set_title(conditions[i], weight='bold', fontsize=30)
        ax.tick_params(axis='both',labelsize=20)
        sns.despine(ax=ax, offset=1, trim=True, bottom=True)
        ax.tick_params(axis='x', length=0)
    
    # Adjust layout and show the plot
    plt.subplots_adjust(right=2)  
    plt.tight_layout()
   

    
    if fig_name:
        plt.rcParams['pdf.fonttype'] = 42
        plt.savefig(fig_name, dpi=300, format='pdf', bbox_inches='tight')
    plt.show()


