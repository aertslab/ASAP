## Topic modeling evaluation: 
print('importing modules...')
import pandas as pd
import argparse
import pickle
import os
import re
import glob
import sys
import matplotlib.pyplot as plt
import plotnine as pn
from plotnine import *
from pycisTopic.cistopic_class import *
from pycisTopic.clust_vis import *
from pycisTopic.lda_models import *
import pyranges
import gzip
print('finished importing modules')

print()

print('functions...')
def make_argument_parser():
    """
    Creates an ArgumentParser to read the options for this script from
    sys.argv
    """
    parser = argparse.ArgumentParser(
        description=" Evaluate topic models.",)
    parser.add_argument('--input_cto', '-i', type=str, required=True,
                        help='Path to cisTopic object pickle file.')
    
    parser.add_argument('--topic_model_dir', '-tmd', type=str, required=True,
                        help='Path topic model dir.')

    parser.add_argument('--selected_model', '-m', type=int, required=True,
                        help='Number of topics.')

    parser.add_argument('--output_dir', '-o', type=str, required=True,
                        help='Path to output dir.')
    
    return parser


def load_and_add_topic_model(model_dir, model, cistopic_obj):
    """
    Load and add topic model
    """

    model_name = f'Topic{model}.pkl.gz'
    print(model_name)
    selected_model = pickle.load(gzip.open(os.path.join(model_dir, model_name),"rb"))
    cistopic_obj.add_LDA_model(selected_model)

def clustering_and_dim_red(cistopic_obj, model):
    """
    Finding clusters and running UMAP
    """
    ntopics = str(model)
    find_clusters(cistopic_obj,
              target  = 'cell',
              k = 10,
              res = [0.6, 1.2, 2.0, 2.5, 3.0, 3.5],
              prefix = f'{ntopics}_pycisTopic_',
              scale = True,
              split_pattern = '-')
    print('Finding clusters finished')
    print()
    
    print('Run UMAP.....')
    run_umap(cistopic_obj, target  = 'cell', scale=True, reduction_name = f'UMAP_{ntopics}')
    print('UMAP done!')
    print('Run tSNE.....')
    run_tsne(cistopic_obj, target  = 'cell', scale=True, reduction_name = f'tSNE_{ntopics}')
    print('tSNE done!')
    print()

#function that adds label at desired location per cell type using a label_dict
def plot_metadata_save(
    cistopic_obj: CistopicObject,
    reduction_name: str,
    variables: list[str],
    target: str = "cell",
    remove_nan: bool = True,
    show_label: bool = True,
    show_legend: bool = False,
    cmap: Optional[Union[str, "matplotlib.cm"]] = cm.viridis,
    dot_size: int = 10,
    text_size: int = 10,
    alpha: float = 1.0,
    seed: int = 555,
    color_dictionary: Optional[Dict[str, str]] = {},
    figsize: tuple[float, float] = (6.4, 4.8),
    num_columns: int = 1,
    selected_features: Optional[List[str]] = None,
    save: Optional[str] = None
):
    """
    Plot categorical and continuous metadata into dimensionality reduction.

    Parameters
    ----------
    cistopic_obj: `class::CistopicObject`
            A cisTopic object with dimensionality reductions in `class::CistopicObject.projections`.
    reduction_name: str
            Name of the dimensionality reduction to use
    variables: list[str]
            List of variables to plot. They should be included in `class::CistopicObject.cell_data` and `class::CistopicObject.region_data`, depending on which
            target is specified.
    target: str
            Whether cells ('cell') or regions ('region') should be used. Default: 'cell'
    remove_nan: bool
            Whether to remove data points for which the variable value is 'nan'. Default: True
    show_label: bool
            For categorical variables, whether to show the label in the plot. Default: True
    show_legend: bool
            For categorical variables, whether to show the legend next to the plot. Default: False
    cmap: str or 'matplotlib.cm'
            For continuous variables, color map to use for the legend color bar. Default: cm.viridis
    dot_size: int
            Dot size in the plot. Default: 10
    text_size: int
            For categorical variables and if show_label is True, size of the labels in the plot. Default: 10
    alpha: float
            Transparency value for the dots in the plot. Default: 1
    seed: int
            Random seed used to select random colors. Default: 555
    color_dictionary: dict, optional
            A dictionary containing an entry per variable, whose values are dictionaries with variable levels as keys and corresponding colors as values.
            Default: None
    figsize: tuple[float, float], optional
            Size of the figure. If num_columns is 1, this is the size for each figure; if num_columns is above 1, this is the overall size of the figure (if keeping
            default, it will be the size of each subplot in the figure). Default: (6.4, 4.8)
    num_columns: int
            For multiplot figures, indicates the number of columns (the number of rows will be automatically determined based on the number of plots). Default: 1
    selected_features: list,[str] optional
            A list with selected features (cells or regions) to plot. This is recommended when working with regions (e.g. selecting
            regions in binarized topics), as working with all regions can be time consuming. Default: None (use all features)
    save: str, optional
            Path to save plot. Default: None.

    """
    if target == "cell":
        data_mat = cistopic_obj.cell_data
    if target == "region":
        data_mat = cistopic_obj.region_data

    embedding = cistopic_obj.projections[target][reduction_name]

    if selected_features is not None:
        data_mat = data_mat.loc[selected_features]
        embedding = embedding.loc[selected_features]

    data_mat = data_mat.loc[embedding.index.to_list()]
    pdf = None
    if (save is not None) and (num_columns == 1):
        pdf = matplotlib.backends.backend_pdf.PdfPages(save)

    if num_columns > 1:
        num_rows = int(np.ceil(len(variables) / num_columns))
        if figsize == (6.4, 4.8):
            figsize = (6.4 * num_columns, 4.8 * num_rows)
        i = 1

    fig = plt.figure(figsize=figsize)

    for var in variables:
        var_data = data_mat.copy().loc[:, var].dropna().to_list()
        if isinstance(var_data[0], str):
            if remove_nan and (data_mat[var].isnull().sum() > 0):
                var_data = data_mat.copy().loc[:, var].dropna().to_list()
                emb_nan = embedding.loc[
                    data_mat.copy().loc[:, var].dropna().index.tolist()
                ]
                label_pd = pd.concat(
                    [emb_nan, data_mat.loc[:, [var]].dropna()], axis=1, sort=False
                )
            else:
                var_data = (
                    data_mat.copy().astype(str).fillna("NA").loc[:, var].to_list()
                )
                label_pd = pd.concat(
                    [embedding, data_mat.astype(str).fillna("NA").loc[:, [var]]],
                    axis=1,
                    sort=False,
                )

            if color_dictionary is None:
                color_dictionary = {}
            categories = set(var_data)

            if var in color_dictionary:
                color_dict = color_dictionary[var]
            else:
                random.seed(seed)
                color = [
                    mcolors.to_rgb("#" + "%06x" % random.randint(0, 0xFFFFFF))
                    for i in range(len(categories))
                ]
                color_dict = dict(zip(categories, color))

            if num_columns > 1:
                plt.subplot(num_rows, num_columns, i)
                i = i + 1

            if remove_nan and (data_mat[var].isnull().sum() > 0):
                plt.scatter(
                    emb_nan.iloc[:, 0],
                    emb_nan.iloc[:, 1],
                    c=data_mat.loc[:, var].dropna().apply(lambda x: color_dict[x]),
                    s=dot_size,
                    alpha=alpha,
                )
                plt.xlabel(emb_nan.columns[0])
                plt.ylabel(emb_nan.columns[1])
            else:
                plt.scatter(
                    embedding.iloc[:, 0],
                    embedding.iloc[:, 1],
                    c=data_mat.astype(str)
                    .fillna("NA")
                    .loc[:, var]
                    .apply(lambda x: color_dict[x]),
                    s=dot_size,
                    alpha=alpha,
                )
                plt.xlabel(embedding.columns[0])
                plt.ylabel(embedding.columns[1])

            if show_label:
                label_pos = label_pd.groupby(var).agg(
                    {label_pd.columns[0]: np.mean, label_pd.columns[1]: np.mean}
                )
                texts = []
                for label in label_pos.index.tolist():
                    texts.append(
                        plt.text(
                            label_pos.loc[label][0],
                            label_pos.loc[label][1],
                            label,
                            horizontalalignment="center",
                            verticalalignment="center",
                            size=text_size,
                            weight="bold",
                            color=color_dict[label],
                            path_effects=[
                                PathEffects.withStroke(linewidth=3, foreground="w")
                            ],
                        )
                    )
                adjust_text(texts)

            plt.title(var)
            patchList = []
            for key in color_dict:
                data_key = mpatches.Patch(color=color_dict[key], label=key)
                patchList.append(data_key)
            if show_legend:
                plt.legend(
                    handles=patchList, bbox_to_anchor=(1.04, 1), loc="upper left"
                )

            if num_columns == 1:
                if save is not None:
                    pdf.savefig(fig, bbox_inches="tight")
                plt.show()
        else:
            var_data = data_mat.copy().loc[:, var].to_list()
            o = np.argsort(var_data)
            if num_columns > 1:
                plt.subplot(num_rows, num_columns, i)
                i = i + 1
            plt.scatter(
                embedding.iloc[o, 0],
                embedding.iloc[o, 1],
                c=subset_list(var_data, o),
                cmap=cmap,
                s=dot_size,
                alpha=alpha,
            )
            plt.xlabel(embedding.columns[0])
            plt.ylabel(embedding.columns[1])
            plt.title(var)
            # setup the colorbar
            normalize = mcolors.Normalize(
                vmin=np.array(var_data).min(), vmax=np.array(var_data).max()
            )
            scalarmappaple = cm.ScalarMappable(norm=normalize, cmap=cmap)
            scalarmappaple.set_array(var_data)
            plt.colorbar(scalarmappaple)
            if num_columns == 1:
                if save is not None:
                    pdf.savefig(fig, bbox_inches="tight")
                plt.show()

    if num_columns > 1:
        plt.tight_layout()
        if save is not None:
            fig.savefig(save, bbox_inches="tight")
        plt.show()
    if (save is not None) & (num_columns == 1):
        pdf = pdf.close()
    return fig

def visualize(cistopic_obj, model, output):
    """
    Plot and save all visualizations
    """
    ntopics = str(model)
    res = [0.6, 1.2, 2.0, 2.5, 3.0, 3.5]
    for i in res:
        i_str = str(i)
        print('plotting for resolution ' + i_str)
        figure = plot_metadata_save(cistopic_obj,
                         reduction_name=f'UMAP_{ntopics}',
                         variables=[f'{ntopics}_pycisTopic_leiden_10_{i_str}'],
                         target='cell',
                         text_size=10,
                         dot_size=0.5,
                         figsize=(10,10))
        figure.savefig(output + f'/{ntopics}_pycisTopic_leiden_10_{i_str}_UMAP.png', format="png", dpi=300)

        figure = plot_metadata_save(cistopic_obj,
                         reduction_name=f'tSNE_{ntopics}',
                         variables=[f'{ntopics}_pycisTopic_leiden_10_{i_str}'],
                         target='cell',
                         text_size=10,
                         dot_size=0.5,
                         figsize=(10,10))
        figure.savefig(output + f'/{ntopics}_pycisTopic_leiden_10_{i_str}_tSNE.png', format="png", dpi=300)
        
    print('plotting QC and sample_id')
    figure = plot_metadata_save(cistopic_obj,
                     reduction_name=f'UMAP_{ntopics}',
                     variables=['log10_unique_fragments_count','sample_id'],
                     target='cell', num_columns=2,
                     text_size=10,
                     dot_size=1,
                     figsize=(10,5))
    figure.savefig(output + f'/{ntopics}_QC_sample_id_UMAP.png', format="png", dpi=300)
    
    figure = plot_metadata_save(cistopic_obj,
                     reduction_name=f'tSNE_{ntopics}',
                     variables=['log10_unique_fragments_count','sample_id'],
                     target='cell', num_columns=2,
                     text_size=10,
                     dot_size=1,
                     figsize=(10,5))
    figure.savefig(output + f'/{ntopics}_QC_sample_id_tSNE.png', format="png", dpi=300)
    
    print('plotting method')
    #plot clusters
    figure = plot_metadata_save(cistopic_obj,
                     reduction_name=f'UMAP_{ntopics}',
                     variables=['method'],
                     target='cell',
                    remove_nan = False,
                     text_size=10,
                     dot_size=0.5,
                     figsize=(10,10))
    figure.savefig(output + f'/{ntopics}_method_UMAP.png', format="png", dpi=300)
    
    figure = plot_metadata_save(cistopic_obj,
                     reduction_name=f'tSNE_{ntopics}',
                     variables=['method'],
                     target='cell',
                    remove_nan = False,
                     text_size=10,
                     dot_size=0.5,
                     figsize=(10,10))
    figure.savefig(output + f'/{ntopics}_method_tSNE.png', format="png", dpi=300)
    
    print(f'plotting annotation')
    variables = ['sample_id','consensus_cell_type_RNA','CISTOPIC_1_final_annot','CISTOPIC_1_final_annot_NA','method','pool_status']
    for var in variables:
        figure = plot_metadata_save(cistopic_obj,
                         reduction_name=f'UMAP_{ntopics}',
                         variables=[var],
                         target='cell',
                         text_size=10,
                         dot_size=0.5,
                         figsize=(10,10))
        figure.savefig(output + f'/{ntopics}_{var}_UMAP.png', format="png", dpi=300)

        figure = plot_metadata_save(cistopic_obj,
                         reduction_name=f'tSNE_{ntopics}',
                         variables=[var],
                         target='cell',
                         text_size=10,
                         dot_size=0.5,
                         figsize=(10,10))
        figure.savefig(output + f'/{ntopics}_{var}_tSNE.png', format="png", dpi=300)
    
    print('plotting done')     

    print('plot also not annotated cells')
    figure = plot_metadata_save(cistopic_obj,
             reduction_name=f'UMAP_{ntopics}',
             variables=['CISTOPIC_1_final_annot'],
             target='cell',
             remove_nan = False,
             text_size=10,
             dot_size=0.5,
             figsize=(10,10))
    figure.savefig(output + f'/{ntopics}_CISTOPIC_1_final_annot_all_UMAP.png', format="png", dpi=300)
        
    figure = plot_metadata_save(cistopic_obj,
             reduction_name=f'tSNE_{ntopics}',
             variables=['CISTOPIC_1_final_annot'],
             target='cell',
             remove_nan = False,
             text_size=10,
             dot_size=0.5,
             figsize=(10,10))
    figure.savefig(output + f'/{ntopics}_CISTOPIC_1_final_annot_all_tSNE.png', format="png", dpi=300)
    print('plotting done')   

print('functions finished')
print()

print('starting main')
print()
def main():
    """
    The main executable function
    """
    print('make argument parser')
    parser = make_argument_parser()
    args = parser.parse_args()

    print('Opening CTO...')
    filename = args.input_cto
    infile = open(filename, 'rb')
    cistopic_obj = pickle.load(infile)
    infile.close()
    print('Input cisTopic_object: ', filename)
    print('Opening done')
    print()
    
    model_dir = args.topic_model_dir
    print('Model dir: ', model_dir)
    print()
    
    model = args.selected_model
    print('Selected topic model: ', str(model))
    print()
    
    output = args.output_dir
    print('Output dir:', output)
    if not os.path.exists(output):
        print('Output dir does not exist, creating...')
        os.makedirs(output)
        print(f"Directory '{output}' created.")
    else:
        print(f"Directory '{output}' already exists.")
    print()
   
    print('Load and select model:')
    load_and_add_topic_model(model_dir, model, cistopic_obj)
    print('model added')
    print()

    print('Clustering and UMAP....')
    clustering_and_dim_red(cistopic_obj, model)
    print()

    print('Opening CTO...')
    ntopics = str(model)

    print('visualization')
    visualize(cistopic_obj, model, output)
    print('visualization done!')

if __name__ == "__main__":
    main()