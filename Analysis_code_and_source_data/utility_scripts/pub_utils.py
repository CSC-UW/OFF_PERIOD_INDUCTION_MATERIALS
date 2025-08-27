import acr
from acr.utils import *
import matplotlib.pyplot as plt
import polars as pl
import pingouin as pg
import os

def get_subject_list(type='som', exp=None):
    import acr.utils as au
    animals = [sub for sub in au.sub_exp_types if au.sub_exp_types[sub] == type]
    if exp is not None:
        animals = [sub for sub in animals if sub in au.sub_exp_manifest[exp]]
        mnfst = au.sub_exp_manifest[exp]
        exps = [mnfst[sub][0] for sub in animals]
        return animals, exps
    return animals

def write_comp_data(df, name):
    path = f'{acr.utils.pub_data_root}/comparison_data/{name}.csv'
    df.to_csv(path, index=False)
    return

def read_comp_data(name):
    path = f'{acr.utils.pub_data_root}/comparison_data/{name}.csv'
    df = pd.read_csv(path)
    return df


def write_quick_data(df, name):
    if not os.path.exists('./quick_temp_data'):
        os.makedirs('./quick_temp_data')
    df.to_csv(f'./quick_temp_data/{name}.csv', index=False)
    return

def read_quick_data(name):
    df = pd.read_csv(f'./quick_temp_data/{name}.csv')
    return df

def save_fig(figure_path, pad=0.2, dpi=600, transparent=True):
    plt.tight_layout(pad=pad)
    plt.savefig(figure_path, dpi=dpi, transparent=transparent)
    return
def write_plot_info(f, ax, figure_path, stats=None):
    txt_path = figure_path.replace('.svg', '_info.txt')
    with open(txt_path, 'w') as f:
        f.write(f'{figure_path}\n')
        if stats is not None:
            f.write(f'p-val: {stats["p-val"]}\n')
        elif stats is None:
            f.write('no stats\n')
        # write the y axis limits and tick values
        f.write(f'ylim: {ax.get_ylim()}\n')
        f.write(f'ytick: {ax.get_yticks()}\n')
        f.write(f'yticklabels: {ax.get_yticklabels()}\n')
        f.write(f'xlim: {ax.get_xlim()}\n')
        f.write(f'xtick: {ax.get_xticks()}\n')
        f.write(f'xticklabels: {ax.get_xticklabels()}\n')
    return 


def add_boxplot(ax, data, positions=[0.5], widths=0.06, color=NNXR_GRAY, means=True, mean_color='#ffaf2b', mean_size=18, mew=4):
    box = ax.boxplot(data, positions=positions, widths=widths, showfliers=False, patch_artist=True, capprops=dict(color='none', linewidth=0), whiskerprops=dict(color='k', linewidth=3), medianprops=dict(color='k', linewidth=3, zorder=101), showmeans=means, meanprops=dict(mfc=mean_color, mec=mean_color, mew=mew, marker='_', markersize=mean_size, zorder=102))
    box['boxes'][0].set_facecolor(color)
    box['boxes'][0].set_alpha(0.8)
    box['boxes'][0].set_linewidth(0)
    box['boxes'][0].set_zorder(100)
    return ax, box


def _box_onesided(nnxr, nnxo, stats=False, color=ACR_BLUE, means=True, mean_color='#ffaf2b', mean_size=18, mew=4, fsize=(2.5, 3)):
    f, ax = plt.subplots(1, 1, figsize=fsize)

    # boxplot - only NNXo
    box_o = ax.boxplot(nnxo, positions=[0.65], widths=0.06, patch_artist=True, capprops=dict(color='none', linewidth=0), whiskerprops=dict(color='black', linewidth=3), medianprops=dict(color='k', linewidth=3, zorder=101), showfliers=False, showmeans=means, meanprops=dict(mfc=mean_color, mec=mean_color, mew=mew, marker='_', markersize=mean_size, zorder=102))
    box_o['boxes'][0].set_facecolor(color)
    box_o['boxes'][0].set_alpha(0.8)
    box_o['boxes'][0].set_linewidth(0)
    box_o['boxes'][0].set_zorder(100)

    # line plots
    for i in range(len(nnxr)):
        ax.plot([0.4, 0.6], [nnxr[i], nnxo[i]], color=color, alpha=0.85, linewidth=3.5, solid_capstyle='round', solid_joinstyle='round')

    # scatter plots for individual points
    for i in range(len(nnxr)):
        ax.scatter(0.4, nnxr[i], color=NNXR_GRAY, alpha=0.7, s=30, zorder=202)
        ax.scatter(0.6, nnxo[i], color=color, alpha=0.7, s=30, zorder=203)
    if stats:
        stats = pg.ttest(nnxr, nnxo, paired=True)
        return f, ax, stats
    else:
        return f, ax
    
    
def swa_boxplot(nnxr, nnxo, stats=False, color=ACR_BLUE, one_sided=False, means=True, mean_color='#ffaf2b', mean_size=18, mew=4, ax=None, fsize=(2.5, 3)):
    
    if one_sided:
        return _box_onesided(nnxr, nnxo, stats, color, means, mean_color, mean_size, mew, fsize)
    
    if ax is None:
        f, ax = plt.subplots(1, 1, figsize=fsize)
    else:
        f = ax.get_figure()
    
    
    # boxplots
    box = ax.boxplot(nnxr, positions=[0.35], widths=0.06, patch_artist=True, capprops=dict(color='none', linewidth=0), whiskerprops=dict(color='k', linewidth=3), medianprops=dict(color='k', linewidth=3, zorder=101), showfliers=False, showmeans=means, meanprops=dict(mfc=mean_color, mec=mean_color, mew=mew, marker='_', markersize=mean_size, zorder=102))
    box['boxes'][0].set_facecolor(NNXR_GRAY)
    box['boxes'][0].set_alpha(0.8)
    box['boxes'][0].set_linewidth(0)
    box['boxes'][0].set_zorder(100)

    box_o = ax.boxplot(nnxo, positions=[0.65], widths=0.06, patch_artist=True, capprops=dict(color='none', linewidth=0), whiskerprops=dict(color='black', linewidth=3), medianprops=dict(color='k', linewidth=3, zorder=101), showfliers=False, showmeans=means, meanprops=dict(mfc=mean_color, mec=mean_color, mew=mew, marker='_', markersize=mean_size, zorder=102))
    box_o['boxes'][0].set_facecolor(color)
    box_o['boxes'][0].set_alpha(0.8)
    box_o['boxes'][0].set_linewidth(0)
    box_o['boxes'][0].set_zorder(100)

    # line plots
    for i in range(len(nnxr)):
        ax.plot([0.40, 0.6], [nnxr[i], nnxo[i]], color=color, alpha=0.85, linewidth=3.5, solid_capstyle='round', solid_joinstyle='round')

    # scatter plots for individual points
    for i in range(len(nnxr)):
        ax.scatter(0.4, nnxr[i], color=NNXR_GRAY, alpha=0.7, s=30, zorder=202)
        ax.scatter(0.6, nnxo[i], color=color, alpha=0.7, s=30, zorder=203)
    if stats:
        stats = pg.ttest(nnxr, nnxo, paired=True)
        return f, ax, stats
    else:
        return f, ax

def _drop_sub_channels_pd(df, sub_channel_dict):
    pldf = pl.DataFrame(df)
    for sub, channels in sub_channel_dict.items():
        pldf = pldf.filter(~((pl.col('subject') == sub) & (pl.col('channel').is_in(channels))))
    return pldf.to_pandas()

def drop_sub_channels(df, sub_channel_dict):
    # drop any rows containing BOTH the subject and channel
    # drop only rows where the subject matches subject AND channel matches channel
    if type(df) == pd.DataFrame:
        return _drop_sub_channels_pd(df, sub_channel_dict)
    
    for sub, channels in sub_channel_dict.items():
        df = df.filter(~((pl.col('subject') == sub) & (pl.col('channel').is_in(channels))))
    return df


def _layers_to_pldf(df, subject):
    """Requires a dataframe with only a single subject

    Parameters
    ----------
    df : _type_
        _description_
    subject : _type_
        _description_

    Returns
    -------
    _type_
        _description_
    """
    if len(df['subject'].unique()) > 1:
        print(f'Multiple subjects in df: {df["subject"].unique()}')
        return df
    if 'layer' not in df.columns:
        df = df.with_columns(pl.lit('no_layer_info').alias('layer'))
    chan_map = acr.info_pipeline.subject_info_section(subject, 'channel_map')
    if chan_map == None:
        print(f'No channel map for {subject}')
        return df
    if len(chan_map) == 0:
        print(f'No channel map for {subject}')
        return df

    stores = df['store'].unique()
    #add layer information:
    for store in stores:
        for chan in df.prb(store)['channel'].unique():
            layer = chan_map[store][str(chan)]['layer']
            df = df.with_columns(pl.when((pl.col('channel') == chan) & (pl.col('store') == store)).then(pl.lit(layer)).otherwise(pl.col('layer')).alias('layer'))
    return df


def add_layer_info_to_df(df, subject):
    
    return_pl = False
    if type(df) == pl.DataFrame:
        return _layers_to_pldf(df, subject)
    if 'layer' not in df.columns:
        df['layer'] = 'no_layer_info'
    chan_map = acr.info_pipeline.subject_info_section(subject, 'channel_map')
    if chan_map == None:
        print(f'No channel map for {subject}')
        return df
    if len(chan_map) == 0:
        print(f'No channel map for {subject}')
        return df

    #check if layer information has already been added:
    stores = df.sbj(subject)['store'].unique() if 'store' in df.columns else df.sbj(subject)['probe'].unique()
    chans = df.sbj(subject).prb(stores[0])['channel'].unique()
    test_chan1 = chans[0]
    test_chan2 = chans[-1]
    layer1 = chan_map[stores[0]][str(test_chan1)]['layer']
    layer2 = chan_map[stores[0]][str(test_chan2)]['layer']
    
    if (df.sbj(subject).prb(stores[0]).chnl(test_chan1)['layer'].values[0] == layer1) & (df.sbj(subject).prb(stores[0]).chnl(test_chan2)['layer'].values[0] == layer2):
        print(f'Layer information already added for {subject}')
        return df

    #add layer information:
    for store in stores:
        for chan in df.prb(store)['channel'].unique():
            layer = chan_map[store][str(chan)]['layer']
            df.loc[(df.subject==subject) & (df.prb(store)['channel']==chan), 'layer'] = layer
    
    return df


def sub_regions_to_df(df):
    locs = acr.utils.sub_probe_locations
    df['region'] = 'none'
    for sub in df['subject'].unique():
        df.loc[df['subject'] == sub, 'region'] = locs[sub]
    return df

def write_stats(stats_df, name):
    if os.path.exists(f'./statistical_results')==False:
        os.makedirs(f'./statistical_results')
    stats_df.to_csv(f'./statistical_results/{name}.csv')
    return

def write_source_data(df, name):
    # make sure the statistical_results folder exists
    if os.path.exists(f'./statistical_source_data')==False:
        os.makedirs(f'./statistical_source_data')
    df.to_csv(f'./statistical_source_data/{name}.csv')
    return

def remove_subject(sub_to_remove, subjects, experiments):
    """
    Remove a subject and its corresponding experiment from both lists.
    
    Parameters:
    -----------
    sub_to_remove : str
        The subject to remove from the lists
    subjects : list
        List of subjects
    experiments : list
        List of experiments corresponding to subjects
        
    Returns:
    --------
    tuple
        Updated (subjects, experiments) lists with the specified subject removed
    """
    try:
        # Find the index of the subject to remove
        index_to_remove = subjects.index(sub_to_remove)
        
        # Remove the subject and corresponding experiment
        updated_subjects = subjects.copy()
        updated_experiments = experiments.copy()
        
        updated_subjects.pop(index_to_remove)
        updated_experiments.pop(index_to_remove)
        
        return updated_subjects, updated_experiments
        
    except ValueError:
        print(f"Subject {sub_to_remove} not found in subjects list")
        return subjects, experiments