try:
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
except ImportError:
    import warnings
    warnings.warn('Plotting not available, matplotlib not installed')

import numpy as np

def get_a_color():
    #TODO
    return

def plot_class_distributions(outfile, value_maps, fuse_buried_core = False, fuse_misc_interaction = False, order = None, hide_small_numbers = False, center_small_sumbers = True, hide_legend = False, title = '', ignore_unmapped = True, ignore_disorder = False, fontsize = 18, small_fontsize = 14, custom_bbox = None, custom_class_order = None):

    color_map = {
        'Core':'red',
        'Buried':'darkred',
        'Surface':'blue',
        'Ligand Interaction':'orange',
        'Protein Interaction':'green',
        'Peptide Interaction':'olive',
        'DNA Interaction':'cyan',
        'Disorder':'purple',
        'Metal Interaction':'silver',
        'Ion Interaction':'black',
        'RNA Interaction':'lightblue',
        'Misc Interaction':'orange'
    }

    fuse_buried_core_map = {'Buried':'Core'}
    fuse_misc_interaction_map = {
        'Ligand Interaction':'Misc Interaction',
        'DNA Interaction':'Misc Interaction',
        'Metal Interaction':'Misc Interaction',
        'Ion Interaction':'Misc Interaction',
        'RNA Interaction':'Misc Interaction',
        'Peptide Interaction' : 'Misc Interaction',}

    classes = ['Core', 'Buried', 'Surface',
               'Ligand Interaction', 'Protein Interaction', 'DNA Interaction', 'Metal Interaction', 'Ion Interaction', 'RNA Interaction', 'Peptide Interaction',
               'Disorder'
    ]

    if custom_class_order is not None:
        classes = custom_class_order

    consistency_map = {
        'metal interaction' : 'Metal Interaction',
        'ion interaction' : 'Ion Interaction',
        'DNA interaction' : 'DNA Interaction',
        'Protein interaction' : 'Protein Interaction',
        'RNA interaction' : 'RNA Interaction',
        'ligand interaction' : 'Ligand Interaction',
        'Peptide interaction' : 'Peptide Interaction',
    }

    consistent_value_map = {}
    for name in value_maps:
        consistent_value_map[name] = {}
        for classification in value_maps[name]:
            if classification in consistency_map:
                consistent_classification = consistency_map[classification]
            else:
                consistent_classification = classification
            if ignore_unmapped and classification == '-':
                continue
            if ignore_disorder and classification == 'Disorder':
                continue
            consistent_value_map[name][consistent_classification] = value_maps[name][classification]
    value_maps = consistent_value_map

    if fuse_buried_core:
        classes.remove('Buried')
        fused_value_maps = {}
        for name in value_maps:
            fused_value_maps[name] = {}
            for classification in value_maps[name]:
                if classification in fuse_buried_core_map:
                    fused_class = fuse_buried_core_map[classification]
                else:
                    fused_class = classification
                if not fused_class in fused_value_maps[name]:
                    fused_value_maps[name][fused_class] = 0
                fused_value_maps[name][fused_class] += value_maps[name][classification]
        value_maps = fused_value_maps

    if fuse_misc_interaction:
        for int_type in fuse_misc_interaction_map:
            classes.remove(int_type)

        classes.remove('Disorder')
        if 'Misc Interaction' not in classes:
            classes.append('Misc Interaction')
        classes.append('Disorder')
        fused_value_maps = {}
        for name in value_maps:
            fused_value_maps[name] = {}
            for classification in value_maps[name]:
                if classification in fuse_misc_interaction_map:
                    fused_class = fuse_misc_interaction_map[classification]
                else:
                    fused_class = classification
                if not fused_class in fused_value_maps[name]:
                    fused_value_maps[name][fused_class] = 0
                fused_value_maps[name][fused_class] += value_maps[name][classification]
        value_maps = fused_value_maps

    if ignore_disorder:
        classes.remove('Disorder')

    x_labels = []

    N = len(value_maps)

    max_values = {}
    percentage_value_maps = {}
    for name in value_maps:
        max_values[name] = 0.
        percentage_value_maps[name] = {}
        value_sum = 0.
        for classification in value_maps[name]:
            value_sum += value_maps[name][classification]

        for classification in value_maps[name]:
            perc_value = value_maps[name][classification]/value_sum
            percentage_value_maps[name][classification] = perc_value
            if perc_value > max_values[name]:
                max_values[name] = perc_value


    ind = np.arange(N)    # the x locations for the groups
    width = 0.8 #1.0/float(N)       # the width of the bars: can also be len(x) sequence

    fig = plt.figure(figsize=(4.5*N,9.))
    ax = plt.subplot(111)

    bars = []
    new_bar = True
    plots = []

    if order is None:
        order = list(value_maps.keys())

    for dataset_pos, name in enumerate(order):
        x_label = name
        x_labels.append(x_label)
        last_bar_small = False
        last_bar_ha = 'center'
        max_value = max_values[name]
        for classification in reversed(classes):
            if classification in percentage_value_maps[name]:
                value = percentage_value_maps[name][classification]
            else:
                value = 0.
                continue

            if classification in color_map:
                col = color_map[classification]
            else:
                col = get_a_color()
                color_map[classification] = col

            if new_bar:
                p = ax.bar([dataset_pos+width*0.5], value, width,color = col)
                cumulated_value = value
                new_bar = False
            else:
                p = ax.bar([dataset_pos+width*0.5], value, width,bottom=cumulated_value,color = col)
                cumulated_value += value
            if value > 0.045:
                plt.text(dataset_pos+width*0.5,cumulated_value-value*0.5,'%.2f%%' % (value*100.0),ha="center", va="center", color="white", fontsize=fontsize)
                last_bar_small = False
                last_bar_ha = 'center'
            elif hide_small_numbers:
                pass
            elif center_small_sumbers:
                plt.text(dataset_pos+width*0.5,cumulated_value-value*0.5,'%.2f%%' % (value*100.0),ha="center", va="center", color="white", fontsize=small_fontsize)
            elif not last_bar_small or last_bar_ha == 'right':
                plt.text(dataset_pos+width*0.5,cumulated_value-value*0.5,'%.2f%%' % (value*100.0),ha="center", va="center", color="white", fontsize=small_fontsize)
                last_bar_small = True
                last_bar_ha = 'center'
            elif last_bar_ha == 'center':
                plt.text(dataset_pos,cumulated_value-value*0.5,'%.2f%%' % (value*100.0),ha="left", va="center", color="white", fontsize=small_fontsize)
                last_bar_small = True
                last_bar_ha = 'left'
            elif last_bar_ha == 'left':
                plt.text(dataset_pos+width,cumulated_value-value*0.5,'%.2f%%' % (value*100.0),ha="right", va="center", color="white", fontsize=small_fontsize)
                last_bar_small = True
                last_bar_ha = 'right'
            plots.append(p)
        new_bar = True

    box = ax.get_position()
    ax.set_position([box.x0, box.y0, box.width * 0.6, box.height])

    plt.ylabel('Structural Classification',fontsize=fontsize)
    plt.title(title,fontsize=fontsize)
    plt.xticks([x+(width/2.0) for x in ind],['\n'.join(x.split(' ')) for x in x_labels],fontsize=fontsize)
    plt.yticks(np.arange(0, 1.01, 0.1),["%s%%" % str(x*10) for x in np.arange(11)],fontsize=fontsize)

    plt.xlim((0,N-1+width))
    plt.ylim((0,1))
    if not hide_legend:
        anchor_map = {1:2.7,
                      2:1.85,
                      3:1.55,
                      4:1.4,
                      5:1.325,
                      6:1.275,
                      9:0.9
        } #Improve by manual labor
        if N in anchor_map:
            x_anchor = anchor_map[N]
        else:
            x_anchor = 2.-N*0.15

        bbox = (x_anchor, 0.5)

        if custom_bbox is not None:
            bbox = custom_bbox

        ax.legend([p[0] for p in reversed(plots)], classes,loc = 7,fontsize=fontsize,bbox_to_anchor = bbox)

    plt.savefig(outfile,bbox_inches='tight')

def makeViolins(violins, outfile, session_name, add=''):
    fs = 24  # fontsize
    matplotlib.rcParams.update({'font.size': fs})
    plt.clf()

    for violin_tag in violins:
        fig, axes = plt.subplots(nrows=1, ncols=1, figsize=(12, 12))
        classes = list(violins[violin_tag].keys())
        pos = list(range(len(classes)))
        data = [violins[violin_tag][cl] for cl in classes]
        axes.violinplot(data, pos, widths=0.5,
                        showmeans=True, showextrema=True, showmedians=True,
                        bw_method='silverman')

        axes.set_xticks(pos)
        axes.set_xticklabels(classes, rotation='vertical')
        # axes[n].set_xlabel(classes, rotation='vertical')

        # for ax in axes.flatten():
        #    ax.set_yticklabels([])

        fig.suptitle("Violin Plots %s - %s value" % (session_name, violin_tag))
        # fig.subplots_adjust(hspace=0.0)
        plt.tight_layout()
        # plt.show()
        violin_file = '%s.violin_plots_%s%s.png' % (outfile, violin_tag, add)
        plt.savefig(violin_file)

def plotScatter(scatter_plots, outfile, session_name):
    """
    x_values = []
    y_values = []

    for pdb in dssp_map:
        for (chain,res) in dssp_map[pdb]:
            (acc,aa) = dssp_map[pdb][(chain,res)]
            measure = measure_map[pdb][chain][res]
            x_values.append(acc)
            if not anti_corr:
                y_values.append(measure)
            else:
                y_values.append(-measure)

    fig, ax = plt.subplots()
    ax.scatter(x_values,y_values,alpha=0.01)
    ax.set_xlabel('RSA', fontsize=15)
    ax.set_ylabel(y_axis_name, fontsize=15)

    plt.savefig(outfile,dpi=300)
    """
    fs = 10  # fontsize
    plt.clf()
    scores = ['LI score', 'CI score', 'SI score', 'MI score', 'LoI score']
    degrees = ['LI degree', 'CI degree', 'SI degree', 'MI degree', 'LoI degree']
    labels = ['Ligand Interaction Score', 'Chain Interaction Score', 'Short Interaction Score', 'Medium Interaction Score', 'Long Interaction Score']
    for scatter_tag in scatter_plots:
        fig, axes = plt.subplots(nrows=2, ncols=3, figsize=(12, 12))

        n = 0
        for score in scores:
            b = n
            a = 0
            if n > 1:
                b = n - 2
                a = 1
            for Class in scatter_plots[scatter_tag][score]:
                if Class == 'Core':
                    axes[a, b].scatter(scatter_plots[scatter_tag][score][Class][0], scatter_plots[scatter_tag][score][Class][1], alpha=0.05, color='red')
                    axes[a, b].set_xlabel(labels[n], fontsize=fs)
                    axes[a, b].set_ylabel(scatter_tag, fontsize=fs)
                """else:
                    axes[a,b].scatter(scatter_plots[scatter_tag][score][Class][0],scatter_plots[scatter_tag][score][Class][1],alpha = 0.05,color='blue')
                    axes[a,b].set_xlabel(labels[n], fontsize=fs)
                    axes[a,b].set_ylabel(scatter_tag, fontsize=fs)"""
            n += 1

        fig.suptitle("Scatter Plots Interaction Scores %s - %s value" % (session_name, scatter_tag))
        # fig.subplots_adjust(hspace=0.0)
        plt.tight_layout()
        # plt.show()
        scatter_file = '%s.scatter_plots_Iscore_%s.png' % (outfile, scatter_tag)

        plt.savefig(scatter_file)

        plt.clf()
        fig, axes = plt.subplots(nrows=2, ncols=3, figsize=(12, 12))
        n = 0
        for score in degrees:

            b = n
            a = 0
            if n > 1:
                b = n - 2
                a = 1
            for Class in scatter_plots[scatter_tag][score]:
                if Class == 'Core':
                    axes[a, b].scatter(scatter_plots[scatter_tag][score][Class][0], scatter_plots[scatter_tag][score][Class][1], alpha=0.05, color='red')
                    axes[a, b].set_xlabel(labels[n], fontsize=fs)
                    axes[a, b].set_ylabel(scatter_tag, fontsize=fs)
                """else:
                    axes[a,b].scatter(scatter_plots[scatter_tag][score][Class][0],scatter_plots[scatter_tag][score][Class][1],alpha = 0.05,color='blue')
                    axes[a,b].set_xlabel(labels[n], fontsize=fs)
                    axes[a,b].set_ylabel(scatter_tag, fontsize=fs)"""
            n += 1

        fig.suptitle("Scatter Plots Interaction Degrees %s - %s value" % (session_name, scatter_tag))
        # fig.subplots_adjust(hspace=0.0)
        plt.tight_layout()
        # plt.show()
        scatter_file = '%s.scatter_plots_Idegrees_%s.png' % (outfile, scatter_tag)

        plt.savefig(scatter_file)
