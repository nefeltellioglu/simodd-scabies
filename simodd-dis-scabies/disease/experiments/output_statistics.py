import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Rectangle

from .plotting_disease import setup_x_axis


mpl.rcParams['xtick.major.pad'] = '8'
mpl.rcParams['ytick.major.pad'] = '8'

try:
    import brewer2mpl as brew

    colors = ['#8dd3c7', '#ffffb3', '#bebada', '#fb8072', '#80b1d3', '#fdb462', '#b3de69']
    #colors = brew.get_map('Pastel2', 'qualitative', 7).mpl_colors
    #colors = ['k'] + brew.get_map('Pastel2', 'qualitative', 7).mpl_colors
    #colors[6] = 'c' # replace invisible yellow
    rb_colors = ['k'] + brew.get_map('YlOrRd', 'sequential', 9).mpl_colors
except ImportError:
    colors = ['k', 'r', 'b', 'g', 'm', 'c', '0.5']

## -------------------------------------------------------------------------------------------
## from http://www.sron.nl/~pault/
## -------------------------------------------------------------------------------------------

# colour table in HTML hex format
hexcols = ['#332288', '#88CCEE', '#44AA99', '#117733', '#999933', '#DDCC77',
           '#CC6677', '#882255', '#AA4499', '#661100', '#6699CC', '#AA4466',
           '#4477AA']

greysafecols = ['#809BC8', '#FF6666', '#FFCC66', '#64C204']

xarr = [[12],
        [12, 6],
        [12, 6, 5],
        [12, 6, 5, 3],
        [0, 1, 3, 5, 6],
        [0, 1, 3, 5, 6, 8],
        [0, 1, 2, 3, 5, 6, 8],
        [0, 1, 2, 3, 4, 5, 6, 8],
        [0, 1, 2, 3, 4, 5, 6, 7, 8],
        [0, 1, 2, 3, 4, 5, 9, 6, 7, 8],
        [0, 10, 1, 2, 3, 4, 5, 9, 6, 7, 8],
        [0, 10, 1, 2, 3, 4, 5, 9, 6, 11, 7, 8]]


# get specified nr of distinct colours in HTML hex format.
# in: nr - number of colours [1..12]
# returns: list of distinct colours in HTML hex
def get_distinct(nr, bw=False, br=False):
    # check if nr is in correct range

    if nr <= 2 and bw:
        return ['#333333', '#AAAAAA']

    if nr <= 2 and br:
        print("BBBBBBB")
        return ['#FFFFFF', '#FF0000']

#    if nr <= 2 and br:
#        return ['#333333', '#CC6677']

    if nr < 1 or nr > 12:
        print("wrong nr of distinct colours!")
        return
    # get list of indices
    lst = xarr[nr - 1]
    # generate colour list by stepping through indices and looking them up
    # in the colour table
    i_col = 0
    col = [0] * nr
    for idx in lst:
        col[i_col] = hexcols[idx]
        i_col += 1
    return col


## -------------------------------------------------------------------------------------------

# ideally this should be expanded, and be settable by main.py
class lp(object):
    lw = 1
    fs = 16 # 16 is very big (John thinks 12 looks better)


def plot_lines(ax, x_values, y_values, x_label, y_label, y_errs=None, series_labels=None, legend_title=None,
               legend_loc='upper left', linestyle='-', lineweight=None, bw=False, br=False, reversed=False, xmin=None):
    """
    Plotting function for lines, with or without error bars.
    """

    # if data values submitted as a single series, wrap in list
    if not (isinstance(y_values[0], list)
            or isinstance(y_values[0], np.ndarray)
            or isinstance(y_values[0], tuple)):
        y_values = [y_values]

    #print len(x_values), len(y_values[0])
    assert len(x_values) == len(y_values[0])

    if len(y_values) <= 12:
        loc_colors = get_distinct(len(y_values), bw, br)
    else:
        # get max number of colours
        loc_colors = get_distinct(12)

    # same with error values
    if y_errs is not None and (not isinstance(y_errs[0], (list, np.ndarray, tuple))):
        y_errs = [y_errs]

    # expand y_errs and series_labels if not present
    if not y_errs:
        y_errs = [None] * len(y_values)
    if not series_labels:
        series_labels = [None] * len(y_values)

    if reversed:
        y_values.reverse()
        loc_colors.reverse()

    # plot lines
    for i, (cur_values, cur_errs, cur_label) in enumerate(zip(y_values, y_errs, series_labels)):
        #print bw, cur_label, loc_colors[i]
        if cur_errs is not None:
            ax.errorbar(x_values, cur_values, yerr=cur_errs, label=cur_label, color=loc_colors[i],
                        lw=lp.lw if not lineweight else lineweight, ls=linestyle)
        else:
            #print len(x_values), len(cur_values)
            ax.plot(x_values, cur_values, label=cur_label, color=loc_colors[i % len(loc_colors)],
                    lw=lp.lw if not lineweight else lineweight, ls=linestyle)

    # set labels and legend
    ax.set_xlabel(x_label, fontsize=lp.fs)
    ax.set_ylabel(y_label, fontsize=lp.fs)

    if xmin is not None:
        ax.set_xlim(xmin=xmin)

    if series_labels[0]:
        leg = ax.legend(scatterpoints=1, numpoints=1, ncol=2,
                        title=legend_title, labelspacing=0.1, columnspacing=1,
                        loc=legend_loc)
        plt.setp(leg.get_texts(), fontsize='small')
        leg.get_frame().set_alpha(0.0)

    plt.setp(ax.get_xticklabels(), fontsize=lp.fs)
    plt.setp(ax.get_yticklabels(), fontsize=lp.fs)


def plot_points_lin_reg(ax, times, y_values, x_label, y_label, series_labels=None, y_min=None, y_max=None,
                        x_offset=None, align=None,
                        legend_title=None, legend_loc='upper left', plot_title=None, bw=False):
    """
    Plot a set of data points together with a line of best fit.
    """
    if not (isinstance(y_values[0], (list, np.ndarray, tuple))):
        y_values = [y_values]

    if not series_labels:
        series_labels = [None] * len(y_values)

    loc_colors = get_distinct(len(y_values), bw)

    if align == 'centre':
        offset = (times['et'] - times['st']) / (len(y_values[0]) + 1) / 2.0
        #print times['st'], times['et'], offset
        x_values = np.linspace(times['st'] + offset, times['et'] - offset, len(y_values[0])) / times['t_per_year']
        #print xvals
    else:
        x_values = np.linspace(times['st'], times['et'], len(y_values[0])) / times['t_per_year']
        #x_values = np.linspace(times['st'], times['et'], len(y_values[0])) / times['t_per_year']

    for i, (cur_values, cur_label) in enumerate(zip(y_values, series_labels)):
        p = np.polyfit(x_values[~np.isnan(cur_values)], cur_values[~np.isnan(cur_values)], 1)
        #print cur_values
        ax.plot(x_values, cur_values, color=loc_colors[i], marker='x', lw=0)
        ax.plot(x_values, p[0] * x_values + p[1], label=cur_label, color=loc_colors[i], linestyle='-', lw=lp.lw)

    if plot_title:
        ax.set_title(plot_title, fontsize=lp.fs)
    ax.set_xlabel(x_label, fontsize=lp.fs)
    ax.set_ylabel(y_label, fontsize=lp.fs)
    ncol = 2 if len(y_values) > 5 else 1
    if series_labels[0]:
        leg = ax.legend(scatterpoints=1, numpoints=1, ncol=ncol, title=legend_title, labelspacing=0.1, columnspacing=1,
                        loc=legend_loc)
        plt.setp(leg.get_texts(), fontsize='small')
        leg.get_frame().set_alpha(0.0)

    if y_min is not None:
        ax.set_ylim(ymin=y_min)
    if y_max is not None:
        ax.set_ylim(ymax=y_max)

    plt.setp(ax.get_xticklabels(), fontsize=lp.fs)
    plt.setp(ax.get_yticklabels(), fontsize=lp.fs)

    setup_x_axis(ax, times, x_offset, set_xlim=False)


def plot_stacked_area(ax, x_values, y_values, x_label, y_label, series_labels=None,
                      legend_title=None, plot_title=None, bw=False):
    if not series_labels:
        series_labels = [None] * len(y_values)

    #print series_labels

    p = []

    colors = get_distinct(len(y_values), bw)

    ax.fill_between(x_values, 0, y_values[0], lw=0, edgecolor='grey', facecolor=colors[0])
    p.append(Rectangle((0, 0), 1, 1, lw=0, edgecolor='grey', facecolor=colors[0]))  # for legend
    for i, (cur_bottom, cur_top) in enumerate(zip(y_values[:-1], y_values[1:])):
        ax.fill_between(x_values, cur_bottom, cur_top, lw=0, edgecolor='grey', facecolor=colors[i + 1])
        p.append(Rectangle((0, 0), 1, 1, lw=0, edgecolor='grey', facecolor=colors[i + 1]))

    ax.set_ylim(ymax=1.0)

    # set labels and legend
    if plot_title:
        ax.set_title(plot_title, fontsize=lp.fs)
    ax.set_xlabel(x_label, fontsize=lp.fs)
    ax.set_ylabel(y_label, fontsize=lp.fs)
    if series_labels[0]:
        leg = ax.legend(p, series_labels, ncol=2, title=legend_title, labelspacing=0.1, columnspacing=1,
                        loc='upper right')
        plt.setp(leg.get_texts(), fontsize='small')
        leg.get_frame().set_alpha(0.8)
        leg.get_frame().set_linewidth(0.0)
        #leg = ax.legend(p, series_labels)

    plt.setp(ax.get_xticklabels(), fontsize=lp.fs)
    plt.setp(ax.get_yticklabels(), fontsize=lp.fs)


def output_bars(ofile, x_values, y_values, x_label, y_label, x_ticklabels=None, y_errs=None, ymin=0.0, ymax=None,
                series_labels=None,
                legend_loc=2, rotate_ticks=False, plot_title=None, bw=False):
    fig = plt.figure()
    ax = fig.add_subplot(111)

    plot_bars(ax, x_values, y_values, x_label, y_label, x_ticklabels, y_errs, ymin, ymax, series_labels,
              legend_loc, rotate_ticks, plot_title, bw)

    fig.tight_layout()
    fig.savefig(ofile)


def plot_bars(ax, x_values, y_values, x_label, y_label, x_ticklabels=None, y_errs=None, ymin=0.0, ymax=None,
              series_labels=None,
              legend_loc=2, rotate_ticks=False, plot_title=None, bw=False):
    # if data values submitted as a single series, wrap in list
    if not (isinstance(y_values[0], list)
            or isinstance(y_values[0], np.ndarray)
            or isinstance(y_values[0], tuple)):
        y_values = [y_values]

    #print len(x_values), len(y_values[0])
    assert len(x_values) == len(y_values[0])

    loc_colors = get_distinct(len(y_values), bw)

    # same with error values
    if y_errs is not None and (not isinstance(y_errs[0], (list, np.ndarray, tuple))):
        y_errs = [y_errs]

    # expand y_errs and series_labels if not present
    if not y_errs:
        y_errs = [None] * len(y_values)
    if not series_labels:
        series_labels = [None] * len(y_values)

    width = 0.8 / len(y_values)

    x_values = np.array(x_values)

    # plot bars
    for i, (cur_values, cur_errs, cur_label) in enumerate(zip(y_values, y_errs, series_labels)):
        ax.bar(x_values + i * width, cur_values, width, lw=0, edgecolor='grey', color=loc_colors[i], ecolor='dimgrey',
               align='center', yerr=cur_errs, label=cur_label)

    if plot_title:
        ax.set_title(plot_title, fontsize=lp.fs)
    ax.set_ylabel(y_label, fontsize=lp.fs)
    ax.set_xlabel(x_label, fontsize=lp.fs)
    if len(y_values) == 1:
        ax.set_xlim(xmin=-1, xmax=len(x_values))
    else:
        ax.set_xlim(xmin=-0.5, xmax=len(x_values))

    ax.set_ylim(ymin=ymin)
    if ymax is not None:
        ax.set_ylim(ymax=ymax)

    if x_ticklabels:
        ax.set_xticks(np.array(x_values) + len(y_values) / 2 * width)
        ax.set_xticklabels(tuple(x_ticklabels), rotation=(20 if rotate_ticks else 0))

    ax.xaxis.set_ticks_position('none')

    if series_labels[0]:
        leg = ax.legend(loc=legend_loc)
        leg.get_frame().set_alpha(0.0)

    plt.setp(ax.get_xticklabels(), fontsize=lp.fs)
    plt.setp(ax.get_yticklabels(), fontsize=lp.fs)


def plot_statistics(ofile, x_values, y_values, x_label, y_label, y_errs=None, series_labels=None,
                    ymin=None, xmax=None, xlogscale=False, ylogscale=False, title=None, bw=False,
                    legend_loc=None, legend_title=None):
    # create figure
    fig = plt.figure()
    ax = fig.add_subplot(111)

    plot_lines(ax, x_values, y_values, x_label, y_label, y_errs, series_labels, bw=bw,
               legend_loc=legend_loc, legend_title=legend_title)

    if xlogscale:
        ax.set_xscale('log')
    if ylogscale:
        ax.set_yscale('log')

    if ymin:
        if isinstance(ymin, (int, float)):
            ax.set_ylim(ymin=ymin)
        elif ymin == 'zero':
            ax.set_ylim(ymin=0)

    if title:
        ax.set_title(title, fontsize=lp.fs)

    fig.tight_layout()
    fig.savefig(ofile)
