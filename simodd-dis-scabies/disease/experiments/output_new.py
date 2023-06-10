import collections

from matplotlib.patches import Rectangle
import numpy as np
import matplotlib as plt


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
def get_distinct(nr):
    # check if nr is in correct range
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


## validation ###

def validate_data(x_values, y_values, y_errs, series_labels):
    # if y_values are submitted as a single series, wrap in list
    if not (isinstance(y_values[0], collections.Iterable)):
        y_values = [y_values]

    assert len(x_values) == len(y_values[0])

    # if y_errs are submitted as a single series, wrap in list
    if y_errs is not None and not isinstance(y_errs[0], collections.Iterable):
        y_errs = [y_errs]

    # expand y_errs and series_labels if not present
    if not y_errs:
        y_errs = [None] * len(y_values)
    if not series_labels:
        series_labels = [None] * len(y_values)

    return x_values, y_values, y_errs, series_labels


## customisation ###

def setup_plot(ax, x_label, y_label, **kwargs):
    fs = kwargs.get('fontsize', 'normal')

    # set title
    if 'plot_title' in kwargs:
        ax.set_title(kwargs['plot_title'], fontsize=fs)

    # set axis labels
    ax.set_xlabel(x_label, fontsize=fs)
    ax.set_ylabel(y_label, fontsize=fs)

    # set axis limits
    if 'y_min' in kwargs:
        ax.set_ylim(ymin=kwargs['y_min'])
    if 'y_max' in kwargs:
        ax.set_ylim(ymax=kwargs['y_max'])

    if 'x_min' in kwargs:
        ax.set_ylim(ymin=kwargs['x_min'])
    if 'x_max' in kwargs:
        ax.set_ylim(ymax=kwargs['x_max'])

    # set tick label sizes
    plt.setp(ax.get_xticklabels(), fontsize=fs)
    plt.setp(ax.get_yticklabels(), fontsize=fs)


def setup_legend(ax, series_labels, **kwargs):
    if series_labels[0]:
        leg = ax.legend(
            scatterpoints=1, numpoints=1, ncol=2,
            labelspacing=0.1, columnspacing=1,
            loc=kwargs.get('legend_loc', 'upper right'),
            title=kwargs.get('legend_title', None))
        leg.get_frame().set_alpha(0.0)


def remove_spines(ax, spines=('top', 'right')):
    # Remove specified axes lines ('spines')
    for spine in spines:
        ax.spines[spine].set_visible(False)
    ax.yaxis.set_ticks_position('none')
    ax.xaxis.set_ticks_position('none')


## helper functions
def get_x_values_from_times(times, y_values, align=None):
    # move alignment of x_values outside of plotting...
    if align == 'centre':
        offset = (times['et'] - times['st']) / (len(y_values[0]) + 1) / 2.0
        x_values = np.linspace(times['st'] + offset, times['et'] - offset, len(y_values[0])) / times['t_per_year']
    else:
        x_values = np.linspace(times['st'], times['et'], len(y_values[0])) / times['t_per_year']
    return x_values


def set_x_axis_from_times(ax, times, x_offset):
    pass


## plotting functions ###

def plot_lines(ax, x_values, y_values, y_errs, series_labels, x_label, y_label, **kwargs):
    x_values, y_values, y_errs, series_labels = validate_data(x_values, y_values, y_errs, series_labels)

    colors = get_distinct(len(y_values))

    for i, (cur_values, cur_errs, cur_label) in enumerate(zip(y_values, y_errs, series_labels)):
        if cur_errs is not None:
            ax.errorbar(x_values, cur_values, yerr=cur_errs, label=cur_label, color=colors[i % len(colors)],
                        lw=kwargs.get('lw', 1))
        else:
            ax.plot(x_values, cur_values, label=cur_label, color=colors[i % len(colors)],
                    lw=kwargs.get('lw', 1))

    setup_plot(x_label, y_label, kwargs)


def plot_bars(ax, x_values, y_values, y_errs, series_labels, x_label, y_label, **kwargs):
    x_values, y_values, y_errs, series_labels = validate_data(x_values, y_values, y_errs, series_labels)

    colors = get_distinct(len(y_values))

    width = 0.8 / len(y_values)
    x_values = np.array(x_values)
    for i, (cur_values, cur_errs, cur_label) in enumerate(zip(y_values, y_errs, series_labels)):
        ax.bar(x_values + i * width, cur_values, width, lw=0, edgecolor='grey', color=colors[i], ecolor='dimgrey',
               align='center', yerr=cur_errs, label=cur_label)

    setup_plot(x_label, y_label, kwargs)

    if 'x_ticklabels' in kwargs:
        ax.set_xticks(np.array(x_values) + len(y_values) / 2 * width)
        ax.set_xticklabels(tuple(kwargs['x_ticklabels']), rotation=kwargs.get('rotate_ticks', 0))


def plot_stacked_area(ax, x_values, y_values, series_labels, x_label, y_label, **kwargs):
    x_values, y_values, y_errs, series_labels = validate_data(x_values, y_values, None, series_labels)

    colors = get_distinct(len(y_values))

    p = []
    ax.fill_between(x_values, 0, y_values[0], lw=0, edgecolor='grey', facecolor=colors[0])
    p.append(Rectangle((0, 0), 1, 1, lw=0, edgecolor='grey', facecolor=colors[0]))  # for legend
    for i, (cur_bottom, cur_top) in enumerate(zip(y_values[:-1], y_values[1:])):
        ax.fill_between(x_values, cur_bottom, cur_top, lw=0, edgecolor='grey', facecolor=colors[i + 1])
        p.append(Rectangle((0, 0), 1, 1, lw=0, edgecolor='grey', facecolor=colors[i + 1]))

    setup_plot(x_label, y_label, kwargs)

    # setup legend
    if series_labels[0]:
        leg = ax.legend(
            p, series_labels,
            ncol=2, title=kwargs.get('legend_title', None),
            labelspacing=0.1, columnspacing=1, loc='upper right')
        plt.setp(leg.get_texts(), fontsize='small')
        leg.get_frame().set_alpha(0.8)
        leg.get_frame().set_linewidth(0.0)


def plot_points_lin_reg(ax, x_values, y_values, series_labels, x_label, y_label, **kwargs):
    x_values, y_values, y_errs, series_labels = validate_data(x_values, y_values, None, series_labels)

    colors = get_distinct(len(y_values))

    for i, (cur_values, cur_label) in enumerate(zip(y_values, series_labels)):
        p = np.polyfit(x_values[~np.isnan(cur_values)], cur_values[~np.isnan(cur_values)], 1)
        ax.plot(x_values, cur_values, color=colors[i], marker='x', lw=0)
        ax.plot(x_values, p[0] * x_values + p[1], label=cur_label, color=colors[i], lw=lp.lw)

    setup_plot(x_label, y_label, kwargs)


def plot_scatter(ax, x_values, y_values, x_label, y_label, **kwargs):
    pass


# nb: differs as multiple series will need unique x values as well as y...


def plot_heatmap(ax, **kwargs):
    pass


def plot_contours(ax, **kwargs):
    pass
