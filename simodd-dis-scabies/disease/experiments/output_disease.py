from collections import Iterable
from matplotlib.colors import LogNorm

from .output_statistics import *



## helper / processing functions # - # - # - # - # - # - # - # - # - # - #

def get_times(p, start_time=0, end_time=0, days_per_tick=50):
    """
    si and ei are start_index and end_index, i.e., start and end times in units of *timesteps*
    st and et are start_time and end_time, i.e., start and end times in units of *days*
    """
    start_year = 0
    final_year = p['years'][-1] if end_time is 0 else p['years'][-2] + end_time / 364.0
    #    print "final_year", final_year
    end_year = final_year - p['years'][-2]
    times = {
        'si': int(p['t_per_year'] * start_year),
        'ei': int(p['t_per_year'] * end_year),
        'st': int((p['burn_in'] + p['epi_burn_in'] + p['years'][-2]) * p['t_per_year']) - start_time,
        'et': int((p['burn_in'] + p['epi_burn_in'] + final_year) * p['t_per_year']) - start_time,
        't_per_year': p['t_per_year'],
        'years_per_tick': 1,
        'days_per_tick': days_per_tick
    }
    #    print p['burn_in'], p['epi_burn_in'], p['years']
    #    print times
    return times


def setup_x_axis_years(ax, times, offset=0, set_xlim=True, years_per_tick=20):
    #print "s", offset
    if not years_per_tick:
        years_per_tick = (times['et'] - times['st']) // (times['t_per_year'] * 5)
    if set_xlim:
        ax.set_xlim(times['st'], times['et'])
        ax.set_xticks(list(range(times['st'],
                            times['et'] + (times['t_per_year'] * years_per_tick),
                            times['t_per_year'] * years_per_tick)))
    ax.set_xticklabels(list(range(int(offset + times['st'] // times['t_per_year']),
                             int(offset + times['et'] // times['t_per_year']) + 1,
                             int(years_per_tick))))
    ax.set_xlabel('Year')


def setup_x_axis_days(ax, times, xmin=0):
    days_per_tick = (times['et'] - times['st']) // (times['t_per_year'] / 364.0 * 6)
    ax.set_xlim(times['st'], times['et'])
    t_per_day = times['t_per_year'] / 364.0
    xtics = np.arange(times['st'], times['et'], t_per_day * days_per_tick, dtype=np.int)
    ax.set_xticks(xtics)
    xticlabels = xtics - times['st'] + (xmin if xmin else 0)  # effectively, handle offset
    xticlabels /= t_per_day
    ax.set_xticklabels(xticlabels)
    ax.set_xlabel('Day')


def setup_x_axis(ax, times, xmin=0, set_xlim=True):
    # if duration is less than two years, show x axis as days
    if times['et'] - times['st'] <= 2 * times['t_per_year']:
        setup_x_axis_days(ax, times, xmin)
    else:
        setup_x_axis_years(ax, times, xmin, set_xlim)


def convert_counts_to_props(bins, cutoffs):
    # calculate sum for each age category, and convert counts to proportions
    sums = np.zeros(len(cutoffs))
    for st in list(bins.values()):
        for i, v in enumerate(st):
            sums[i] += v
    props = {}
    for k in bins:
        props[k] = bins[k] / sums * 100
    return props


def output_timeseries(ofile, times, series, ylabel, y_errs=None, labels=None,
                      logscale=False, xmin=None, ymin=None, ymax=None, ylines=None, x_offset=0, legend_title=None):
    fig = plt.figure()
    ax = fig.add_subplot(111)

    plot_timeseries_ax(ax, times, series, ylabel, y_errs, labels, logscale,
                       xmin, ymin, ymax, ylines, x_offset=x_offset, legend_title=legend_title)

    fig.tight_layout()
    fig.savefig(ofile)


def plot_timeseries_ax(ax, times, series, ylabel, y_errs=None, labels=None,
                       logscale=False, xmin=None, ymin=None, ymax=None, ylines=None, align=None,
                       x_offset=0, legend_title=None, legend_loc='upper left', plot_title=None,
                       linestyle='-', lineweight=None, bw=False, br=False, reversed=False):
    # wrap series if only a single list passed
    if not isinstance(series[0], Iterable):
        series = [series]

    print(('xoff', x_offset))
    # compute xvalues (possibly also allow to passed explicitly?)
    if align == 'centre':
        offset = (times['et'] - times['st']) / (len(series[0]) + 1) / 2.0
        #print times['st'], times['et'], offset
        xvals = np.linspace(times['st'] + offset, times['et'] - offset, len(series[0]))
        #print xvals
    else:
        xvals = np.linspace(times['st'], times['et'], len(series[0]))

    plot_lines(ax, xvals, series, '', ylabel, y_errs, labels, legend_title, legend_loc, linestyle, lineweight, bw, br, reversed)

    if ymax:
        ax.set_ylim(ymax=ymax)
    if ymin is not None:
        ax.set_ylim(ymin=ymin)
    if xmin:
        times['st'] += xmin

    if ylines:
        for cur_y in ylines:
            ax.hlines(cur_y, ax.get_xlim()[0], ax.get_xlim()[1], colors='k', linestyle='dotted')

    #print x_offset
    setup_x_axis(ax, times, x_offset)
    if plot_title:
        ax.set_title(plot_title, fontsize=16)

    if logscale:
        ax.set_yscale('log')


def output_stacked_bars_timeseries(ofile, times, series, ylabel, labels=None):
    fig = plt.figure()
    ax = fig.add_subplot(111)

    plot_stacked_bars_timeseries_ax(ax, times, series, ylabel, labels=labels)

    fig.savefig(ofile)


def plot_stacked_bars_timeseries_ax(ax, times, series, ylabel, x_offset=0, labels=None, legend_title=None,
                                    plot_title=None):
    xvals = np.linspace(times['st'], times['et'], len(series[0]))

    y = np.row_stack(series)
    y_stack = np.cumsum(y, axis=0)

    plot_stacked_area(ax, xvals, y_stack, '', ylabel, labels, legend_title, plot_title)

    setup_x_axis(ax, times, x_offset)


def output_general(ofile, xvals, yvals, xlabel, ylabel, y_errs=None, labels=None,
                   logscale=False, title=None, xmin=None, xmax=None, ymin=None, ymax=None, ylines=None):
    """
    General plotting function almost identical to output_timeseries except can plot over any xvals so
    doesn't need the special 'times' passed in as input
    """

    fig = plt.figure()
    ax = fig.add_subplot(111)

    plot_lines(ax, xvals, yvals, xlabel, ylabel, y_errs, labels)

    if xmin is not None:
        ax.set_xlim(xmin=xmin)
    if xmax:
        ax.set_xlim(xmax=xmax)
    if ymin is not None:
        ax.set_ylim(ymin=ymin)
    if ymax:
        ax.set_ylim(ymax=ymax)

    if ylines:
        for cur_y in ylines:
            ax.hlines(cur_y, ax.get_xlim()[0], ax.get_xlim()[1], colors='k', linestyle='dotted')

    if logscale:
        ax.set_yscale('log')

    if title:
        ax.set_title(title)

    fig.savefig(ofile)


def output_scatter(ofile, xvals, yvals, xlabel, ylabel,
                   logscale=False, title=None, xmin=None, xmax=None, ymin=None, ymax=None, ylines=None):
    """
    Function to output a scatter plot of xvals and yvals which are both lists of equal length. Doesn't
    include ability to plot error bars (as for single points)

    Nic might have a function that does the same thing
    """

    fig = plt.figure()
    ax = fig.add_subplot(111)

    ax.scatter(xvals, yvals)

    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)

    if xmin is not None:
        ax.set_xlim(xmin=xmin)
    if xmax:
        ax.set_xlim(xmax=xmax)
    if ymin is not None:
        ax.set_ylim(ymin=ymin)
    if ymax:
        ax.set_ylim(ymax=ymax)

    if ylines:
        for cur_y in ylines:
            ax.hlines(cur_y, ax.get_xlim()[0], ax.get_xlim()[1], colors='k', linestyle='dotted')

    if logscale:
        ax.set_yscale('log')

    if title:
        ax.set_title(title)

    fig.savefig(ofile)


def output_contour(data, x_param, y_param, ofile, zlog=False, contour_values=None, tgt_value=None, plot_values=None,
                   title=None):
    """
    Produce a contour map and, incidentally, calculate the y values corresponding to a target
    contour line for each value of x_param.

    TODO: These two functionalities should be separated.
    """

    x_values = x_param['values']
    y_values = y_param['values']

    fig = plt.figure()
    ax = fig.add_subplot(111)
    norm = LogNorm() if zlog else None

    if not contour_values: contour_values = 15

    img = ax.contourf(np.array(x_values), np.array(y_values), np.array(data), contour_values, cmap=plt.cm.rainbow,
                      vmax=abs(data).max(), vmin=-abs(data).max(), norm=norm)
    fig.colorbar(img)

    img = ax.contour(np.array(x_values), np.array(y_values), np.array(data), contour_values, linewidths=0.5, colors='k')
    ax.set_xlabel(x_param['name'])
    ax.set_ylabel(y_param['name'])

    fig.savefig(ofile)


def output_prop_detected(data, times, ofile, mark_times=None, labels=None, y_max=None):
    fig = plt.figure()
    ax = fig.add_subplot(111)
    for i, (props, marks, cur_label) in enumerate(zip(data, mark_times, labels)):
        ax.plot(list(range(times['st'], times['et'])), props, color=colors[i + 1], lw=lp.lw, label='.')
        #print mark_times
        #print [data[x-times['st']] for x in mark_times]
        if mark_times:
            ax.plot(marks, [props[x - times['st']] for x in marks], 'o', color=colors[i + 1], ms=5)  # , mew=0)
        for mt in marks:
            l = ax.axvline(mt, ls='--', ymax=props[mt - times['st']] * (1 / y_max if y_max else 1),
                           color=colors[i + 1], alpha=0.8)
            l.set_dashes([1, 1])
    ax.set_ylim((0, 1))
    if y_max:
        ax.set_ylim(ymax=y_max)

    setup_x_axis_days(ax, times)
    ax.set_ylabel('Cumulative proportion of cases detected')
    if labels:
        h, l = ax.get_legend_handles_labels()
        leg = ax.legend(h, labels, loc='lower right')
        leg.get_frame().set_alpha(0.0)
    fig.savefig(ofile)


def output_prop_detected_case(data, ofile):
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.plot(list(range(1, len(data) + 1)), data)
    ax.set_xlabel('Detected cases')
    ax.set_ylabel('Cumulative fraction of cases detected')
    fig.savefig(ofile)
