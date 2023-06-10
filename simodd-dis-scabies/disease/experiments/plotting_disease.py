## NB: this module is deprecated and will be removed soon

# all necessary (available) output functions are in output_disease and output_statistics


import numpy as np
import matplotlib.pyplot as plt

try:
    import brewer2mpl as brew

    colors = ['k'] + brew.get_map('Dark2', 'qualitative', 8).mpl_colors  # Set1
except ImportError:
    colors = ['b', 'r', 'g', 'm', 'c', 'k', '0.5']


class lp(object):
    lw = 1


## general / helper functions # - # - # - # - # - # - # - # - # - # - # - # 

def setup_x_axis_years(ax, times, offset=0, set_xlim=True):
    #print "s", offset
    years_per_tick = (times['et'] - times['st']) / (times['t_per_year'] * 5)
    # why is this divided by 10?, it produces an error if want a range of less than 10 years [can't work out why]
    if set_xlim:
        ax.set_xlim(times['st'], times['et'])
        ax.set_xticks(list(range(times['st'],
                            times['et'] + (times['t_per_year'] * years_per_tick),
                            times['t_per_year'] * years_per_tick)))
    ax.set_xticklabels(list(range(offset + times['st'] / times['t_per_year'],
                             offset + times['et'] / times['t_per_year'] + 1,
                             years_per_tick)))
    ax.set_xlabel('Year')


def setup_x_axis_days(ax, times, xmin=0):
    days_per_tick = (times['et'] - times['st']) / (times['t_per_year'] / 364.0 * 6)
    ax.set_xlim(times['st'], times['et'])
    t_per_day = times['t_per_year'] / 364.0
    xtics = np.arange(times['st'], times['et'], t_per_day * days_per_tick, dtype=np.int)
    ax.set_xticks(xtics)
    xticlabels = xtics - times['st'] + (xmin if xmin else 0)  # effectively, handle offset
    xticlabels /= t_per_day
    ax.set_xticklabels(xticlabels)
    ax.set_xlabel('Day')


def setup_x_axis(ax, times, xmin=0, set_xlim=True):
    # if duration is less than one year, show x axis as days
    if times['et'] - times['st'] <= times['t_per_year']:
        setup_x_axis_days(ax, times, xmin)
    else:
        setup_x_axis_years(ax, times, xmin, set_xlim)


## OLD functions # - # - # - # - # - # - # - # - # - # - # - # - # - # - # - #

def plot_infection(ax, series, y_errs, start_index, end_index,
                   start_time, end_time, ylabel="Proportion infected", interval=1):
    print("************* SUPERCEDED BY output_timeseries")

    xvals = np.linspace(start_time, end_time, len(series[0]))
    for i, cur_counts in enumerate(series):
        print((len(yvals)))
        if len(yvals) < len(xvals):
            yvals.extend([0] * (len(xvals) - len(yvals)))
        ax.plot(xvals, yvals, color=colors[i], lw=lp.lw)
    ax.set_ylabel(ylabel)


def plot_infection_multi(ax, inc_mean, inc_std, start_index, end_index,
                         start_time, end_time, t_dur, label="", color='k'):
    print("************* SUPERCEDED BY output_timeseries")

    yvals = list(inc_mean[start_index:end_index])
    xvals = list(range(start_time, end_time, t_dur))
    yerr = list(inc_std[start_index:end_index])
    if len(yvals) < len(xvals):
        yvals.extend([0] * (len(xvals) - len(yvals)))
    if len(yerr) < len(xvals):
        yerr.extend([0] * (len(xvals) - len(yerr)))
    ax.errorbar(xvals, yvals, label=label, yerr=yerr, color=color)
    ax.set_ylim(ymin=0.0)
    ax.set_ylabel('Incidence')


def plot_vaccine_coverage(ax, data, start_time, end_time, t_per_year, start_index,
                          logscale=False):
    print("************* SUPERCEDED BY output_timeseries")

    #    print "%%% IN plot_vaccine_coverage %%%"
    #    print start_index
    #print start_time, end_time 
    xvals = list(range(start_time, end_time, t_per_year))
    #print xvals
    #print len(xvals)
    for cur_cov in data:
        #print cur_cov[0]
        yvals = cur_cov[1][start_index:]
        #print len(yvals)
        ax.plot(xvals, yvals, label=cur_cov[0])
    if logscale: ax.set_yscale('log')
    ax.set_ylabel('Coverage')
    ax.set_ylim((0, 1))
    leg = ax.legend(loc=2)
    leg.get_frame().set_alpha(0.0)


def plot_states(ax, props, times, logscale=False, marker=-1, pop_size=-1):
    xvals = np.arange(times['st'], times['et'])

    print("************* SUPERCEDED BY output_timeseries")

    if marker > 0:
        ax.vlines(marker, 0.0005, 1, color='g')
    leg_labels = []
    for i, p in enumerate(props):
        yvals = list(p[2][times['si']:times['ei']])
        if len(yvals) < len(xvals):
            yvals.extend([0] * (len(xvals) - len(yvals)))
        ax.plot(xvals, yvals, color=colors[i])
        leg_labels.append(p[0])
    if logscale:
        ax.set_yscale('log')
        ymin = 1.0 / pop_size if pop_size > 0 else 0.0001
        ax.set_ylim((ymin, 1.0))
    leg = ax.legend(leg_labels, bbox_to_anchor=(0., 1.02, 1., .102), loc=3,
                    ncol=len(props), mode="expand", borderaxespad=0.)

    ax.set_ylabel('Fraction of population')


def plot_introductions(ax, times, start_time, end_time):
    colours = ['b', 'g', 'r', 'c', 'm', 'y', '0.25']
    for cur_size in range(1, 8):
        xvals = [x[1] for x in times if start_time < x[1] < end_time \
                                            and x[0] == cur_size or (cur_size >= 7 and x[0] >= cur_size)]
        if len(xvals) > 0:
            ax.scatter(xvals, [0.001] * len(xvals), c=colours[cur_size - 1],
                       marker='d', edgecolors='none')


def plot_hh_case_props(ax, props, start_index, end_index,
                       start_time, end_time, t_dur):
    xvals = [x for x in zip(*props)[0] if x > start_time]
    i_begin = len(props) - len(xvals)

    tprops = list(zip(*props))
    colours = ['b', 'g', 'r', 'c', 'm', 'y', '0.25']
    b = np.zeros(len(tprops[1][i_begin:]))
    for i in range(1, len(tprops)):
        ax.bar(xvals, tprops[i][i_begin:], color=colours[i - 1], bottom=b,
               align='center', linewidth=0, width=16)
        b += tprops[i][i_begin:]

    #    ax.plot(xvals, zip(*props)[1][i_begin:], c='b')#, edgecolors='none')
    #    ax.plot(xvals, zip(*props)[2][i_begin:], c='g')#, edgecolors='none')
    #    ax.plot(xvals, zip(*props)[3][i_begin:], c='r')#, edgecolors='none')
    #    ax.plot(xvals, zip(*props)[4][i_begin:], c='c')#, edgecolors='none')
    #    ax.plot(xvals, zip(*props)[5][i_begin:], c='m')#, edgecolors='none')
    #    ax.plot(xvals, zip(*props)[6][i_begin:], c='y')#, edgecolors='none')

    ax.set_ylabel('Number of cases by household size')


def plot_immunity_by_hh_size(ax, rates, start_index, end_index, start_time,
                             end_time, t_dur):
    print("************* SUPERCEDED BY output_timeseries")

    for i, cur_size in enumerate(zip(*rates)):
        xvals = list(range(start_time, end_time + 1, t_dur))
        yvals = list(cur_size[start_index:end_index])
        print(t_dur)
        print(xvals)
        print((len(xvals), len(yvals)))
        if len(yvals) < len(xvals):
            yvals.extend([yvals[-1]] * (len(xvals) - len(yvals)))
        ax.plot(xvals, yvals, color=colors[i])
    leg_labels = ['%d' % (x + 1) for x in range(len(list(zip(*rates))))]
    leg_labels[-1] += '+'
    leg = ax.legend(leg_labels, loc=3, ncol=2, title='Household size',
                    labelspacing=0.2, columnspacing=0.4, handletextpad=0.2)
    leg.get_frame().set_alpha(0.0)
    plt.setp(leg.get_texts(), fontsize='small')


def plot_inc_over_time(ax, data, age_bins, times, ylabel, ymax=-1, units='years'):
    leg_labels = ['%d-%d %s' % (x, y, units) \
                  for x, y in zip(age_bins[:-1], age_bins[1:])]
    for i, x in enumerate(zip(*data)):
        x_val = list(range(times['st'], times['et'], (times['et'] - times['st']) / len(x)))
        print((times['st'], times['et'], len(x_val), len(x)))
        print(x_val)
        ax.plot(x_val, x, color=colors[i % len(colors)])
    if ymax > 0:
        ax.set_ylim(ymax=ymax)
    leg = ax.legend(leg_labels)  #, prop={'size':'large'})
    leg.get_frame().set_alpha(0.0)
    ax.set_ylabel(ylabel)


def plot_hh_fracs(ax, data, c='b', t_step=1):
    #    print data
    d = list(zip(*data))
    x = np.array(d[0]) * t_step
    y = np.array(d[1])
    #    print np.mean(y)
    ax.plot(x, y, color=c)
    ax.set_ylabel("Household fraction")
    #    ax.set_ylim(ymax=0.5)
    ax.set_ylim(ymin=0)


def plot_susc_prop(ax, data, c, times):
    # plots a single output series
    print(times)
    xvals = list(range(times['st'], times['et'] + 1, 1))  #364/times['t_per_year'])
    yvals = data
    print((len(xvals), len(yvals)))
    if len(yvals) < len(xvals):
        yvals.extend([yvals[-1]] * (len(xvals) - len(yvals)))
    ax.plot(xvals, yvals, color=c, marker='o', markeredgewidth=0, markersize=3)
    ax.set_yscale('log')
    ax.set_ylim(ymin=1e-4, ymax=1)


def plot_susc_props(ax, data, times):
    # plots a single output figure
    for i in range(0, max(data.keys()) + 1):
        plot_susc_prop(ax, data[i], colors[i], times)
    leg = ax.legend(list(range(0, max(data.keys()) + 1)), title='Number Susceptible')  #, prop={'size':'large'})
    leg.get_frame().set_alpha(0.0)


def plot_hh_fracs_by_hh_size(ax, data, t_step=1):
    for i in range(1, 8):
        if data[i]:
            plot_hh_fracs(ax, data[i], colors[i - 1], t_step=t_step)


def plot_hh_fracs_multi(ax, hhf_mean, hhf_sd, label="", color='b'):
    ax.set_xlabel("Days")
    ax.set_ylabel("Household fraction")
    xvals = list(range(len(hhf_mean)))
    if hhf_sd:
        ax.errorbar(xvals, hhf_mean, yerr=hhf_sd, label=label, color=color)
    else:
        ax.plot(xvals, hhf_mean, label=label, color=color)
    ax.set_ylim(ymin=0.0)


def plot_infection_immunity(ax, counts, ab, start_index, start_time=0,
                            points_per_year=52, years_per_tick=5):
    plot_infection(ax, counts, start_index, start_time,
                   points_per_year, years_per_tick)
    ax2 = ax.twinx()
    ax2.plot(list(range(start_index + start_time, start_time + len(ab))),
             ab[start_index:])
    ax2.set_ylabel('Mean antibody level')
    ax.set_xlim(xmin=start_index + start_time, xmax=start_time + len(counts))
    ax.set_xticks(list(range(start_index + start_time, start_time + len(counts),
                        points_per_year * years_per_tick)))
    ax.set_xticklabels(list(range(0, len(counts) / points_per_year, years_per_tick)))


def plot_cases(ax, counts, bins, xlabel):
    ax.bar(bins[:-1] - 0.4, counts, width=0.8)
    ax.set_xlabel(xlabel)
    ax.set_ylabel('Cases')
    ax.set_xlim(xmin=-0.5)


def plot_incidence_line(ax, counts, errs, xlabel, label=None, color='b'):
    ax.plot(list(range(len(counts))), counts, label=label, lw=1.5, color=color)
    if errs is not None:
        ax.errorbar(list(range(len(counts))), counts, yerr=errs, color=color, capsize=0, alpha=0.4)
    ax.set_xlabel(xlabel)
    ax.set_ylabel('Incidence')
    ax.set_ylim(ymin=0.0)
    ax.set_xlim(xmax=100)


def plot_incidence(ax, counts, state_labels, bins, xlabel,
                   labels=None, y_max=-1, x_max=-1):
    b = np.zeros(len(counts[state_labels[0]]))
    for i, label in enumerate(state_labels):
        if label not in counts: continue
        ax.bar(bins[:-1], counts[label], color=colors[i], bottom=b,
               align='center', linewidth=0, width=0.8)
        b += counts[label]

    #    ax.bar(bins[:-1]-0.4, counts, width=0.8)
    ax.set_xlabel(xlabel)
    ax.set_ylabel('Incidence')
    ax.set_xlim(xmin=bins[0] - 0.5)
    if x_max > 0:
        ax.set_xlim(xmax=x_max)
    else:
        ax.set_xlim(xmax=bins[-1] - 0.5)

    if labels:
        ax.set_xticks(list(range(len(labels))))
        ax.set_xticklabels(labels)

    if y_max > 0:
        ax.set_ylim(ymax=y_max)


def plot_age_inc_diffs(ax, counts):
    ax.bar(list(range(len(counts))), counts, align='center', linewidth=0, width=0.8)


def plot_inc_total_sweep(ax, means, stdevs, param_values):
    xvals = list(range(len(param_values)))
    ax.bar(xvals, means, yerr=stdevs, align='center')
    ax.set_xticks(xvals)
    ax.set_xticklabels([str(x) for x in param_values])
    ax.set_ylabel('Cases')


def plot_ab_snapshot(ax, snapshot, binsize=5):
    binned = []
    cur_bin = []
    for i, x in list(snapshot.items()):
        #        print i, x
        cur_bin.extend(x)
        if i % binsize == 0:
            binned.append(cur_bin)
            cur_bin = []
        #    print binned

    ax.boxplot(binned, positions=list(range(0, len(binned))))  #, widths=3)
    ax.set_xlabel('Age')
    ax.set_ylabel('Antibody level')


def bin_ab_level(snapshot, cutoffs=(5, 62.5, 125),
                 bins=(1, 2, 3, 4, 5, 6, 7, 8, 9, 14, 19, 24, 34, 44, 59, 74)):
    cutoffs = list(cutoffs) + [10000]
    bins = list(bins) + [2000]
    counts = np.zeros((len(bins), len(cutoffs)))
    for cur_age, ab_levels in list(snapshot.items()):
        i = 0
        while cur_age >= bins[i]:
            i += 1
        for ab in ab_levels:
            j = 0
            while ab >= cutoffs[j]:
                j += 1
            counts[i][j] += 1

    return counts


def build_ab_legend_labels(cutoffs):
    """
    Create a list of series labels based upon cutoffs as follows:
    ['<c[0]', 'c[0]--c[1]-1', ..., 'c[n-1]--c[n]-1', 'c[n]+']
    """

    labels = ['<%g EU/mL' % cutoffs[0]]
    labels += ['%g to <%g EU/mL' % (x, y) for x, y in zip(cutoffs[:-1], cutoffs[1:])]
    labels += [''.join(('\\u2265', '%g EU/mL' % cutoffs[-1]))]
    return labels


def build_tick_labels(bins):
    labels = ['%d' % bins[0]]
    for x, y in zip(bins[:-1], bins[1:]):
        if y - x == 1:
            labels.append('%d' % y)
        else:
            labels.append('%d-%d' % (x + 1, y))
    labels.append('%d+' % (bins[-1] + 1))
    return labels


def plot_state_snapshot(ax, props, time, state_labels, age_cutoffs):
    x = np.arange(len(age_cutoffs) + 1)
    plots = []
    b = np.zeros(len(x))
    for i, label in enumerate(state_labels):
        plots.append(ax.bar(x, props[label], color=colors[i],
                            bottom=b, align='center', lw=0))
        b += props[label]
    #    legend([p[0] for p in plots], state_labels, bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
    #    leg = ax.legend([p[0] for p in plots], state_labels)
    ax.set_xlim(-0.4, len(x) - 0.6)
    ax.set_xticks(list(range(len(x))))
    ax.xaxis.set_ticks_position('none')
    ax.set_xticklabels(build_tick_labels(age_cutoffs), size='x-small')
    ax.set_xlabel('Age group (years)')
    ax.set_ylim((0, 100))
    ax.set_ylabel('Percentage')
    ax.set_yticklabels(['%d%%' % z for z in range(0, 120, 20)], size='x-small')
    ax.set_title('year %03d day %03d' % (time / 365, time % 365))


def plot_ab_snapshot_bars(ax, counts, time, cutoffs=(5, 62.5, 125),
                          bins=(1, 2, 3, 4, 5, 6, 7, 8, 9, 14, 19, 24, 34, 44, 59, 74),
                          colours=('#ff0000', '#ffa500', '#ffff00', '#9acd32')):
    x = np.arange(len(counts))
    plots = []
    b = np.zeros((len(counts)))
    for i in range(len(cutoffs) + 1):
        plots.append(ax.bar(x, counts[:, i], color=colours[i], bottom=b, align='center'))
        b += counts[:, i]
    leg = ax.legend([p[0] for p in plots], build_ab_legend_labels(cutoffs), mode='expand', ncol=4,
                    bbox_to_anchor=(0., -.16, 1., .102), loc=8, handletextpad=0.3, markerscale=0.5)
    leg.get_frame().set_alpha(0.0)
    plt.setp(leg.get_texts(), fontsize='x-small')
    ax.set_xticks(list(range(len(bins) + 1)))
    ax.xaxis.set_ticks_position('none')
    ax.set_xticklabels(build_tick_labels(bins), size='x-small')
    ax.set_xlim((-0.4, len(bins) + 0.4))
    ax.set_xlabel('Age group (years)')
    ax.set_yticklabels(['%d%%' % x for x in range(0, 120, 20)], size='x-small')
    ax.set_ylim((0, 100))
    ax.set_ylabel('Percentage positive')
    ax.set_title('Antibody levels by age group (day %d)' % time)


def plot_hh_risk(ax, data, binsize=5):
    """
    Plot distribution of household sizes (using boxplot) by age of household.
    """

    # bin data into five year intervals
    binned = []
    cur_bin = []
    for i, x in list(data.items()):
        cur_bin.extend(x)
        if i % binsize == 0:
            binned.append(cur_bin)
            cur_bin = []

    r = ax.boxplot(binned, positions=list(range(0, len(data), binsize)), widths=3)
    for k in list(r.keys()):
        plt.setp(r[k], color='black')
    plt.setp(r['medians'], color='red')
    #    ax.plot(range(len(data)), [np.mean(x) for x in data], color='k', lw=lp.lw)
    ax.set_xlabel('Household age')
    ax.set_ylim(ymin=0, ymax=6)
    ax.set_xticks(list(range(0, len(data), 5)))
    ax.set_ylabel('Number at risk')


def plot_hh_risk_series(ax, data):
    for cur_hh in list(data.values()):
        ax.plot(list(range(len(cur_hh))), cur_hh)  #, color=colors[max(cur_hh)])
