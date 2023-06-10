## BELOW HERE STORED FOR ARCHIVAL PURPOSES...
## snapshot # - # - # - # - # - # - # - # - # - # - # - # - # - # - # - # - #

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


def output_state_snapshot(data,
                          start_index, end_index, start_time, end_time,
                          t_dur, years_per_tick, logscale,
                          snapshot, time, state_labels, age_cutoffs, ofile):
    """
    Plots four graphs:
    1. disease states by age (bars)
    2. disease states over time (current time marked by line)
    3. disease states of household containing children
    4. disease states of households not containing children

    Can be stitched together to make movie of outbreak.
    """


    # for each snapshot, bin individiuals by age and state
    cutoffs = list(age_cutoffs) + [2000]
    all_bins = {}
    kids_bins = {}
    nokids_bins = {}
    for label in state_labels:
        all_bins[label] = np.zeros(len(cutoffs))
        kids_bins[label] = np.zeros(len(cutoffs))
        nokids_bins[label] = np.zeros(len(cutoffs))
    for cur_ind in snapshot:
        index = 0
        while cur_ind['age'] >= cutoffs[index]:
            index += 1
        all_bins[cur_ind['state'].label][index] += 1
        if cur_ind['hh_type'].endswith('kids'):
            kids_bins[cur_ind['state'].label][index] += 1
        else:
            nokids_bins[cur_ind['state'].label][index] += 1

    all_props = convert_counts_to_props(all_bins, cutoffs)
    kids_props = convert_counts_to_props(kids_bins, cutoffs)
    nokids_props = convert_counts_to_props(nokids_bins, cutoffs)

    # plot figure
    fig = plt.figure()
    ax = fig.add_subplot(221)
    plot_state_snapshot(ax, all_props, time, state_labels, age_cutoffs)
    ax = fig.add_subplot(222)
    plot_states(ax, data, start_index, end_index, start_time, end_time,
                t_dur, logscale, marker=time)
    #    setup_x_axis(ax, start_time, end_time, years_per_tick)
    ax = fig.add_subplot(223)
    plot_state_snapshot(ax, kids_props, time, state_labels, age_cutoffs)
    ax = fig.add_subplot(224)
    plot_state_snapshot(ax, nokids_props, time, state_labels, age_cutoffs)
    fig.savefig(ofile)


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


def plot_ab_snapshot_bars(ax, counts, time, cutoffs=(5, 62.5, 125),
                          bins=(1, 2, 3, 4, 5, 6, 7, 8, 9, 14, 19, 24, 34, 44, 59, 74),
                          colours=('#ff0000', '#ffa500', '#ffff00', '#9acd32')):
    x = np.arange(len(counts))
    plots = []
    b = np.zeros((len(counts)))
    for i in range(len(cutoffs) + 1):
        plots.append(ax.bar(x, counts[:, i], color=colours[i], bottom=b, align='center'))
        b += counts[:, i]
    leg = ax.legend([p[0] for p in plots], build_ab_legend_labels(cutoffs), mode='expand',
                    ncol=4, bbox_to_anchor=(0., -.16, 1., .102), loc=8, handletextpad=0.3, markerscale=0.5)
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


def output_ab_snapshot(ss, ofile, cutoffs=(0.2, 0.5, 0.8)):
    """
    Plot a snapshot of antibody levels (by age?)
    """

    counts = bin_ab_level(ss[1], cutoffs)
    sums = np.array([[float(sum(x))] for x in counts])
    props = counts / sums * 100

    fig = plt.figure()
    ax = fig.add_subplot(111)
    plot_ab_snapshot_bars(ax, props, ss[0], cutoffs)
    fig.savefig(ofile)


def output_immunity_ratio(snapshots, state_labels, age_cutoffs, ofile):
    cutoffs = list(age_cutoffs) + [2000]

    kids_bins = {}
    nokids_bins = {}

    ratios = []

    for snapshot in snapshots:
        for label in state_labels:
            kids_bins[label] = np.zeros(len(cutoffs))
            nokids_bins[label] = np.zeros(len(cutoffs))

        for cur_ind in snapshot[1]:
            index = 0
            while cur_ind['age'] >= cutoffs[index]:
                index += 1
            if cur_ind['hh_type'].endswith('kids'):
                kids_bins[cur_ind['state'].label][index] += 1
            else:
                nokids_bins[cur_ind['state'].label][index] += 1

        kids_props = convert_counts_to_props(kids_bins, cutoffs)
        nokids_props = convert_counts_to_props(nokids_bins, cutoffs)

        cur_ratios = [x / y for x, y in zip(kids_props['R'], nokids_props['R'])]

        ratios.append(cur_ratios)

    of = file(ofile, 'w')
    for row in ratios:
        of.write(','.join([str(x) for x in row]))
        of.write('\n')
    of.close()


def output_heatmap(data, x_param, y_param, ofile, zlog=False):
    # compute boundaries of each heatmap cell such that actual parameter
    # value occupies the centre of the cell
    # NB: may not be the most legitimate way to display this info... (?)
    x_values = [np.mean([x, y]) \
                for x, y in zip(x_param['values'][:-1], x_param['values'][1:])]
    x_values.insert(0, x_values[0] - (x_values[1] - x_values[0]))
    x_values.append(x_values[-1] + (x_values[-1] - x_values[-2]))
    y_values = [np.mean([x, y]) \
                for x, y in zip(y_param['values'][:-1], y_param['values'][1:])]
    y_values.insert(0, y_values[0] - (y_values[1] - y_values[0]))
    y_values.append(y_values[-1] + (y_values[-1] - y_values[-2]))

    X, Y = np.meshgrid(np.array(x_values), np.array(y_values))

    fig = plt.figure()
    ax = fig.add_subplot(111)
    norm = LogNorm() if zlog else None
    img = ax.pcolor(X, Y, np.array(data), cmap='jet', norm=norm)
    ax.set_xlabel(x_param['name'])
    ax.set_ylabel(y_param['name'])

    ax.set_xlim((x_values[0], x_values[-1]))
    ax.set_ylim((y_values[0], y_values[-1]))

    fig.colorbar(img)
    fig.savefig(ofile)


def get_y_err(contour, cur):
    pass


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

    y_fits = {}

    # if contour values and target value have been specified, find corresponding y
    # values for each x value
    contour_values = img.levels
    print((type(contour_values)))
    print(contour_values)
    print(tgt_value)
    if tgt_value in contour_values:
        i_tgt = [contour_values.tolist().index(tgt_value)]

        if not plot_values: plot_values = x_values

        for i, cur_x in enumerate(plot_values):
            # fn returns the distance between the current x, y location and the contour
            # line corresponding to the target R0 value
            fn = lambda cur_y: img.find_nearest_contour(cur_x, cur_y, i_tgt, pixel=False)[-1]
            best_y = fmin(fn, [0.1])
            y_fits[cur_x] = best_y[0]

            # plot intersections
            ax.plot((cur_x, cur_x), (0.005, best_y), 'w')
            ax.plot((0, cur_x), (best_y, best_y), 'w')

    if title:
        ax.set_title(title)

    fig.savefig(ofile)

    return y_fits


def output_outbreak_sizes(data, R0_values, ofile):
    # plot the distribution of epidemic outbreak sizes for various values of R0

    fig = plt.figure()
    ax = fig.add_subplot(111)

    for cur_R0_value, cur_data in zip(R0_values, data):
        print(cur_data)
        hist, bin_edges = np.histogram(cur_data, bins=101, range=(0.0, 1.0))
        hist = np.array(hist, dtype=np.float32) / float(sum(hist))
        print(hist)
        print(bin_edges)
        ax.plot([(x + y) / 2.0 for x, y in zip(bin_edges[:-1], bin_edges[1:])],
                hist, label='R0=%g' % cur_R0_value)

    leg = ax.legend()
    ax.set_xlabel('Outbreak size (proportion of population)')
    ax.set_ylabel('Proportion of runs')
    leg.get_frame().set_alpha(0.0)
    fig.savefig(ofile)
