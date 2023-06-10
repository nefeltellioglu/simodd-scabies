from scipy import stats

from disease.experiments.output_disease import *
from disease.experiments.estimate_R import *
from disease.experiments.param_sweeps import go_single


def go_final_attack_rate(obj):
    """
    An example of how an 'immediate data gathering' function should interpose between sweep and go_single
    """

    print("@@@_go_final_attack_rate")
    times = get_times(obj['p'])

    go_single(obj['p'], obj['disease_type'], obj['cmatrix'], obj['seed'], p_string=obj['p_string'])
    disease = obj['disease_type'](obj['p'], obj['cmatrix'], Random(obj['seed']),
                                  os.path.join(obj['p']['prefix'], 'disease.hd5'), mode='r')
    s_list = disease.observers['disease'].get_counts_by_state('S')
    pop_size = disease.observers['disease'].get_pop_sizes()[-1]

    output = {}
    output['final_attack_rate'] = (s_list[0] - s_list[-1]) / float(pop_size)
    output['final_shap'] = disease.observers['cases'].get_shap(times['st'], times['et'])
    disease.done()
    return output


def go_susc_cluster(obj):
    print("@@@_go_susc_cluster")
    go_single(obj['p'], obj['disease_type'], obj['cmatrix'], obj['seed'], p_string=obj['p_string'])
    bin_years = 5
    disease = obj['disease_type'](obj['p'], obj['cmatrix'], Random(obj['seed']),
                                  os.path.join(obj['p']['prefix'], 'disease.hd5'), mode='r')

    times = get_times(obj['p'])

    hh_at_risk = disease.observers['cases'].get_hh_at_risk_clustering(times['st'], times['et'],
                                                                      times['t_per_year'] * bin_years)
    series_by_hh_size = list(zip(*hh_at_risk))[2:]  # transform to be time series, skip hh of size 0 and 1
    years = (times['et'] - times['st']) / times['t_per_year']
    slopes = [stats.linregress(list(range(0, years, bin_years)), x)[0] for x in series_by_hh_size]
    print(("slopes:", slopes))
    output = {}
    output['susc_clust_slopes'] = slopes
    output['susc_clust_begin'] = hh_at_risk[0][2:]  # avg number of susceptibles during first period
    output['susc_clust_end'] = hh_at_risk[-1][2:]  # avg number of susceptibles during final period
    output['susc_clust_inc'] = [(y - x) / x for x, y in zip(hh_at_risk[0][2:], hh_at_risk[-1][2:])]

    disease.done()
    return output




