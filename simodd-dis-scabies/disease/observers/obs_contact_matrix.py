"""
Observer object for contact matrices
"""
import tables as tb
import os

import numpy as np
import matplotlib.pyplot as plt

from disease.observers.obs_base import Observer


class ContactMatrixObserver(Observer):
    def __init__(self, h5file, max_hh_size=10, age_classes=None, t_interval_basic=1, t_interval_derived=1,
                 sample_size=100):
        self.max_hh_size = max_hh_size
        self.max_age = 101
        self.age_classes = age_classes
        self.t_interval_basic = t_interval_basic
        self.t_interval_derived = t_interval_derived
        self.sample_size = sample_size

        desc = {
            't': tb.UInt32Col(),
            # basic matrices
            'hh': tb.Float32Col(shape=(self.max_age, self.max_age)),
            'ext': tb.Float32Col(shape=(self.max_age, self.max_age)),
            'comm': tb.Float32Col(shape=(len(self.age_classes), len(self.age_classes))),
            # derived matrices
            'hh_by_age': tb.Float32Col(shape=(self.max_hh_size, self.max_age)),
            'age_by_hh': tb.Float32Col(shape=(self.max_age, self.max_hh_size)),
            'hh_by_hh': tb.Float32Col(shape=(self.max_hh_size, self.max_hh_size)),
            'hh_by_hh_scale': tb.Float32Col(shape=(self.max_hh_size, self.max_hh_size))
        }
        super(ContactMatrixObserver, self).__init__(h5file, 'contact', desc, 'Contact Matrix Observer')

    def update(self, t, disease, pop, rng, **kwargs):
        if t % self.t_interval_basic > 0:
            return
        self.row['t'] = t
        self.row['hh'] = self.get_hh_matrix(pop, age_scale=False)
        self.row['ext'] = self.get_extended_family_contact_matrix(pop)
        self.row['comm'] = disease.cmatrix.C

        if t % self.t_interval_derived > 0:
            self.row.append()
            self.h5file.flush()
            return

        hh_cm, hh_cm_scale, hh_age_cm, age_hh_cm = self.get_derived_contact_matrices(pop, disease.cmatrix, rng)
        self.row['hh_by_age'] = hh_age_cm
        self.row['age_by_hh'] = age_hh_cm
        self.row['hh_by_hh'] = hh_cm
        self.row['hh_by_hh_scale'] = hh_cm_scale
        self.row.append()
        self.h5file.flush()

    # matrix storage functions # - # - # - # - # - # - # - # - # - # - # - # - # - # - #   

    def get_derived_contact_matrices(self, P, cmatrix, rng):
        """
        Create a "hh size by hh size" matrix of "between household" contact rates.
        Create a "hh by age" and "age by hh" matrices of "between household" contact rates.
        """

        hh_cm = np.zeros((self.max_hh_size, self.max_hh_size), dtype=np.float)
        hh_cm_scale = np.zeros((self.max_hh_size, self.max_hh_size), dtype=np.float)
        hh_age_cm = np.zeros((self.max_hh_size, self.max_age), dtype=np.float)
        age_hh_cm = np.zeros((self.max_age, self.max_hh_size), dtype=np.float)
        hh_size_dist = np.zeros(self.max_hh_size)
        age_dist = np.zeros(self.max_age, dtype=np.float)

        sample = rng.sample(list(P.I.values()), self.sample_size) if self.sample_size > 0 else list(P.I.values())

        for ind_A in sample:
            hh_size_dist[min(P.hh_size(ind_A), self.max_hh_size) - 1] += 1
            age_dist[min(ind_A.age, self.max_age - 1)] += 1
            for ind_B in list(P.I.values()):
                # skip self and household members
                if ind_A.ID == ind_B.ID or ind_B in P.housemates(ind_A):
                    continue
                cur_contact = cmatrix.C[cmatrix.age_map[ind_A.age]][cmatrix.age_map[ind_B.age]]
                hh_cm[min(P.hh_size(ind_A), self.max_hh_size) - 1][
                    min(P.hh_size(ind_B), self.max_hh_size) - 1] += cur_contact
                hh_cm_scale[min(P.hh_size(ind_A), self.max_hh_size) - 1][
                    min(P.hh_size(ind_B), self.max_hh_size) - 1] += cur_contact
                hh_age_cm[min(P.hh_size(ind_B), self.max_hh_size) - 1][ind_A.age] += cur_contact
                age_hh_cm[min(ind_B.age, self.max_age - 1)][min(P.hh_size(ind_A), self.max_hh_size) - 1] += cur_contact

        #print cm.shape
        #print hh_size_dist       
        hh_cm_scale /= hh_size_dist
        hh_cm_scale /= hh_size_dist[None].T

        # normalise hh age matrix by age distribution of samples (to get contact patterns per individual)
        nonzero = age_dist > 0
        hh_age_cm.T[nonzero] /= age_dist[None].T[nonzero]

        #print age_hh_cm
        #print hh_size_dist
        nonzero = hh_size_dist > 0
        age_hh_cm.T[nonzero] /= hh_size_dist[None].T[nonzero]
        #print age_hh_cm

        return hh_cm, hh_cm_scale, hh_age_cm, age_hh_cm

    def get_hh_matrix(self, P, age_scale=False):
        """
        Basically, for each person, we want to know what aged people they 
        encounter in household contexts.  
    
        So, start with a matrix of all zeros, from age 0...100 x 0...100.
    
        Iterate over individuals in population (age = i)
    
        For each person in their household, other than themselves (age = j), 
        add 1 to matrix(i, j)
    
        At end, divide all entries by 2 (as we will have counted each 
        "encounter" twice.
    
        This doesn't precisely give a "rate" as such, because there is no 
        time unit, but it does give a relative level of contact between 
        individuals of different ages.
        """

        cm = np.zeros((self.max_age, self.max_age))

        for ind in list(P.I.values()):
            for cur_hm in P.housemates(ind):
                cm[ind.age][cur_hm.age] += 1

        cm /= 2

        if age_scale:
            age_dist = np.array(P.age_dist(self.max_age, self.max_age, norm=False)[0], dtype=np.float32)
            #print cm
            #print age_dist
            cm /= age_dist
            #print cm
            #print np.max(cm)
            #cm /= age_dist[None].T

        return cm

    def get_extended_family_contact_matrix(self, P, age_scale=False):
        cm = np.zeros((self.max_age, self.max_age))

        #print len(P.I.values())
        efm_counts = []
        for ind in list(P.I.values()):
            efms = P.ancestors(ind, None, [])
            efm_counts.append(len(efms))
            for cur_efm in efms:
                cm[ind.age][cur_efm.age] += 1
                cm[cur_efm.age][ind.age] += 1

        #print len(P.I.values()), np.mean(efm_counts)

        cm /= 2

        if age_scale:
            age_dist = np.array(P.age_dist(self.max_age, self.max_age, norm=False)[0], dtype=np.float32)
            nonzero = age_dist > 0
            #print cm.shape
            #print age_dist[None].T.shape

            #cm /= age_dist
            #cm /= age_dist[None].T
            cm.T[nonzero] /= age_dist[None].T[nonzero]
            #        cm.T[nonzero] /= age_dist[None].T[nonzero][None].T

        return cm

    #    def create_age_map(cutoffs):
    #
    #        age_map = []
    #        prev_cutoff = 0
    #        for index, cur_cutoff in enumerate(cutoffs):
    #            age_map.extend([index] * (cur_cutoff - prev_cutoff))
    #            prev_cutoff = cur_cutoff
    #        return age_map
    #
    #
    #    def get_hh_contact_matrix_cutoff(P, cutoffs, age_scale=False):
    #
    #        age_map = create_age_map(cutoffs)
    #
    #        cm = np.zeros((len(cutoffs), len(cutoffs)))
    #
    #        for ind in P.I.values():
    #            for cur_hm in P.housemates(ind):
    #                cm[age_map[ind.age]][age_map[cur_hm.age]] += 1
    #
    #        cm /= 2
    #
    #        if age_scale:
    #    #        print "XXX"
    #
    #            age_dist = np.array(P.age_dist(
    #                [0]+cutoffs if cutoffs[0]>0 else cutoffs+[101],
    #                cutoffs[-1], norm=False)[0],
    #                    dtype=np.float32)
    #
    #            ### NB: this [0]+ is problematic if cutoffs already begins with 0
    #
    #    #        print len(age_dist), age_dist
    #    #        print cm.size
    #            cm /= age_dist
    #    #        cm /= age_dist[None].T
    #    #        cm /= (age_dist * age_dist[None].T)
    #
    #        return cm

    ## matrix output functions # - # - # - # - # - # - # - # - # - # - # - # - # - # - #

    def output_all(self, p, times):
        for row in self.data:
            self.output_contact_matrix_base(row['hh'],
                                            os.path.join(p['prefix'], 'within_hh_cm_%08d' % row['t']))  # , vmax=0.165)
            self.output_contact_matrix_base(row['ext'], os.path.join(p['prefix'], 'ext_fam_cm_%08d' % row['t']),
                                            vmax=70)
            self.output_contact_matrix_base(row['hh'] + 0.5 * row['ext'],
                                            os.path.join(p['prefix'], 'combo_cm_%08d' % row['t']), vmax=70)

            self.output_contact_matrix_base(row['comm'], os.path.join(p['prefix'], 'comm_cm_%08d' % row['t']),
                                            vmax=0.006)

            self.output_contact_matrix_base(row['hh_by_hh'], os.path.join(p['prefix'], 'between_hh_cm_%08d' % row['t']),
                                            xlabel='Household size', ylabel='Household size')

            cur_hh_scale = row['hh_by_hh_scale']
            for i in range(cur_hh_scale.shape[0]):
                for j in range(cur_hh_scale.shape[1]):
                    cur_hh_scale[i][j] /= ((i + 1) * (j + 1))

            self.output_contact_matrix_base(cur_hh_scale,
                                            os.path.join(p['prefix'], 'between_hh_cm_scaled_%08d' % row['t']),
                                            xlabel='Household size', ylabel='Household size')

            self.output_contact_matrix_base(row['hh_by_age'], os.path.join(p['prefix'], 'hh_by_age_cm_%08d' % row['t']),
                                            xlabel='Age', ylabel='Household size', vmax=3.8)
            self.output_contact_matrix_base(row['age_by_hh'], os.path.join(p['prefix'], 'age_by_hh_cm_%08d' % row['t']),
                                            xlabel='Household size', ylabel='Age', vmax=0.55)

    def output_contact_matrix_base(self, cm_data, ofile, xlabel='Age', ylabel='Age', vmax=-1.0):
        fig = plt.figure()
        ax = fig.add_subplot(111)
        img = self.output_contact_matrix_base_ax(cm_data, ax, xlabel, ylabel, vmax)
        fig.colorbar(img)
        fig.savefig(ofile + '.png')

    def output_contact_matrix_base_ax(self, cm_data, ax, xlabel='Age', ylabel='Age of contact', vmax=-1.0, plot_title=None):
        #fig = plt.figure()
        #ax = fig.add_subplot(111)
        nan_locs = np.isnan(cm_data)
        cm_data[nan_locs] = 0
        if vmax > 0:
            img = ax.pcolor(cm_data, vmin=0.0, vmax=vmax, cmap='jet')
        else:
            img = ax.pcolor(cm_data, vmin=0.0, cmap='jet')
        ax.set_xlim(xmax=cm_data.shape[1])
        ax.set_ylim(ymax=cm_data.shape[0])
        #ax.set_aspect(1)
        if xlabel == 'Household size':
            ax.set_xticks(np.arange(cm_data.shape[1]) + 0.5, minor=False)
            ax.set_xticklabels(np.arange(cm_data.shape[1]) + 1)
        if ylabel == 'Household size':
            ax.set_yticks(np.arange(cm_data.shape[0]) + 0.5, minor=False)
            ax.set_yticklabels(np.arange(cm_data.shape[0]) + 1)
        ax.tick_params(length=0)
        if plot_title:
            ax.set_title(plot_title, fontsize=16)
        ax.set_xlabel(xlabel, fontsize=16)
        ax.set_ylabel(ylabel, fontsize=16)

        plt.setp(ax.get_xticklabels(), fontsize=16)
        plt.setp(ax.get_yticklabels(), fontsize=16)

        return img
        #fig.colorbar(img)
        #fig.savefig(ofile+'.png')
        #np.savetxt(ofile+'.csv', cm_data, fmt='%g', delimiter=',')

    def output_contact_matrix(self, cm_data, ofile, size=10, x_label=None, y_label=None):
        data = cm_data  # np.reshape(cm_data, cm_data.size, order='F').reshape((cm_data.shape))

        fig = plt.figure()
        ax = fig.add_subplot(111)
        size_x = len(cm_data[0])
        size_y = len(cm_data)

        cmap = cm.jet
        X, Y = np.meshgrid(np.array(list(range(size_x + 1))), np.array(list(range(size_y + 1))))
        img = ax.pcolormesh(X, Y, data, cmap=cmap)
        if size_x == size_y:
            ax.set_aspect(1)
        if x_label:
            ax.set_xlabel(x_label)
        if y_label:
            ax.set_ylabel(y_label)
        fig.colorbar(img)
        fig.savefig(ofile)

    ## Accepting that this is a mess for now; to clean

    def output_comm_contact_matrix(self, cmatrix, ofile, max_age=100,
                                   vmin=None, vmax=None,
                                   age_scale=False, log_scale=False, phi=None, colbar=True):
        #print "&&&&&&&&&&&&&&&&&"
        cm_data = np.reshape(cmatrix.C, cmatrix.C.size, order='F').reshape((cmatrix.C.shape))

        fig = plt.figure()
        ax = fig.add_subplot(111)
        if cmatrix.epsilon > -1:
            cm_data[abs(cm_data) < 1e-6] = np.NaN
        else:
            cm_data[abs(cm_data) < 1e-6] = 1e-4
        cm_data_m = np.ma.array(cm_data, mask=np.isnan(cm_data))

        if age_scale:
            cm_data_m = np.array(cm_data_m) * np.array([ \
                cmatrix.get_class_size(x) \
                for x in cmatrix.age_classes])
        X, Y = np.meshgrid(np.array(cmatrix.age_classes + [max_age]),
                           np.array(cmatrix.age_classes + [max_age]))
        cmap = cm.jet
        cmap.set_bad('k', 1.)
        if vmax:
            if log_scale:
                img = ax.pcolormesh(X, Y, cm_data_m,
                                    vmin=vmin, vmax=vmax,
                                    norm=LogNorm(vmin=vmin, vmax=vmax), cmap=cmap)
            else:
                img = ax.pcolormesh(X, Y, cm_data_m,
                                    vmin=vmin, vmax=vmax, cmap=cmap)
        else:
            if log_scale:
                img = ax.pcolormesh(X, Y, cm_data_m,
                                    vmin=vmin, norm=LogNorm(), cmap=cmap)
            else:
                img = ax.pcolormesh(X, Y, cm_data_m,
                                    vmin=vmin, cmap=cmap)

        ax.set_aspect(1)
        ax.set_xlabel('Individual Age')
        ax.set_ylabel('Contact Age')
        if colbar == True:
            fig.colorbar(img)  #, ticks=[
            #        1e-4, 2e-4, 3e-4, 4e-4, 5e-4, 6e-4, 7e-4, 8e-4, 9e-4,
            #        1e-3, 2e-3, 3e-3, 4e-3, 5e-3, 6e-3, 7e-3, 8e-3, 9e-3,
            #        1e-2, 2e-2, 3e-2, 4e-2])
            #    fig.savefig(ofile, dpi=600)
        fig.savefig(ofile)

    #    np.savetxt(ofile+'.csv', cm_data_m, fmt='%g', delimiter=',')

    #from matplotlib.colors import LogNorm

    #def output_extended_family_contact_matrix(self, P, ofile, age_scale=False):
    #    cm_data = get_extended_family_contact_matrix(P, age_scale)
    #    output_contact_matrix_base(cm_data, ofile, age_scale)

    def output_combined_family_contact_matrix(self, P, ofile, age_scale=False):
        hh_data = get_hh_contact_matrix(P, age_scale)
        ext_data = get_extended_family_contact_matrix(P, age_scale)
        for ext_wt in [0.0, 0.2, 0.4, 0.6, 0.8, 1.0]:
            combo_data = hh_data + ext_wt * ext_data
            output_contact_matrix_base(combo_data, ofile + '_%.1f' % ext_wt)


#    def output_hh_contact_matrix(self, P, ofile, cutoffs=None, age_scale=False, vmax=None):
#        fig = plt.figure()
#        ax = fig.add_subplot(111)
#        if cutoffs:
#            cm_data = get_hh_contact_matrix_cutoff(P, cutoffs, age_scale)
#            print cutoffs
#            X, Y = np.meshgrid(np.array([0]+cutoffs), np.array([0]+cutoffs))
#            print X
#            print Y
#            print cm_data
#            if vmax:
#                img = ax.pcolor(X, Y, cm_data, vmin=0.0, vmax=vmax,cmap='jet')
#            else:
#                img = ax.pcolor(X, Y, cm_data, vmin=0.0, cmap='jet')
#        else:
#            cm_data = get_hh_contact_matrix(P, age_scale)
#            img = ax.pcolor(cm_data[:100,:100], vmin=0.0, cmap='jet')    
#    
#    
#        ax.set_xlim(xmax=100)
#        ax.set_ylim(ymax=100)
#        ax.set_aspect(1)
#        ax.set_xlabel('Age')
#        ax.set_ylabel('Age')
#        fig.colorbar(img)
#        fig.savefig(ofile+'.png')
#        np.savetxt(ofile+'.csv', cm_data, fmt='%g', delimiter=',')
#
