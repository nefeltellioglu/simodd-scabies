from math import pi, exp, sqrt

import numpy as np


mossong_age_classes = list(range(0, 80, 5))

#based on lib contact matrix
mossong_a_levels=[12.4470138097, 23.635755183, 22.6722320486, 25.128843179660006, 15.609865529899999, 12.585700007789999, 10.573382443429999, 10.121678666649998, 9.34681309277, 7.24893763342, 7.691463691499999, 5.6530755933300005, 3.914026770499999, 0.78918120543, 0.6068943621999999, 0.35479825455999997]

#givencmatrix
import pandas as pd
#liberia_matrix=pd.read_csv('/Users/ntellioglu/repos/simodd-dis-development/example/data/Liberia_contact_matrix.csv',header=None).to_numpy()

"""
# by-age activity levels from POLYMOD
# (these values are modified by c2000 household size dist)

mossong_a_levels = [7.346, 11.7465074, 15.0756694, 14.8946084,
                    11.66090728, 11.71655332, 12.036, 11.69172403,
                    11.478, 11.83833327, 10.4181811, 11.02727324,
                    8.0980591, 8.2170918, 6.08298147, 6.32693737]

# by-age activity levels, modified to include a separate 0-1y age group
mossong_age_classes_lo = [0, 1] + list(range(5, 80, 5))
mossong_a_levels_lo = [2.346, 7.346, 11.7465074, 15.0756694, 14.8946084,
                       11.66090728, 11.71655332, 12.036, 11.69172403,
                       11.478, 11.83833327, 10.4181811, 11.02727324,
                       8.0980591, 8.2170918, 6.08298147, 6.32693737]

# (these are raw TOTAL activity levels by age)
mossong_a_levels_unmod = [10.21,  14.81,  18.22,  17.58,
                          13.57, 13.57, 14.14, 14.14,
                          13.83, 13.83, 12.3, 12.3,
                          9.21, 9.21, 6.89, 6.89]

#givencmatrix
#import pandas as pd
#cmatrix=pd.read_csv('/Users/ntellioglu/repos/simodd-dis-development/example/data/Liberia_contact_matrix.csv',header=None).to_numpy()
"""



class KnownContactMatrix(object):
    def __init__(self, givenmatrix=None,age_classes=mossong_age_classes):
        self.C = givenmatrix
        self.EC = None
        self.max_age = 100
        self.age_classes = age_classes
        self.age_class_ids = {}
        self._create_age_maps(age_classes)

    
    def _create_age_maps(self, cutoffs):
        """
        Create maps from :
            (a) age (years) to age_classes (index); and
            (b) age (years) to age_class_ids (left edge of age class)
        """

        # an age map that maps ages to age class indices
        self.age_map = []
        prev_cutoff = 0
        for index, cur_cutoff in enumerate(cutoffs[1:] + [self.max_age + 2]):
            self.age_map.extend([index] * (cur_cutoff - prev_cutoff))
            prev_cutoff = cur_cutoff

        # an age map that maps ages to left edge of matching age class
        self.age_ids_map = []
        prev_cutoff = 0
        for index, cur_cutoff in enumerate(cutoffs[1:] + [self.max_age + 2]):
            self.age_ids_map.extend([prev_cutoff] * (cur_cutoff - prev_cutoff))
            prev_cutoff = cur_cutoff
            
    def init_a_levels(self, P):
        pass

    def init_age_classes(self, P):
        for age, next_age in zip(mossong_age_classes, mossong_age_classes[1:] + [self.max_age + 1]):
            self.age_class_ids[age] = P.ind_ids_by_age(age, next_age - 1)

    def init_contact_matrix_gaussian(self, epsilon, sigma_2):
        #self.C = givenmatrix
        self._expand_contact_matrix()
    
    def init_contact_matrix(self, epsilon):
        """
        Initialise the contact matrix.
        
        Epsilon in the range [0,1] is a convex combination parameter:
        epsilon = 1 ---> proportionate mixing
        epsilon = 0 ---> preferred mixing (i.e., on the diagonal only)
        """
        #self.C = givenmatrix
        #init_mat=self.C
        #self.C = np.zeros((self.max_age + 1, self.max_age + 1,))
        #for i in range(self.max_age + 1):
        #    ci = self.age_map[i]
        #    for j in range(len(self.age_classes)):
        #        cj = self.age_map[j]
        #        self.C[i][j] = init_mat[ci][cj]/5
        self._expand_contact_matrix()
        #self.EC=self.C

    def _expand_contact_matrix(self):
        """
        Create an expanded version of the contact matrix, with a row for each individual
        age.  This speeds up calculation of FOI slightly, with only a very small cost in
        storage redundancy.
        """
        
        self.EC = np.zeros((self.max_age + 1, len(self.age_classes)))
        for i in range(self.max_age + 1):
            ci = self.age_map[i]
            for j in range(len(self.age_classes)):
                #cj = self.age_map[j]
                self.EC[i][j] = self.C[ci][j]

    def update_age_classes(self, births, deaths, immigrants, birthdays):
        pass


class ContactMatrixConstant(object):
    def __init__(self):
        self.C = None
        self.EC = None
        self.max_age = 100
        self.age_classes = [0]
        self.age_class_ids = {}
        self._create_age_maps()

    def _create_age_maps(self):
        """
        Create maps from :
            (a) age (years) to age_classes (index); and
            (b) age (years) to age_class_ids (left edge of age class)
        """
        # an age map that maps ages to age class indices
        self.age_map = [0] * (self.max_age + 1)

        # an age map that maps ages to left edge of matching age class
        self.age_ids_map = [0] * (self.max_age + 1)

    def init_a_levels(self, P):
        pass

    def init_age_classes(self, P):
        for age, next_age in zip(mossong_age_classes, mossong_age_classes[1:] + [self.max_age + 1]):
            self.age_class_ids[age] = P.ind_ids_by_age(age, next_age - 1)

    def init_contact_matrix_gaussian(self, epsilon, sigma_2):

        self.C = np.zeros((1, 1))

        self.C[0][0] = 1.0/10000

        self._expand_contact_matrix()

    def _expand_contact_matrix(self):
        """
        Create an expanded version of the contact matrix, with a row for each individual
        age.  This speeds up calculation of FOI slightly, with only a very small cost in
        storage redundancy.
        """

        self.EC = np.zeros((self.max_age + 1, len(self.age_classes)))
        for i in range(self.max_age + 1):
            ci = self.age_map[i]
            for j in range(len(self.age_classes)):
                #cj = self.age_map[j]
                self.EC[i][j] = self.C[ci][j]

    def update_age_classes(self, births, deaths, immigrants, birthdays):
        pass


class ContactMatrixFlat(object):
    def __init__(self):
        self.C = None
        self.EC = None
        self.max_age = 100
        self.age_classes = [0]
        self.age_class_ids = {}
        self._create_age_maps()

    def _create_age_maps(self):
        """
        Create maps from :
            (a) age (years) to age_classes (index); and
            (b) age (years) to age_class_ids (left edge of age class)
        """
        # an age map that maps ages to age class indices
        self.age_map = [0] * (self.max_age + 1)

        # an age map that maps ages to left edge of matching age class
        self.age_ids_map = [0] * (self.max_age + 1)

    def init_a_levels(self, P):
        pass

    def init_age_classes(self, P):
        for age, next_age in zip(mossong_age_classes, mossong_age_classes[1:] + [self.max_age + 1]):
            self.age_class_ids[age] = P.ind_ids_by_age(age, next_age - 1)


    def init_contact_matrix_gaussian(self, epsilon, sigma_2):

        # calculate total number of contacts
        D = sum([aj * len(self.age_class_ids[j]) for aj, j in zip(mossong_a_levels, mossong_age_classes)])

        # average across population size
        a_level = D / sum([len(self.age_class_ids[j]) for j in mossong_age_classes])

        self.C = np.zeros((1, 1))

        self.C[0][0] = (a_level * a_level) / D

        self._expand_contact_matrix()

    def _expand_contact_matrix(self):
        """
        Create an expanded version of the contact matrix, with a row for each individual
        age.  This speeds up calculation of FOI slightly, with only a very small cost in
        storage redundancy.
        """

        self.EC = np.zeros((self.max_age + 1, len(self.age_classes)))
        for i in range(self.max_age + 1):
            ci = self.age_map[i]
            for j in range(len(self.age_classes)):
                #cj = self.age_map[j]
                self.EC[i][j] = self.C[ci][j]

    def update_age_classes(self, births, deaths, immigrants, birthdays):
        pass


class ContactMatrix(object):
    def __init__(self, age_classes=mossong_age_classes,
                 a_levels=mossong_a_levels, max_age=100, smooth=True):

        assert len(age_classes) == len(a_levels)

        # the contact matrix and expanded contact matrix
        self.C = None
        self.EC = None

        # the maximum age in the contact matrix 
        self.max_age = max_age

        # age_classes stores a list of the LEFT edge of each age class
        self.age_classes = age_classes

        # prev_age_class maps each age to the preceding age class
        # (used for updating as people age)
        self.prev_age_class = {}

        # age_class_ids stores a current list of IDs keyed by age class
        self.age_class_ids = {}

        # activity levels (avg # contacts per day) for each age class
        self.base_a_levels = np.array(a_levels)  # total daily contacts
        self.comm_a_levels = self.base_a_levels  # daily community contacts adjusting for expected household size

        # create a map from age to age_class
        self.age_map = []
        self.age_ids_map = []
        self._create_age_maps(age_classes)
        
        if smooth:
            self._smooth_matrix()

    def _smooth_matrix(self):
        new_age_classes = list(range(max(self.age_classes)))
        polyf = np.poly1d(np.polyfit(self.age_classes, self.base_a_levels, 9))
        self.base_a_levels = [polyf(x) for x in new_age_classes]
        self.age_classes = new_age_classes

        #self.age_classes.extend([x for x in range(new_age_classes[-1], self.max_age)])
        #final_a = self.base_a_levels[-1]
        #len_diff = len(self.age_classes)-len(self.base_a_levels)
        #self.base_a_levels.extend([final_a]*len_diff)

        self.comm_a_levels = self.base_a_levels
        self._create_age_maps(self.age_classes)
        
    def _create_age_maps(self, cutoffs):
        """
        Create maps from :
            (a) age (years) to age_classes (index); and
            (b) age (years) to age_class_ids (left edge of age class)
        """

        # an age map that maps ages to age class indices
        self.age_map = []
        prev_cutoff = 0
        for index, cur_cutoff in enumerate(cutoffs[1:] + [self.max_age + 2]):
            self.age_map.extend([index] * (cur_cutoff - prev_cutoff))
            prev_cutoff = cur_cutoff

        # an age map that maps ages to left edge of matching age class
        self.age_ids_map = []
        prev_cutoff = 0
        for index, cur_cutoff in enumerate(cutoffs[1:] + [self.max_age + 2]):
            self.age_ids_map.extend([prev_cutoff] * (cur_cutoff - prev_cutoff))
            prev_cutoff = cur_cutoff

    def init_age_classes(self, P):
        """Set up the age_class_ids and prev_age_class structures."""

        for age, next_age in zip(self.age_classes,
                                 self.age_classes[1:] + [self.max_age + 1]):
            self.age_class_ids[age] = P.ind_ids_by_age(age, next_age - 1)
            self.prev_age_class[next_age] = age

    def init_a_levels(self, P):
        """
        adjust activity levels used to specify community mixing matrices according to
        current expected number of household contacts

        NB: If using this, base_a_levels should be the TOTAL number of contacts per day
        """
        #mean_housemate_counts = P.get_mean_housemates_by_age(cutoffs=mossong_age_classes)
        mean_housemate_counts = P.get_mean_housemates_by_age(cutoffs=list(range(75)))

        self.comm_a_levels = self.base_a_levels - np.array(mean_housemate_counts)


    def init_contact_matrix_gaussian(self, epsilon, sigma_2):
        """
        Initialise the contact matrix with a Gaussian kernel.

        epsilon is the degree of age assortativity in the contact matrix
        sigma_2 is the variance of the Gaussian distribution.
        exp_hh_contacts is the expected number of household contacts to account for
            (ie, subtract from activity levels)
        """

        # calculate total number of contacts
        D = sum([aj * len(self.age_class_ids[j]) for aj, j in
                 zip(self.comm_a_levels, self.age_classes)])

        # create the Gaussian kernel
        kernel = np.zeros((self.max_age + 1, self.max_age + 1))

        sigma_2 /= 1.0

        for i in range(self.max_age + 1):
            for j in range(self.max_age + 1):
                kernel[i][j] = 1 / sqrt(2 * pi * sigma_2) * exp(-(i - j) ** 2 / (2.0 * sigma_2))

        # normalise rows (to sum to 1.0)
        row_sums = np.sum(kernel, axis=1)
        kernel /= row_sums[None].T

        #        np.savetxt('/home/ngeard/Desktop/kernel.csv', kernel, fmt='%g', delimiter=',')

        kernel_ref = np.zeros((len(self.age_classes), len(self.age_classes)))

        # aggregate kernel cells to reduce to dimensions: len(age_classes) x len(age_classes)
        for i, (li, ri) in enumerate(zip(self.age_classes, self.age_classes[1:] + [self.max_age])):
            for j, (lj, rj) in enumerate(zip(self.age_classes, self.age_classes[1:] + [self.max_age])):
                kernel_ref[i][j] = np.mean(kernel[li:ri, lj:rj]) * (rj - lj)

                #        np.savetxt('/home/ngeard/Desktop/kernel_ref.csv', kernel_ref, fmt='%g', delimiter=',')

        self.C = np.zeros((len(self.age_classes), len(self.age_classes)))

        # make sure we are working with floats...
        epsilon = float(epsilon)

        for i, ai in enumerate(self.comm_a_levels):
            for j, aj in enumerate(self.comm_a_levels):
                pref = float(ai) / len(self.age_class_ids[self.age_classes[j]]) if len(
                    self.age_class_ids[self.age_classes[j]]) > 0 else 0
                self.C[i][j] = epsilon * ai * aj / D + (1.0 - epsilon) * pref * kernel_ref[i][j]

        self._expand_contact_matrix()
        #self.EC = self.C



    def init_contact_matrix(self, epsilon):
        """
        Initialise the contact matrix.
        
        Epsilon in the range [0,1] is a convex combination parameter:
        epsilon = 1 ---> proportionate mixing
        epsilon = 0 ---> preferred mixing (i.e., on the diagonal only)
        """

        # calculate total number of contacts
        D = sum([aj * len(self.age_class_ids[j]) for aj, j in
                 zip(self.comm_a_levels, self.age_classes)])

        self.C = np.zeros((len(self.age_classes), len(self.age_classes)))

        # make sure we are working with floats...
        epsilon = float(epsilon)

        for i, ai in enumerate(self.comm_a_levels):
            for j, aj in enumerate(self.comm_a_levels):
                pref = float(ai) / len(self.age_class_ids[self.age_classes[j]]) \
                    if i == j and self.age_class_ids[self.age_classes[j]] else 0
                self.C[i][j] = \
                    epsilon * ai * aj / D + (1.0 - epsilon) * pref

        self._expand_contact_matrix()

    def _expand_contact_matrix(self):
        """
        Create an expanded version of the contact matrix, with a row for each individual
        age.  This speeds up calculation of FOI slightly, with only a very small cost in
        storage redundancy.
        """

        self.EC = np.zeros((self.max_age + 1, len(self.age_classes)))
        for i in range(self.max_age + 1):
            ci = self.age_map[i]
            for j in range(len(self.age_classes)):
                #cj = self.age_map[j]
                self.EC[i][j] = self.C[ci][j]

    def get_class_size(self, age_class_index):
        """Return the number of people in the specified age class."""

        return len(self.age_class_ids[age_class_index])

    def update_age_classes(self, births, deaths, immigrants, birthdays):
        """
        Update age_class_ids on the basis of demographic events.
        
        NOTE: this won't actually modify the contact matrix; currently,
        init_contact_matrix would have to be called again for this to occur.
        """

        # birthdays - move people from one age class to next
        for cur_age_class in self.age_classes:
            for cur_ind in birthdays[cur_age_class - 1]:
                self.age_class_ids[
                    self.prev_age_class[cur_age_class]].remove(cur_ind)
                self.age_class_ids[cur_age_class].append(cur_ind)

        # births - add people to first age class
        self.age_class_ids[0].extend([x.ID for x in births])

        # deaths - remove people from appropriate age class
        for dth in deaths:
            self.age_class_ids[self.age_ids_map[dth.age]].remove(dth.ID)

        # immigrants - add people to appropriate age class
        for imm in immigrants:
            self.age_class_ids[self.age_ids_map[imm.age]].append(imm.ID)

#    def get_contact_rate(self, ind_age, contact_age_class):
#        """
#        Return the contact rate between an individual and contact of
#        specified ages.
#
#        NOTE: this will probably need to be updated.
#        """
#
#        return self.EC[ind_age][contact_age_class]

#    def set_contact_rates(self, new_matrix):
#        self.C = new_matrix
#        self.expand_contact_matrix()

    def get_normalised_matrix(self):
        """
        Return a copy of the (expanded) matrix with row sums normalised to 1.0
        """
        n_matrix = np.copy(self.EC)
        row_sums = np.sum(self.EC, axis=1)
        n_matrix /= row_sums[None].T

        return n_matrix

    def export_contact_matrix(self, ofile):
        np.savetxt(ofile, self.C, fmt='%g', delimiter=',')

      
class ContactMatrix2(object):
    def __init__(self, age_classes=mossong_age_classes,
                 a_levels=mossong_a_levels, 
                 given_matrix=None,
                 max_age=100, smooth=True):

        assert len(age_classes) == len(a_levels)

        # the contact matrix and expanded contact matrix
        self.C = None
        self.EC = None
        self.given_matrix=given_matrix
        # the maximum age in the contact matrix 
        self.max_age = max_age

        # age_classes stores a list of the LEFT edge of each age class
        self.age_classes = age_classes

        # prev_age_class maps each age to the preceding age class
        # (used for updating as people age)
        self.prev_age_class = {}

        # age_class_ids stores a current list of IDs keyed by age class
        self.age_class_ids = {}

        # activity levels (avg # contacts per day) for each age class
        self.base_a_levels = np.array(a_levels)  # total daily contacts
        self.comm_a_levels = self.base_a_levels  # daily community contacts adjusting for expected household size

        # create a map from age to age_class
        self.age_map = []
        self.age_ids_map = []
        self._create_age_maps(age_classes)

        if smooth:
            self._smooth_matrix()

    def _smooth_matrix(self):
        new_age_classes = list(range(max(self.age_classes)))
        inter_matrix= np.zeros((len(self.age_classes), len(new_age_classes)))
        new_matrix= np.zeros((len(new_age_classes), len(new_age_classes)))
        polyf = np.poly1d(np.polyfit(self.age_classes, self.base_a_levels, 9))
        self.base_a_levels = [polyf(x) for x in new_age_classes]
        
        #self.age_classes.extend([x for x in range(new_age_classes[-1], self.max_age)])
        #final_a = self.base_a_levels[-1]
        #len_diff = len(self.age_classes)-len(self.base_a_levels)
        #self.base_a_levels.extend([final_a]*len_diff)
        for i in range(len(self.given_matrix)):
            polyf = np.poly1d(np.polyfit(self.age_classes, self.given_matrix[i], 1))
            inter_matrix[i] = [polyf(x) for x in new_age_classes]
            #display(self.given_matrix[i])
            #display(inter_matrix[i])
        #display(inter_matrix)
        for j in range(len(new_age_classes)):
            polyf = np.poly1d(np.polyfit(self.age_classes, [el[j] for el in inter_matrix], 1))
            new_column=[polyf(x) for x in new_age_classes]
            new_matrix[:, j] =new_column
            #display(new_matrix[:, j])
        self.given_matrix=  new_matrix 
        #display(sum(self.given_matrix<0))
        #display(self.given_matrix)
        
        
        self.age_classes = new_age_classes
        self.comm_a_levels = self.base_a_levels
        self._create_age_maps(self.age_classes)

    def _create_age_maps(self, cutoffs):
        """
        Create maps from :
            (a) age (years) to age_classes (index); and
            (b) age (years) to age_class_ids (left edge of age class)
        """

        # an age map that maps ages to age class indices
        self.age_map = []
        prev_cutoff = 0
        for index, cur_cutoff in enumerate(cutoffs[1:] + [self.max_age + 2]):
            self.age_map.extend([index] * (cur_cutoff - prev_cutoff))
            prev_cutoff = cur_cutoff

        # an age map that maps ages to left edge of matching age class
        self.age_ids_map = []
        prev_cutoff = 0
        for index, cur_cutoff in enumerate(cutoffs[1:] + [self.max_age + 2]):
            self.age_ids_map.extend([prev_cutoff] * (cur_cutoff - prev_cutoff))
            prev_cutoff = cur_cutoff

    def init_age_classes(self, P):
        """Set up the age_class_ids and prev_age_class structures."""

        for age, next_age in zip(self.age_classes,
                                 self.age_classes[1:] + [self.max_age + 1]):
            self.age_class_ids[age] = P.ind_ids_by_age(age, next_age - 1)
            self.prev_age_class[next_age] = age

    def init_a_levels(self, P):
        """
        adjust activity levels used to specify community mixing matrices according to
        current expected number of household contacts

        NB: If using this, base_a_levels should be the TOTAL number of contacts per day
        """
        #mean_housemate_counts = P.get_mean_housemates_by_age(cutoffs=mossong_age_classes)
        mean_housemate_counts = P.get_mean_housemates_by_age(cutoffs=list(range(75)))

        self.comm_a_levels = self.base_a_levels - np.array(mean_housemate_counts)


    def init_contact_matrix_gaussian(self, epsilon, sigma_2):
        """
        Initialise the contact matrix with a Gaussian kernel.

        epsilon is the degree of age assortativity in the contact matrix
        sigma_2 is the variance of the Gaussian distribution.
        exp_hh_contacts is the expected number of household contacts to account for
            (ie, subtract from activity levels)
        """

        # calculate total number of contacts
        D = sum([aj * len(self.age_class_ids[j]) for aj, j in
                 zip(self.comm_a_levels, self.age_classes)])

        # create the Gaussian kernel
        kernel = np.zeros((self.max_age + 1, self.max_age + 1))

        sigma_2 /= 1.0

        for i in range(self.max_age + 1):
            for j in range(self.max_age + 1):
                kernel[i][j] = 1 / sqrt(2 * pi * sigma_2) * exp(-(i - j) ** 2 / (2.0 * sigma_2))

        # normalise rows (to sum to 1.0)
        row_sums = np.sum(kernel, axis=1)
        kernel /= row_sums[None].T

        #        np.savetxt('/home/ngeard/Desktop/kernel.csv', kernel, fmt='%g', delimiter=',')

        kernel_ref = np.zeros((len(self.age_classes), len(self.age_classes)))

        # aggregate kernel cells to reduce to dimensions: len(age_classes) x len(age_classes)
        for i, (li, ri) in enumerate(zip(self.age_classes, self.age_classes[1:] + [self.max_age])):
            for j, (lj, rj) in enumerate(zip(self.age_classes, self.age_classes[1:] + [self.max_age])):
                kernel_ref[i][j] = np.mean(kernel[li:ri, lj:rj]) * (rj - lj)

                #        np.savetxt('/home/ngeard/Desktop/kernel_ref.csv', kernel_ref, fmt='%g', delimiter=',')

        self.C = np.zeros((len(self.age_classes), len(self.age_classes)))

        # make sure we are working with floats...
        epsilon = float(epsilon)

        for i, ai in enumerate(self.comm_a_levels):
            for j, aj in enumerate(self.comm_a_levels):
                pref = float(ai) / len(self.age_class_ids[self.age_classes[j]]) if len(
                    self.age_class_ids[self.age_classes[j]]) > 0 else 0
                self.C[i][j] = epsilon * ai * aj / D + (1.0 - epsilon) * pref * kernel_ref[i][j]

        self._expand_contact_matrix()
        #self.EC = self.C



    def init_contact_matrix(self, epsilon):
        """
        Initialise the contact matrix.
        
        Epsilon in the range [0,1] is a convex combination parameter:
        epsilon = 1 ---> proportionate mixing
        epsilon = 0 ---> preferred mixing (i.e., on the diagonal only)
        """

        # calculate total number of contacts
        D = sum([aj * len(self.age_class_ids[j]) for aj, j in
                 zip(self.comm_a_levels, self.age_classes)])

        self.C = np.zeros((len(self.age_classes), len(self.age_classes)))

        # make sure we are working with floats...
        epsilon = float(epsilon)
        #display(self.comm_a_levels)
        for i, ai in enumerate(self.comm_a_levels):
            for j, aj in enumerate(self.comm_a_levels):
                #pref = float(ai) / len(self.age_class_ids[self.age_classes[j]]) \
                #    if i == j and self.age_class_ids[self.age_classes[j]] else 0
                
                #display("i_%s     j_%s"%(i,j))
                #display(self.given_matrix[i][j])
                self.C[i][j] = self.given_matrix[i][j] / len(self.age_class_ids[self.age_classes[j]])#D 

        self._expand_contact_matrix()

    def _expand_contact_matrix(self):
        """
        Create an expanded version of the contact matrix, with a row for each individual
        age.  This speeds up calculation of FOI slightly, with only a very small cost in
        storage redundancy.
        """

        self.EC = np.zeros((self.max_age + 1, len(self.age_classes)))
        for i in range(self.max_age + 1):
            ci = self.age_map[i]
            for j in range(len(self.age_classes)):
                #cj = self.age_map[j]
                self.EC[i][j] = self.C[ci][j]

    def get_class_size(self, age_class_index):
        """Return the number of people in the specified age class."""

        return len(self.age_class_ids[age_class_index])

    def update_age_classes(self, births, deaths, immigrants, birthdays):
        """
        Update age_class_ids on the basis of demographic events.
        
        NOTE: this won't actually modify the contact matrix; currently,
        init_contact_matrix would have to be called again for this to occur.
        """

        # birthdays - move people from one age class to next
        for cur_age_class in self.age_classes:
            for cur_ind in birthdays[cur_age_class - 1]:
                self.age_class_ids[
                    self.prev_age_class[cur_age_class]].remove(cur_ind)
                self.age_class_ids[cur_age_class].append(cur_ind)

        # births - add people to first age class
        self.age_class_ids[0].extend([x.ID for x in births])

        # deaths - remove people from appropriate age class
        for dth in deaths:
            self.age_class_ids[self.age_ids_map[dth.age]].remove(dth.ID)

        # immigrants - add people to appropriate age class
        for imm in immigrants:
            self.age_class_ids[self.age_ids_map[imm.age]].append(imm.ID)

#    def get_contact_rate(self, ind_age, contact_age_class):
#        """
#        Return the contact rate between an individual and contact of
#        specified ages.
#
#        NOTE: this will probably need to be updated.
#        """
#
#        return self.EC[ind_age][contact_age_class]

#    def set_contact_rates(self, new_matrix):
#        self.C = new_matrix
#        self.expand_contact_matrix()

    def get_normalised_matrix(self):
        """
        Return a copy of the (expanded) matrix with row sums normalised to 1.0
        """
        n_matrix = np.copy(self.EC)
        row_sums = np.sum(self.EC, axis=1)
        n_matrix /= row_sums[None].T

        return n_matrix

    def export_contact_matrix(self, ofile):
        np.savetxt(ofile, self.C, fmt='%g', delimiter=',')
