def cov_gen_interrupt(time_points, cov_levels):
    """
    Helper function for generating heterogeneous patterns of disease coverage.
    Given, eg, time_points=[2,5,3], cov_levels=[0.9, 0.6, 0.8], the coverage pattern:
        [0.9, 0.9, 0.6, 0.6, 0.6, 0.6, 0.6, 0.8, 0.8, 0.8]
    will be generated (ie 0.9 for the first 2 years, followed by 0.6 for five years,
    then 0.8 for three years, and beyond)
    """
    cov = []
    for tp, lev in zip(time_points, cov_levels):
        cov.extend([lev] * tp)
    return cov


class VaccineBase(object):
    """
    Base vaccine class
    """

    def __init__(self, label, min_age,max_age, age_days, coverage, t_per_year, offset, noncompliance_percentage,
                 precond, fail_prob=0.0, v_state='S',vaccination_per_year=1,immediately_effective=True):

        self.label = label  # vaccine id
        self.immediately_effective=immediately_effective #May42021
        self.vaccination_per_year=vaccination_per_year  #May42021
        self.min_age = min_age  # minimum age when vaccine is given (year component)
        self.max_age = max_age  # maximum age when vaccine is given (year component)
        self.age_days = age_days  # age when vaccine is given (day component)
        self.t_per_year = t_per_year  # store locally for convenience
        self.v_factor = 364 // self.t_per_year
        self.coverage = coverage  # the coverage for each vaccination: % of population in given age interval got vaccinated
        self.v_state = v_state  # the disease state that being vaccinated moves you into
        self.offset = offset  # offset due to burn in durations:give vaccine after year 10 for ex.
        self.precond = precond  # vaccine preconditions (ie, primary course needed before booster)
        self.fail_prob = fail_prob  # probability of vaccine failure
        self.noncompliance_percentage = noncompliance_percentage
        #elf.rng = rng
        self.count = 0  # number of times this vaccine has been provided
        self.count_ever = 0
        
    def get_t_years(self, t):
        """
        work out where we are in the coverage schedule
        if t_years < 0, vaccination hasn't started yet
        if we are out of the scheduled vaccination period, vaccines are not given anymore.
        """
#        print "VACCINE:", self.label, "t_years:", (t / self.t_per_year - self.offset), self.coverage[(t / self.t_per_year - self.offset)]
        #step_in_schedule=(t * self.vaccination_per_year) // (self.t_per_year ) - (self.offset * self.vaccination_per_year)
        step_in_schedule = ((t - (self.offset * self.t_per_year)) * self.vaccination_per_year) / (self.t_per_year ) 
        if isinstance(step_in_schedule, int) or step_in_schedule.is_integer():
            if step_in_schedule > (len(self.coverage) - 1):
                step_in_schedule=-1
        #return min((t * self.vaccination_per_year) // (self.t_per_year ) - (self.offset * self.vaccination_per_year), len(self.coverage) - 1)
        return int(step_in_schedule)
    def get_candidates(self, t, P):
        """ Individuals in the given age interval got vaccinated. """
#        print "CAND:", self.label, ":", self.age, self.age_days, self.v_factor, len(P.individuals_by_age(self.age))
        #return [x for x in P.individuals_by_age(self.age) if
        #        x.age_days // self.v_factor == self.age_days // self.v_factor]
        #return [x for x in P.individuals_by_age(self.age, 150) if
        #        x.age_days // self.v_factor == self.age_days // self.v_factor]
        return [x for x in P.individuals_by_age(self.min_age, self.max_age)]
    def check_vaccination_period(self,t):
        """ Checks if there is a vaccination scheduled at a given time point"""
        #vaccine_time_counter_in_a_year= (t * self.vaccination_per_year) % self.t_per_year
        vaccine_time_counter_in_a_year= (t - (self.offset * self.t_per_year)) % self.t_per_year
        
        if self.get_t_years(t) >= 0 and vaccine_time_counter_in_a_year == 0: 
            #if we are in the vaccination period and this is the time where vaccinations are applied
            return True
        else:
            return False
           
    def vaccinate_pop(self, t, P, states, rng):
        #check if we are at the time of the year where vaccinations are scheduled
        #if self.check_vaccination_period(t): 
        # generate a pool of candidates who are of the target age for the vaccine
        candidates = self.get_candidates(t, P)
    #   print "VACCINE:", self.label, "vaccinate_pop: candidates:", len(candidates)
        for ind in candidates:
            if self.check_coverage(t, ind, disease=self, rng=rng):
                self.apply_vaccine(ind, states)
                    


    def check_basic_conditions_met(self, t_years, ind):
        """
        Returns true if basic conditions are not met:
        -- vaccine is operating
        -- individual is of appropriate age to receive vaccine
        -- individual has not previously received vaccine
        """
        # fail if basic vaccine conditions are not met
        fail = t_years < 0  # or (self.label in ind.vaccines_received)

        # fail if vaccine has unmet preconditions (i.e., booster dependant on primary)
        if self.precond:
            if not set(self.precond).issubset(ind.vaccines_received):
                fail = True

        return not fail

    def check_success(self, ind,rng):
        """
        Returns true if vaccination is successful (ie does not fail)
        :param ind:
        :return:
        """
        # check for vaccine failure
        if self.fail_prob > 0.0:

            # if precondition vaccines were successful, this one is too
            if self.precond:
                p_index = ind.vaccines_received.index(self.precond[0])
                if ind.vaccines_effective[p_index]:
                    ind.vaccines_effective.append(True)
                    return True

            if rng.random() <= self.fail_prob:
                ind.vaccines_effective.append(False)
                return False
            else:
                ind.vaccines_effective.append(True)
                return True
        else:
            ind.vaccines_effective.append(True)
            return True

    def check_coverage(self, t, ind, **kwargs):
        """
        Check if individual will receive vaccine.  Must be defined by subclasses.
        """
        pass

    def give_vaccine(self, t, ind):
        """
        Administer vaccine to individual.
        """
#        print "Vaccine", self.label, "given!"
        self.count += 1
        self.count_ever += 1
        ind.vaccines_received.append(self.label)
        ind.vaccine_times.append(t)
    def apply_vaccine(self, ind, states):
        """
        Implement the effects of the vaccine (as may be ineffectvie)
        """
        ind.next_state = states[self.v_state]
        
        

class VaccineSingle(VaccineBase):
    """
    Basic single vaccination, decision made on an individual basis.
    """

    def __init__(self, label,  min_age,max_age, age_days, coverage, t_per_year, offset,
                 noncompliance_percentage, precond=None, fail_prob=0.0,v_state='S',vaccination_per_year=1,immediately_effective=True):
        super(VaccineSingle, self).__init__(label,  min_age,max_age, age_days, 
                 coverage, t_per_year, offset, noncompliance_percentage,
                 precond, fail_prob, v_state,vaccination_per_year,immediately_effective)
    
        self.hh_compliance = {}

    def check_compliance(self, ind, rng,t_years, hh_compliance):
        hh_id = ind.groups['household']

        # if this compliance decision for a given hh 
        #has been made previously, stick with that
        if hh_id in hh_compliance:
            compliance = hh_compliance[hh_id]
        # otherwise work out whether compliance will occur
        else:
            compliance = (rng.random() < (1 - self.noncompliance_percentage))

        # store for future
        hh_compliance[hh_id] = compliance

        return compliance

    def check_coverage(self, t, ind, rng, disease, **kwargs):
        # if passed end of set coverage values, use final value
        t_years = self.get_t_years(t)

        if not self.check_basic_conditions_met(t_years, ind):
            return False

        if self.coverage[t_years] == 0:
            return False

        # check if individual receives vaccine
        if self.check_compliance(ind, rng, t_years, self.hh_compliance) and\
            rng.random() <= ((self.coverage[t_years])/(1 -  self.noncompliance_percentage)):
            self.give_vaccine(t, ind)
        else:
            return False

        return self.check_success(ind,rng)


class VaccineHH(VaccineBase):
    """
    Vaccination on the basis of households.  Decision is made for firstborn and applies for all
    subsequent children born to that household.
    """

    def __init__(self, label,  min_age,max_age, age_days, coverage, t_per_year, 
                 offset, noncompliance_percentage,precond=None, fail_prob=0.0,           v_state='S',vaccination_per_year=1,immediately_effective=True):
        super(VaccineHH, self).__init__(label,  min_age,max_age, age_days, coverage, 
                                        t_per_year, offset, noncompliance_percentage,
                                        precond, fail_prob, v_state,vaccination_per_year,immediately_effective)
        self.hh_vacc = {}
        self.hh_compliance = {}

    def check_compliance(self, ind, rng,t_years, hh_compliance):
        hh_id = ind.groups['household']

        # if this compliance decision for a given hh 
        #has been made previously, stick with that
        if hh_id in hh_compliance:
            compliance = hh_compliance[hh_id]
        # otherwise work out whether compliance will occur
        else:
            compliance = (rng.random() < (1 - self.noncompliance_percentage))

        # store for future
        hh_compliance[hh_id] = compliance

        return compliance
    
    def check_household(self, ind, rng,t_years, hh_vacc):
        hh_id = ind.groups['household']

        # if this hh has made a previous vaccine decision, stick with that
        if hh_id in hh_vacc:
            vaccine_given = hh_vacc[hh_id]
        # otherwise work out whether vaccination will occur
        #assuming that hh compliance has been checked before
        else:
            vaccine_given = (rng.random() < ((self.coverage[t_years])/(1 - self.noncompliance_percentage)))

        # store for future
        hh_vacc[hh_id] = vaccine_given

        return vaccine_given

    def check_coverage(self, t, ind, rng, disease, **kwargs):

        t_years = self.get_t_years(t)

        if not self.check_basic_conditions_met(t_years, ind):
            return False

        if self.coverage[t_years] == 0:
            return False

        # check if individual receives vaccine
        if self.check_compliance(ind, rng, t_years, self.hh_compliance) and\
            self.check_household(ind, rng, t_years, self.hh_vacc):
            self.give_vaccine(t, ind)
        else:
            return False

        return self.check_success(ind,rng)
