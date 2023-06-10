"""
Basic disease states.
"""

from math import exp

from disease.general.state_base import State, StateTimed


class Susceptible(State):
    """Susceptible state"""

    def __init__(self, order):
        super(Susceptible, self).__init__(order, 'S', 'k')
        self.at_risk = True
        self.infectious = False

    def test_exposure(self, states, ind, foi, rng):
        if rng.random() < 1.0 - exp(-foi):
            if 'E' in states:
                ind.next_state = states['E']
            else:
                ind.next_state = states['I']
            return "infection"
        else:
            return False


###############################################################################

class Exposed(StateTimed):
    """Exposed state"""

    def __init__(self, order, duration):
        super(Exposed, self).__init__(order, duration, 'E', 'o')
        self.at_risk = False
        self.infectious = False

    def update(self, t, ind, states):
        if super(Exposed, self).update(t, ind, states):
            ind.next_state = states['I']


###############################################################################

class Infected(StateTimed):
    """Infected state"""

    def __init__(self, order, duration):
        super(Infected, self).__init__(order, duration, 'I', 'r')
        self.at_risk = False
        self.infectious = True
        self.current = set()

    def enter(self, t, ind):
        super(Infected, self).enter(t, ind)
        self.current.add(ind.ID)

    def exit(self, t, ind):
        super(Infected, self).exit(t, ind)
        self.current.remove(ind.ID)

    def update(self, t, ind, states):
        if super(Infected, self).update(t, ind, states):
    #        ind.next_state = states['R']
            ind.next_state = states['S']

    

###############################################################################

class Removed(State):
    """Removed state"""

    def __init__(self, order, label='R', color='b'):
        super(Removed, self).__init__(order, label, color)
        self.at_risk = False
        self.infectious = False


class RemovedTemp(StateTimed):
    """Removed state with timing -- (ie, non-permanent immunity)"""

    def __init__(self, order, duration, label='R', color='b'):
        super(RemovedTemp, self).__init__(order, duration, label, color)
        self.at_risk = False
        self.infectious = False

    def update(self, t, ind, states):
        if super(RemovedTemp, self).update(t, ind, states):
            ind.next_state = states['S']


            ###############################################################################


class Maternal(StateTimed):
    """Maternal immunity state"""

    def __init__(self, order, duration, label='M', color='b'):  # same color as removed state ?
        super(Maternal, self).__init__(order, duration, label, color)
        self.at_risk = False
        self.infectious = False

    def update(self, t, ind, states):
        if super(Maternal, self).update(t, ind, states):
            ind.next_state = states['S']


#class Susceptible_Maternal(Susceptible):
# NB: I think this should be deprecated with new implementation of MSIR model.
#
#    def __init__(self, successors, order, duration, maternal_dur):
#        super(Susceptible_Maternal, self).__init__(successors, order, duration)
#        self.maternal_dur = maternal_dur
#        
#    def enter(self, t, ind):
#        super(Susceptible_Maternal, self).enter(t, ind)
#        if ind.age==0 and ind.age_days==0:
#            mother = ind.parents[0]
#            # transfer a newborn into the maternal class if their mother's most recent
#            # exposure to infection was less than 'maternal_dur' timesteps ago
#            if len(mother.infections) > 0 and (t - mother.infections[-1] < self.maternal_dur):
#                ind.next_state = self.successors['M']
