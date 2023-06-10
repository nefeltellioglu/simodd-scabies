"""
Base class for disease state.
"""


class State(object):
    def __init__(self, order, label, color=None):
        self.order = order
        self.count = 0
        self.label = label
        self.color = color
        self.at_risk = False
        self.infectious = False

    def enter(self, t, ind):
        """Enter the state.  Perform any necessary initialisation."""
        self.count += 1
        ind.state = self

    def exit(self, t, ind):
        """Exit the state.  Perform any necessary clean-up"""
        if ind.state:
            ind.state.count -= 1
            ind.prev_state = ind.state

    def update(self, t, ind, states):
        """Update the state of this individual according to internal logic."""
        # decrement all counters and remove counters that have 'finished'
        labels = [key for key in ind.counters]
        for key in labels:
            if ind.counters[key] < 1:
                del ind.counters[key]
            else:
                ind.counters[key] -= 1


class StateTimed(State):
    """
    A state that an individual only belongs to for a specified (randomly sampled) duration.
    """

    def __init__(self, order, duration, label, color=None):
        super(StateTimed, self).__init__(order, label, color)
        self.duration = duration

    def enter(self, t, ind):
        super(StateTimed, self).enter(t, ind)
        ind.durations[self.label] = self.duration.get_duration()
        ind.counters[self.label] = ind.durations[self.label]

    def update(self, t, ind, states):
        super(StateTimed, self).update(t, ind, states)
        # output True if time in state has 'finished'
        if self.label not in ind.counters:
            return True
        else:
            return False
