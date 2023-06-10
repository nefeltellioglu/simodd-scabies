"""
.. module:: household
.. moduleauthor:: Nic Geard <ngeard@unimelb.edu.au>
"""


class Household:
    """
    A household object, only used for logging purposes.  That is, it has no functional role, and
    can safely be omitted.

    :param time: time household was founded
    :type time: int
    :param code: a code for filtering events.
    :type code: string
    :param msg: a text message describing this event.
    :type msg: string
    :param size: current number of people in household.
    :type size: int
    :param adam: `True` if individual is part of the bootstrapped population.
    :type adam: bool

    """

    def __init__(self, time=0, code=None, msg=None, size=None, adam=False):
        self.founded = time
        self.log = []
        self.snapshots = []
        self.adam = adam
        if code is not None:
            self.add_log(time, code, msg, size)

    def add_log(self, time, code, msg, size=None):
        """
        Add a new log item.

        :param code: a code for filtering events.
        :type code: string
        :param msg: a text message describing this event.
        :type msg: string
        :param size: current number of people in household.
        :type size: int

        Possible log codes for households:

        - f: formation
        - cb: birth of child
        - cd: death of child (dependent)
        - c+: gain of child (relocated)
        - gd: death of adult (guardian)
        - m: merging (due to remarriage)
        - s: separation

        """

        if size is None:
            size = self.log[-1]['size']
        self.log.append({
            'age': time - self.founded,
            'code': code,
            'msg': msg,
            'time': time,
            'size': size
        })

    def size_at_time(self, t):
        """
        On the basis of stored logs, retrieve the household size at the 
        specified time.

        :param t: time to retrieve size for.
        :type t: int
        :returns: the size of the household at the specified time.
        """

        # return -1 if no logs yet, or t occurs before creation of hh
        if not self.log or self.log[0]['time'] > t:
            return -1
        size = self.log[0]['size']
        for cur_item in self.log:
            if cur_item['time'] < t:
                size = cur_item['size']
            else:
                return size
        return self.log[-1]['size']

    def count_events_in_range(self, code, begin, end):
        """
        On the basis of stored logs, retrieve the number of events of type
        code that occurred between the specified beginning and end times.

        :param code: a code for filtering events.
        :type code: string
        :param begin: the start of the measurement time period.
        :type begin: int
        :param end: the end of the measurement time period.
        :type end: int
        :returns: the number of events occurring.
        """
        return len([
            x for x in self.log
            if x['code'] == code and begin < x['time'] < end])
