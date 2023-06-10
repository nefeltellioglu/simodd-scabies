import tables as tb

"""
Base class for disease observers.
"""

class Observer(object):

    def __init__(self, h5file, label, description, title):
        self.label = label
        self.h5file = h5file
        self.exists = False
        self.data = None
        self.row = None
        if self.h5file.mode in ['w', 'a']:
            self.create_storage(description, title)
        else:
            self.load_storage()

    def create_storage(self, description, title):
        if '/%s' % self.label not in self.h5file:
            group = self.h5file.create_group('/', self.label, title)
        filters = tb.Filters(complevel=9)   # TODO: investigate BLOSC (?)

        if description:
            self.data = self.h5file.create_table(group, 'base_data', description, filters=filters)
            self.row = self.data.row
        self.exists = True
#    __create_storage = create_storage

    def load_storage(self):
        if self.h5file.__contains__('/%s' % self.label):
            if 'base_data' in self.h5file.get_node(self.h5file.root, self.label):
                self.data = self.h5file.get_node(self.h5file.root, self.label).base_data
                self.row = self.data.row
            self.exists = True
        else:
            self.exists = False
#    __load_storage = load_storage

    def update(self, t, **kwargs):
        pass

