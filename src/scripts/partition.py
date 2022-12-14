class Partition:
    def __init__(  # pylint:  disable=too-many-arguments
        self, alignment, by_type, reads=None, coords=None, start=None
    ):

        self.alignment = alignment
        self.reads = reads if reads else []
        self.by_type = by_type
        self.coords = coords if coords else []
        self.start = start
        self.number = None
        self.parent = None

    def get_read_names(self):
        return self.alignment.get_read_names(self.reads)

    def get_coords(self, bases_freq_threshold):
        return self.alignment.get_coords(bases_freq_threshold, self)

    def get_sub_arr(self, coords):
        sub_arr = self.alignment.arr[self.reads]
        sub_arr = sub_arr[:, coords]
        return sub_arr

    def get_reference_size(self):
        return self.alignment.get_reference_size()
