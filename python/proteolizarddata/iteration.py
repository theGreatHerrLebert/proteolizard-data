import numpy as np


class SliceIterator:
    def __init__(self, dh, slice_length=2.5, load_precursors=True, load_fragments=False):
        """
        An iterator for easy iteration of a timsTOF Experiment, going over the data in slices at a time
        :param dh: a data handle to connect to a timsTOF Experiment
        :param slice_length: length of slices of data to be loaded in MINUTES
        :param load_precursors: if true, precursor frames will be loaded into slice
        :param load_fragments: if true, fragment frames will be loaded into slice
        """
        self.dh = dh
        self.rts = self.dh.meta_data.Time.values / 60
        self.ids = self.dh.meta_data.Id.values

        self.precursor_ids = set(self.dh.precursor_frames)
        self.fragment_ids = set(self.dh.fragment_frames)

        self.batch_dict = self.group_batches(slice_length)
        self.batch_counter = 0

        self.__load_precursors = load_precursors
        self.__load_fragments = load_fragments

        self.slice = self.get_slice(precursors=self.__load_precursors, fragments=self.__load_fragments)

    def get_slice(self, precursors=True, fragments=False):
        """

        :param precursors:
        :param fragments:
        :return:
        """
        assert precursors | fragments, "Need to at least load precursors or fragmets."
        batch = self.batch_dict[self.batch_counter]

        precursor_ids, fragment_ids = [], []

        if precursors:
            precursor_ids = self.filter_batch(batch, return_precursors=True)

        if fragments:
            fragment_ids = self.filter_batch(batch, return_precursors=True)

        return self.dh.get_slice(precursor_ids, fragment_ids)

    def __iter__(self):
        return self.slice

    def __next__(self):

        slic = self.slice

        if self.batch_counter < len(self.batch_dict) - 1:
            self.batch_counter += 1
            self.slice = self.get_slice(precursors=self.__load_precursors, fragments=self.__load_fragments)
            return slic

        raise StopIteration

    def group_batches(self, slice_length):
        """

        :param slice_length:
        :return:
        """
        batch_ids = [int(np.floor(x / slice_length)) for x in self.rts]
        batch_dict = {}

        for b, i in zip(batch_ids, self.ids):
            if b in batch_dict:
                batch_dict[b].append(i)

            else:
                batch_dict[b] = [i]

        return batch_dict

    def filter_batch(self, id_batch, return_precursors=True):
        """

        :param id_batch:
        :param return_precursors:
        :return:
        """
        if return_precursors:
            return np.sort(np.array(list(filter(lambda x: x in self.precursor_ids, id_batch))))
        return np.sort(np.array(list(filter(lambda x: x in self.fragment_ids, id_batch))))
