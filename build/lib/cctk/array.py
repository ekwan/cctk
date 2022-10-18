import numpy as np
import copy

class OneIndexedArray(np.ndarray):
    """
    Wrapper for ``np.ndarray`` that's indexed from one, not zero, to store atomic numbers and geometries.
    This only works on 1D or 2D arrays. Additionally, only the first index of a 2D array will be 1-indexed.

    Note that ``enumerate(one_indexed_array)`` will throw ``IndexError`` -- instead, use ``enumerate(one_indexed_array, start=1)``.
    """

    def __new__(cls, obj, **kwargs):
        new = np.array(obj, **kwargs).view(cls)
        return new

    def __getitem__(self, index):
        index = copy.deepcopy(index)
        if isinstance(index, slice):
            if index.start is None:
                start = 0
            else:
                start = index.start - 1
            if index.stop is None:
                stop = -1
            else:
                stop = index.stop - 1
            new_index = slice(start, stop, index.step)
            return super().__getitem__(new_index)
        elif isinstance(index, int):
            if index > 0:
                return super().__getitem__(index-1)
            elif index == 0:
                raise IndexError("this is a 1-indexed array: no element 0!")
            elif index < 0:
                return super().__getitem__(index)
        elif (isinstance(index, tuple)) and (len(index) == 2):
            if index[0] is None:
                return super().__getitem__((index[0], index[1]))
            elif index[0] > 0:
                return super().__getitem__((index[0]-1, index[1]))
            elif index[0] == 0:
                raise IndexError("this is a 1-indexed array: no element 0!")
            elif index[0] < 0:
                return super().__getitem__((index[0], index[1]))
        elif (isinstance(index, tuple)) and (len(index) == 1):
            return self.__getitem__(index[0])
        elif isinstance(index, np.ndarray):
            if index.dtype == bool:
                return super().__getitem__(index)
            elif index.ndim == 1:
                index[index >= 1] += -1
                return super().__getitem__(index)
            else:
                index[0][index >= 1] += -1
                return super().__getitem__(index)
        elif isinstance(index, list):
            if isinstance(index[0], bool):
                return super().__getitem__(index)
            elif isinstance(index[0], list):
                if isinstance(index[0][0], bool):
                    return super().__getitem__(index)
                for i, v in enumerate(index[0]):
                    if v >= 1:
                        index[i] += -1
                return super().__getitem__(index)
            else:
                for i, v in enumerate(index):
                    if v >= 1:
                        index[i] += -1
                return super().__getitem__(index)
        else:
            return super().__getitem__(index)

    def __setitem__(self, index, value):
        index = copy.deepcopy(index)
        if isinstance(index, int):
            if index > 0:
                if self.ndim == 1:
                    super().__setitem__(index-1, value)
                elif self.ndim == 2:
                    super().__setitem__(index, value)
                else:
                    raise TypeError("this datatype is only defined for 1D and 2D ndarrays")
            elif index == 0:
                raise IndexError("this is a 1-indexed array: no element 0!")
            elif index < 0:
                super().__setitem__(index, value)
        elif (isinstance(index, tuple)) and (len(index) == 2):
            if index[0] is None:
                super().__setitem__((index[0], index[1]), value)
            elif index[0] > 0:
                super().__setitem__((index[0]-1, index[1]), value)
            elif index[0] == 0:
                raise IndexError("this is a 1-indexed array: no element 0!")
            elif index[0] < 0:
                super().__setitem__((index[0], index[1]), value)
        elif (isinstance(index, tuple)) and (len(index) == 1):
            return self.__setitem__(index[0], value)
        elif isinstance(index, np.ndarray):
            if index.dtype == bool:
                super().__setitem__(index, value)
            elif index.ndim == 1:
                index[index >= 1] += -1
                super().__setitem__(index, value)
            else:
                index[0][index >= 1] += -1
                super().__setitem__(index, value)
        elif isinstance(index, list):
            if isinstance(index[0], bool):
                super().__setitem__(index, value)
            elif isinstance(index[0], list):
                if isinstance(index[0][0], bool):
                    super().__setitem__(index, value)
                for i, v in enumerate(index[0]):
                    if v >= 1:
                        index[i] += -1
                super().__setitem__(index, value)
            else:
                for i, v in enumerate(index):
                    if v >= 1:
                        index[i] += -1
                super().__setitem__(index, value)
        else:
            super().__setitem__(index, value)
#            raise IndexError(f"invalid index {index} for OneIndexedArray")

    def __iter__(self):
        for idx in range(1,len(self)+1):
            yield self.__getitem__(idx)

    def __hash__(self):
        return hash(self.data.tobytes())

    def __str__(self):
        return self.view(np.ndarray).__str__()
