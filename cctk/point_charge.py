import numpy as np

class PointCharge():
    """
    Represents a point charge.

    Attributes:
        coordinates (np.ndarray): 3-element ndarray
        charge (float): charge
    """

    def __init__(self, coordinates, charge):
        assert isinstance(coordinates, (np.ndarray, list)), "coordinates must be list or ndarray!"
        assert len(coordinates) == 3, "coordinates must have len 3!"
        self.coordinates = np.array(coordinates)

        assert isinstance(charge, (float, int)), "charge must be numeric"
        self.charge = float(charge)
