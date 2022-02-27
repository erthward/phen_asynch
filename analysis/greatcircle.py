"""
DETH: 05-01-21: taken from:
    https://python.plainenglish.io/calculating-great-circle-distances-
                                                        in-python-cf98f64c1ea0
"""
import math

DEGREES_IN_RADIAN = 57.29577951
MEAN_EARTH_RADIUS_KM = 6371
KILOMETRES_IN_MILE = 1.60934


class GreatCircle(object):

    """
    This class has attributes for the names and locations of a pair of cities,
    and for their distances.
    After the city attributes are set, call calculate to set the distance attributes.
    Validation is carried out and the valid attributes set. This should be checked before using
    calculated attributes.
    """

    def __init__(self):

        """
        Create a set of attributes with default values.
        """

        self.name1 = None
        self.latitude1_degrees = 0
        self.longitude1_degrees = 0
        self.latitude1_radians = 0
        self.longitude1_radians = 0

        self.name2 = None
        self.latitude2_degrees = 0
        self.longitude2_degrees = 0
        self.latitude2_radians = 0
        self.longitude2_radians = 0

        self.central_angle_radians = 0
        self.central_angle_degrees = 0
        self.distance_kilometres = 0
        self.distance_miles = 0
        self.valid = False

    def calculate(self):

        """
        Central method to set calculated attributes, which it
        does by calling other private functions.
        """

        self.__validate_degrees()

        if self.valid:
            self.__calculate_radians()
            self.__calculate_central_angle()
            self.__calculate_distance()

    def __validate_degrees(self):

        """
        Check latitudes and longitudes are within valid ranges,
        setting the valid attribute accordingly.
        """

        self.valid = True

        if self.latitude1_degrees < -90.0 or self.latitude1_degrees > 90.0:
            self.valid = False

        if self.longitude1_degrees < -180.0 or self.longitude1_degrees > 180.0:
            self.valid = False

        if self.latitude2_degrees < -90.0 or self.latitude2_degrees > 90.0:
            self.valid = False

        if self.longitude2_degrees < -180.0 or self.longitude2_degrees > 180.0:
            self.valid = False

    def __calculate_radians(self):

        """
        Calculate radians from degrees by dividing by constant.
        """

        self.latitude1_radians = self.latitude1_degrees / DEGREES_IN_RADIAN
        self.longitude1_radians = self.longitude1_degrees / DEGREES_IN_RADIAN

        self.latitude2_radians = self.latitude2_degrees / DEGREES_IN_RADIAN
        self.longitude2_radians = self.longitude2_degrees / DEGREES_IN_RADIAN

    def __calculate_central_angle(self):

        """
        Slightly complex formula for calculating the central angle
        between two points on the surface of a sphere.
        """

        if self.longitude1_radians > self.longitude2_radians:
            longitudes_abs_diff = self.longitude1_radians - self.longitude2_radians
        else:
            longitudes_abs_diff = self.longitude2_radians - self.longitude1_radians

        self.central_angle_radians = math.acos( math.sin(self.latitude1_radians)
                                         * math.sin(self.latitude2_radians)
                                         + math.cos(self.latitude1_radians)
                                         * math.cos(self.latitude2_radians)
                                         * math.cos(longitudes_abs_diff))

        self.central_angle_degrees = self.central_angle_radians * DEGREES_IN_RADIAN

    def __calculate_distance(self):

        """
        Because we are using radians, this is a simple formula multiplying the radius
        by the angle, the actual units used being irrelevant.
        Also the distance in miles is calculated from kilometres.
        """

        self.distance_kilometres = MEAN_EARTH_RADIUS_KM * self.central_angle_radians

        self.distance_miles = self.distance_kilometres / KILOMETRES_IN_MILE 
