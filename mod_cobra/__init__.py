__author__ = 'anna'


ZERO_THRESHOLD = 1e-4
ROUND_COEFFICIENT = 4


def round_value(value):
    return round(float(value), ROUND_COEFFICIENT)
