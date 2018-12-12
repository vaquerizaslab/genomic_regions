"""
Helper functions for working with genomic regions.
"""

import re
import numpy as np


def human_format(num, precision=0):
    """
    Format a number as a string, suffixing letter for 1000 (K), 100000 (M), ...

    :param num: any number larger than zero
    :param precision: number of positions after decimal point
    :return: string representing the number
    """
    magnitude = 0
    while abs(num) >= 1000:
        magnitude += 1
        num /= 1000.0
    # add more suffixes if you need them
    return '{:.{prec}f}{}'.format(num, ['', 'k', 'M', 'G', 'T', 'P'][magnitude], prec=precision)


def str_to_int(num_string, decimal_separator='.', thousand_separator=','):
    """
    Convert a string denoting a genomic location to int (base pairs).

    :param num_string: input :class:`~str`
    :param decimal_separator: Decimal separator used for float conversion
    :param thousand_separator: Thousand separator (to be ignored)
    :return: int (base pairs)
    """
    try:
        num_string = num_string.replace(thousand_separator, '').lower()
    except AttributeError:
        pass

    try:
        return int(num_string)
    except ValueError:
        i = 0
        while i < len(num_string) and (num_string[i].isdigit() or num_string[i] == decimal_separator):
            i += 1

        try:
            number = float(num_string[:i])
            suffix = num_string[i:]

            multipliers = {
                'gb': 1000000000,
                'mb': 1000000,
                'kb': 1000,
                'bp': 1,
                'g': 1000000000,
                'm': 1000000,
                'k': 1000,
                'b': 1,
            }

            return int(number * multipliers[suffix])
        except (KeyError, ValueError):
            raise ValueError("Cannot convert '{}' to integer!".format(num_string))


def _convert(text):
    return int(text) if text.isdigit() else text.lower()


def natural_sort(l):
    """
    Sort list of chromosome names "naturally"

    :param l: :class:`~list` of strings
    :return: naturally sorted :class:`~list` of strings
    """
    return sorted(l, key=lambda key: [_convert(c) for c in re.split('([0-9]+)', key)])


def apply_sliding_func(a, window, func=np.ma.mean):
    """
    Apply function on a sliding window over an array, ignoring Numpy NaN values.

    :param a: Numpy array on which function is applied
    :param window: The sliding window is i - window:i + window + 1
                   so total window is twice this parameter.
    :param func: Function to apply
    """
    out = np.empty(a.shape)
    for i in range(len(a)):
        window_start = max(0, i - window)
        window_end = min(len(a), i + window + 1)
        cur_window = a[window_start:window_end]
        out[i] = func(cur_window[~np.isnan(cur_window)])
    return out


def intervals_weighted_mean(intervals):
    """
    Weighted mean for merging multiple interval scores into one.

    Takes into account the size of each interval as score weight.

    :param intervals: list of :class:`~tuple` (start, end, score)
    :return: float score
    """
    intervals = np.array(intervals)
    if len(intervals) == 0:
        return np.nan
    mask = np.isfinite(intervals[:, 2])
    valid = intervals[mask]
    if len(valid) == 0:
        return np.nan
    weights = (valid[:, 1] - valid[:, 0])
    weights += 1
    # safety
    weights = [weight if weight > 0 else 1 for weight in weights]
    return np.average(valid[:, 2], weights=weights)
