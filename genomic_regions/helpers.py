import os
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
    num_string = num_string.replace(thousand_separator, '').lower()
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


def natural_sort(l):
    def convert(text):
        return int(text) if text.isdigit() else text.lower()
    return sorted(l, key=lambda key: [convert(c) for c in re.split('([0-9]+)', key)])


def natural_cmp(pa, pb):
    i = 0
    j = 0

    if pa[-2] == '/' == pb[-2]:
        pa = pa[:-2]
        pb = pb[:-2]

    try:
        while i < len(pa) and j < len(pb):
            if pa[i].isdigit() and pb[j].isdigit():
                # digit comparison
                while pa[i] == '0':
                    i += 1
                while pb[j] == '0':
                    j += 1
                while pa[i].isdigit() and pb[j].isdigit() and pa[i] == pb[j]:
                    i += 1
                    j += 1
                if pa[i].isdigit() and pb[j].isdigit():
                    k = 0
                    try:
                        while pa[i+k].isdigit() and pb[j+k].isdigit():
                            k += 1
                    except IndexError:
                        if i+k < len(pa):
                            return 1
                        if j+k < len(pb):
                            return -1
                        # both ran over index
                        return int(pa[i]) - int(pb[j])
                    return 1 if pa[i+k].isdigit() else -1 if pb[j+k].isdigit() else int(pa[i]) - int(pb[j])
                elif pa[i].isdigit():
                    return 1
                elif pb[j].isdigit():
                    return -1
                elif i != j:
                    return 1 if i < j else -1
            else:
                # string comparison
                if pa[i] != pb[j]:
                    return -1 if pa[i] < pb[j] else 1
                i += 1
                j += 1
    except IndexError:
        pass
    return 1 if i < len(pa) else -1 if j < len(pb) else 0


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


def which(program):
    """
    Check if executable exists in PATH
    :param program: executable name or path to executable
    :return: full path if found, None if not
    """
    def is_exe(fpath):
        return os.path.isfile(fpath) and os.access(fpath, os.X_OK)

    fpath, fname = os.path.split(program)
    if fpath:
        if is_exe(program):
            return program
    else:
        for path in os.environ["PATH"].split(os.pathsep):
            path = path.strip('"')
            exe_file = os.path.join(path, program)
            if is_exe(exe_file):
                return exe_file

    return None
