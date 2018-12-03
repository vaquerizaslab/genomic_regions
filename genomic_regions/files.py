import logging
import os.path
import pyBigWig
import sys
from collections import defaultdict

import numpy as np

try:
    from itertools import izip as zip
except ImportError:
    pass
import os
import errno

# configure logging
logger = logging.getLogger(__name__)


def mkdir(dir_name):
    dir_name = os.path.expanduser(dir_name)

    try:
        os.makedirs(dir_name)
    except OSError as exception:
        if exception.errno != errno.EEXIST:
            raise

    if not dir_name.endswith('/'):
        dir_name += '/'

    return dir_name


def which(program):
    """
    Check if executable exists in PATH
    :param program: executable name or path to executable
    :return: full path if found, None if not
    """
    import os

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


def write_bed(file_name, regions, mode='w', **kwargs):
    if file_name == '-':
        bed_file = sys.stdout
        must_close = False
    elif hasattr(file_name, 'write'):
        must_close = False
        bed_file = file_name
    else:
        bed_file = open(file_name, mode)
        must_close = True

    try:
        for region in regions:
            bed_file.write(region.as_bed_line(**kwargs) + '\n')
    finally:
        if must_close:
            bed_file.close()
        else:
            bed_file.flush()

    return file_name


def write_gff(file_name, regions, mode='w', **kwargs):
    if file_name == '-':
        gff_file = sys.stdout
        must_close = False
    elif hasattr(file_name, 'write'):
        must_close = False
        gff_file = file_name
    else:
        gff_file = open(file_name, mode)
        must_close = True

    try:
        for region in regions:
            gff_file.write(region.as_gff_line(**kwargs) + '\n')
    finally:
        if must_close:
            gff_file.close()
        else:
            gff_file.flush()

    return file_name


def write_bigwig(file_name, regions, mode='w', score_field='score'):
    logger.debug("Writing output...")
    bw = pyBigWig.open(file_name, mode)
    # write header

    chromosomes = []
    chromosome_lengths = defaultdict(int)
    interval_chromosomes = []
    interval_starts = []
    interval_ends = []
    interval_values = []
    for region in regions:
        if not isinstance(region.chromosome, str):
            chromosome = region.chromosome.decode() if isinstance(region.chromosome, bytes) \
                else region.chromosome.encode('ascii', 'ignore')
        else:
            chromosome = region.chromosome

        if chromosome not in chromosome_lengths:
            chromosomes.append(chromosome)
        chromosome_lengths[chromosome] = region.end

        interval_chromosomes.append(chromosome)
        interval_starts.append(region.start - 1)
        interval_ends.append(region.end)
        try:
            score = float(getattr(region, score_field))
        except AttributeError:
            score = np.nan
        interval_values.append(score)

    header = []
    for chromosome in chromosomes:
        chromosome_length = chromosome_lengths[chromosome]
        header.append((chromosome, chromosome_length))
    print(header)
    bw.addHeader(header)

    bw.addEntries(interval_chromosomes, interval_starts, ends=interval_ends, values=interval_values)

    bw.close()
    return file_name
