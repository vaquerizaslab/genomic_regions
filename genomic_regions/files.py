"""
Module providing functions for working with files
revolving around :class:`~GenomicRegion` objects.
"""

import logging
import os.path
import pyBigWig
import sys
from collections import defaultdict
import tempfile
import shutil
import pysam
import subprocess

import numpy as np

try:
    from itertools import izip as zip
except ImportError:
    pass
import os

# configure logging
logger = logging.getLogger(__name__)


__all__ = ['write_bed', 'write_bigwig', 'write_gff']


def create_temporary_copy(src_file_name, preserve_extension=True):
    """
    Copies the source file into a temporary file.
    Returns a _TemporaryFileWrapper, whose destructor deletes the temp file
    (i.e. the temp file is deleted when the object goes out of scope).
    """
    src_file_name = os.path.expanduser(src_file_name)
    tf_suffix = ''
    if preserve_extension:
        _, tf_suffix = os.path.splitext(src_file_name)
    tf = tempfile.NamedTemporaryFile(suffix=tf_suffix, delete=False)
    shutil.copy2(src_file_name, tf.name)

    logger.debug("Working from temporary output file: {}. "
                 "Original file: {}".format(tf.name, src_file_name))

    return tf.name


def create_temporary_output(original_file_name, preserve_extension=True):
    tf_suffix = ''
    if preserve_extension:
        _, tf_suffix = os.path.splitext(original_file_name)
    tmp_file = tempfile.NamedTemporaryFile(delete=False, suffix=tf_suffix)
    tmp_file.close()
    logger.debug("Working from temporary output file: {}. "
                 "Original file: {}".format(tmp_file.name, original_file_name))
    return tmp_file.name


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


def read_chromosome_sizes(file_name):
    """
    Obtain chromosome sizes from <genome>.chrom.sizes file (UCSC).

    :param file_name: Path to <genome>.chrom.sizes file
    :return: :class:`~dict` {chromosome: size, ...}
    """
    chrom_sizes = {}
    with open(os.path.expanduser(file_name), 'r') as chrom_sizes_file:
        for line in chrom_sizes_file:
            line = line.rstrip()
            if line != '':
                chromosome, chromosome_length = line.split("\t")
                chrom_sizes[chromosome] = int(chromosome_length)
    return chrom_sizes


def write_bed(file_name, regions, mode='w', **kwargs):
    """
    Write :class:`~GenomicRegion` objects to BED file.

    :param file_name: Path to output BED file
    :param regions: :class:`~list` of :class:`~GenomicRegion` objects
    :param mode: File mode. Default: 'w'
    :param kwargs: Keyword arguments passed to GenomicRegion.as_bed_line
    :return: Path to output BED file
    """
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
    """
    Write :class:`~GenomicRegion` objects to GFF file.

    :param file_name: Path to output GFF file
    :param regions: :class:`~list` of :class:`~GenomicRegion` objects
    :param mode: File mode. Default: 'w'
    :param kwargs: Keyword arguments passed to GenomicRegion.as_gff_line
    :return: Path to output GFF file
    """
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


def write_bigwig(file_name, regions, mode='w', score_field='score',
                 chromosome_sizes=None):
    """
    Write :class:`~GenomicRegion` objects to BigWig file.

    :param file_name: Path to output BigWig file
    :param regions: :class:`~list` of :class:`~GenomicRegion` objects
    :param mode: File mode. Default: 'w'
    :param score_field: Attribute from each object to use as interval score.
    :param chromosome_sizes: Dictionary of the form {chromosome: size}
    :return: Path to output BigWig file
    """
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

    if chromosome_sizes is None:
        chromosome_sizes = chromosome_lengths

    header = []
    for chromosome in chromosomes:
        chromosome_length = chromosome_sizes[chromosome]
        header.append((chromosome, chromosome_length))
    bw.addHeader(header)

    bw.addEntries(interval_chromosomes, interval_starts, ends=interval_ends, values=interval_values)

    bw.close()
    return file_name


def sort_natural_sam(sam_file, output_file=None, sambamba=True, _sambamba_path='sambamba'):
    if which(_sambamba_path) is None and sambamba:
        logger.info('Cannot find {} on this machine, falling back to samtools sort. '
                    'This is not a problem, but if you want to speed up your SAM/BAM '
                    'file sorting, install sambamba and ensure it is in your PATH!'.format(_sambamba_path))
        sambamba = False

    basename, extension = os.path.splitext(sam_file)

    replace_input = False
    if output_file is None:
        with tempfile.NamedTemporaryFile(delete=False, prefix='gr_', suffix=extension) as f:
            output_file = f.name
        replace_input = True

    if sambamba:
        sambamba_command = [_sambamba_path, 'sort', '-N', '-o', output_file, sam_file]
        ret = subprocess.call(sambamba_command)
        if ret != 0:
            sambamba = False
            logger.warning("{} failed, falling back to pysam/samtools".format(_sambamba_path))
    if not sambamba:
        pysam.sort('-n', '-o', output_file, sam_file)

    if replace_input:
        logger.info("Replacing input SAM/BAM file {} with sorted version...".format(sam_file))
        shutil.copy(output_file, sam_file)
        os.remove(output_file)
        output_file = sam_file

    return output_file
