"""
This module provides functions and classes to work with genomic regions
(also referred to as genomic intervals).

Its main classes are:

* :class:`~GenomicRegion`: A class that represents a genomic region/interval
* :class:`~RegionBased`: The base class for collections of genomic regions

The aim of this module, besides providing an intuitive set of tools operating
on genomic regions and collections thereof, is to supply a unified interface
for the different representations of genomic region sets. Specifically, it
gives the user access to the same methods with identical syntax regardless of
what type of genomic regions file the user currently works with (BED, GFF,
BigWig, Tabix, ...).

Most of the time, is is enough to open a file with :func:`load` - the module
will figure out the underlying file type automatically. Please refer to the
documentation for further details.
"""

from __future__ import division, print_function

import os.path

import numpy as np
import pybedtools
import pysam

from .files import write_bigwig, write_bed, write_gff
from .helpers import natural_sort, str_to_int, apply_sliding_func, intervals_weighted_mean

try:
    from itertools import izip as zip
except ImportError:
    pass
from collections import defaultdict
import copy
import re
import shlex
import warnings
from bisect import bisect_left
from future.utils import string_types
from builtins import object
import pandas
import pyBigWig
import intervaltree
import logging

logger = logging.getLogger(__name__)

__all__ = ['GenomicRegion', 'Bed', 'Bedpe', 'RegionWrapper',
           'RegionBased', 'GenomicDataFrame', 'Tabix',
           'BigWig', 'as_region', 'load', 'merge_overlapping_regions']


def as_region(region):
    """
    Convert string to :class:`~GenomicRegion`.

    This function attempts to convert any string passed to it
    to a :class:`~GenomicRegion`. Strings are expected to be
    of the form <chromosome>[:<start>-<end>[:[strand]], e.g.
    chr1:1-1000, 2:2mb-5mb:-, chrX:1.5kb-3mb, ...

    Numbers can be abbreviated as '12k', '1.5Mb', etc.

    When fed a :class:`~GenomicRegion`, it will simply be
    returned, making the use of this function as an
    "if-necessary" converter possible.

    :param region: str or :class:`~GenomicRegion`
    :return: :class:`~GenomicRegion`
    """
    if isinstance(region, string_types):
        return GenomicRegion.from_string(region)
    elif isinstance(region, GenomicRegion):
        return region
    raise ValueError("region parameter cannot be converted to GenomicRegion!")


def merge_overlapping_regions(regions):
    """
    Merge overlapping regions in list.

    Provided with a list of :class:`~GenomicRegion` objects,
    this function will determine overlapping regions and merge
    them. The output is a list of non-overlapping regions.

    :param regions: :class:`list` of :class:`~GenomicRegion` objects
    :return: :class:`list` of merged :class:`~GenomicRegion` objects
    """
    sorted_regions = sorted(regions, key=lambda r: (r.chromosome, r.start))

    merged_regions = []
    current_regions = []
    last_end = None
    for region in sorted_regions:
        if len(current_regions) == 0:
            current_regions.append(region)
            last_end = region.end
        elif region.chromosome == current_regions[0].chromosome and region.start < last_end:
            current_regions.append(region)
            last_end = max(last_end, region.end)
        else:
            merged_region = GenomicRegion(chromosome=current_regions[0].chromosome,
                                          start=current_regions[0].start, end=last_end,
                                          strand=current_regions[0].strand)
            merged_regions.append(merged_region)
            current_regions = [region]
            last_end = region.end

    merged_region = GenomicRegion(chromosome=current_regions[0].chromosome,
                                  start=current_regions[0].start, end=last_end,
                                  strand=current_regions[0].strand)
    merged_regions.append(merged_region)

    return merged_regions


def load(file_name, *args, **kwargs):
    """
    Open file containing genomic regions as :class:`~RegionBased` object.

    'Magic' function that wraps a file containing genomic regions in a
    :class:`~RegionBased` interface, thus providing the same methods
    to different types of genomic data formats.

    Compatible formats include: BED, GFF/GTF, BigWig, and Tabix (i.e. BED,
    GFF, and compatible files index with tabix).

    SAM files are also detected, but opened using pysam.AlignmentFile,
    so they are not compatiable with the :class:`~RegionBased` interface
    (yet).

    :param file_name: Path to genomic regions file
    :param args: Additional arguments passed to downstream class
    :param kwargs: Additional keyword arguments passed to downstream class
    :return: :class:`~RegionBased`
    :raises: ValueError if file type not supported
    """
    import os
    file_name = os.path.expanduser(file_name)

    # SAM/BAM
    import pysam
    try:
        sb = pysam.AlignmentFile(file_name, 'rb')
        if kwargs.get('mode', 'r') != 'rb':
            sb.close()
            sb = pysam.AlignmentFile(file_name, *args, **kwargs)
        return sb
    except (ValueError, IOError):
        pass

    # Tabix
    try:
        f = Tabix(file_name, *args, **kwargs)
        return f
    except (IOError, OSError, ValueError, TypeError):
        pass

    # BEDPE
    if file_name.endswith('.bedpe'):
        try:
            f = Bedpe(file_name, *args, **kwargs)
            _ = f.regions[0]
            return f
        except (ValueError, TypeError):
            pass

    import pybedtools
    f = Bed(file_name, *args, **kwargs)
    try:
        ft = f.file_type
        if ft != 'empty':
            return f
    except (IndexError, pybedtools.MalformedBedLineError, UnicodeDecodeError):
        pass

    try:
        import pyBigWig
        f = pyBigWig.open(file_name, 'r')
        if kwargs.get('mode', 'r') != 'r':
            f.close()
            f = pyBigWig.open(file_name, *args, **kwargs)

        return BigWig(f)
    except (ImportError, RuntimeError):
        raise ValueError("File type not recognised ({}).".format(file_name))


class GenomicRegion(object):
    """
    Class representing a genomic region.

    .. attribute:: chromosome

        Name of the chromosome this region is located on

    .. attribute:: start

        Start position of the region in base pairs

    .. attribute:: end

        End position of the region in base pairs

    .. attribute:: strand

        Strand this region is on. Can be a str ('+', '-', '.'),
        None, or an int (+1, -1)

    .. attribute:: ix

        Index of the region in the context of a set of
         genomic regions.

    """

    def __init__(self, chromosome=None, start=None, end=None, strand=None, ix=None, **kwargs):
        """
        Initialize this object.

        :param chromosome: Name of the chromosome this region is located on
        :param start: Start position of the region in base pairs
        :param end: End position of the region in base pairs
        :param strand: Strand this region is on. Can be a str ('+', '-', '.'),
                       None, or an int (+1, -1)
        :param ix: Index of the region in the context of a set of genomic
                   regions.
        """
        self.start = start
        if end is None:
            end = start
        self.end = end
        if strand == "+":
            strand = 1
        elif strand == "-":
            strand = -1
        elif strand == 0 or strand == "0" or strand == ".":
            strand = None
        self.strand = strand
        self.chromosome = chromosome.decode() if isinstance(chromosome, bytes) else chromosome
        self.ix = ix

        for name, value in kwargs.items():
            self.set_attribute(name, value)

    def set_attribute(self, attribute, value):
        """
        Safely set an attribute on the :class:`~GenomicRegion` object.

        This automatically decodes bytes objects into UTF-8 strings.
        If you do not care about this, you can also use
        region.<attribute> = <value> directly.

        :param attribute: Name of the attribute to be set
        :param value: Value of the attribute to be set
        """
        setattr(self, attribute.decode() if isinstance(attribute, bytes) else attribute,
                value.decode() if isinstance(value, bytes) else value)

    @property
    def attributes(self):
        """
        Return all visible attributes of this :class:`~GenomicRegion`.

        Returns all attribute names that do not start with an underscore.
        :return: list of attribute names
        """
        return [name for name in self.__dict__.keys() if not name.startswith('_')]

    def update(self, **kwargs):
        for key, value in kwargs.items():
            self.set_attribute(key, value)

    @classmethod
    def from_string(cls, region_string):
        """
        Convert a string into a :class:`~GenomicRegion`.

        This is a very useful convenience function to quickly
        define a :class:`~GenomicRegion` object from a descriptor
        string. Numbers can be abbreviated as '12k', '1.5M', etc.

        :param region_string: A string of the form
                              <chromosome>[:<start>-<end>[:<strand>]]
                              (with square brackets indicating optional
                              parts of the string). If any optional
                              part of the string is omitted, intuitive
                              defaults will be chosen.
        :return: :class:`~GenomicRegion`
        """
        chromosome = None
        start = None
        end = None
        strand = None

        # strip whitespace
        no_space_region_string = "".join(region_string.split())
        fields = no_space_region_string.split(':')

        if len(fields) > 3:
            raise ValueError("Genomic range string must be of the form "
                             "<chromosome>[:<start>-<end>:[<strand>]]")

        # there is chromosome information
        if len(fields) > 0:
            chromosome = fields[0]

        # there is range information
        if len(fields) > 1 and fields[1] != '':
            start_end_bp = fields[1].split('-')
            if len(start_end_bp) > 0:
                start = str_to_int(start_end_bp[0])

            if len(start_end_bp) > 1:
                end = str_to_int(start_end_bp[1])

                if not end >= start:
                    raise ValueError("The end coordinate must be bigger than the start.")

        # there is strand information
        if len(fields) > 2:
            if fields[2] == '+' or fields[2] == '+1' or fields[2] == '1':
                strand = 1
            elif fields[2] == '-' or fields[2] == '-1':
                strand = -1
            else:
                raise ValueError("Strand only can be one of '+', '-', '+1', '-1', and '1'")
        return cls(chromosome=chromosome, start=start, end=end, strand=strand)

    def to_string(self):
        """
        Convert this :class:`~GenomicRegion` to its string representation.

        :return: str
        """
        region_string = ''
        if self.chromosome is not None:
            region_string += '%s' % self.chromosome

            if self.start is not None:
                region_string += ':%d' % self.start

                if self.end is not None:
                    region_string += '-%d' % self.end

                if self.strand is not None:
                    if self.strand == 1:
                        region_string += ':+'
                    else:
                        region_string += ':-'
        return region_string

    def __repr__(self):
        return self.to_string()

    def overlaps(self, region):
        """
        Check if this region overlaps with the specified region.

        :param region: :class:`~GenomicRegion` object or string
        """
        region = as_region(region)

        if region.chromosome != self.chromosome:
            return False

        if self.end is None or region.start is None or region.start <= self.end:
            if self.start is None or region.end is None or region.end >= self.start:
                return True
        return False

    def overlap(self, region):
        """
        Return the overlap in base pairs between this
        region and another region.

        :param region: :class:`~GenomicRegion` to find overlap for
        :return: overlap as int in base pairs
        """
        region = as_region(region)

        if region.chromosome != self.chromosome:
            return 0

        return max(0, min(self.end, region.end) - max(self.start, region.start))

    def contains(self, region):
        """
        Check if the specified region is completely contained in this region.

        :param region: :class:`~GenomicRegion` object or string
        """
        region = as_region(region)

        if region.chromosome != self.chromosome:
            return False

        if region.start >= self.start and region.end <= self.end:
            return True
        return False

    def _equals(self, region):
        """
        Is this region in the same location as another region?

        Checks if chromosome, start, end, and strand are equal
        between this object and another :class:`~GenomicRegion`

        :param region: :class:`~GenomicRegion`
        :return: True if equal, False if not
        """

        region = as_region(region)

        if region.chromosome != self.chromosome:
            return False
        if region.start != self.start:
            return False
        if region.end != self.end:
            return False
        if region.strand != self.strand:
            return False
        return True

    def is_reverse(self):
        """
        Return True if this region is on the reverse strand
        of the reference genome.

        :return: True if on '-' strand, False otherwise.
        """
        if self.strand == -1 or self.strand == '-':
            return True
        return False

    def is_forward(self):
        """
        Return True if this region is on the forward strand
        of the reference genome.

        :return: True if on '+' strand, False otherwise.
        """
        if self.strand == 1 or self.strand == '+' or self.strand == 0:
            return True
        return False

    @property
    def strand_string(self):
        """
        Return the 'strand' attribute as string.

        :return: strand as str ('+', '-', or '.')
        """
        if self.is_forward():
            return '+'
        if self.is_reverse():
            return '-'
        return '.'

    @property
    def center(self):
        """
        Return the center coordinate of the :class:`~GenomicRegion`.

        :return: float
        """
        return self.start + (self.end - self.start) / 2

    @property
    def five_prime(self):
        """
        Return the position of the 5' end of this :class:`~GenomicRegion`
        on the reference.

        :return: int
        """
        return self.end if self.is_reverse() else self.start

    @property
    def three_prime(self):
        """
        Return the position of the 3' end of this :class:`~GenomicRegion`
        on the reference.

        :return: int
        """
        return self.start if self.is_reverse() else self.end

    def copy(self):
        """
        Return a (shallow) copy of this :class:`~GenomicRegion`

        :return: :class:`~GenomicRegion`
        """
        d = {attribute: getattr(self, attribute) for attribute in self.attributes}
        return GenomicRegion(**d)

    def __eq__(self, other):
        return self._equals(other)

    def __ne__(self, other):
        return not self._equals(other)

    def __len__(self):
        return self.end - self.start

    def as_bed_line(self, score_field='score', name_field='name'):
        """
        Return a representation of this object as line in a BED file.

        :param score_field: name of the attribute to be used in the
                            'score' field of the BED line
        :param name_field: name of the attribute to be used in the
                           'name' field of the BED line
        :return: str
        """
        try:
            score = getattr(self, score_field)
        except AttributeError:
            score = '.'

        try:
            name = getattr(self, name_field)
        except AttributeError:
            name = '.'

        return "{}\t{}\t{}\t{}\t{}\t{}".format(self.chromosome, self.start, self.end,
                                               name, score, self.strand_string)

    def as_gff_line(self, source_field='source', feature_field='feature', score_field='score',
                    frame_field='frame', float_format='.2e'):
        """
        Return a representation of this object as line in a GFF file.

        :param source_field: name of the attribute to be used in the
                            'source' field of the GFF line
        :param feature_field: name of the attribute to be used in the
                              'feature' field of the GFF line
        :param score_field: name of the attribute to be used in the
                            'score' field of the GFF line
        :param frame_field: name of the attribute to be used in the
                            'frame' field of the GFF line
        :param float_format: Formatting string for the float fields

        :return: str
        """
        try:
            source = getattr(self, source_field)
        except AttributeError:
            source = '.'

        try:
            feature = getattr(self, feature_field)
        except AttributeError:
            feature = '.'

        try:
            score = "{:{float_format}}".format(getattr(self, score_field), float_format=float_format)
        except AttributeError:
            score = '.'

        try:
            frame = getattr(self, frame_field)
        except AttributeError:
            frame = '.'

        no_group_items = {'start', 'end', 'chromosome', 'source', 'feature', 'score',
                          'frame', 'ix', 'strand', 'fields'}
        group = ''
        for attribute in self.attributes:
            if attribute not in no_group_items:
                a = getattr(self, attribute)
                if isinstance(a, float):
                    a = "{:{float_format}}".format(a, float_format=float_format)
                elif isinstance(a, string_types):
                    a = '"{}"'.format(a)
                group += '{} {}; '.format(attribute, a)

        return "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}".format(self.chromosome, source,
                                                           feature, self.start + 1,
                                                           self.end, score,
                                                           self.strand_string, frame,
                                                           group)

    def expand(self,
               absolute=None, relative=None,
               absolute_left=0, absolute_right=0,
               relative_left=0.0, relative_right=0.0,
               copy=True, from_center=False):
        """
        Expand this region by a relative or an absolute amount.

        :param absolute: Absolute amount in base pairs by which to
                         expand the region represented by this
                         :class:`~GenomicRegion` object on both
                         sides. New region start will be
                         <old start - absolute>, new region end
                         will be <old end + absolute>
        :param relative: Relative amount as fraction of region by which to
                         expand the region represented by this
                         :class:`~GenomicRegion` object on both
                         sides. New region start will be
                         <old start - relative*len(self)>,
                         new region end will be
                         <old end + relative*(len(self)>
        :param absolute_left: Absolute amount in base pairs by which to
                              expand the region represented by this
                              :class:`~GenomicRegion` object on the
                              left side
        :param absolute_right: Absolute amount in base pairs by which to
                               expand the region represented by this
                               :class:`~GenomicRegion` object on the
                               right side
        :param relative_left: Relative amount in base pairs by which to
                              expand the region represented by this
                              :class:`~GenomicRegion` object on the
                              left side
        :param relative_right: Relative amount in base pairs by which to
                               expand the region represented by this
                               :class:`~GenomicRegion` object on the
                               right side
        :param copy: If True, return a copy of the original region,
                     if False will modify the existing region in place
        :param from_center: If True measures distance from center rather
                            than start and end of the old region
        :return: :class:`~GenomicRegion`
        """
        if absolute is not None:
            absolute_left, absolute_right = absolute, absolute
        if relative is not None:
            relative_left, relative_right = relative, relative

        extend_left_bp = str_to_int(absolute_left) + int(relative_left * len(self))
        extend_right_bp = str_to_int(absolute_right) + int(relative_right * len(self))

        new_region = self.copy() if copy else self
        if from_center:
            center = self.center
            new_region.start = max(0, int(center) - extend_left_bp)
            new_region.end = int(center) + extend_right_bp
        else:
            new_region.start = max(0, int(self.start) - extend_left_bp)
            new_region.end = int(self.end) + extend_right_bp
        return new_region

    def __add__(self, distance):
        new_region = self.copy()
        new_region.start += str_to_int(distance)
        new_region.end += str_to_int(distance)
        return new_region

    def __sub__(self, distance):
        return self.__add__(-str_to_int(distance))

    def fix_chromosome(self, copy=False):
        """
        Change chromosome representation from chr<NN> to <NN> or vice versa.

        :param copy: If True, make copy of region, otherwise
                     will modify existing region in place.
        :return: :class:`~GenomicRegion`
        """
        region = self.copy() if copy else self
        if region.chromosome.startswith('chr'):
            region.chromosome = region.chromosome[3:]
        else:
            region.chromosome = 'chr' + region.chromosome
        return region


class RegionBased(object):
    """
    Base class for working with genomic regions.

    Guide for inheriting classes which functions to override:

    MUST (basic functionality):
        _region_iter
        _get_regions

    SHOULD (works if above are implemented, but is highly inefficient):
        _region_subset
        _region_intervals

    CAN (override for potential speed benefits or added functionality):
        _region_len
        chromosomes
        chromosome_lengths
        region_bins
    """

    def __init__(self):
        pass

    @property
    def _estimate_region_bounds(self):
        return True

    @property
    def file_type(self):
        return 'region'

    def _region_iter(self, *args, **kwargs):
        raise NotImplementedError("Function not implemented")

    def _get_regions(self, item, *args, **kwargs):
        raise NotImplementedError("Function not implemented")

    def _add_region(self, region, *args, **kwargs):
        raise NotImplementedError("Function not implemented")

    def _region_subset(self, region, *args, **kwargs):
        return self.regions[self.region_bins(region)]

    def _region_intervals(self, region, score_field='score', *args, **kwargs):
        kwargs.setdefault('lazy', True)
        intervals = []
        for region in self.regions(region, *args, **kwargs):
            intervals.append((region.start, region.end, getattr(region, score_field)))
        return intervals

    def _region_len(self):
        return sum(1 for _ in self.regions)

    def __len__(self):
        return self._region_len()

    def __getitem__(self, item):
        return self._get_regions(item)

    def __iter__(self):
        return self.regions()

    @property
    def regions(self):
        """
        Iterate over genomic regions in this object.

        Will return a :class:`~GenomicRegion` object in every iteration.
        Can also be used to get the number of regions by calling
        len() on the object returned by this method.

        :return: RegionIter
        """

        class RegionIter(object):
            def __init__(self, region_based):
                self._region_based = region_based

            def __len__(self):
                return self._region_based._region_len()

            def __iter__(self):
                return self()

            def _fix_chromosome(self, regions):
                for r in regions:
                    r.fix_chromosome(copy=True)

            def __call__(self, key=None, *args, **kwargs):
                fix_chromosome = kwargs.pop('fix_chromosome', False)

                if key is None:
                    iterator = self._region_based._region_iter(*args, **kwargs)
                else:
                    if isinstance(key, string_types) or isinstance(key, GenomicRegion):
                        iterator = self._region_based.region_subset(key, *args, **kwargs)
                    else:
                        iterator = self._region_based._get_regions(key, *args, **kwargs)

                if fix_chromosome:
                    return self._fix_chromosome(iterator)
                else:
                    return iterator

            def __getitem__(self, item):
                if isinstance(item, string_types) or isinstance(item, GenomicRegion):
                    return self._region_based.region_subset(item)
                return self._region_based._get_regions(item)

        return RegionIter(self)

    def _convert_region(self, region):
        """
        Take any object that can be interpreted as a region and return a :class:`GenomicRegion`.

        :param region: Any object interpretable as genomic region (string, :class:`GenomicRegion`)
        :return: :class:`GenomicRegion`
        """
        if isinstance(region, string_types):
            region = GenomicRegion.from_string(region)

        if isinstance(region, GenomicRegion):
            if region.start is None and self._estimate_region_bounds:
                region.start = 1

            if region.end is None and self._estimate_region_bounds:
                chromosome_lengths = self.chromosome_lengths
                if region.chromosome in chromosome_lengths:
                    region.end = chromosome_lengths[region.chromosome]
        return region

    def chromosomes(self):
        """
        Get a list of chromosome names.
        """
        chromosomes_set = set()
        chromosomes = []
        for region in self.regions:
            if region.chromosome not in chromosomes_set:
                chromosomes_set.add(region.chromosome)
                chromosomes.append(region.chromosome)
        return chromosomes

    @property
    def chromosome_lengths(self):
        """
        Returns a dictionary of chromosomes and their length
        in bp.
        """
        chr_lens = {}
        for r in self.regions(lazy=True):
            if chr_lens.get(r.chromosome) is None:
                chr_lens[r.chromosome] = r.end
                continue
            if r.end > chr_lens[r.chromosome]:
                chr_lens[r.chromosome] = r.end
        return chr_lens

    def region_bins(self, region):
        """
        Takes a genomic region and returns a slice of the bin
        indices that are covered by the region.

        :param region: String or class:`~GenomicRegion`
                       object for which covered bins will
                       be returned.
        :return: slice
        """
        region = self._convert_region(region)

        start_ix = None
        end_ix = None
        for i, r in enumerate(self.regions):
            ix = r.ix if hasattr(r, 'ix') and r.ix is not None else i
            if not (r.chromosome == region.chromosome and r.start <= region.end and r.end >= region.start):
                continue
            if start_ix is None:
                start_ix = ix
                end_ix = ix + 1
                continue
            end_ix = ix + 1
        return slice(start_ix, end_ix)

    def find_region(self, query_regions, _regions_dict=None, _region_ends=None, _chromosomes=None):
        """
        Find the region that is at the center of a region.

        :param query_regions: Region selector string, :class:~GenomicRegion, or
                              list of the former
        :return: index (or list of indexes) of the region at the center of the
                 query region
        """
        is_single = False
        if isinstance(query_regions, string_types):
            is_single = True
            query_regions = [GenomicRegion.from_string(query_regions)]

        if isinstance(query_regions, GenomicRegion):
            is_single = True
            query_regions = [query_regions]

        if _regions_dict is None or _region_ends is None or _chromosomes is None:
            regions_dict = defaultdict(list)
            region_ends = defaultdict(list)
            chromosomes = set()

            for region in self.regions:
                regions_dict[region.chromosome].append(region)
                region_ends[region.chromosome].append(region.end)
                chromosomes.add(region.chromosome)
        else:
            regions_dict = _regions_dict
            region_ends = _region_ends
            chromosomes = _chromosomes

        hit_regions = []
        for query_region in query_regions:
            if isinstance(query_region, string_types):
                query_region = GenomicRegion.from_string(query_region)

            if query_region.chromosome not in chromosomes:
                hit_regions.append(None)
                continue

            center = query_region.start + (query_region.end - query_region.start) / 2
            ix = bisect_left(region_ends[query_region.chromosome], center)
            try:
                hit_regions.append(regions_dict[query_region.chromosome][ix])
            except IndexError:
                hit_regions.append(None)
        if is_single:
            return hit_regions[0]
        return hit_regions

    def subset(self, region, *args, **kwargs):
        return self.region_subset(region, *args, **kwargs)

    def region_subset(self, region, *args, **kwargs):
        """
        Takes a class:`~GenomicRegion` and returns all regions that
        overlap with the supplied region.

        :param region: String or class:`~GenomicRegion`
                       object for which covered bins will
                       be returned.
        """
        region = self._convert_region(region)
        return self._region_subset(region, *args, **kwargs)

    def region_intervals(self, region, bins=None, bin_size=None, smoothing_window=None,
                         nan_replacement=None, zero_to_nan=False, score_field='score',
                         *args, **kwargs):
        """
        Return equally-sized genomic intervals and associated scores.

        Use either bins or bin_size argument to control binning.

        :param region:  String or class:`~GenomicRegion`
                        object denoting the region to be binned
        :param bins: Number of bins to divide the region into
        :param bin_size: Size of each bin (alternative to bins argument)
        :param smoothing_window: Size of window (in bins) to smooth scores
                                 over
        :param nan_replacement: NaN values in the scores will be replaced
                                with this value
        :param zero_to_nan: If True, will convert bins with score 0 to NaN
        :param args: Arguments passed to _region_intervals
        :param kwargs: Keyword arguments passed to _region_intervals
        :return: iterator of tuples: (start, end, score)
        """
        region = self._convert_region(region)
        if not isinstance(region, GenomicRegion):
            raise ValueError("Region must be a GenomicRegion object or equivalent string!")

        raw_intervals = self._region_intervals(region, score_field=score_field, *args, **kwargs)

        if raw_intervals is None:
            raw_intervals = []

        if bins is None and bin_size is None:
            return raw_intervals

        if bins is not None:
            return RegionBased.bin_intervals(raw_intervals, bins,
                                             interval_range=[region.start, region.end],
                                             smoothing_window=smoothing_window,
                                             nan_replacement=nan_replacement,
                                             zero_to_nan=zero_to_nan)

        if bin_size is not None:
            return RegionBased.bin_intervals_equidistant(raw_intervals, bin_size,
                                                         interval_range=[region.start, region.end],
                                                         smoothing_window=smoothing_window,
                                                         nan_replacement=nan_replacement,
                                                         zero_to_nan=zero_to_nan)

    def intervals(self, *args, **kwargs):
        """
        Alias for region_intervals.
        """
        return self.region_intervals(*args, **kwargs)

    def binned_regions(self, region=None, bins=None, bin_size=None, smoothing_window=None,
                       nan_replacement=None, zero_to_nan=False, *args, **kwargs):
        """
        Same as region_intervals, but returns :class:`~GenomicRegion`
        objects instead of tuples.

        :param region:  String or class:`~GenomicRegion`
                        object denoting the region to be binned
        :param bins: Number of bins to divide the region into
        :param bin_size: Size of each bin (alternative to bins argument)
        :param smoothing_window: Size of window (in bins) to smooth scores
                                 over
        :param nan_replacement: NaN values in the scores will be replaced
                                with this value
        :param zero_to_nan: If True, will convert bins with score 0 to NaN
        :param args: Arguments passed to _region_intervals
        :param kwargs: Keyword arguments passed to _region_intervals
        :return: iterator of :class:`~GenomicRegion` objects
        """
        region = self._convert_region(region)
        br = []
        if region is None:
            for chromosome in self.chromosomes():
                interval_bins = self.region_intervals(chromosome, bins=bins, bin_size=bin_size,
                                                      smoothing_window=smoothing_window,
                                                      nan_replacement=nan_replacement,
                                                      zero_to_nan=zero_to_nan, *args, **kwargs)
                br += [GenomicRegion(chromosome=chromosome, start=interval_bin[0],
                                     end=interval_bin[1], score=interval_bin[2])
                       for interval_bin in interval_bins]
        else:
            interval_bins = self.region_intervals(region, bins=bins, bin_size=bin_size,
                                                  smoothing_window=smoothing_window,
                                                  nan_replacement=nan_replacement,
                                                  zero_to_nan=zero_to_nan, *args, **kwargs)
            br += [GenomicRegion(chromosome=region.chromosome, start=interval_bin[0],
                                 end=interval_bin[1], score=interval_bin[2])
                   for interval_bin in interval_bins]
        return br

    @staticmethod
    def bin_intervals(intervals, bins, interval_range=None, smoothing_window=None,
                      nan_replacement=None, zero_to_nan=False):
        """
        Bin a given set of intervals into a fixed number of bins.

        :param intervals: iterator of tuples (start, end, score)
        :param bins: Number of bins to divide the region into
        :param interval_range: Optional. Tuple (start, end) in base pairs
                               of range of interval to be binned. Useful if
                               intervals argument does not cover to exact
                               genomic range to be binned.
        :param smoothing_window: Size of window (in bins) to smooth scores
                                 over
        :param nan_replacement: NaN values in the scores will be replaced
                                with this value
        :param zero_to_nan: If True, will convert bins with score 0 to NaN
        :return: iterator of tuples: (start, end, score)
        """
        if intervals is None:
            return []

        intervals = np.array(list(intervals))

        if interval_range is None:
            try:
                interval_range = (min(intervals[:, 0]), max(intervals[:, 1]))
            except (IndexError, TypeError):
                raise ValueError("intervals cannot be None or length 0 if not providing interval_range!")

        bin_size = (interval_range[1] - interval_range[0] + 1) / bins
        logger.debug("Bin size: {}".format(bin_size))

        return RegionBased._bin_intervals_equidist(intervals, bin_size, interval_range, bins=bins,
                                                   smoothing_window=smoothing_window,
                                                   nan_replacement=nan_replacement,
                                                   zero_to_nan=zero_to_nan)

    @staticmethod
    def bin_intervals_equidistant(intervals, bin_size, interval_range=None, smoothing_window=None,
                                  nan_replacement=None, zero_to_nan=False):
        """
        Bin a given set of intervals into bins with a fixed size.

        :param intervals: iterator of tuples (start, end, score)
        :param bin_size: Size of each bin in base pairs
        :param interval_range: Optional. Tuple (start, end) in base pairs
                               of range of interval to be binned. Useful if
                               intervals argument does not cover to exact
                               genomic range to be binned.
        :param smoothing_window: Size of window (in bins) to smooth scores
                                 over
        :param nan_replacement: NaN values in the scores will be replaced
                                with this value
        :param zero_to_nan: If True, will convert bins with score 0 to NaN
        :return: iterator of tuples: (start, end, score)
        """
        if intervals is None:
            return []

        intervals = np.array(list(intervals))

        if interval_range is None:
            try:
                interval_range = (min(intervals[:, 0]), max(intervals[:, 1]))
            except (IndexError, TypeError):
                raise ValueError("intervals cannot be None or length 0 if not providing interval_range!")

        if isinstance(interval_range, GenomicRegion):
            interval_range = (interval_range.start, interval_range.end)

        return RegionBased._bin_intervals_equidist(intervals, bin_size, interval_range,
                                                   smoothing_window=smoothing_window,
                                                   nan_replacement=nan_replacement,
                                                   zero_to_nan=zero_to_nan)

    @staticmethod
    def _bin_intervals_equidist(intervals, bin_size, interval_range, bins=None, smoothing_window=None,
                                nan_replacement=None, zero_to_nan=False):
        bin_size = str_to_int(bin_size)
        if bins is None:
            bins = int((interval_range[1] - interval_range[0] + 1) / bin_size + .5)

        current_interval = 0
        bin_coordinates = []
        bin_weighted_sum = [0.0] * bins
        bin_weighted_count = [0.0] * bins
        bin_start = interval_range[0]
        last_start = intervals[0][0]
        for bin_counter in range(bins):
            bin_end = int(interval_range[0] + bin_size + (bin_size * bin_counter) + 0.5) - 1
            bin_coordinates.append((bin_start, bin_end))

            if current_interval < len(intervals):
                interval = intervals[current_interval]
                if last_start > interval[0]:
                    raise ValueError("Intervals / regions must be sorted by start coordinate for binning!")
                last_start = interval[0]
            else:
                interval = None

            # add all successive, fully-contained intervals to bin
            while interval is not None and (interval[0] <= interval[1] <= bin_end and interval[1] >= bin_start):
                value = interval[2]
                if zero_to_nan and value < 10e-8:
                    value = np.nan

                if not np.isfinite(value):
                    if nan_replacement is not None:
                        value = nan_replacement
                    else:
                        value = None

                if value is not None:
                    f = (interval[1] + 1 - interval[0]) / bin_size
                    bin_weighted_sum[bin_counter] += f * value
                    bin_weighted_count[bin_counter] += f

                current_interval += 1
                if current_interval < len(intervals):
                    interval = intervals[current_interval]
                else:
                    interval = None

            # add partially-contained interval to bin
            if interval is not None and (interval[0] <= bin_end and interval[1] >= bin_start):
                value = interval[2]
                if zero_to_nan and value < 10e-8:
                    value = np.nan

                if not np.isfinite(value):
                    if nan_replacement is not None:
                        value = nan_replacement
                    else:
                        value = None

                if value is not None:
                    f = (min(bin_end, interval[1] + 1) - max(bin_start, interval[0])) / bin_size
                    bin_weighted_sum[bin_counter] += f * value
                    bin_weighted_count[bin_counter] += f

            bin_start = bin_end + 1

        with np.errstate(divide='ignore', invalid='ignore'):
            result = np.true_divide(bin_weighted_sum, bin_weighted_count)

            if nan_replacement is not None:
                result[~ np.isfinite(result)] = nan_replacement  # -inf inf NaN
            else:
                result[~ np.isfinite(result)] = np.nan  # -inf inf NaN

        if smoothing_window is not None:
            with warnings.catch_warnings():
                warnings.simplefilter('ignore')
                result = apply_sliding_func(result, smoothing_window)

        return tuple((bin_coordinates[i][0], bin_coordinates[i][1], result[i]) for i in range(len(result)))

    def add_region(self, region, *args, **kwargs):
        """
        Add a genomic region to this object.

        This method offers some flexibility in the types of objects
        that can be loaded. See parameters for details.

        :param region: Can be a :class:`~GenomicRegion`, a str in the form
                       '<chromosome>:<start>-<end>[:<strand>], a dict with
                       at least the fields 'chromosome', 'start', and
                       'end', optionally 'ix', or a list of length 3
                       (chromosome, start, end) or 4 (ix, chromosome,
                       start, end).
        """
        ix = -1

        if isinstance(region, GenomicRegion):
            return self._add_region(region.copy(), *args, **kwargs)
        elif isinstance(region, string_types):
            return self._add_region(GenomicRegion.from_string(region), *args, **kwargs)
        elif type(region) is dict:
            return self._add_region(GenomicRegion(**copy.copy(region)), *args, **kwargs)
        else:
            try:
                offset = 0
                if len(region) == 4:
                    ix = region[0]
                    offset += 1
                chromosome = region[offset]
                start = region[offset + 1]
                end = region[offset + 2]
                strand = 1
            except TypeError:
                raise ValueError("Node parameter has to be GenomicRegion, dict, or list")

        new_region = GenomicRegion(chromosome=chromosome, start=start, end=end,
                                   strand=strand, ix=ix)
        return self._add_region(new_region, *args, **kwargs)

    @property
    def regions_dict(self):
        """
        Return a dictionary with region index as keys
        and regions as values.

        :return: dict {region.ix: region, ...}
        """
        regions_dict = dict()
        for i, r in enumerate(self.regions):
            regions_dict[getattr(r, 'ix', i)] = r
        return regions_dict

    def to_bed(self, file_name, subset=None, **kwargs):
        """
        Export regions as BED file

        :param file_name: Path of file to write regions to
        :param subset: optional :class:`~GenomicRegion` or str to
                       write only regions overlapping this region
        :param kwargs: Passed to :func:`write_bed`
        """
        write_bed(file_name, self.regions(subset, lazy=True), **kwargs)

    def to_gff(self, file_name, subset=None, **kwargs):
        """
        Export regions as GFF file

        :param file_name: Path of file to write regions to
        :param subset: optional :class:`~GenomicRegion` or str to
                       write only regions overlapping this region
        :param kwargs: Passed to :func:`write_gff`
        """
        write_gff(file_name, self.regions(subset, lazy=True), **kwargs)

    def to_bigwig(self, file_name, subset=None, **kwargs):
        """
        Export regions as BigWig file.

        :param file_name: Path of file to write regions to
        :param subset: optional :class:`~GenomicRegion` or str to
                       write only regions overlapping this region
        :param kwargs: Passed to :func:`write_bigwig`
        """
        write_bigwig(file_name, self.regions(subset, lazy=True), mode='w', **kwargs)

    def contains(self, region, full=False):
        region = as_region(region)
        hits = 0
        if full:
            for r in self.regions(region, lazy=True):
                if r.contains(region):
                    hits += 1
        else:
            for _ in self.regions(region, lazy=True):
                hits += 1

        return hits > 0

    def __contains__(self, item):
        return self.contains(item)


class RegionWrapper(RegionBased):
    """
    Provide :class:`~RegionBased` functionality to any list of regions.

    This class uses interval trees internally to provide fast region subsetting.
    On initialisation these trees will be generated, which might take some time.
    """
    def __init__(self, regions):
        super(RegionWrapper, self).__init__()

        self._regions = []
        for i, region in enumerate(regions):
            self._regions.append(region)

        self._region_trees = {}
        self._update_trees()

    def _update_trees(self):
        region_intervals = defaultdict(list)

        for i, region in enumerate(self._regions):
            interval = intervaltree.Interval(region.start - 1, region.end, data=(i, region))
            region_intervals[region.chromosome].append(interval)

        self._region_trees = {}
        for chromosome, intervals in region_intervals.items():
            self._region_trees[chromosome] = intervaltree.IntervalTree(intervals)

    def _get_regions(self, item, *args, **kwargs):
        return self._regions[item]

    def _region_iter(self, *args, **kwargs):
        for region in self._regions:
            yield region

    def _region_subset(self, region, *args, **kwargs):
        sort = kwargs.get("sort", False)
        try:
            tree = self._region_trees[region.chromosome]
            if sort:
                intervals = sorted(tree[region.start - 1:region.end], key=lambda r: r.begin)
            else:
                intervals = sorted(tree[region.start - 1:region.end], key=lambda r: r.data[0])
            for interval in intervals:
                yield interval.data[1]
        except KeyError:
            pass

    def _region_len(self):
        return len(self._regions)

    def _add_region(self, region, *args, **kwargs):
        self._regions.append(region)
        if kwargs.get('flush', True):
            self._update_trees()

    def chromosomes(self):
        return list(self._region_trees.keys())


class PbtRegion(GenomicRegion):
    def __init__(self, interval):
        self._interval = interval

    def __getattr__(self, item):
        try:
            return self._interval.attrs[item]
        except KeyError:
            return getattr(self._interval, item)

    @property
    def chromosome(self):
        return self._interval.chrom

    @property
    def strand(self):
        if self._interval.strand == '+':
            return 1
        elif self._interval.strand == '-':
            return -1
        return None

    @property
    def score(self):
        try:
            return float(self._interval.score)
        except (TypeError, ValueError):
            pass

        if len(self._interval.fields) == 4:  # likely bedGraph!
            try:
                return float(self._interval.fields[3])
            except ValueError:
                pass

        return np.nan

    @property
    def name(self):
        try:
            return self._interval.name
        except (TypeError, ValueError):
            warnings.warn("Pybedtools could not retrieve interval name. Continuing anyways.")
        return None

    @property
    def source(self):
        if self._interval.file_type == 'gff':
            return self._interval.fields[1]
        return None

    @property
    def feature(self):
        if self._interval.file_type == 'gff':
            return self._interval.fields[2]
        return None

    @property
    def frame(self):
        if self._interval.file_type == 'gff':
            return self._interval.fields[7] if len(self._interval.fields) > 7 else '.'

    @property
    def attributes(self):
        a = set(self.attrs.keys())
        a.add('chromosome')
        a.add('start')
        a.add('end')
        a.add('strand')
        a.add('score')
        a.add('name')
        if self._interval.file_type == 'gff':
            a.add('source')
            a.add('feature')
            a.add('frame')
        return list(a)

    def expand(self, *args, **kwargs):
        kwargs['copy'] = True
        return GenomicRegion.expand(self, *args, **kwargs)


class Bed(pybedtools.BedTool, RegionBased):
    """
    Data type representing a BED file.

    Extends :class:`~pybedtools.BedTool` and therefore
    provides all the methods of the original class, such
    as intersect, etc.
    """

    def __init__(self, *args, **kwargs):
        pybedtools.BedTool.__init__(self, *args, **kwargs)

    def __exit__(self, exec_type, exec_val, exec_tb):
        pass

    def __enter__(self):
        return self

    def _region_iter(self, lazy=False, *args, **kwargs):
        if lazy:
            lazy_region = PbtRegion(None)
        else:
            lazy_region = None

        for interval in self.intervals:
            yield self._interval_to_region(interval, lazy_region=lazy_region)

    def _get_regions(self, item, *args, **kwargs):
        if isinstance(item, string_types):
            item = GenomicRegion.from_string(item)

        if not isinstance(item, GenomicRegion):
            intervals = pybedtools.BedTool.__getitem__(self, item)
            if isinstance(intervals, pybedtools.Interval):
                return self._interval_to_region(intervals)
            elif isinstance(intervals, GenomicRegion):
                return intervals
            else:
                regions = []
                for interval in intervals:
                    regions.append(self._interval_to_region(interval))
                return regions

        start = item.start if item.start is not None else 1

        query_interval = pybedtools.cbedtools.Interval(chrom=item.chromosome,
                                                       start=start,
                                                       end=item.end)

        regions = []
        for interval in self.all_hits(query_interval):
            region = self._interval_to_region(interval)
            regions.append(region)
        return regions

    def _region_subset(self, region, lazy=False, *args, **kwargs):
        if lazy:
            lazy_region = PbtRegion(None)
        else:
            lazy_region = None
        for interval in self.filter(lambda i: i.chrom == region.chromosome
                                              and i.start <= region.end
                                              and i.end >= region.start):
            yield self._interval_to_region(interval, lazy_region=lazy_region)

    def _region_len(self):
        return sum(1 for _ in self.intervals)

    def _interval_to_region(self, interval, lazy_region=None):
        if lazy_region is None:
            return PbtRegion(interval)
        else:
            lazy_region._interval = interval
            return lazy_region

    def merge_overlapping(self, stat=intervals_weighted_mean, sort=True):
        """
        Merge overlapping BED intervals.

        :param stat: Function to use for scoring the merged interval.
        :param sort: Sort bed file intervals by position before merging.
        :return: iterator of merged intervals
        """
        if sort:
            bed = self
        else:
            bed = self.sort()

        current_intervals = []
        for interval in bed:
            if len(current_intervals) == 0 or (current_intervals[-1].start < interval.end and
                                               current_intervals[-1].end > interval.start and
                                               current_intervals[-1].chrom == interval.chrom):
                current_intervals.append(interval)
            else:
                # merge
                intervals = np.array([(current.start, current.end,
                                       float(current.score) if current.score != '.' else np.nan)
                                      for current in current_intervals])
                merged_score = "{:0.6f}".format(stat(intervals))
                merged_strand = current_intervals[0].strand
                merged_start = min(intervals[:, 0])
                merged_end = max(intervals[:, 1])
                merged_chrom = current_intervals[0].chrom
                merged_name = current_intervals[0].name
                merged_interval = pybedtools.Interval(merged_chrom, merged_start, merged_end, name=merged_name,
                                                      score=merged_score, strand=merged_strand)
                current_intervals = [interval]
                yield merged_interval


class Bedpe(Bed):
    """
    Represents a BEDPE file (genomic region pairs).

    Access each region of the pair with
    chromosome<1|2>
    start<1|2>
    end<1|2>
    strand<1|2>
    """
    def __init__(self, *args, **kwargs):
        Bed.__init__(self, *args, **kwargs)

    @property
    def file_type(self):
        return 'bedpe'

    def _interval_to_region(self, interval, lazy_region=None):
        fields = interval.fields

        if len(fields) < 6:
            raise ValueError("File does not appear to be BEDPE (columns: {})".format(len(fields)))

        try:
            score = float(fields[7])
        except (IndexError, TypeError, ValueError):
            score = np.nan

        try:
            name = fields[6]
        except IndexError:
            name = '.'

        try:
            strand1 = fields[8]
        except IndexError:
            strand1 = '.'

        try:
            strand2 = fields[9]
        except IndexError:
            strand2 = '.'

        if lazy_region is None:
            return GenomicRegion(chromosome=fields[0], start=int(fields[1]), end=int(fields[2]),
                                 chromosome1=fields[0], start1=int(fields[1]), end1=int(fields[2]),
                                 chromosome2=fields[3], start2=int(fields[4]), end2=int(fields[5]),
                                 strand=strand1, strand1=strand1, strand2=strand2,
                                 score=score, fields=fields,
                                 name=name)
        lazy_region.update(chromosome=fields[0], start=int(fields[1]), end=int(fields[2]),
                           chromosome1=fields[0], start1=int(fields[1]), end1=int(fields[2]),
                           chromosome2=fields[3], start2=int(fields[4]), end2=int(fields[5]),
                           strand=strand1, strand1=strand1, strand2=strand2,
                           score=score, fields=fields,
                           name=name)
        return lazy_region


class BigWig(RegionBased):
    """
    Represents a BigWig file.

    Forwards function and property calls that do not belong to RegionBased
    to :class:`~pyBigWig.BigWig`.
    """
    def __init__(self, bw):
        RegionBased.__init__(self)
        if isinstance(bw, string_types):
            bw = pyBigWig.open(bw)
        self.bw = bw
        self._intervals = None

    @property
    def file_type(self):
        return 'bw'

    def __exit__(self, exec_type, exec_val, exec_tb):
        pass

    def __enter__(self):
        return self

    def _region_iter(self, lazy=False, *args, **kwargs):
        if lazy:
            lazy_region = GenomicRegion(chromosome=None, start=0, end=0)
        else:
            lazy_region = None

        chromosome_lengths = self.chromosome_lengths
        chromosomes = self.chromosomes()
        for chromosome in chromosomes:
            for start, end, score in self.bw.intervals(chromosome, 1, chromosome_lengths[chromosome]):
                if lazy_region is None:
                    yield GenomicRegion(chromosome=chromosome, start=start + 1, end=end, score=score)
                else:
                    lazy_region.update(chromosome=chromosome, start=start + 1, end=end, score=score)
                    yield lazy_region

    def _get_regions(self, item, *args, **kwargs):
        if isinstance(item, string_types):
            item = GenomicRegion.from_string(item)

        if not isinstance(item, GenomicRegion):
            return self.bw[item]

        return self.subset(item)

    def _region_subset(self, region, lazy=False, *args, **kwargs):
        if isinstance(region, GenomicRegion):
            regions = [region]
        else:
            regions = region

        if lazy:
            lazy_region = GenomicRegion(chromosome=None, start=0, end=0)
        else:
            lazy_region = None

        for r in regions:
            for start, end, score in self.region_intervals(r):
                if lazy_region is None:
                    yield GenomicRegion(chromosome=r.chromosome, start=start, end=end, score=score)
                else:
                    lazy_region.update(chromosome=r.chromosome, start=start, end=end, score=score)
                    yield lazy_region

    def _region_intervals(self, region, *args, **kwargs):
        if self._intervals is None:
            with warnings.catch_warnings():
                warnings.simplefilter('ignore')
                try:
                    intervals = self.bw.intervals(region.chromosome, region.start, region.end)
                except RuntimeError:
                    logger.debug("Invalid interval bounds? {}".format(region))
                    raise
        else:
            intervals = self._memory_intervals(region)

        interval_list = []
        if intervals is not None:
            for interval in intervals:
                interval_list.append((interval[0] + 1, interval[1], interval[2]))
        return interval_list

    def _region_len(self):
        return sum(1 for _ in self.regions)

    def chromosomes(self):
        return natural_sort(list(self.chromosome_lengths.keys()))

    @property
    def chromosome_lengths(self):
        return self.bw.chroms()

    def __getattr__(self, name):
        try:
            func = getattr(self.__dict__['bw'], name)
            return func
        except AttributeError:
            if name == '__enter__':
                return BigWig.__enter__
            elif name == '__exit__':
                return BigWig.__exit__
            raise

    def __getitem__(self, item):
        return self._get_regions(item)

    def load_intervals_into_memory(self):
        """
        Load entire BigWig file into memory.

        May speed up interval search over slow file systems or connections.
        """
        self._intervals = dict()
        for chromosome in self.bw.chroms().keys():
            chromosome_intervals = []
            for start, end, score in self.bw.intervals(chromosome):
                interval = intervaltree.Interval(start, end, data=score)
                chromosome_intervals.append(interval)

            self._intervals[chromosome] = intervaltree.IntervalTree(chromosome_intervals)

    def _memory_intervals(self, region):
        chromosome_intervals = self._intervals[region.chromosome]
        return [(interval.begin, interval.end, interval.data)
                for interval in sorted(chromosome_intervals[region.start - 1:region.end],
                                       key=lambda r: r[0])]

    def region_stats(self, region, bins=1, stat='mean'):
        """
        BigWig.stats with region query.

        :param region: :class:`~GenomicRegion`
        :param bins: Number of bins with stats to return
        :param stat: name of statistic to use (default: mean)
        :return: interval stats
        """
        if isinstance(region, string_types):
            region = GenomicRegion.from_string(region)

        chroms = self.bw.chroms()
        r_start = region.start - 1 if region.start is not None else 0
        r_end = region.end if region.end is not None else chroms[region.chromosome]

        return self.stats(region.chromosome, r_start, r_end, type=stat, nBins=bins)

    def intervals(self, region, bins=None, bin_size=None, smoothing_window=None,
                  nan_replacement=None, zero_to_nan=False, *args, **kwargs):
        return self.region_intervals(region, bins=bins, bin_size=bin_size, smoothing_window=smoothing_window,
                                     nan_replacement=nan_replacement, zero_to_nan=zero_to_nan,
                                     *args, **kwargs)


class BedRegion(GenomicRegion):
    """
    Represents a BED file line.
    """
    def __init__(self, bed_line, ix=None):
        try:
            self.fields = bed_line.split("\t")
        except AttributeError:
            self.fields = bed_line

        self._attributes = ('chromosome', 'start', 'end', 'name', 'score', 'strand',
                            'thick_start', 'thick_end', 'item_rgb', 'block_count',
                            'block_sizes', 'block_starts')[:len(self.fields)]
        self.ix = ix

    @property
    def attributes(self):
        return self._attributes

    @property
    def chromosome(self):
        return self.fields[0]

    @property
    def start(self):
        return int(self.fields[1])

    @property
    def end(self):
        return int(self.fields[2])

    @property
    def name(self):
        try:
            return self.fields[3]
        except IndexError:
            return None

    @property
    def strand(self):
        try:
            s = self.fields[5]
            return -1 if s == '-' else 1
        except IndexError:
            return 1

    @property
    def score(self):
        try:
            return float(self.fields[4])
        except (IndexError, ValueError):
            return np.nan

    @property
    def thick_start(self):
        try:
            return int(self.fields[6])
        except IndexError:
            return None

    @property
    def thick_end(self):
        try:
            return int(self.fields[7])
        except IndexError:
            return None

    @property
    def item_rgb(self):
        try:
            return self.fields[8].split(',')
        except IndexError:
            return None

    @property
    def block_count(self):
        try:
            return int(self.fields[9])
        except IndexError:
            return None

    @property
    def block_sizes(self):
        try:
            return [int(s) for s in self.fields[10].split(',')]
        except IndexError:
            return None

    @property
    def block_starts(self):
        try:
            return [int(s) for s in self.fields[11].split(',')]
        except IndexError:
            return None


class GffRegion(GenomicRegion):
    """
    Represents a GFF file line.
    """
    def __init__(self, gff_line, ix=None):
        try:
            self.fields = gff_line.split("\t")
        except AttributeError:
            self.fields = gff_line

        self.ix = ix
        self._attribute_dict = None

    @property
    def attributes(self):
        a = ['ix', 'chromosome', 'source', 'feature', 'start', 'end',
             'score', 'strand', 'frame']

        return a + list(self.attribute_dict.keys())

    @property
    def attribute_dict(self):
        if self._attribute_dict is None:
            self._attribute_dict = dict()

            attribute_fields = re.split(";\s*", self.fields[8])
            for field in attribute_fields:
                try:
                    key, value = shlex.split(field)
                except ValueError:
                    try:
                        key, value = re.split('=', field)
                    except ValueError:
                        continue
                self._attribute_dict[key] = value
        return self._attribute_dict

    def __getattr__(self, item):
        try:
            return self.attribute_dict[item]
        except (IndexError, KeyError):
            raise AttributeError("Attribute {} cannot be found".format(item))

    @property
    def seqname(self):
        return self.fields[0]

    @property
    def source(self):
        return self.fields[1]

    @property
    def feature(self):
        return self.fields[2]

    @property
    def chromosome(self):
        return self.seqname

    @property
    def start(self):
        return int(self.fields[3])

    @property
    def end(self):
        return int(self.fields[4])

    @property
    def name(self):
        return None

    @property
    def strand(self):
        try:
            s = self.fields[6]
            return -1 if s == '-' else 1
        except IndexError:
            return 1

    @property
    def score(self):
        try:
            return float(self.fields[4])
        except (IndexError, ValueError):
            return np.nan

    @property
    def frame(self):
        return self.fields[7]


class Tabix(RegionBased):
    """
    Represents a Tabix file.

    Tabix-indexed files offer large speed improvements over
    regular BED/VCF/GFF files.
    """
    def __init__(self, file_name, preset=None):
        self._file_name = file_name
        self._file = pysam.TabixFile(file_name, parser=pysam.asTuple())
        RegionBased.__init__(self)

        self._file_type = self._get_file_extension()
        if preset is None:
            preset = self.file_type

        if isinstance(preset, string_types):
            if preset == 'gff' or preset == 'gtf':
                self._region_object = GffRegion
            elif preset == 'bed' or preset == 'bdg':
                self._region_object = BedRegion
            elif preset == 'vcf':
                self._region_object = BedRegion
            else:
                raise ValueError("Preset {} not valid".format(preset))
        else:
            self._region_object = preset

    @property
    def _estimate_region_bounds(self):
        return False

    @property
    def file_type(self):
        return self._file_type

    def _get_file_extension(self):
        fn = self._file_name
        if fn.endswith('.gz') or fn.endswith('.gzip'):
            fn = os.path.splitext(fn)[0]
        extension = os.path.splitext(fn)[1]
        return extension[1:]

    def _region_iter(self, *args, **kwargs):
        for chromosome in self.chromosomes():
            for region in self.region_subset(chromosome):
                yield region

    def _region_subset(self, region, *args, **kwargs):
        try:
            for fields in self._file.fetch(region.chromosome, region.start, region.end):
                yield self._region_object(fields)
        except ValueError:
            if region.chromosome not in self.chromosomes():
                warnings.warn('{} not in list of contigs'.format(region.chromosome))
            else:
                raise

    def chromosomes(self):
        return self._file.contigs

    @staticmethod
    def to_tabix(file_name, preset=None, _tabix_path='tabix'):
        tabix_command = [_tabix_path]
        if preset is not None:
            tabix_command += ['-p', preset]
        tabix_command += file_name


class GenomicDataFrame(pandas.DataFrame, RegionBased):
    """
    Represents :class:`~pandas.DataFrame` as RegionBased object.

    For full functionality, must contains the columns:
    chromosome
    start
    end

    """
    def __init__(self, *args, **kwargs):
        RegionBased.__init__(self)
        pandas.DataFrame.__init__(self, *args, **kwargs)

    @property
    def _estimate_region_bounds(self):
        return True

    @property
    def file_type(self):
        return "dataframe"

    def _region_iter(self, lazy=False, *args, **kwargs):
        if lazy:
            lazy_region = GenomicRegion(chromosome=None, start=0, end=0)
        else:
            lazy_region = None

        for ix, (index, row) in enumerate(self.iterrows()):
            yield self._row_to_region(row, ix=ix, lazy_region=lazy_region)

    def _region_subset(self, region, lazy=False, *args, **kwargs):
        if lazy:
            lazy_region = GenomicRegion(chromosome=None, start=0, end=0)
        else:
            lazy_region = None

        df_sub = self.query('chromosome == "{}" and start < {} and end > {}'.format(
            region.chromosome, region.end, region.start
        ))
        for index, row in df_sub.iterrows():
            yield self._row_to_region(row, ix=index, lazy_region=lazy_region)

    def _region_len(self):
        return self.shape[0]

    def _get_regions(self, item, *args, **kwargs):
        df_sub = self.iloc[item]

        # single row selected
        if isinstance(df_sub, pandas.Series):
            return self._row_to_region(df_sub)

        # multiple rows selected
        return [self._row_to_region(row, ix=index)
                for index, row in df_sub.iterrows()]

    def chromosomes(self):
        try:
            return self.chromosome.unique()
        except AttributeError:
            raise AttributeError("To support chromosome queries, DataFrame "
                                 "MUST contain 'chromosome' column")

    def _sub_rows(self, region):
        if isinstance(region, string_types):
            region = GenomicRegion.from_string(region)

        if isinstance(region, GenomicRegion):
            regions = [region]
        else:
            regions = region

        for r in regions:
            if isinstance(r, string_types):
                r = GenomicRegion.from_string(r)

            query = ''
            if r.chromosome is not None:
                query += 'chromosome == "{}" and '.format(r.chromosome)
            if r.start is not None:
                query += 'start >= {} and '.format(r.start)
            if r.end is not None:
                query += 'end <= {} and '.format(r.end)
            query = query[:-5]

            sub_df = self.query(query)
            for ix, row in sub_df.iterrows():
                yield ix, row

    def _row_to_region(self, row, ix=None, lazy_region=None):
        attributes = {'ix': ix}
        for key, value in row.items():
            attributes[key] = value
        if lazy_region is None:
            return GenomicRegion(**attributes)
        lazy_region.update(**attributes)
        return lazy_region

    @classmethod
    def read_table(cls, file_name, **kwargs):
        return cls(pandas.read_table(file_name, **kwargs))
