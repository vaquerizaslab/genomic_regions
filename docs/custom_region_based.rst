.. _custom-region-based:

#############################
Custom RegionBased subclasses
#############################

Are you working with a file format that is not natively supported by
:code:`genomic_regions`? This guide will help you subclass :class:`~RegionBased`
yourself, so you can make use of the :class:`~RegionBased` functionality using
any data format.


To subclass :class:`~RegionBased`, you need to override a couple of
methods that form the basis of all other :class:`~RegionBased` methods:

* :func:`__init__`: Use this to make data-type specific initialisations
* :func:`_region_iter`: Provides basic functionality to iterate over regions
* :func:`_get_regions`: Provides basic functionality to specifically select regions

With the above methods you get all basic :class:`~RegionBased` functionality,
but for additional speed benefits you should also override:

* :func:`_region_subset`: Speeds up region selection by interval
* :func:`_region_len`: Return the number of regions in the object

In addition, you may override any of the other methods to speed them up, such
as the chromosome list :code:`chromosomes`.

In the following, we will use the simple :class:`~RegionWrapper` as implemented
in this module for illustration:


********
__init__
********

:class:`~RegionWrapper` uses a simple list to store :class:`~GenomicRegion`
objects, and interval trees from the :code:`intervaltree` module
to allow the region subsetting.

.. code:: python

    from genomic_regions import RegionBased
    from collections import defaultdict
    import intervaltree

These are set up in the :code:`__init__` method. Each chromosome gets a
separate interval tree, which is stored in a :code:`dict`:

.. code:: python

    class RegionWrapper(RegionBased):
        def __init__(self, regions):
            super(RegionWrapper, self).__init__()

            region_intervals = defaultdict(list)  # temporary variable used to hold intervals

            self._regions = []  # internal list of regions
            for i, region in enumerate(regions):
                self._regions.append(region)

                # in the "data" argument, we store both the
                # region and its original position in the list
                interval = intervaltree.Interval(region.start - 1, region.end, data=(i, region))
                region_intervals[region.chromosome].append(interval)

            self.region_trees = {}
            for chromosome, intervals in region_intervals.items():
                self.region_trees[chromosome] = intervaltree.IntervalTree(intervals)

************
_region_iter
************

To iterate over the regions, we simply iterate over the regions list.
:code:`_region_iter` should return an iterator:

.. code:: python

    def _region_iter(self, *args, **kwargs):
        for region in self._regions:
            yield region

************
_get_regions
************

To select specific regions, we can also use basic list subsetting:

.. code:: python

    def _get_regions(self, item, *args, **kwargs):
        return self._regions[item]

**************
_region_subset
**************

Due to the use of interval trees, region subsetting is also not very complicated:

.. code:: python

    def _region_subset(self, region, *args, **kwargs):
        # select the intervaltree by chromosome
        tree = self.region_trees[region.chromosome]

        # we sort by the region position int he list here, because that information
        # is lost in intervaltree
        intervals = sorted(tree[region.start:region.end], key=lambda r: r.data[0])
        for interval in intervals:
            yield interval.data[1]  # iterate over the overlapping regions

**************
_other methods
**************

Finally, we override two additional methods: :code:`_region_len` and :code:`chromosomes`.
Both of these would normally be calculated by iterating over all regions to obtain
the necessary information, but we can speed this up greatly by relying on the
internal data structure we chose for :class:`~RegionWrapper`.

.. code:: python

    def _region_len(self):
        return len(self._regions)

    def chromosomes(self):
        return list(self.region_trees.keys())

And that is all you need to subclass