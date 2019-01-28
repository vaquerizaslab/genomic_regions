import pytest
import os
from genomic_regions.files import write_bed, write_bigwig, which
from genomic_regions import *


class TestGenomicRegion:
    def test_from_string(self):
        region1 = GenomicRegion.from_string('chr1')
        assert region1.chromosome == 'chr1'
        assert region1.start is None
        assert region1.end is None
        assert region1.strand is None

        region2 = GenomicRegion.from_string('chr1:0')
        assert region2.chromosome == 'chr1'
        assert region2.start == 0
        assert region2.end == 0
        assert region2.strand is None

        region3 = GenomicRegion.from_string('chr1:0-4956')
        assert region3.chromosome == 'chr1'
        assert region3.start == 0
        assert region3.end == 4956
        assert region3.strand is None

        region4 = GenomicRegion.from_string('chr1:0-4956:-')
        assert region4.chromosome == 'chr1'
        assert region4.start == 0
        assert region4.end == 4956
        assert region4.strand == -1

        region5 = GenomicRegion.from_string('chr1:0-4956:+1')
        assert region5.chromosome == 'chr1'
        assert region5.start == 0
        assert region5.end == 4956
        assert region5.strand == 1

        with pytest.raises(ValueError):
            # invalid start
            GenomicRegion.from_string('chr1:x-4956:-')
        with pytest.raises(ValueError):
            # too many fields
            GenomicRegion.from_string('chr1:0:4956:-')
        with pytest.raises(ValueError):
            # invalid strand
            GenomicRegion.from_string('chr1:0-4956:0')


class RegionBasedTestFactory:
    def setup_method(self, method):
        self.regions = None
        self.empty_regions = None

    def test_get_item(self):
        region = self.regions.regions[0]
        assert isinstance(region, GenomicRegion)
        assert region.chromosome == 'chr1'
        assert region.start == 1
        assert region.end == 1000
        assert region.strand is None

    def test_len(self):
        assert len(self.regions.regions) == 29

    @pytest.mark.parametrize("lazy", [True, False])
    def test_iter(self, lazy):
        region_iter = self.regions.regions(lazy=lazy)

        for i, region in enumerate(region_iter):
            start = 1 + i * 1000
            chromosome = 'chr1'
            if i > 22:
                start -= 23000
                chromosome = 'chr3'
            elif i > 8:
                start -= 9000
                chromosome = 'chr2'

            assert region.chromosome == chromosome
            assert region.start == start

    @pytest.mark.parametrize("lazy", [True, False])
    def test_region_subset(self, lazy):
        region_iter = self.regions.regions('chr1', lazy=lazy)

        for i, region in enumerate(region_iter):
            start = 1 + i * 1000

            assert region.chromosome == 'chr1'
            assert region.start == start

    def test_region_bins(self):
        bins = self.regions.region_bins(GenomicRegion(chromosome='chr1', start=3400, end=8100))
        assert bins.start == 3
        assert bins.stop == 9

        bins = self.regions.region_bins('chr2:1-5000')
        assert bins.start == 9
        assert bins.stop == 14

        bins = self.regions.region_bins('chr2:1-5001')
        assert bins.start == 9
        assert bins.stop == 15

    @pytest.mark.parametrize("lazy", [True, False])
    def test_subset(self, lazy):
        # this is essentially the same as region_bins
        intersect = self.regions.subset(GenomicRegion(chromosome='chr1', start=3400, end=8100),
                                        lazy=lazy)
        assert len(list(intersect)) == 6


class TestRegionWrapper(RegionBasedTestFactory):
    def setup_method(self, method):
        chromosomes = [
            {'name': 'chr1', 'end': 10000},
            {'name': 'chr2', 'end': 15000},
            {'name': 'chr3', 'end': 7000}
        ]

        regions = []
        for chromosome in chromosomes:
            for start in range(1, chromosome["end"] - 1000, 1000):
                regions.append(GenomicRegion(start=start, end=start + 999, chromosome=chromosome["name"]))
        self.regions = RegionWrapper(regions)
        self.empty_regions = RegionWrapper([])

    def test_add_region(self):
        # GenomicRegion
        self.empty_regions.add_region(GenomicRegion(start=1, end=1000, chromosome='chr1'))
        assert self.empty_regions[0].start == 1
        assert self.empty_regions[0].end == 1000
        assert self.empty_regions[0].chromosome == 'chr1'

        # dict
        self.empty_regions.add_region({'start': 1001, 'end': 2000, 'chromosome': 'chr1'})
        assert self.empty_regions[1].start == 1001
        assert self.empty_regions[1].end == 2000
        assert self.empty_regions[1].chromosome == 'chr1'

        # list
        self.empty_regions.add_region(['chr1', 2001, 3000])
        assert self.empty_regions[2].start == 2001
        assert self.empty_regions[2].end == 3000
        assert self.empty_regions[2].chromosome == 'chr1'


class TestBed(RegionBasedTestFactory):
    @pytest.fixture(autouse=True)
    def setup_method(self, tmpdir):
        chromosomes = [
            {'name': 'chr1', 'end': 10000},
            {'name': 'chr2', 'end': 15000},
            {'name': 'chr3', 'end': 7000}
        ]

        regions = []
        for chromosome in chromosomes:
            for start in range(1, chromosome["end"] - 1000, 1000):
                regions.append(GenomicRegion(start=start, end=start + 999, chromosome=chromosome["name"]))

        bed_file = os.path.join(str(tmpdir), 'test.bed')
        write_bed(bed_file, regions)

        self.regions = Bed(bed_file)


class TestBigWig(RegionBasedTestFactory):
    @pytest.fixture(autouse=True)
    def setup_method(self, tmpdir):
        chromosomes = [
            {'name': 'chr1', 'end': 10000},
            {'name': 'chr2', 'end': 15000},
            {'name': 'chr3', 'end': 7000}
        ]

        regions = []
        for chromosome in chromosomes:
            for start in range(1, chromosome["end"] - 1000, 1000):
                regions.append(GenomicRegion(start=start, end=start + 999, chromosome=chromosome["name"]))

        bw_file = os.path.join(str(tmpdir), 'test.bw')
        write_bigwig(bw_file, regions)

        self.regions = BigWig(bw_file)

    def test_get_item(self):
        pass


class TestBigWigMemory(RegionBasedTestFactory):
    @pytest.fixture(autouse=True)
    def setup_method(self, tmpdir):
        chromosomes = [
            {'name': 'chr1', 'end': 10000},
            {'name': 'chr2', 'end': 15000},
            {'name': 'chr3', 'end': 7000}
        ]

        regions = []
        for chromosome in chromosomes:
            for start in range(1, chromosome["end"] - 1000, 1000):
                regions.append(GenomicRegion(start=start, end=start + 999, chromosome=chromosome["name"]))

        bw_file = os.path.join(str(tmpdir), 'test.bw')
        write_bigwig(bw_file, regions)

        self.regions = BigWig(bw_file)
        self.regions.load_intervals_into_memory()

    def test_get_item(self):
        pass


@pytest.mark.skipif(which('tabix') is None, reason="Cannot find tabix in PATH")
class TestTabix(RegionBasedTestFactory):
    @pytest.fixture(autouse=True)
    def setup_method(self, tmpdir):
        chromosomes = [
            {'name': 'chr1', 'end': 10000},
            {'name': 'chr2', 'end': 15000},
            {'name': 'chr3', 'end': 7000}
        ]

        regions = []
        for chromosome in chromosomes:
            for start in range(1, chromosome["end"] - 1000, 1000):
                regions.append(GenomicRegion(start=start, end=start + 999, chromosome=chromosome["name"]))

        bed_file = os.path.join(str(tmpdir), 'test.bed')
        write_bed(bed_file, regions)
        bed = Bed(bed_file)
        bed.tabix()

        tabix_file = os.path.join(str(tmpdir), 'test.bed.gz')

        self.regions = Tabix(tabix_file)

    def test_get_item(self):
        pass


class TestGenomicDataFrame(RegionBasedTestFactory):
    @pytest.fixture(autouse=True)
    def setup_method(self, tmpdir):
        chromosomes = [
            {'name': 'chr1', 'end': 10000},
            {'name': 'chr2', 'end': 15000},
            {'name': 'chr3', 'end': 7000}
        ]

        regions = []
        for chromosome in chromosomes:
            for start in range(1, chromosome["end"] - 1000, 1000):
                regions.append(GenomicRegion(start=start, end=start + 999, chromosome=chromosome["name"]))

        chromosome_list = [r.chromosome for r in regions]
        start_list = [r.start for r in regions]
        end_list = [r.end for r in regions]

        data = {
            'chromosome': chromosome_list,
            'start': start_list,
            'end': end_list,
        }

        self.regions = GenomicDataFrame(data)
