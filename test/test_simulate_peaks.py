'''
Created on Aug 6, 2012

@author: gabrielp
'''
import unittest
from src.simulate_peaks import *
from numpy import * #array, ones, mean, zeros, concatenate
from numpy.testing import *
import pkg_resources
import pybedtools
from collections import namedtuple
import filecmp

class Test(unittest.TestCase):

    """
    
    Tests simulate peaks module
    
    """
    
    def test_Weighted_interval_add_peaks(self):
        
        """
        
        Tests add peaks function in weighted interval
        
        """
        
        interval = Weighted_interval(pybedtools.Interval("chr1", 10, 100))
        
        result = interval.add_peak((90, 110))
        self.assertFalse(result, "failled overlapping edge add at end")
        
        result = interval.add_peak((5, 60))
        self.assertFalse(result, "failled overlapping edge add at start")
        
        result = interval.add_peak((5, 110))
        self.assertFalse(result, "failled overlapping edge add at both ends")
        
        result = interval.add_peak((15, 25))
        self.assertTrue(result, "failled standard adding function")
        
        result = interval.add_peak((20,30))
        self.assertFalse(result, "failled adding over another peak")
        

    def test_create_genome(self):
        
        """
        
        Tests the create genome function
        
        """
        
        #Basic test of function
        handle = pkg_resources.resource_filename(__name__, "../test/test.bed")

        true_intervals = pybedtools.BedTool(handle)
        result = create_genome(handle)
        
        for true_interval, test_interval in zip(true_intervals, result):
            assert len(true_interval) == len(test_interval.weights)
            
        #handle null input
        self.assertRaises(TypeError, create_genome, None)
        
        #handle input with nothing in it
        handle = pkg_resources.resource_filename(__name__, "../test/empty.bed")
        result = create_genome(handle)
        assert result == []
        
        #handle = pkg_resources.resource_filename(__name__, "../test/not_bed.bed")
        #result = create_genome(handle)

    def test_assign_random_peaks(self):
        genome = [Weighted_interval(pybedtools.Interval("chr1", 1, 100000)),
                   Weighted_interval(pybedtools.Interval("chr1", 1, 100000))]
        
        peak_size = 50
        num_peaks = 10
        result = assign_random_peaks(genome, peak_size, num_peaks)
        
        
        total_peaks = 0
        for location in result:
            if location.peaks is None:
                continue
            
            print location.peaks
            total_peaks += len(location.peaks)
            
            #test peak size
            for start, stop in location.peaks:
                self.assertGreater(start, 0, "start is less than zero")
                self.assertGreater(stop, 0, "stop is less than zero")
                assert peak_size == stop - start
        
        #test proper number of peaks
        self.assertEqual(total_peaks, num_peaks)
    
    def test_assign_random_peaks_loopy(self):
        
        """
        
        Tests to make sure genomes of too small size still work.  
        Simply prints out less than the requested number of peaks
        
        No test, just make sure we don't go into an infanate loop here
        
        """
        genome = [Weighted_interval(pybedtools.Interval("chr1", 1, 100)),
                   Weighted_interval(pybedtools.Interval("chr2", 1, 100))]
        
        peak_size = 50
        num_peaks = 10
        result = assign_random_peaks(genome, peak_size, num_peaks)
        
        
        total_peaks = 0
        for location in result:
            if location.peaks is None:
                continue
            
            print location.peaks
            total_peaks += len(location.peaks)
            
            #test peak size
            for start, stop in location.peaks:
                self.assertGreater(start, 0, "start is less than zero")
                self.assertGreater(stop, 0, "stop is less than zero")
                assert peak_size == stop - start
        
        
    def test_assign_static_peaks(self):
        
        """
        
        Tests assigning static peaks from a file
        
        """
        
        genome = [Weighted_interval(pybedtools.Interval("chr1", 1, 100000)),
                   Weighted_interval(pybedtools.Interval("chr2", 1, 100000))]
        
        #tests error mode
        peaks = pkg_resources.resource_filename(__name__, "../test/peaks_test_err.bed")
        self.assertRaises(ValueError, assign_static_peaks, genome, peaks)
        
        #tests correect mode
        genome = [Weighted_interval(pybedtools.Interval("chr1", 1, 100000)),
                   Weighted_interval(pybedtools.Interval("chr2", 1, 100000))]
        peaks = pkg_resources.resource_filename(__name__, "../test/peaks_test.bed")
        genome = assign_static_peaks(genome, peaks)
        self.assertSetEqual(genome[0].peaks, set([(50,100)]))
        self.assertSetEqual(genome[1].peaks, set([(50,100)]))
        
        #should test value error 
    def test_distribute_background_weights(self):
        
        """
        
        Tests distrbiute background weights, currently a gama distribution, so the test 
        is to verify that the mean is the mean of a gama
        
        """
        
        genome = [Weighted_interval(pybedtools.Interval("chr1", 1, 100000)),
                   Weighted_interval(pybedtools.Interval("chr1", 1, 100000))]
                
        result = distribute_background_weights(genome, 5, 5)
        
        for interval in result:
            self.assertAlmostEqual(mean(interval.weights), 5 * 5, delta=1)
            
        #should also test for null input, but we'll skip that for now
    
    def test_distribute_peak_weights(self):
        
        """
        
        tests distrbiute peaks weights, does this by assuming a uinformly 1 background distriubtion and verifying that
        the attachment model increases the weights up to about 5
        TODO: increase should be exactly 5, need to figure out whats going on
        
        """
        genome = [Weighted_interval(pybedtools.Interval("chr1", 1, 100000)),
                  Weighted_interval(pybedtools.Interval("chr1", 1, 100000))]
        
        for interval in genome:
            interval.weights = ones(len(interval))
            interval.peaks = set([(1, 100), (200, 300)])
            
        enrichment_coeff = 5
        genome = distribute_peak_weights(genome, enrichment_coeff)
        
        
        #this is hacky, but for some reason numpy doesn't work
        total_count = 0
        n_items = 0 
        for interval in genome:
            for start, stop in interval.peaks:
                total_count += sum(interval.weights[start:stop])
                n_items += stop - start
        
        
        #this should work because the background is uniformly one in this test
        self.assertAlmostEqual(total_count / n_items, enrichment_coeff, delta=.5) 
        
    def test_distribute_reads(self):
        
        """
        
        Tests distribute reads function
        
        """
        
        genome = [Weighted_interval(pybedtools.Interval("chr1", 1, 100000))]
        
        
        genome[0].weights = array([1, 1, 1, 1])
        num_reads = 4
        result = distribute_reads(genome, num_reads)
        assert_array_equal(result[0].start_sites, array([1, 1, 1, 1]))
        
        genome[0].weights = array([2, 2, 2, 2])
        num_reads = 4
        result = distribute_reads(genome, num_reads)
        assert_array_equal(result[0].start_sites, array([1, 1, 1, 1]))
        
        genome[0].weights = array([2, 2, 0, 0])
        num_reads = 4
        result = distribute_reads(genome, num_reads)
        assert_array_equal(result[0].start_sites, array([2, 2, 0, 0]))
        
        genome[0].weights = array([1, 2, 0, 0])
        num_reads = 4
        result = distribute_reads(genome, num_reads)
        assert_array_equal(result[0].start_sites, array([1, 3, 0, 0]))
        
        #This test sort of breaks, need to assign exactly the number of reads
        #Need to figure out how to assign exactly correct numner of reads
        genome[0].weights = array([1, 2, 1, 1])
        num_reads = 4
        result = distribute_reads(genome, num_reads)
        assert_array_equal(result[0].start_sites, array([1, 2, 1, 1]))
        
        #test two regions
        genome = [Weighted_interval(pybedtools.Interval("chr1", 1, 100000)),
                  Weighted_interval(pybedtools.Interval("chr1", 1, 100000))]
        genome[0].weights = array([1, 1, 1, 1])
        genome[1].weights = array([1, 1, 1, 1])
        num_reads = 8
        result = distribute_reads(genome, num_reads)
        assert_array_equal(result[0].start_sites, array([1, 1, 1, 1]))
        assert_array_equal(result[1].start_sites, array([1, 1, 1, 1]))
        
    def test_output_bam(self):
        
        """
        
        tests output function, a bit difficult to test at the moment because I'm not sure what the 
        output should look like
        
        """
        
        genome = [Weighted_interval(pybedtools.Interval("chr1", 1, 100000)),
                  Weighted_interval(pybedtools.Interval("chr2", 1, 100000))]
        
        genome[0].start_sites = ones(4)
        genome[1].start_sites = ones(4)
        header_file = pkg_resources.resource_filename(__name__, "../test/header_test.txt")

        read_length = 50
        output_name = "foo.bam"
        output_bam(genome, header_file, read_length, output_name)
        true = pkg_resources.resource_filename(__name__, "../test/bam_test.bam")
        assert filecmp.cmp('foo.bam', true )

        
        
    def test_output_bed(self):
        
        """
        
        Tests that a bed file is outputed
        
        """
        
        genome = [Weighted_interval(pybedtools.Interval("chr1", 1, 100000)),
                  Weighted_interval(pybedtools.Interval("chr2", 1, 100000))]
        
        for interval in genome:
            interval.weights = ones(len(interval))
            interval.peaks = set([(1, 100), (200, 300)])

        output_name = "foo.bed"
        output_bed(genome, output_name)
        
        true = open(pkg_resources.resource_filename(__name__, "../test/bed_test.bed"))
        test = open("foo.bed")
        for true_line, test_line in zip(true, test):
            assert true_line == test_line
    
    def test_make_header(self):
        
        """
        
        Tests the make header function
        
        """
        
        header_file = pkg_resources.resource_filename(__name__, "../test/header_test.txt")
        
        header = make_header(header_file)
        
        self.assertDictEqual(header, {'HD': {'VN': '1.0'},
                                      'SQ' : [{'SN' : 'chr1', 'LN' : 50},
                                              {'SN' : 'chr2', 'LN' : 100}]
                                      }) 
        
if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()
    