'''
Created on Jul 25, 2012

@author: gabrielp
'''
import unittest
from peak_caller.src.call_peak import *
import peaks
import scipy
from numpy import *
import scipy
from scipy import interpolate
class Test(unittest.TestCase):


    """
    
    Performs unit tests on get_FDR_cutoff_mode function
    Function is currently deperacated, tests not done
    
    """
    def test_get_FDR_cutoff_mode(self):
        pass
    
    """
    
    Performs unit tests on get_FDR_cutoff_mean function
    
    Difficult to test because of random sampling
    
    TODO: Figure out how to manually calculate FDR cutoff, right now 
    I just took the old result and am using that.  Math appears to be correct,
    but this is a really bad practice
    
    """
    def test_get_FDR_cutoff_mean(self):
        
        
        #Case: Not enough reads, expected result: returns the passed min cutoff 
        result = get_FDR_cutoff_mean([], 100)
        assert result == 2
        
        #Case: Not enough reads with different min cutoff value
        
        result = get_FDR_cutoff_mean([], 100, mincut = 10)
        assert result == 10
        
        #Case, enough reads, mean is higher than minimum cutoff, expected result: mean is returned
        #setup, create read lengths 
        
        #mean can also be estimated by emperically (N * L) / G
        #N is number of reads (20)
        #L is number of nucleotides per read (30)
        #G is genome size (in this case an RNA (100)
        #Read depth should be 6
        
        read_lengths = [30] * 100
        result = get_FDR_cutoff_mean(read_lengths, 1000)
        self.assertEqual(result, 7)
        
        #Second similar case 
        read_lengths = [30] * 20
        result = get_FDR_cutoff_mean(read_lengths, 100)
        assert result == 11 
       
        
        #Case: enough reads, mean is lower than minimum cutoff, expected result: minimum cutoff is returned
        result = get_FDR_cutoff_mean(read_lengths, 100, mincut = 20)
        assert result == 20
        
        """
    
    Quality control function, not important to main running of program
    Not tested
    
    """
    def test_plotSpline(self):
        pass
        #plotSections([5] * 10, ["1|5", "7|9"], 3)
        #assert 1 == 0
    
    """
    
    tests plotSections function, appears to have computer specific issues
    
    """
    def test_plotSections(self):
        """

        Plots each section individually, I think
        Wiggle is a list representing a wiggle track
        sections is a list of strings of format "start|stop" where start and stop are both integers
        threshold is an integer 

        """     
        pass
            
    """
    
    Performs unit testing on find_univariateSpline
    As this is mostly a wrapper for a scipy function I will not test spline calling
    beyond basics.
    
    These tests will verify that the resid logic works, and eventually these 
    will be factored into two functions.  Its bad form to have to toggle your return
    value in the function
    
    """
    def test_find_univariateSpline(self):
        
        #Case null inputs, expected:  Everything goes to hell, not testing
        
        #Case resid is false, expected: returns univariateSpline spline, just verifies that it is the same as 
        #the scipy result
        
        #setup 
        x1 = range(10)
        x2 = range(10)
        x2.reverse()
        xvals = range(20)
        data = x1 + x2
        smoothing = 5 
        #expected
        expected = scipy.interpolate.UnivariateSpline(xvals, data, k=3, s=smoothing)
        
        #test
        result = find_univariateSpline(smoothing, xvals, data, 3, resid=False)
        
        #hacky test, but they should be about the same
        self.assertAlmostEqual(expected.get_residual(), result.get_residual()) 
        
        #Case resid is true. Expected: returns residual sum of squared error between the 
        #spline and the actual curve 
        #TODO: write better tests, this only tests very basic case
        residual = find_univariateSpline(smoothing, xvals, data, 3, resid=True)
        self.assertAlmostEqual(expected.get_residual(), result.get_residual()) 
        
        #this works because the number of turns is 1 so the residual * turns is the residual
        self.assertAlmostEqual(residual, sqrt(expected.get_residual()))
        
        #Case: calculation of univarate spline breaks
        #There should be more error handling this is a general catch, we should be more sensetavie
        #to different types of errors
        assert Inf == find_univariateSpline(None, None, None, None, resid=True)
        
    """
    
    Performs unit testing on poissonP
    
    Will not test math behind calling poisson fuction more than nessessary as it has already been tested by scipy
    
    """ 
    def test_poissonP(self):
        #Case: fewer than 3 reads fall within a peak region: Expected result should return as though there are 3 expected reads
        result = poissonP(50, 3, 50, 2)
        self.assertAlmostEqual(result, (1 - 0.64723188878223115)) #manually calculated on personal computer stats.poisson.cdf(3, 3)
        
        #Case: More than 3 reads fall within a peak region. Expected: Result should return as normal
        result = poissonP(50, 3, 50, 4)
        self.assertAlmostEqual(result, (1 - 0.4334701203667089)) #manually calculated stats.poisson.cdf(3, 4)
    
        #Case: poissonp throws an exception. Expected: Result should return 1
        result = poissonP(None, None, None, None)    
        assert 1 == result
    
    """
    
    Performs unit testing on call_peaks
    
    will not test peaks_from_info here, just the error handling of call peaks
    Need to create a dummy call peaks function or stub or something to keep everything from
    getting called
    
    The way this is currently written it is difficult to test logic other than error checking
    
    """ 
    
    def test_call_peaks(self):
        #peaks_from_info = poissonP
        """call_peaks("chr1|foo|5|10|-", 
                   50, 
                   bam_fileobj=None, 
                   bam_file=None)"""
        assert 1 == 1
        #Case: Bam file is passed
        
        #Case: Bam object is passed
        
        #Case: Both bam file and object are passed
    
    """
    
    Performs unit testing on peaks_from_info
    
    Testing of this is currently handled in the allup test
    This method is currently to complicated to test as is
    
    """
    def test_peaks_from_info(self):
        pass

    """
    
    tests generating start and stops function, this is kind of a badly written function,
    need better tests and edge cases 
    
    """
    def test_generate_starts_and_stops(self):
        cutoff = 71.5879289615 
        xvals = array([ 0,  1,  2,  3,  4,  5,  6,  7,  8,  9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22,
                  23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45,
                  46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68])
        data = array([2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0]) 
        degree = 3 
        weights = None
        threshold = 3
        sect_length = 69 #might not be nessessary
        spline = find_univariateSpline(cutoff, xvals, data, degree, weights, resid=False)

        starts_and_stops, starts, stops = generate_starts_and_stops(threshold, sect_length, xvals, spline)
        
if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()