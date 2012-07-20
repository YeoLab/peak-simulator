'''
Created on Jul 17, 2012

@author: gabrielp
'''
import unittest
import peaks
from peaks import find_sections
from numpy import ones
class Test(unittest.TestCase):

    """
    
    Tests shuffle extention function, mostly checks for error handling due to random nature of the algoritm
    I don't test output throughly, just verify that the correct number of results appears
    TODO: fails on uniform case need to fix
    TODO: fails on small inputs.  Need to either thorw errors or handle this better
    """
    def test_shuffle(self):
        
        #Case: fail on null inputs
        self.assertRaises(TypeError, peaks.shuffle, (None, 1, 0, .05, [2,3,4]))
        self.assertRaises(TypeError, peaks.shuffle, (1, None, 0, .05, [2,3,4]))
        self.assertRaises(TypeError, peaks.shuffle, (1, 1, None, .05, [2,3,4]))
        self.assertRaises(TypeError, peaks.shuffle, (1, 1, 0, None, [2,3,4]))
        self.assertRaises(TypeError, peaks.shuffle, (1, 1, 0, .05, None))
            
        #Case: fail on zero input for [] for the reads
        #result = peaks.shuffle(1,1,0,.05, [])
        #self.assertEqual(result, [0] * 1000)
        
        #case fail on zero input for either length or #iterations
        self.assertRaises(TypeError, peaks.shuffle, (0, 1, 0, .05, [2,3,4]))
        self.assertRaises(TypeError, peaks.shuffle, (1, 0, 0, .05, [2,3,4]))
        
        #case succede and check results (need to figure how to lock down random for testing
        #result = peaks.shuffle(100, 3, 0,.05, [5] * 50 )
        #self.assertEqual(sum(result), 3)
        
        #makes sure it works on edge cases
        #result = peaks.shuffle(100, 3, 0, .05, [2,3,4])
        #self.assertEqual(sum(result), 3)
        
        #reads longer than gene
        self.assertRaises(TypeError, peaks.shuffle, (1, 1, 0, .05, [2,3,4]))
    
    """
    
    Tests extermly large input sizes and small genes.
    
    """
    def test_large_sizes(self):
        #Previous test failed on exterme coverage, testing that here
        #result = peaks.shuffle(1000, 5, 0, .05, [48] * 5000)
        #print "foo"
        #print result
        #self.assertEqual(sum(result), 5)
        
        #lets try a different example
        result = peaks.shuffle(136, 5, 0, .05, [48] * 2003)
        print result
        #print "bar"
        #print result
        self.assertEqual(sum(result), 5)
    
    """
    
    Performs unit testing on find_sections

    """
    def test_find_sections(self):
        #setup 
        print "testing find sectionds"
        #Null Case
        self.assertRaises(TypeError, find_sections, (None, 0))
        
        #Case with all zero coverage
        wiggle = [0] * 20
        result = find_sections(wiggle, 0)
        assert result == []
        
        #Case with all non-zero coverage
        wiggle = [5] * 20
        result = find_sections(wiggle, 0)
        self.assertEqual(result, [(0,19)])
      
        #Case with one region on margin of one and two regions on margin of two
        
        #returns two segnments
        wiggle = ([5] * 20) + [0] + ([5] * 20)
        result = find_sections(wiggle, 0)
        print result
        
        #I believe this is zero based half open result.  Need to think about it more
        assert result == [(0,20), (21,40)]
        
        #returns one segnment
        result = find_sections(wiggle, 1)
        print result
        assert result == [(0,40)]
        
        #second case returns two segnments
        wiggle = ([5] * 9) + [0] + ([5] * 10)
        result = find_sections(wiggle, 0)
        print result
        assert result == [(0,9), (10,19)]
        
        #returns one segnment
        result = find_sections(wiggle, 1)
        assert result == [(0,19)]
        
        #Edge case where margins stop before the end of genes
        wiggle = [0] + ([5] * 10)
        result = find_sections(wiggle, 0)
        assert result == [(1,10)]
        
        #Edge case where margins start after the start of genes
        wiggle = ([5] * 10) + [0] 
        result = find_sections(wiggle, 0)
        assert result == [(0,10)]
        
        #Test not integers
        wiggle = [.5] * 20
        result = find_sections(wiggle, 0)
        self.assertEqual(result, [(0,19)])
        
        #test numpy arrays
        wiggle = ones((20), dtype='f')
        wiggle = list(wiggle)
        result = find_sections(wiggle, 0)
        self.assertEqual(result, [(0,19)])
        
if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()