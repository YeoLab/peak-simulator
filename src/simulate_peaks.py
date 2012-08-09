'''
Created on Aug 2, 2012

@author: gabrielp
'''


#Things to include
#gene expression levels -> coverage for any given gene

import random
import numpy
from numpy import round, arange, cumsum, where, zeros
import pysam
import pybedtools
from optparse import OptionParser
#from collections import namedtuple

class Weighted_interval(pybedtools.Interval):
    """
    
    Wrapper for interval that adds weights to the object
    I've figured out a hacky solution to wrapping this object, but I need a better method
    
    
    """

    def __init__(self, interval):
        
        """
        
        wrapper for interval that includes a weights object, still need to figure out how to engineer this better
        
        """
        
        if len(interval.fields) < 6:
            otherfields = None
        else:
            otherfields = interval.fields[6:]
        
        if interval.name is None:
            name = "."
        else:
            name = interval.name
        
        if interval.score is None:
            score = "."
        else:
            score = interval.score
            
        if interval.strand is None:
            strand = "."
        else:
            strand = interval.strand
        
        pybedtools.Interval.__init__(self, interval.chrom, interval.start, interval.end, name, score, strand, otherfields)
     
        self.weights = zeros(len(self))
        self.peaks = set([])

        #should define accessor method for weights that makes it fail if its refefined to something other than the size of the interval
        
def create_genome(bed):
    
    """
    
    Returns a dictionary (might convert to a named tuple 
    
    bed: the location of a bed file to read from
    
    
    """
    
    if bed is None:
        raise TypeError("bed file is of type None")
    
    #Think about modeling both positive and negative strand here (not for now)
   
    genome_list = []
    tool = pybedtools.BedTool(bed)
    for feature in tool:
        genome_list.append(Weighted_interval(feature))
 
    return genome_list        

def assign_peaks(genome, peak_size, num_peaks):
    
    """

    Returns a list of locations for peaks. Peaks is a list of tuples interval (a pointer), start, stop (locations in the interval)
    randomly assigned to locations in the genome
    
    Input
    genome: genome (will decide on format later)
    peak_size: size of peak to generate
        
    """
     
    #this logic should uniformly distribute peaks across the entire transcriptome
    genome_sizes = map(len, genome)
    total_genome_size = sum(genome_sizes)
    
    #cumulative genome locations so we can assign peaks
    genome_locations = cumsum(genome_sizes)
    
    for count in range(num_peaks):
    
        #initially we don't know what it is so its a duplicate
        duplicate = True
        
        #Small potental for infinate loop here
        while duplicate:
            duplicate = False
            genomic_start = random.randint(0, total_genome_size - 1)
            interval = where((genome_locations > genomic_start) == True)[0][0]
            interval_start = genomic_start - interval
            interval_stop  = interval_start + peak_size
            #make sure the included peak has no overlaps, 
            #if it is re-generate start and stop sites
            
            for cur_start, cur_stop in genome[interval].peaks:
                if ((cur_start < interval_start and cur_stop > interval_start) or 
                    (cur_start < interval_stop and cur_stop  > interval_stop)):
                    duplicate = True
                    break
                    
        
        #adds a peak to the genomic interval
        genome[interval].peaks.add((interval_start, interval_stop))
    return genome



def distribute_background_weights(genome, shape, scale):
    
    """

    Distributes background weights across the entire genome, returns array of weights
    
    genome - array, or list of arrays to represent genome, still optimizing
    
    shape - the shape of the gamma distribution
    
    scale - the scale of the gamma distribution
    
    """
    
    #Figure out the more pythonic way of doing this, I feel like this is hacky  
    for item in genome:
        item.weights = numpy.random.gamma(shape, scale, len(item.weights))
    
    return genome


def distribute_peak_weights(genome, enrichment_coeff):

    """
    PRECONDITION: Must have called assign_peaks and distribute_background_weights on a genome list before calling distribute_peak_weights
    Returns a genome with both peak weights and background weights calculated as an array
    
    peaks - list of Weighted_interval [(start, stop)] that already has background distribution called and peaks defined
    enrichment_coeff - the ratio above background to enrich peaks
    
    """
    
   
    #calculate the average background weight
    average_background_weight = numpy.mean(numpy.concatenate(map(lambda x: x.weights, genome)))
    average_peak_weight = average_background_weight * enrichment_coeff 
    
    #convert the peaks into a list of tuples so I can index them 
    #all together better
    peak_list = []
    for interval in genome:
        for start, stop in interval.peaks:
            peak_list.append((interval, start, stop))
    
    
    #gets the total size of all peaks and the total background weights 
    total_size = 0
    peak_weights = zeros(len(peak_list))
    for i, locs in enumerate(peak_list):
        interval, start, stop = locs
        total_size += stop - start
        peak_weights[i] = sum(interval.weights[start:stop])
    
    #figures the total weight left to assign to peaks
    total_peak_weight = int((total_size * average_peak_weight) - sum(peak_weights))
     
    #for our purposes weight will be discrete
    #randomly distributes all remaining weights in a powerlaw 
    #form to all peak objects
    for count in range(total_peak_weight):
        peak_probs = cumsum(peak_weights * 1. / sum(peak_weights))
        peak = where((peak_probs > random.random()) == True)[0][0]
        peak_weights[peak] += 1
        interval, start, stop = peak_list[peak]
        interval.weights[random.randint(start, stop)] += 1

    #TODO distribute weights normally instead of unifromily 
    return genome
    
        
    
    
        
    #pick a peak based on its total binding weight
    #update weight at peak
    #iterate until weight is done being added
    
    #redistribute peaks as a bionimal distribution 
def distribute_reads(genome, total_reads):
    """
    
    Distributes reads along the weights array returns an array of number of reads at a specific start location
    
    weight - an array of weights to distribute reads along
    num_reads - total number of reads to distribute
    read_length - length of reads to distribute
    
    TODO: print out 
    
    """
    
    total_weight = sum(map(lambda x: sum(x.weights), genome))
    
    for interval in genome:
        interval.start_sites = numpy.zeros(len(interval.weights))
    
    #distribute read starts
    
    #this is broken and where I Will come in tomorrow and fix
    for interval in genome:
        for read_location, weight in enumerate(interval.weights):
            fractional_weight = float(weight) / float(total_weight)
            num_reads = fractional_weight * total_reads
            interval.start_sites[read_location] = round(num_reads)
        
    return genome 

def make_header(genome_lengths):
    
    """
    
    Makes the header of a bam file in python given a file that points to a list of chromosomes and their lengths
    
    """
    header = { 'HD': {'VN': '1.0'} }
    
    chrs = []
    
    #creates all chrs in the genome lengths file first line is 
    #chromosome second line is length of chromosome
    for line in open(genome_lengths):
        line = line.split()
        chrs.append({'SN' : line[0], 'LN' : int(line[1])})
    
    header['SQ'] = chrs
    
    return header

def output_bam(genome, genome_lengths, read_length, outfile_name):
    
    """
    
    Prints a list of reads to standard out
    
    Inputs
    lst: a list of reads to output as a bam file
    
    """
    
    header = make_header(genome_lengths)

    outfile = pysam.Samfile(outfile_name, "wb", header = header)
    
    #prints out all reads in reads
    read_count = 0
    for interval in genome:
        for read_location, num_read in enumerate(interval.start_sites):
            for read in arange(num_read):
                cur_read = pysam.AlignedRead()
                cur_read.qname = "read_%i" % (read_count)
                cur_read.seq = "A" * read_length
                cur_read.flag = 16
                cur_read.rname = 0 
                cur_read.pos = read_location
                cur_read.mapq = 255
                cur_read.cigar = ( (0, read_length - 1), (1, 1) )
                cur_read.mrnm = 0 
                cur_read.isize = 0
                cur_read.qual = "<" * read_length
                outfile.write(cur_read)
                
                read_count += 1
    outfile.close()
    
    
def output_bed(genome, outfile):
    
    """
    
    accepts a peaks list and the name of an outfile and prints out all peaks in the list to the outfile in bed format
    
    """
    #hacky way to make peaks into a bedfile, need to change from chr1 at some point
    bedstring = ""
    for interval in genome:
        for start, stop in interval.peaks:
            bedstring += "%s\t%s\t%s\t\n" % (interval.chrom, start, stop) 

    tool = pybedtools.BedTool(bedstring, from_string= True)
    tool.saveas(outfile)
    
def run():
    usage = """\npython simulate_peaks.py more to come once I figure things out"""
    description = """Simulate peaks.  Gabriel Pratt 2012. 
                     A peak simulator that is used to validate various CLIP-seq peak finding algorithms
                     Refer to: https://github.com/YeoLab/simulate_peaks/wiki for instructions. 
                     Questions should be directed to gpratt@ucsd.edu."""  
    parser = OptionParser(usage=usage, description=description)
    parser.add_option("--locations", "-l", dest="bed", help="A bed file defining where in the genome to distribute reads", type="string", metavar="FILE.bam")
    parser.add_option("--genome", "-g", dest="genome", help="A bed file defining the genome of intrest chromosome <tab> length", type="string", metavar="FILE.bam")
    parser.add_option("--reads", "-r", type="int", default = 10000, dest="reads", help="The aproximate number of reads to assign default:%default")
    parser.add_option("--peak_size", "-p", type="int", default = 50, dest="peak_size", help="The size of peaks to create default:%default")
    parser.add_option("--num_peaks", "-n", type="int", default = 100, dest="num_peaks", help="The number of peaks to create")
    parser.add_option("--peak_coeff", "-c", type="int", default = 5, dest="peak_coeff", help="Coefficent that increases the weight of peak locations, higher value = more obvious peaks default:%default")
    parser.add_option("--gamma_shape", type="int", default = 5, dest="gamma_shape", help="Coefficent that increases the weight of peak locations, higher value = more obvious peaks")
    parser.add_option("--gamma_scale", type="int", default = 10, dest="gamma_scale", help="Coefficent that increases the weight of peak locations, higher value = more obvious peaks")
    parser.add_option("--outfile", "-o", dest="outfile", default="out", help="a bed and bam file root for outputting results, default:%default")
    
    #Things to add
    #Distribution, 
    
    (options,args) = parser.parse_args()

    genome = create_genome(options.bed)

    genome = assign_peaks(genome, options.peak_size, options.num_peaks)

    genome = distribute_background_weights(genome, options.gamma_shape, options.gamma_scale)

    genome = distribute_peak_weights(genome, options.peak_coeff)

    genome = distribute_reads(genome, options.reads)
 
    output_bam(genome, options.genome, 50, options.outfile + ".bam")
    output_bed(genome, options.outfile + ".bed")
    #
    #pass

if __name__ == "__main__":
    run()