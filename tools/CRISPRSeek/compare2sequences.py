"""
    This tool recursively map sequenced files to any number of index files in given order. 
 
    usage: compare2sequences.py [options]
    --inputfile1: Sequence input file path 1 that contains sequences to be searched for potential gRNAs
    --inputfile2: Sequence input file path 2 that contains sequences to be searched for potential gRNAs
    --outputdir: output directory
    --output: output file
    --user_email: user email address 
    --maxmismatch: Maximum mismatch allowed in off target search, default 3. Warning: will be considerably slower if set > 3    
    --searchDirection: Both, 1to2 and 2to1
    --findPairedgRNAOnly: Choose whether to only search for paired gRNAs in such an orientation that the first one is on minus strand called reverse gRNA and the second one is on plus strand called forward gRNA. TRUE or FALSE, default FALSE
    --mingap: Min Gap, default 0
    --maxgap: Max Gap, default 20
    --PAMsize: PAM length, default 3
    --PAM: PAM sequence after the gRNA, default NGG
    --PAMPattern: Regular expression of protospacer-adjacent motif (PAM), default N[A|G]G$
    --gRNAsize: The size of the gRNA, default 20
    --weights: numeric vector size of gRNA length, default 0, 0, 0.014, 0, 0, 0.395, 0.317, 0, 0.389, 0.079, 0.445, 0.508, 0.613, 0.851, 0.732, 0.828, 0.615, 0.804, 0.685, 0.583 which is used in Hsu et al., 2013 cited in the reference section  

    usage: compare2sequences.py --inputfile --outputdir --user_email --output --bsgenomename
"""
# imports
import os, re, string, sys
from sys import argv, exit, stderr
from optparse import OptionParser
from os.path import dirname, exists, join
from os import system
from subprocess import Popen, PIPE

bin_dir = dirname(argv[0])

# error
def stop_err( msg ):
    sys.stderr.write( "%s\n" % msg )
    sys.exit()

# main
def main():
    # define options
    parser = OptionParser()
    parser.add_option("", "--inputfile1", dest="inputfile1")
    parser.add_option("", "--inputfile2", dest="inputfile2")
    parser.add_option("-o", "--outputdir", dest="outputdir")
    parser.add_option("-f", "--output", dest="output")
    parser.add_option("-u", "--user_email", dest="useremail")
    parser.add_option("-m", "--maxmismatch", dest="maxmismatch")
    parser.add_option("", "--searchDirection", dest="searchDirection")
    parser.add_option("", "--findPairedgRNAOnly", dest="findPairedgRNAOnly")
    parser.add_option("", "--mingap", dest="mingap")
    parser.add_option("", "--maxgap", dest="maxgap")
    parser.add_option("", "--PAMsize", dest="PAMsize")
    parser.add_option("", "--PAM", dest="PAM")
    parser.add_option("", "--PAMPattern", dest="PAMPattern")
    parser.add_option("", "--gRNAsize", dest="gRNAsize")
    parser.add_option("", "--weights", dest="weights")
    parser.add_option("-j", "--jobsubmit", dest="jobsubmit")
    parser.add_option("-s", "--servicename", dest="servicename")

    # parse
    options, args = parser.parse_args()
    
    try:
        # retrieve options
        inputfile1                = options.inputfile1
        inputfile2                = options.inputfile2
        outputdir                 = options.outputdir
        output                    = options.output
        username                  = options.useremail
        maxmismatch               = options.maxmismatch
        searchDirection           = options.searchDirection
        findPairedgRNAOnly        = options.findPairedgRNAOnly
        mingap                    = options.mingap  	
        maxgap                    = options.maxgap  	
        PAMsize                   = options.PAMsize 	
        PAMPattern                = options.PAMPattern 	
        PAM                       = options.PAM  	
        gRNAsize                  = options.gRNAsize 	
        weights                   = options.weights	
 
  
	if not inputfile1 :
		print >>stderr, 'Error: Input file 1 is not defined.'
		exit( 128 )

	if not inputfile2 :
		print >>stderr, 'Error: Input file 2 is not defined.'
		exit( 128 )

	if not outputdir :
		print >>stderr, 'Error: Output dir is not defined.'
		exit( 128 )

        print outputdir 
        
    
    except Exception, ex:
        stop_err('Error running compare2sequences.py\n' + str(ex))

    runCRISPRSeek( inputfile1, inputfile2, outputdir, output, username, maxmismatch, 
                 searchDirection, findPairedgRNAOnly, mingap, maxgap, PAMsize, PAMPattern, 
                 PAM, gRNAsize, weights)
    sys.exit(0)

def runCRISPRSeek( inputfile1, inputfile2, outputdir, output, username, maxmismatch, 
                 searchDirection, findPairedgRNAOnly, mingap, maxgap, PAMsize, PAMPattern, 
                 PAM, gRNAsize, weights):
        if (inputfile1.startswith('__gt__')):
            fp = open ( outputdir+"/inputfile1.fa", 'w') 
	    print >>fp, '%s'%parse_content(inputfile1)
            fp.close()
            inputfile1=outputdir+"/inputfile1.fa"

        if (inputfile2.startswith('__gt__')):
            fp = open ( outputdir+"/inputfile2.fa", 'w') 
	    print >>fp, '%s'%parse_content(inputfile2)
            fp.close()
            inputfile2=outputdir+"/inputfile2.fa"
            
        
        print  outputdir+"/compare2sequences.R"

        fp = open ( outputdir+"/compare2sequences.R", 'w') 
  
        Rprog="library(CRISPRseek)\n"
	Rprog+="setwd('"+outputdir+"')\n"
        Rprog+="outputDir <- file.path('"+outputdir+"','CRISPRSeek' )\n"
        Rprog+="compare2Sequences('%s','%s',\n"%(inputfile1, inputfile2)
        Rprog+="   REpatternFile = system.file('extdata', 'NEBenzymes.fa', package = 'CRISPRseek'),\n"
        if (maxmismatch > 0):
           Rprog+="   max.mismatch=%s,\n"%str(maxmismatch)
        if (searchDirection):
           Rprog+="   searchDirection='%s',\n"%str(searchDirection)
        if ( findPairedgRNAOnly == "Yes"):
           Rprog+="   findPairedgRNAOnly=TRUE,\n"
           Rprog+="   min.gap=%s,\n"%str(mingap)
           Rprog+="   max.gap=%s,\n"%str(maxgap)
        if (PAMsize and PAMsize.find('@')==-1):
           Rprog+="   PAM.size=%s,\n"%str(PAMsize)
        if (PAMPattern and PAMPattern.find('@')==-1):
           rep={'__ob__' : '[', '__cb__' : ']', '_p_' : '|', '_d_' : '$'}
           for k,v in rep.iteritems():
              PAMPattern = PAMPattern.replace(k, v)
           Rprog+="   PAM.pattern='%s',\n"%str(PAMPattern)
        if (PAM and PAM.find('@')==-1):
           Rprog+="   PAM='%s',\n"%str(PAM)
        if (gRNAsize and gRNAsize.find('@')==-1):
           Rprog+="   gRNA.size=%s,\n"%str(gRNAsize)
        if (weights and weights.find('@')==-1):
           Rprog+="   weight=c(%s),\n"%weights
        Rprog+="   outputDir = outputDir, overwrite = TRUE)\n\n"

        print Rprog

	print >>fp, '%s'%Rprog
        system( 'mkdir -p '+outputdir+'/CRISPRSeek')     
	fp.close()
	command_line = 'module load R/3.1.0;Rscript --no-save %s/compare2sequences.R'%outputdir;

	print command_line
	exit_status = system( command_line )

	exit( exit_status ) 

def write_workflow( file ):
        fp = open ( file, 'w') 
        sep='\t'

        stepline=step % locals()
        print >>fp, '%s'%stepline
        fp.close()

class NoContentParsedError( Exception ) :
	pass
class ContentFormatError( Exception ) :
	pass

def replace_space(content) :
        content = re.sub('[\s\t,]+', '_', content)
        return content

def parse_content( content, ncols=8, base64=False, verbose=0 ) :
	'''
	This is a function that parse the inputparam content and 
	returns the base64 encoded string if base64 is True otherwise 
	returns the concatenated string with by the conversion to make
	('\t' -> ',', ' ' -> ',', '\n'->':').
	'''

	content = content.replace( '__gt__', '>' ) 
	content = re.sub( '[:,]', '\n', content ) 
        return content

if __name__ == "__main__":
    main()
