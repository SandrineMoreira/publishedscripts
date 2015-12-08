#!/usr/bin/python
# vim: tabstop=8 expandtab shiftwidth=4 softtabstop=4

'''
Created on 11 mars 2013

@author: sam
'''

# IMPORT
import os , sys , argparse, re
import glob
from subprocess import Popen, PIPE
from time import time
from datetime import datetime

# ARGUMENTS
ARGS = {}       # Command line arguments dictionary

def parseoptions( ):
        """ Docstring 
            .... """

        print " ".join( sys.argv )

        parser = argparse.ArgumentParser( description="Merge 2 SAM file. Whenever 1 read map on the same contig in both files, keep the one with the best mapping quality" )
        parser.add_argument( '-s1',  '--sam1',  help="First sam file",   required=True )
        parser.add_argument( '-s2',  '--sam2',  help="Second sam file",   required=True )
        parser.add_argument( '-q',  '--mapq',  default=20, type=int, help="Minimum MAPQ value, otherwise skip",   required=False )
        parser.add_argument( '-o',  '--output', help="Output SAM file" )

        #parser.add_argument( '-b',  '--blast',  default=0, type=int,   help="bla bla",   required=True|False )
        # Noter que les tirets dans les noms d'arguments sont transformes en souligne
        global ARGS             # Modify global variable ARGS
        ARGS = parser.parse_args()
        #globals().update(vars(args))        # Makes variables be seen globally

# MAIN
def main() :
    t1    = time()
    print "BEGIN " + str( datetime.now() )

    # ARGS                                                                           
    parseoptions( )     # Parse sys.argv if you want quick and dirty script
                        # sys.argv[ 0 ] : name of the pgm
                        # ARGS.output to access "output" argument value
    
    # Output file
    fileout = ARGS.output
    if not ARGS.output :
        # splittext removes the last extension (even if more than one '.' in the file name
        fileout = os.path.splitext(ARGS.sam1)[0] + "_+_" + os.path.splitext(ARGS.sam2)[0] + ".sam"
    filereport = os.path.splitext( fileout )[0] + ".report.out"

    # Read and parse sam file 
    dh1         = {}
    ds1         = {}
    lref        = []    # List of reference sequences
    readSAM( ARGS.sam1, dh1, ds1, lref )
    dh2         = {}
    ds2         = {}
    readSAM( ARGS.sam2, dh2, ds2, lref )

    mergeSAM( lref, ARGS.sam1, dh1, ds1, ARGS.sam2, dh2, ds2, fileout, filereport)


    dt = time()-t1
    print "END " + str( datetime.now() ) + " ---- time(s) = " + str(dt)
# dheader : 
# dsam[ qname ][ contig ] = mapq     Store the BEST mapq

def readSAM( samfile, dheader, dsam, lref) :
    print "\n# Parsing SAM file %s"%(samfile,) 
    # SAM file
    n       = 0     # Nb line processed
    #@SQ     SN:Dp_cox1-m8_casA+     LN:244
    pathead = re.compile("^@SQ\s+SN:(\S+)\s+")
    f       = open( samfile , "r" )
    for line in f :
        if line[0] == '@' :     # header
            patmat = re.match( pathead, line )
            if patmat :
                contig = patmat.group(1)
                dheader[ contig ] = line
                if not contig in lref :
                    lref.append( contig )
            continue
        if not line.strip() :   # Empty lines             
            continue   

        # if n > 1000 :
        #     break

        # List of fields
        l = line.split()
        qname   = l[0] 
        contig  = l[2]                  # Reference sequence
        mapq    = int( l[4] ) 

        if not qname in dsam :
            dsam[ qname ] = {}

        if not contig in dsam[ qname ] :
            dsam[ qname ][ contig ] = mapq
        elif dsam[ qname ][ contig ] < mapq :
            dsam[ qname ][ contig ] = mapq

        n += 1 
        # printandflush("%i\tlines"%(n,))
    f.close()

def mergeSAM( lref, samfile1, dh1, ds1, samfile2, dh2, ds2, fileout, filereport ) :
    print "\n# Merge files. \n\tOutput %s\n\tReport %s\n\tMinimum mapq value (-q) %i"%( fileout, filereport, ARGS.mapq) 
    fout    = open( fileout , "w" )    # a : append ; w : overwrite
    frep    = open( filereport, "w")

    n1      = 0
    r1      = 0     # Removed from file 1
    q1      = 0     # Removed from file 1 because of quality
    n2      = 0
    r2      = 0     # Removed from file 2
    q2      = 0     # Removed from file 2 because of quality

    for contig in lref :
        if contig in dh1 :
            fout.write( dh1[ contig ] )
        elif contig in dh2 :
            fout.write( dh2[ contig ] )
    
    frep.write( "File\tContig\tread\tmapq\t< best mapq\n" )
    f       = open( samfile1 , "r" )
    for line in f :
        mapq1   = 0
        mapq2   = 0
        n1      += 1
        if line[0] == '@' :     # header
            continue
        if not line.strip() :   # Empty lines             
            continue   

        l = line.split()
        qname   = l[0] 
        contig  = l[2]                  # Reference sequence
        mapq    = int( l[4] ) 

        mapq1   = ds1[ qname ][ contig ]
        if qname in ds2 and contig in ds2[ qname ] :
            mapq2 = ds2[ qname ][ contig ]

        if mapq >= ARGS.mapq and mapq >= mapq1 and mapq >= mapq2 :
            fout.write( line )
        else :
            r1      += 1
            if mapq < ARGS.mapq :
                q1  += 1
            mapqbest = max(mapq1,mapq2,mapq,ARGS.mapq)
            frep.write( "%s\t%s\t%s\t%i\t< %i\n"%(samfile1, contig,qname,mapq,mapqbest))
    f.close()

    f       = open( samfile2 , "r" )
    for line in f :
        mapq1   = 0
        mapq2   = 0
        n2      += 1
        if line[0] == '@' :     # header
            continue
        if not line.strip() :   # Empty lines             
            continue   

        l = line.split()
        qname   = l[0] 
        contig  = l[2]                  # Reference sequence
        mapq    = int( l[4] ) 

        mapq2   = ds2[ qname ][ contig ]
        if qname in ds1 and contig in ds1[ qname ] :
            mapq1 = ds1[ qname ][ contig ]

        if mapq >= ARGS.mapq and mapq >= mapq2 and mapq > mapq1 :
            fout.write( line )
        else :
            r2      += 1
            if mapq < ARGS.mapq :
                q2  += 1
            mapqbest = max(mapq1,mapq2,mapq,ARGS.mapq)
            frep.write( "%s\t%s\t%s\t%i\t< %i\n"%(samfile2, contig,qname,mapq,mapqbest))
    f.close()

    print "\nFile %s\t%i reads\t%i discarded (%i with mapq < %i)"%( samfile1, n1, r1, q1, ARGS.mapq)
    print "File %s\t%i reads\t%i discarded (%i with mapq < %i)*"%( samfile2, n2, r2, q2, ARGS.mapq)
    print "*If reads have same mapq value, the read from first file is conserved (-> read form the second file discarded).\n"


    fout.close()

def write2samORreport( samfile_curr, ds_curr, ds_other, fout, frep, first=True ) :
    f       = open( samfile_curr , "r" )
    for line in f :
        mapq_curr   = 0
        mapq_oth   = 0
        if line[0] == '@' :     # header
            continue
        if not line.strip() :   # Empty lines             
            continue   

        l = line.split()
        qname   = l[0] 
        contig  = l[2]                  # Reference sequence
        mapq    = int( l[4] ) 

        mapq_curr   = ds_curr[ qname ][ contig ]
        if qname in ds_other and contig in ds_other[ qname ] :
            mapq_oth = ds_other[ qname ][ contig ]


        if mapq >= ARGS.mapq and mapq >= mapq_curr and mapq >= mapq2 :
            fout.write( line )
        else :
            mapqbest = max(mapq1,mapq2,mapq,ARGS.mapq)
            frep.write( "%s\t%s\t%s\t%i\t< %i\n"%(samfile1, contig,qname,mapq,mapqbest))
    f.close()



def printandflush( mess ) :
    sys.stdout.write('\r')
    sys.stdout.write( mess )
    sys.stdout.flush()

def runcmd() :
    filelist = glob.glob( fastq ) 
    for filename in filelist :
        print filename
        fileout = filename.split(".")[0] + "SUFFIX"
        cmd = "xxxxxx " 
        print cmd
        p=Popen( cmd , shell=True)
        #p.wait() 

def readfile( filename ) :
    f = open( filename , "r" )
    # Line by line ; If all lines at once lines = f.readlines()
    for line in f :
        # do some stuff
        print line
    f.close()

    # All lines at once
    #with open( filename , 'r' ) as f :
    #    alist = [ line.strip() for line in f ]

def writeinfile( filename , content ) :
    f  = open( filename , "a" )    # a : append ; w : overwrite
    if isinstance( content , str ) :
        f.write( content + "\n" )
    else :
        f.write( "\n".join(content) )
    f.close()

def removefile( filename ) :
    if os.path.exists( filename ) :
                os.remove( filename )

def outfilename( infile , suffix ) :
    fpat    = re.compile( r"(.*)\.(.*)" )
    match   = fpat.match( infile )
    fout    = match.group(1) + suffix
    return fout

# http://biopython.org/DIST/docs/api/Bio.SeqIO._convert-module.html
def fastq_illumina_to_fastq_sanger_Dico() :
    mapping = "".join([chr(0) for ascii in range(0, 64)] 
                + [chr(33 + q) for q in range(0, 62 + 1)]
                + [chr(0) for ascii in range(127, 256)])
    assert len(mapping) == 256 
    return mapping







if __name__ == "__main__":
    try:
        main()
    
    except KeyboardInterrupt:
        print >>sys.stderr, "Program canceled by user..."

