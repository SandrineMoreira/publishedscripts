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
from Bio import SeqIO
import operator

# ARGUMENTS
ARGS = {}       # Command line arguments dictionary
FLANKING = 10

def parseoptions( ):
        """ Docstring 
            .... """

        print " ".join( sys.argv )

        parser = argparse.ArgumentParser( description="Search the different categories of edited reads (completely, partially)" )
        parser.add_argument( '-f',  '--fasta',   help="Reference fasta",   required=True )
        parser.add_argument( '-g',  '--gff',   help="gff file for the list of SNP. Only SNP present in the gff file are considered",   required=True )
        parser.add_argument( '-s',  '--sam',   help="Alignment sam file",   required=True )
        parser.add_argument( '-n',  '--nbmin', default=1,  type=int, help="Minimal number of reads in a category to be considered" ) 
        parser.add_argument( '-l',  '--flanking', default=10,  type=int, help="Flanking regions length" )
        parser.add_argument( '--nopartial', help="Keep only reads that cover completely the edited region", action='store_true', default=False)
        parser.add_argument( '--sortbyscore', help="Sort results by score rather than by number", action='store_true', default=False)        
        parser.add_argument( '-o',  '--output', help="Output result file" )

        # parser.add_argument( '-m',  '--pmin', default=1,  type=int, help="Periodicity min" ) 
        # parser.add_argument( '-M',  '--pmax', default=10,  type=int, help="Periodicity max" ) 
        # parser.add_argument( '-l',  '--lmin', default=2,  type=int, help="Motif length min" ) 
        # parser.add_argument( '-L',  '--lmax', default=5,  type=int, help="Motif length max" ) 
        # parser.add_argument( '-p',  '--position', default='-',  help="Look before '-' or after '+' edited position" )
        # parser.add_argument( '-e',  '--edited', default=0,  type=int, help="Look on non-edited sequence '0' or edited '1'" )
        # parser.add_argument( '-d',  '--direction', default='3',  type=int, help="Edition proceed 3'->5' '3' or 5'->3' '5' " )


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
    
    if ARGS.output :
        fileout = ARGS.output
    else :
        # splittext removes the last extension (even if more than one '.' in the file name
        fileout = os.path.basename(ARGS.sam).split(".")[0] + "_Edit.txt"
    removefile( fileout )

    # Flnking region length
    FLANKING    = ARGS.flanking

    # Get Fasta seq
    dfasta        = readFasta(ARGS.fasta)

    # Get list of SNPs
    dsnp        = readgff( ARGS.gff, 'snp' )

    # Get list of edited regions
    dreged      = getRegEd( dfasta, dsnp )

    # Find pos edited
    parseSam( ARGS.sam, dreged )
    
    lregedsort  = sorted(dreged.values(), key=lambda reged: reged.name, reverse=False)
    for reg in lregedsort :
        reg.printreged(fileout)

    print "\n#[OUT] Output file : %s"%(fileout,)
    print "\tScore\tsimilarity with the edited sequence = nb_match * 2 + nb_'-' * 1 "
    print "\t'TOT'\tsum of 'A','T','G','C' ; '?' and '-' are not counted"
    print "\t'-'\tPosition not mapped (usually, end of the read)"
    print "\t'.'\tPosition like in edited sequence. When the nucleotide does not correspond, the actual nucleotide is displayed"
    print "\t'?'\tNucleotide are different on the 2 reads of a pair"

    dt = time()-t1
    print "END " + str( datetime.now() ) + " ---- time(s) = " + str(dt)


def parseSam( samfile, dreged ) :
    print "\n#[IN] Parsing SAM file  %s"%(samfile,) 
    f   = open( samfile , "r" )
    n   = 0     # Nb line processed

    for line in f :
        if line[0] == '@' : # header
            continue
        # List of fields
        l = line.split()
        n += 1

        # if n > 1000 :
        #     break
        if n % 500 == 0 :
            printandflush( "\t ... %i \tsam lines "%(n,))


        contig          = l[2]                  # Reference sequence
        contigprefix    = contig
        # if '_ED' in contig :            # In case '_ED' appended at the end 
        #     contigprefix = contig[0:-3] # get prefix contig name

        # read map on a sequence edited    
        if contigprefix in dreged :
            qname   = l[0]          # read name
            qprefix = qname.split('_')[0]
            # flag    = int( l[1] )   
            pos     = int( l[3] )   # position of the leftmost base
            # mapq    = int( l[4] )  
            cigar   = l[5].strip()
            pnext   = int( l[7] )
            tlen    = int( l[8] )   # Mapping length
            seq     = l[9]
            seqmap  = getSeqMapped( seq, cigar) # Padded and/or trimmed sequence that mapped on the reference 

            # Is edited region covered by read
            readbegin   = pos           # if paired read not mapped
            readend     = pos + len(seqmap)
            keep        = 0             # keep this read if overlap an edited site
            re          = dreged[ contigprefix ]  # re: region edited
            for p in re.lpos :      
                if p in range( readbegin-1 , readend ) : # readbegin-1 because read positions begin at 1 and pos in lpos begin at 0
                    keep += 1


            # Record read info
            if keep > 0 :
                seqedit     = [ "-" for i in range (0,len(re.lpos)) ]       # Init edited sequence with unknown value
                if not qprefix in re.dreads :           # If no read in dico for this pair
                    re.dreads[ qprefix ] = seqedit      # New entry in dico
                else :
                    seqedit = list(re.dreads[ qprefix ])      # ... or take the seqedit from previous paired read

                # For all positions edited
                # !!!! This one is complicated. Remember that readbegin/readend are from SAM (begins at 1) and lpos begins at 0
                for i in range(0,len(re.lpos)) :
                    p = (re.lpos[i] + 1) - readbegin + 1        # Find the corresponding pos on the read sequence (!!!!begins at 1)
                    if p in range(1,len(seqmap)+1) :      # if the position is on the read
                        if seqedit[i] == '-' :          # if the nucleotide is unknown, 
                            seqedit[i] = seqmap[p-1]      # replace by the new one
                        elif not seqedit[i] == seqmap[p-1] :  # if nucleotide are different on both pair
                            seqedit[i] = '?'                # Put a '?'
                re.dreads[ qprefix ] = ''.join(seqedit)              # Update the new edited sequence  

    #For all sequence edited
    for re in dreged.values() :
        ledit   = []                # Temporary list of edited sequences
        for seqedit in re.dreads.values() :     # For all sequence of edited positions

            if '-' in seqedit and ARGS.nopartial == True :      #If we don't want partial sequences
                continue

            score   = 0
            # count the number of distinct edited seq
            if not seqedit in ledit :     
                # re.deditseq[seqedit] = 1
                ledit.append(seqedit)
                score  = editscore( seqedit, re.regED, re.lposED )                 # Distance of editions
                d = { "seq" : seqedit, "nb" : 1, "score": score}
                re.lregED.append(d)
            else :
                i = ledit.index(seqedit)
                if re.lregED[i]["seq"] != seqedit :
                    print "OUILLE"
                re.lregED[i]["nb"] += 1
                score   = re.lregED[i]["score"]
                # re.deditseq[seqedit] += 1

            # For each position in edit sequence, update stat dico
            # ltmp = list( seqedit ) 
            if score > (re.maxscore() / 2 ) :           # If the sequence is edited
                for i in range(0, len(seqedit)) :       # For all positions edited 
                    c = seqedit[i]                      # Get the character at this position
                    if not c in re.dregEDstat :           # If character not in dico of different characters, add it
                        re.dregEDstat[c] = [0 for i in range(0,len(re.regED)) ]
                    re.dregEDstat[c][i] += 1              # Increment the nb of this character for this position 


        # print re.name
        # for seqedit in re.deditseq.keys() :
        #     print "%s [%i]"%(seqedit,re.deditseq[seqedit]) 

# distance between seq and the region edited
# Both sequences should have the same length
def editscore( seq, regED, lposED) : 
    if len(seq) != len(regED) :
        return -1
    score   = 0
    for i in lposED :
        if seq[i] == regED[i] :
            score += 2
        elif seq[i] == '-' :
            score += 1
    return score

        
# Get padded sequence that map on reference
def getSeqMapped( seq, cigar ) :

    # dCigar          = { 'M':0, 'D':0, 'N':0, 'X':0, 'P':0, 'I':0 }
    cigarpattern    = re.compile( r"(\d+)([MIDNSHPX])" )    # Cigar element
    cigarmatch      = re.findall( cigarpattern , cigar )
    seqmap          = ""
    pos             = 0     # position on sequence 
    for c in cigarmatch : 
        length  = int( c[0] )
        type    = c[1]
        
        if type in ['M','=','X'] :
            seqmap  += seq[pos:pos+length]
            pos     += length
        elif type in ['D','N'] :
            seqmap  += length*'-'
        elif type in ['I','P','H','S'] :
            pos     += length
    return seqmap

def cigardico( cigar ) :
    dCigar = { 'M':0, 'D':0, 'N':0, 'X':0, 'P':0, 'I':0 }
    cigarpattern = re.compile( r"(\d+)([MIDNSHPX])" )    # Cigar element
    cigarmatch  = re.findall( cigarpattern , cigar )
    
    for c in cigarmatch : 
        length  = int( c[0] )
        type    = c[1]
        
        if not type in dCigar :
            dCigar[ type ] = 0
        dCigar[ type ] = dCigar[ type ] + length
    
    return dCigar

def getRegEd( dfasta, dsnp ) :
    dreged      = {}        # dico of edited regions ; key = seqname , value = reged
    global FLANKING

    for seqname in dfasta.keys() :

        # if this sequence has snp 
        if seqname in dsnp :
            # print seqname
            r   = reged(seqname)
            # Non edited sequence 
            r.seq     = dfasta[seqname].seq

            # Edited sequence
            lseqED    = list(r.seq)

            # Dico of gff snp pos
            dpos = {}       # key = pos, value = gff
            for g in dsnp[ seqname ] :
                dpos[ g.start - 1 ] = g  # In gff, pos begins at 1, but we want to begin at 0

            # List of positions in edited region 
            lposED      = sorted(dpos.keys())
            minpos      = max(0,lposED[0]- FLANKING)
            maxpos      = lposED[-1] + FLANKING
            r.lpos      = [ x for x in range( minpos, maxpos+1)]

            # For all SNP in the sequence
            for i in range( 0, len(r.lpos) ) :
                p = r.lpos[i]
                if p in dpos :
                    r.lposED.append( i )
                    g = dpos[ p ]
                    if 'alt' in g.dattributes :     # SNP sequence (if any, other keep ref seq) 
                        lseqED[ p ] = g.dattributes[ 'alt' ]
                r.regED += lseqED[ p ]    # sequence in region edited

            # for g in dsnp[ seqname ] :
            #     r.lpos.append( g.start )        # Position of SNP
            #     if 'alt' in g.dattributes :     # SNP sequence (if any, other keep ref seq)
            #         lseqED[ g.start ] = g.dattributes[ 'alt' ]
            #     r.regED += lseqED[ g.start ]    # SNP sequence in region edited
            r.seqED   = ''.join( lseqED )       # Final edited sequence

            # Stat for all nucleotides and each position
            for k in r.dregEDstat.keys() :
                r.dregEDstat[k] = [0 for i in range(0,len(r.regED)) ]

            # Add to dico
            # r.printreged()
            dreged[seqname]     = r 
    return dreged

def readFasta( fastafile ) :
    print "\n#[IN] Parsing fasta file %s"%( fastafile, ) 
    ffasta  = open( fastafile, "r" ) 
    dfasta  = {} 
    nfasta  = 0  
    for rec in SeqIO.parse(ffasta, "fasta") :
        afasta = fasta( rec.id, rec.description, "%s"%(rec.seq,) )
        dfasta[rec.id] = afasta
    ffasta.close()
    return dfasta

class fasta : 
    def __init__( self, name, desc="", seq="") :   
        self.name = name
        self.desc = desc
        self.seq  = seq

# Edited region
class reged :
    def __init__( self, name ) :   
        self.name       = name
        self.seq        = ""
        self.seqED      = ""
        self.lpos       = []    # positions of regED (with reference to seq)
        self.regED      = ""    # Edited region of the fully edited sequence
        self.lposED     = []    # Only the edited positions (with reference to regED i.e. if lposED == x then regED[x] is a position edited)

        # self.deditseq   = {}    # key = seq of pos edited, value = nb of reads having this seq of pos edited
        self.dreads     = {}    # temporary dico of reads
                                # { key = read prefix : value = seq of pos edited    ]
        self.lregED     = []    # list of regions edited in the sam file. List of dico { "seq" : , "nb" : , "score" : }
        self.dregEDstat   = {"A":[], "T":[], "C":[], "G":[], "-":[], "?" : []}   #  For EDITED sequences. dico { "A": list of pos , "T": list of pos, "C": list of pos, "G": list of pos, "-", "?" }

    def maxscore( self ) :
        return len( self.lposED ) * 2

    # Print edited positions
    def printreged( self, fileout=None ) :

        # minpos  = max(self.lpos[0] - 5 , 0)
        # maxpos  = min(self.lpos[-1] + 5 , len(self.seq))
        minpos  = self.lpos[0] 
        maxpos  = self.lpos[-1] + 1

        space   = ""
        mess    = ""
        mess    += ">" + self.name + " [%i-%i]\n"%(minpos+1,maxpos+1)

        mess    += self.seq[minpos:maxpos] + "\n" 
        #mess    += self.seqED[minpos:maxpos] + "\n" 
        # Mark edited pos
        mess    += self.editionseq (self.seq[minpos:maxpos], self.seqED[minpos:maxpos],  match = ' ' , mismatch='*' ) 
        mess    += "\tmax_score=%i\n" %(self.maxscore(),)
        mess    += self.seqED[minpos:maxpos] + "\n" 
        # mess    += self.regED + "\n" 

        #mess    += self.seq[minpos:maxpos] + "\n" 

        # Print edited sequence
        nED         = 0         # Nb of edited sequence (having at least maxscore/2)
        ntot        = 0         # Total number of sequences
        lED         = []        # List of edited sequences
        if ARGS.sortbyscore :
            lregEDsort  = sorted(self.lregED, key=operator.itemgetter('score','nb'), reverse=True)
        else :
            lregEDsort  = sorted(self.lregED, key=operator.itemgetter('nb','score'), reverse=True)
        for ED in lregEDsort : # sort by quantity
            if ED['nb'] >= ARGS.nbmin :
                mess    += self.editionseq(self.regED, ED['seq'], match = '.' , mismatch='') 
                mess    += "\t %i score=%i\n"%( ED['nb'],ED['score'])
                ntot    += ED['nb']
                if ED['score'] > (self.maxscore()/2) :
                    nED += ED['nb']
                    lED.append(ED)
            # else :
            # mess += ED['seq'] + "--\n"

        # print stat dico for all sequences 
        mess += "\nSTAT for EDITED sequences\n"
        mess += "Edited sequences = sequences having a score > max_score/2 \n"
        mess += "Nb edited seq = %i / %i \n"%( nED, ntot)
        # space   = "     " 
        # mess += "   " + space + space.join( self.regED ) + "\n"
        # # Each nucleotide
        # for c in sorted(self.dregEDstat.keys() ) :
        #     lval = self.dregEDstat[c]     # Number of this nucleotide along the sequence
        #     mess += c + "  "            
        #     for v in lval :             # Write all the values 
        #         mess += "%6i"%(v)
        #     mess += "\n"
        # # Total
        # mess    += "tot"
        # for i in range(0, len(self.regED)) :
        #     tot     = 0
        #     for c in self.dregEDstat.keys() :
        #         if c in ['A','T','G','C'] :
        #             tot += self.dregEDstat[c][i]
        #     mess += "%6i"%(tot,)
        # mess += "\n\n"

        ### Vertical
        mess += "Pos\tNuc\t"
        for c in sorted( self.dregEDstat.keys() ) :
            mess    +=  c + "\t"
        mess += "TOT\n"
        # for i in range(0, len(self.regED)) :
        for i in range(0, len(self.regED)) :
            tot     = 0

            isED    = ''
            if i in self.lposED :
                isED = '*'

            mess    += "%i\t%s%s\t"%(self.lpos[i], self.regED[i], isED)
            for c in sorted( self.dregEDstat.keys() ) :
                mess    += "%i\t"%(self.dregEDstat[c][i],)
                if c in ['A','T','G','C'] :
                    tot += self.dregEDstat[c][i]
            mess += "%i\n"%(tot,)
        mess += "\n\n"




        #
        nEDA    = 0                     # Total nb of A non-edited
        nEDC    = 0                     # Total nb of C non-edited
        nseqA   = 0                     # Number of sequences having non-edited A
        nseqC   = 0                     # Number of sequences having non-edited C
        nG      = self.regED.count('G') # Number of A->G edited positions
        nT      = self.regED.count('T') # Number of C->T edited positions
        
        # For all edited sequences
        for ED in lED :
            foundA = False
            foundC = False

            # For each edited positions
            # for i in range(0, len(self.regED)) :
            for i in self.lposED : 
                # if the position is non-Edited (in an edited sequence)

                if self.regED[i] == 'G' and ED['seq'][i] == 'A' :
                    nEDA  += ED['nb']
                    foundA  = True

                if self.regED[i] == 'T' and ED['seq'][i] == 'C' :
                    nEDC  += ED['nb']
                    foundC = True
            if foundA :
                nseqA += ED['nb']
            if foundC :
                nseqC += ED['nb']

            # if self.regED[i] == 'G' :
            #     nA      = float(self.dregEDstat['A'][i])
            # elif self.regED[i] == 'T' :
            #     nC      = float(self.dregEDstat['C'][i])
        mess += "Ratio of non-edited nucleotide by pos and by seq: \nEx: Ratio non-edited A = Nb_A/(Nb_posA->G * Nb_seqEdited)\n"

        mess += "A->G: Nb_A=%i\tNb_posA->G=%i\tNb_seq_with_non-edA=%i\tNb_seqEdited=%i\t"%(nEDA, nG, nseqA, nED )
        if nED > 0 and nG > 0 :
            mess += "Ratio= %.3f %%\n"%(float(nEDA)*100/(nG*nED), )
        else :
            mess += "Ratio= -/-\n"
        mess += "C->T: Nb_C=%i\tNb_posC->T=%i\tNb_seq_with_non-edC=%i\tNb_seqEdited=%i\t"%(nEDC, nT, nseqC, nED )
        if nED > 0 and nT > 0 :
            mess += "Ratio= %.3f %%\n"%(float(nEDC)*100/(nT*nED), )
        else :
            mess += "Ratio= -/-\n"

        mess += "- "*40 + "\n\n"

        #re.dregEDstat[c][i] += 1

        if not fileout is None :
            f = open(fileout, "a")
            f.write(mess)
            f.close()
            #print mess
        else :
            print mess

    # spacemismatch : put a space when thereis a mismatch, otherwise, put the seq char
    def editionseq(self, seqRef, seq, match='.', mismatch=' ' ) :
        edition = ""
        # Sequences do not have same length => pbm
        if len( seqRef ) != len( seq ) :
            # return edition
            return "Ref:%i - seq:%i"%(len(seqRef), len(seq))
        # For all positions
        for i in range( 0, len( seqRef ) ) :
            if seqRef[ i ] == seq[ i ] :
                edition += match
            elif len( mismatch ) > 0 :
                edition += mismatch
            else :
                edition += seq[ i ]
        return edition

    def editionseqOLD( self, seq, minpos ) :
        edition = ""
        last    = minpos - 1
        lseq    = list(seq)
        lregED  = list(self.regED)
        # For all edited positions
        #for i in range( 0, len(lpos)) :
        for i in range( 0, len(self.regED)) :
            pos     = self.lpos[i]
            edition += (pos-1 - last ) * ' '
            if lseq[i] == lregED[i] :
                edition += '.'
            else :
                edition += lseq[i]
            last    = pos
        return edition


#######  GFF
def readgff( fgff, feature ) :
    print "\n#[IN] Read GFF file %s"%( fgff )

    f           = open( fgff, "r" )     
    dgff        = {}    # key = contig(sequence), value = gff
    n           = 0
    
    for line in f :         
        if line[0] == '#' :     # comment             
            continue         
        if not line.strip() :   # Empty lines             
            continue         
        l   = line.split("\t")         
        if len(l) < 9 :         # Not a gff line  
            print "\tNot a gff line: %s"%( line.strip(), )           
            continue         
        g   = GFF(l)
        if g.feature == feature.lower() :
            n += 1
            if not g.contig in dgff :
                dgff[ g.contig ] = []
            dgff[ g.contig ].append( g )
            
    f.close()
    print "\t%i line(s) with feature %s"%( n, feature )
    return dgff

class GFF : 
    def __init__(self, lvalues):
        if len(lvalues) < 9 :
            return None 
        self.contig      = lvalues[0].strip() 
        self.source     = lvalues[1].strip() 
        self.feature    = lvalues[2].strip().lower()
        self.start      = int(lvalues[3])
        self.end        = int(lvalues[4])
        self.score      = lvalues[5].strip()
        self.strand     = lvalues[6].strip()
        # if self.strand not in [ '.', '-', '+' ] :
        #     self.strand = '.'
        self.frame      = lvalues[7].strip()
        # if self.frame not in [ '.', '0', '1', '2' ] :
        #     self.frame = '.'
        self.attribute  = lvalues[8].strip()
        self.dattributes = {}
        self.id         = ""
        self.parent = []

        # Parse attributes
        # To find attribute ID : dattributes["id"]
        # !!! attribute keys in LOWER CASE
        latt    = self.attribute.split(';')
        if len( latt ) > 0 :
            for a in latt :
                if '=' in a :
                    k,v = a.split("=")
                    if not k is None and not v is None :
                        self.dattributes[ k.strip().lower() ] = v.strip()
    
        # ID is the non-Edited name
        if "id" in self.dattributes :
            self.id     = self.dattributes["id"]
        else :
            self.id     = "%s_%s_%i-%i"%( self.contig, self.feature, self.start, self.end )
        # Parent
        if "parent" in self.dattributes :
            lparent = self.dattributes["parent"].split(',')
            for p in lparent :
                self.parent.append( p )
        
    def filtergff(self, feature ) :
        return self.feature == feature.lower() 

    def gff2text(self) :
        line = "%s\t%s\t%s\t%i\t%i\t%s\t%s\t%s\t%s"%(self.contig, self.source, self.feature, self.start, self.end, 
                                                self.score, self.strand, self.frame, self.attribute)  
        return line









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


def samFlagPosToInt( samflag, position ) :
    if (( samflag >> position ) & 1 ) == 1 :
        return True
    else :
        return False

def isReadPaired( samflag ) :
    return samFlagPosToInt(samflag , 0 ) 

def isBothMatesMapped( samflag ) :
    return samFlagPosToInt(samflag , 1 ) 

def isReadMapped( samflag ) :
    return not samFlagPosToInt(samflag , 2) 

def isMateMapped( samflag ):
    return not samFlagPosToInt(samflag , 3)

def isReadReverseStrand( samflag ) :
    return samFlagPosToInt(samflag , 4 ) 

def isReadSenseStrand( samflag ) :
    return not samFlagPosToInt(samflag , 4 ) 
def isMateReverseStrand( samflag ) :
    return samFlagPosToInt(samflag , 5 ) 

def isReadFirstInPair( samflag ) :
    return samFlagPosToInt( samflag , 6 ) 

def isReadSecondInPair( samflag ) :
    return samFlagPosToInt(samflag , 7 ) 

def isPrimaryAlignments( samflag ) :
    return not samFlagPosToInt(samflag , 8 )






if __name__ == "__main__":
    try:
        main()
    
    except KeyboardInterrupt:
        print >>sys.stderr, "Program canceled by user..."

