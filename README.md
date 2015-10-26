# DownstreamOfSangerSequencing
#Compare blood and tumor DNA to spot change in translated amino acids resulting from tumor mutations
seq_blood='GCCTGCTGAAAATGACTGAATATAAACTTGTGGTAGTTGGAGCTGGTGGCGTAGGCAAGAGTGCCTTGACGATACAGCTAATTCAGAATCATTTTGTGGACGAATATGATCCAACAATAGAG'
seq_tumor = 'ATGACTGAATATAAACTTGTGGTAGTTGGAGCTGCTGGCGTAGGCAAGAGTGCCTTGACGATACAGCTAATTCAGAATCATTTTGTGGACGAATATGATCCAACAATAGAG'

CancerInvestigation='colorectal'


gencode = {
      'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
      'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
      'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
      'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
      'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
      'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
      'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
      'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
      'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
      'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
      'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
      'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
      'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
      'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
      'TAC':'Y', 'TAT':'Y', 'TAA':'_', 'TAG':'_',
      'TGC':'C', 'TGT':'C', 'TGA':'_', 'TGG':'W'}

def translate_frameshifted( sequence ):
      translate = ''.join([gencode.get(sequence[3*i:3*i+3],'X') for i in range(len(sequence)//3)])
      return translate


blood1=translate_frameshifted(seq_blood[0:])
tumor1=translate_frameshifted(seq_tumor[0:])
blood2=translate_frameshifted(seq_blood[1:])
tumor2=translate_frameshifted(seq_tumor[1:])
blood3=translate_frameshifted(seq_blood[2:])
tumor3=translate_frameshifted(seq_tumor[2:])


ExonSeq='MTEYKLVVVGAGGVGKSALTIQLIQNHFVDEYDPTIE'



def calcSeqIdentity(blood1, ExonSeq):
 numPlaces = min(len(blood1), len(ExonSeq))
 score = 0.0
 for i in range(numPlaces):
  if blood1[i] == ExonSeq[i]:
    score += 1.0
 return 100.0 * score/numPlaces
print "similarity score for blood reading frame1 and Exon=",calcSeqIdentity(blood1, ExonSeq)
if calcSeqIdentity(blood1, ExonSeq)>30:
    print "reading frame1 is correct reading frame"
    print calcSeqIdentity(blood1, ExonSeq)
    from Bio import pairwise2
    alignments = pairwise2.align.globalxx(blood1, ExonSeq)
    from Bio.pairwise2 import format_alignment
    for a in pairwise2.align.globalxx(blood1, ExonSeq):
      print (format_alignment(*a)),'Reading Frame 1 alignment with Exon'
else:
    print "Reading frame 1 isn't the correct reading frame"
if calcSeqIdentity(blood1, ExonSeq)>30:
    for a in pairwise2.align.globalxx(blood1, tumor1):
      print (format_alignment(*a)),'Reading Frame 1 alignment of blood over tumor'
    for a in pairwise2.align.globalxx(blood1, tumor2):
      print (format_alignment(*a)),'Reading Frame 1 alignment of blood over tumor'
    for a in pairwise2.align.globalxx(blood1, tumor3):
      print (format_alignment(*a)),'Reading Frame 1 alignment of blood over tumor'





def calcSeqIdentity(blood2, ExonSeq):
 numPlaces = min(len(blood2), len(ExonSeq))
 score = 0.0
 for i in range(numPlaces):
  if blood2[i] == ExonSeq[i]:
    score += 1.0
 return 100.0 * score/numPlaces
print "similarity score for blood reading frame2 and Exon=",calcSeqIdentity(blood2, ExonSeq)
if calcSeqIdentity(blood2, ExonSeq)>30:
    print "reading frame2 is correct reading frame"
    print calcSeqIdentity(blood2, ExonSeq)
    from Bio import pairwise2
    alignments = pairwise2.align.globalxx(blood2, ExonSeq)
    from Bio.pairwise2 import format_alignment
    for a in pairwise2.align.globalxx(blood2, ExonSeq):
      print (format_alignment(*a)),'Reading Frame 2 alignment with Exon'
else:
    print "Reading frame 2 isn't the correct reading frame"
if calcSeqIdentity(blood2, ExonSeq)>30:
    for a in pairwise2.align.globalxx(blood2, tumor1):
      print (format_alignment(*a)),'Reading Frame 3 alignment of blood over tumor'
    for a in pairwise2.align.globalxx(blood2, tumor2):
      print (format_alignment(*a)),'Reading Frame 1 alignment of blood over tumor'
    for a in pairwise2.align.globalxx(blood2, tumor3):
      print (format_alignment(*a)),'Reading Frame 1 alignment of blood over tumor'

def calcSeqIdentity(blood3, ExonSeq):
 numPlaces = min(len(blood3), len(ExonSeq))
 score = 0.0
 for i in range(numPlaces):
  if blood3[i] == ExonSeq[i]:
    score += 1.0
 return 100.0 * score/numPlaces
print "similarity score for blood reading frame3 and Exon=",calcSeqIdentity(blood3, ExonSeq)
if calcSeqIdentity(blood3, ExonSeq)>30:
    print "reading frame3 is correct reading frame"
    from Bio import pairwise2
    alignments = pairwise2.align.globalxx(blood3, ExonSeq)
    from Bio.pairwise2 import format_alignment
    for a in pairwise2.align.globalxx(blood3, ExonSeq):
      print (format_alignment(*a)),'Reading Frame 3 alignment with Exon'
else:
    print "Reading frame 3 isn't the correct reading frame"

if calcSeqIdentity(blood3, ExonSeq)>30:
    if calcSeqIdentity(blood3, tumor1)>30:
        for a in pairwise2.align.globalxx(blood3, tumor1):
            print (format_alignment(*a)),'Reading Frame 3 alignment of blood over tumor1'
    else:
        print ""

    if calcSeqIdentity(blood3, tumor2)>30:
        for a in pairwise2.align.globalxx(blood3, tumor2):
            print (format_alignment(*a)),'Reading Frame 3 alignment of blood over tumor2'
    else:
        print ""

    if calcSeqIdentity(blood3, tumor3)>30:
        for a in pairwise2.align.globalxx(blood3, tumor3):
            print (format_alignment(*a)),'Reading Frame 3 alignment of blood over tumor3'
    else:
        print ""

if calcSeqIdentity(blood1, ExonSeq)>30:
    if calcSeqIdentity(blood1, tumor1)>30:
        count=0
        for k in range(len(blood1)):
            if blood1[k]!=tumor1[k]:
                count=count+1
                print 'total number of point mutations between blood and tumor=',count
        for k in range(len(blood1)):
            if blood1[k]!=tumor1[k]:
                print blood1[k],"in blood mutated to",tumor1[k],"in tumor at Amino Acid position",k+1
    if calcSeqIdentity(blood1, tumor2)>30:
        count=0
        for k in range(len(blood1)):
            if blood1[k]!=tumor2[k]:
                count=count+1
                print 'total number of point mutations between blood and tumor=',count
        for k in range(len(blood1)):
            if blood1[k]!=tumor2[k]:
                print blood1[k],"in blood mutated to",tumor2[k],"in tumor at Amino Acid position",k+1
    if calcSeqIdentity(blood1, tumor3)>30:
        count=0
        for k in range(len(blood1)):
            if blood1[k]!=tumor3[k]:
                count=count+1
                print 'total number of point mutations between blood and tumor=',count
        for k in range(len(blood1)):
            if blood1[k]!=tumor3[k]:
                print blood1[k],"in blood mutated to",tumor3[k],"in tumor at Amino Acid position",k+1

if calcSeqIdentity(blood2, ExonSeq)>30:
    if calcSeqIdentity(blood2, tumor1)>30:
        count=0
        for k in range(len(blood2)):
            if blood2[k]!=tumor1[k]:
                count=count+1
                print 'total number of point mutations between blood and tumor=',count
        for k in range(len(blood2)):
            if blood2[k]!=tumor1[k]:
                print blood2[k],"in blood mutated to",tumor1[k],"in tumor at Amino Acid position",k+1
    if calcSeqIdentity(blood2, tumor2)>30:
        count=0
        for k in range(len(blood2)):
            if blood2[k]!=tumor2[k]:
                count=count+1
                print 'total number of point mutations between blood and tumor=',count
        for k in range(len(blood2)):
            if blood2[k]!=tumor2[k]:
                print blood2[k],"in blood mutated to",tumor2[k],"in tumor at Amino Acid position",k+1
    if calcSeqIdentity(blood2, tumor3)>30:
        count=0
        for k in range(len(blood2)):
            if blood2[k]!=tumor3[k]:
                count=count+1
                print 'total number of point mutations between blood and tumor=',count
        for k in range(len(blood2)):
            if blood2[k]!=tumor3[k]:
                print blood2[k],"in blood mutated to",tumor3[k],"in tumor at Amino Acid position",k+1

if calcSeqIdentity(blood3, ExonSeq)>30:
    if calcSeqIdentity(blood3, tumor1)>30:
        count=0
        for k in range(len(blood3)):
            if blood3[k]!=tumor1[k]:
                count=count+1
                print 'total number of point mutations between blood and tumor=',count
        for k in range(len(blood3)):
            if blood3[k]!=tumor1[k]:
                print blood3[k],"in blood mutated to",tumor1[k],"in tumor at Amino Acid position",k+1
            if blood3[11]!=tumor1[11]:
                print "KRAS mutation confirmed"

    if calcSeqIdentity(blood3, tumor2)>30:
        count=0
        for k in range(len(blood3)):
            if blood3[k]!=tumor2[k]:
                count=count+1
                print 'total number of point mutations between blood and tumor=',count
        for k in range(len(blood3)):
            if blood3[k]!=tumor2[k]:
                print blood3[k],"in blood mutated to",tumor2[k],"in tumor at Amino Acid position",k+1
    if calcSeqIdentity(blood3, tumor3)>30:
        count=0
        for k in range(len(blood3)):
            if blood3[k]!=tumor3[k]:
                count=count+1
                print 'total number of point mutations between blood and tumor=',count
        for k in range(len(blood3)):
            if blood3[k]!=tumor3[k]:
                print blood3[k],"in blood mutated to",tumor3[k],"in tumor at Amino Acid position",k+1
