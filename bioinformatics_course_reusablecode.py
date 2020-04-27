'''
# Frequency mapper for BioInformatics Course
def FrequencyMap(Text, k):
    freq = {}
    n = len(Text)
    for i in range(n-k+1):
        Pattern = Text[i:i+k]
        freq[Pattern]= 0
        for i in range(n-k+1):
            if Text[i:i+k] == Pattern:
                freq[Pattern] += 1
    return freq



### ### ### ###

Text = ""
k = ""

def FrequentWords(Text, k):
    words = []
    freq = FrequencyMap(Text, k)
    m = max(freq.values())
    for key in freq:
        if freq[key] == m: 
        	pattern = key
        	words.append(pattern)
    return words

# Copy your FrequencyMap() function here.
def FrequencyMap(Text, k):
    freq = {}
    n = len(Text)
    for i in range(n-k+1):
        Pattern = Text[i:i+k]
        freq[Pattern]= 0
        for i in range(n-k+1):
            if Text[i:i+k] == Pattern:
                freq[Pattern] = freq[Pattern] + 1
    return freq


## OR ##

## FInd how many time a pattern is in a text
def PatternCount(Text, Pattern):
    count = 0
    for i in range(len(Text)-len(Pattern)+1):
        if Text[i:i+len(Pattern)] == Pattern:
            count = count+1
    return count 

def FrequencyMap(Text, k):
    freq = {}
    n = len(Text)
    for i in range(n-k+1):
        Pattern = Text[i:i+k]
        freq[Pattern]= PatternCount(Text,Pattern)
        
    return freq

def FrequentWords(Text, k):
    words = []
    freq = FrequencyMap(Text, k)
    m = max(freq.values())
    for key in freq:
        if freq[key] == m:
            words.append(key)
            words.sort()
    return words

FrequentWords("ACGTTGCATGTCGCATGATGCATGAGAGCT",4)

### ### ### ### ###


### ReverseComplement Function ###

# Highest Level Function
def ReverseComplement(Pattern):
    print("Initial Pattern: " + Pattern)
    Pattern = Reverse(Pattern) # reverse all letters in a string
    print("Reversed Pattern: " + Pattern)
    Pattern = Complement(Pattern) # complement each letter in a string
    print("Complement Pattern: " + Pattern)
    return Pattern

def Reverse(Pattern):
    PatternReversed = ""
    for i in range(len(Pattern)):
        PatternReversed = Pattern[i] + PatternReversed
        
    Pattern = PatternReversed
    return PatternReversed

def Complement(Pattern):
    list_pattern = list(Pattern)
    for i in range(len(list_pattern)):
        if list_pattern[i] == "A":
            list_pattern[i] = "T"
        elif list_pattern[i] == "C":
            list_pattern[i] = "G"
        elif list_pattern[i] == "G":
            list_pattern[i] = "C"
        elif list_pattern[i] == "T":
            list_pattern[i] = "A"
        
    Pattern = "".join(list_pattern)
    return Pattern

print(ReverseComplement("ACTGAGTC"))


### Pattern Matching Function ###

def PatternMatching(Pattern, Genome):

    positions = []
    def find_all(Pattern, Genome):
        start = 0
        while True:
            start = Genome.find(Pattern, start)
            if start == -1: return
            yield start
            start += 1 # use start += len(Pattern) to remove overlapping matches
    positions = list(find_all(Pattern, Genome))
    
    return positions

print(PatternMatching('ATAT', 'GATATATGCATATACTT')) 

#### Find frequency of a certain nucleotide in a window
#### of half the genome's size

def PatternCount(Pattern, Text):
    count = 0
    for i in range(len(Text)-len(Pattern)+1):
        if Text[i:i+len(Pattern)] == Pattern:
            count = count+1
    return count 

def SymbolArray(Genome, symbol):
    array = {}
    n = len(Genome)
    ExtendedGenome = Genome + Genome[0:n//2]
    for i in range(n):
        array[i] = PatternCount(symbol, ExtendedGenome[i:i+(n//2)])
    return array

# Efficient version of the code on top

def PatternCount(Pattern, Text):
    count = 0
    for i in range(len(Text)-len(Pattern)+1):
        if Text[i:i+len(Pattern)] == Pattern:
            count = count+1
    return count 

def FasterSymbolArray(Genome, symbol):
    array = {}
    n = len(Genome)
    ExtendedGenome = Genome + Genome[0:n//2]

    # look at the first half of Genome to compute first array value
    array[0] = PatternCount(symbol, Genome[0:n//2])

    for i in range(1, n):
        # start by setting the current array value equal to the previous array value
        array[i] = array[i-1]

        # the current array value can differ from the previous array value by at most 1
        if ExtendedGenome[i-1] == symbol:
            array[i] = array[i]-1
        if ExtendedGenome[i+(n//2)-1] == symbol:
            array[i] = array[i]+1
    return array

#### Code to show a plot for a nucleotide's occurence in the genome with a window of length n//2

def PatternCount(Pattern, Text):
    count = 0
    for i in range(len(Text)-len(Pattern)+1):
        if Text[i:i+len(Pattern)] == Pattern:
            count = count+1
    return count 

def FasterSymbolArray(Genome, symbol):
    array = {}
    n = len(Genome)
    ExtendedGenome = Genome + Genome[0:n//2]

    # look at the first half of Genome to compute first array value
    array[0] = PatternCount(symbol, Genome[0:n//2])

    for i in range(1, n):
        # start by setting the current array value equal to the previous array value
        array[i] = array[i-1]

        # the current array value can differ from the previous array value by at most 1
        if ExtendedGenome[i-1] == symbol:
            array[i] = array[i]-1
        if ExtendedGenome[i+(n//2)-1] == symbol:
            array[i] = array[i]+1
    return array

## This is to read an externally saved genome
with open('e_coli.txt') as file:
    e_coli = file.read();

array = FasterSymbolArray(e_coli, "C")

import matplotlib.pyplot as plt
plt.plot(*zip(*sorted(array.items())))
plt.show()

### SkewArray Function 

def SkewArray(Genome):
    skew = [0]
    score = {"A":0, "T":0, "C":-1, "G":1}
    for i in range(1,len(Genome)+1):
            skew.append(score[Genome[i-1]] + skew[i-1])
    return skew


## Finding the ORI Region with the minimum values outputed

def SkewArray(Genome):
    array = [0]
    Skew = 0
    for i in Genome:
        if i == 'A' or i == 'T':
            Skew += 0
            array.append(Skew)
        if i == 'C':
            Skew -= 1
            array.append(Skew)
        if i == 'G':
            Skew += 1
            array.append(Skew)
    return array

def MinimumSkew(Genome):
    array = SkewArray(Genome)
    positions = []
    count = 0
    minarray = min(array)
    for i in array:
        if i == minarray:
            positions.append(count)
        count +=1
    return positions

### Find distance between two strings

def HammingDistance (p, q):
    count = 0
    for i, j in zip(p, q):
       if i != j:
           count += 1
    return count


### AproximatePatternMatching

def ApproximatePatternMatching(Text, Pattern, d):
    positions = [] # initializing list of positions
    for i in range(len(Text)-len(Pattern)+1):
        if HammingDistance(Text[i:i+len(Pattern)], Pattern) <= d:
            positions.append(i)
    return positions
#range is modified like this because if is it just length of the bigger text, pattern will keep sliding along with empty letters, adding more to the list of positions

def HammingDistance(p, q):
    count = 0
    for i in range(len(p)):
        if p[i] != q[i]:
            count += 1
    return count

### Approximate Pattern Count

def ApproximatePatternCount(Pattern, Text, d):
    count = 0 # initialize count variable
    for i in range(len(Text)-len(Pattern)+1):
        if HammingDistance(Pattern, Text[i:i+len(Pattern)]) <= d:
            count += 1
    return count
#same thing as the pattern matching before, but it is replaced by count

def HammingDistance(p, q):
    count = 0
    for i in range(len(p)):
        if p[i] != q[i]:
            count += 1
    return count

### DO NOT MODIFY THE CODE BELOW THIS LINE ###
import sys
lines = sys.stdin.read().splitlines()
print(ApproximatePatternCount(lines[0],lines[1],int(lines[2])))

#here is my pattern matching for reference. range is set like that because it will continue matching with the empty letters at the end of the long text

def ApproximatePatternMatching(Text, Pattern, d):
    positions = [] # initializing list of positions
    for i in range(len(Text)-len(Pattern)+1):
        if HammingDistance(Text[i:i+len(Pattern)], Pattern) <= d:
            positions.append(i)
    return positions

### Frequent Words with Mismatches Problem
    Possible Solution:
        Use the sliding window technique and slipt each sliding window into single characters,
        then she how often the certain characters are in each window.


##### Motifs.py #####
### Count
def Count(Motifs):
    count = {}
    ## k - Length of kmer/window
    k = len(Motifs[0])
    for symbol in "ACGT":
        count[symbol] = []
        for j in range(k):
            count[symbol].append(0)
    
    t = len(Motifs)
    for i in range(t):
        for j in range(k):
            symbol = Motifs[i][j]
            count[symbol][j] +=1

    return count

### Profile (requires Count function)

def Count(Motifs):
    count = {}
    ## k - Length of kmer/window
    k = len(Motifs[0])
    for symbol in "ACGT":
        count[symbol] = []
        for j in range(k):
            count[symbol].append(0)
    
    t = len(Motifs)
    for i in range(t):
        for j in range(k):
            symbol = Motifs[i][j]
            count[symbol][j] +=1

    return count


def Profile(Motifs):
    t = len(Motifs)
    k = len(Motifs[0])
    profile = Count(Motifs)
    
    for key,v in profile.items():
        v[:] = [x / t for x in v]
    return profile

### Consensus
def Count(Motifs):
    count = {}
    ## k - Length of kmer/window
    k = len(Motifs[0])
    for symbol in "ACGT":
        count[symbol] = []
        for j in range(k):
            count[symbol].append(0)
    
    t = len(Motifs)
    for i in range(t):
        for j in range(k):
            symbol = Motifs[i][j]
            count[symbol][j] +=1

    return count

def Consensus(Motifs):
    k = len(Motifs[0])
    count = Count(Motifs)

    consensus = ''

    for j in range(k):
        m = 0
        frequentSymbol = ''
        for symbol in 'ACGT':
            if count[symbol][j] > m:
                m = count[symbol][j]
                frequentSymbol = symbol
        consensus += frequentSymbol
    return consensus

### Score
def Count(Motifs):
    count = {}
    ## k - Length of kmer/window
    k = len(Motifs[0])
    for symbol in "ACGT":
        count[symbol] = []
        for j in range(k):
            count[symbol].append(0)
    
    t = len(Motifs)
    for i in range(t):
        for j in range(k):
            symbol = Motifs[i][j]
            count[symbol][j] +=1

    return count

def Consensus(Motifs):
    k = len(Motifs[0])
    count = Count(Motifs)

    consensus = ''

    for j in range(k):
        m = 0
        frequentSymbol = ''
        for symbol in 'ACGT':
            if count[symbol][j] > m:
                m = count[symbol][j]
                frequentSymbol = symbol
        consensus += frequentSymbol
    return consensus

def Score(Motifs):
    consensus = Consensus(Motifs)
    count = 0
    for motif in Motifs:
        for index in range(len(motif)):
            if motif[index] != consensus[index]:
                count += 1
    return count

### Probability

def Pr(Text, Profile):
    product = 1
    for index, nucleotide in enumerate(Text):
        product *= Profile[nucleotide][index]
    return product

### 
def ProfileMostProbableKmer(text, k, profile):
    n = len(text)
    pr = {}
    most_likely_kmer = []
    for i in range(n-k+1):
        k_mer = text[i:i+k]
        probability = Pr(k_mer, profile)
        pr[k_mer] = probability
    m = max(pr.values())
    for key, value in pr.items():
        if pr[key] == m:
            most_likely_kmer.append(key)
    return most_likely_kmer[0]

def Pr(Text, Profile):
    product = 1
    for index, nucleotide in enumerate(Text):
        product *= Profile[nucleotide][index]
    return product

### Greedy Motif Search
def Count(Motifs):
    k = len(Motifs[0])
    count = {}
    for symbol in "ACGT":
        count[symbol] = []
        for j in range(k):
            count[symbol].append(0)
    t = len(Motifs)
    for i in range(t):
        for j in range(k):
            symbol = Motifs[i][j]
            count[symbol][j] += 1
    return count


def Consensus(Motifs):
    count = Count(Motifs)
    k = len(Motifs[0])
    consensus = ""
    for j in range(k):
        m = 0
        frequentSymbol = ""
        for symbol in "ACGT":
            if count[symbol][j] > m:
                m = count[symbol][j]
                frequentSymbol = symbol
        consensus += frequentSymbol
    return consensus


def Profile(Motifs):
    t = len(Motifs)
    k = len(Motifs[0])
    profile = Count(Motifs)
    for i in 'ACTG':
        for j in range(k):
            profile[i][j] = profile[i][j]/t  
    return profile

def Score(Motifs):
    k = len(Motifs[0])
    score = 0
    count = Count(Motifs)
    max_symbol = Consensus(Motifs)
    sum1 = 0
    for j in range(k):
        m = 0
        for symbol in "ATCG":
            if count[symbol][j] > m:
                sum1 += count[symbol][j]
    for j in range(k):
        m = 0
        for symbol in "AGTC":
            if count[symbol][j] > m:
                m = count[symbol][j]
        score += m  
    return sum1-score




def Pr(Text, Profile):
    k = len(Profile["A"])
    p=1
    for i in range(len(Text)):
        p=p*Profile[Text[i]][i]
    return p

def ProfileMostProbablePattern(text,k,profile):
    p=-1
    result=text[0:k]
    for i in range(len(text)-k+1):
        seq=text[i:i+k]
        pr=Pr(seq,profile)
        if pr>p:
            p=pr
            result=seq
    return result

def GreedyMotifSearch(Dna,k,t):
    BestMotifs = []
    for i in range(0, t):
        BestMotifs.append(Dna[i][0:k])
    n = len(Dna[0])
    for m in range(n-k+1):
        Motifs = []
        Motifs.append(Dna[0][m:m+k])
        for j in range(1, t):
            P = Profile(Motifs[0:j])
            Motifs.append(ProfileMostProbablePattern(Dna[j], k, P))
        if Score(Motifs) < Score(BestMotifs):
            BestMotifs = Motifs
    return BestMotifs

### Big DNA set

def Count(Motifs):

    count = {}

    k = len(Motifs[0])

    for symbol in "ACGT":

        count[symbol] = []

        for j in range(k):

            count[symbol].append(0)

    t = len(Motifs)

    for i in range(t):

        for j in range(k):

            symbol = Motifs[i][j]

            count[symbol][j] += 1

    return count

 

def Consensus(Motifs):

  

    k = len(Motifs[0])

    count = Count(Motifs)

    consensus = ""

    for j in range(k):

        m = 0

        frequentSymbol = ""

        for symbol in "ACGT":

            if count[symbol][j] > m:

                m = count[symbol][j]

                frequentSymbol = symbol

        consensus += frequentSymbol

    return consensus

 

def Profile(Motifs):

    t = len(Motifs)

    k = len(Motifs[0])

    profile = Count(Motifs)

    for i in 'ACTG':

        for j in range(k):

            profile[i][j] = profile[i][j]/t  

    return profile

 

def Score(Motifs):

    # Insert code here

    score = 0

    k = len(Motifs[0])

    count = Count(Motifs)

    max_symbol = Consensus(Motifs)

    sum1 = 0

    for j in range(k):

        m = 0

        for symbol in "ATCG":

            if count[symbol][j] > m:

                sum1 += count[symbol][j]

    for j in range(k):

        m = 0

        for symbol in "AGTC":

            if count[symbol][j] > m:

                m = count[symbol][j]

        score += m  

    return sum1-score

 

def Pr(Text, Profile):

    p=1

    k = len(Profile["A"])

    for i in range(len(Text)):

        p=p*Profile[Text[i]][i]

    return p

 

#Finally solved it

def ProfileMostProbablePattern(text,k,profile):

    p=-1

    result=text[0:k]

    for i in range(len(text)-k+1):

        seq=text[i:i+k]

        pr=Pr(seq,profile)

        if pr>p:

            p=pr

            result=seq

    return result

 

def GreedyMotifSearch(Dna,k,t):

    BestMotifs = []

    for i in range(0, t):

        BestMotifs.append(Dna[i][0:k])

    n = len(Dna[0])

    for m in range(n-k+1):

        Motifs = []

        Motifs.append(Dna[0][m:m+k])

        for j in range(1, t):

            P = Profile(Motifs[0:j])

            Motifs.append(ProfileMostProbablePattern(Dna[j], k, P))

        if Score(Motifs) < Score(BestMotifs):

            BestMotifs = Motifs

    return BestMotifs


### Count own often a letter appears in their "vertical position" + 1
def CountWithPseudocounts(Motifs):
    count = {}
    aa = []
    cc = []
    gg = []
    tt = []
    z = 0 
    while z <= (len(Motifs[0])-1):
        somaa = 0
        for i in Motifs:
            if i[z] == 'A':
                somaa += 1
        somaa = somaa + 1
        aa.append(somaa)
    
        somac = 0
        for i in Motifs:
            if i[z] == 'C':
                somac += 1
        somac = somac + 1
        cc.append(somac)

        somag = 0
        for i in Motifs:
            if i[z] == 'G':
                somag += 1
        somag = somag + 1
        gg.append(somag)

        somat = 0
        for i in Motifs:
            if i[z] == 'T':
                somat += 1
        somat = somat + 1
        tt.append(somat)

    
        z += 1
    count['A'] = aa       
    count['C'] = cc
    count['G'] = gg
    count['T'] = tt
    return count

#ProfileWithPseudoCounts

def ProfileWithPseudocounts(Motifs):
    t = len(Motifs) + 4
    k = len(Motifs[0])
    profile = CountWithPseudocounts(Motifs)
    
    for key,v in profile.items():
        v[:] = [x / t for x in v]
    return profile


# Input:  A set of kmers Motifs
# Output: CountWithPseudocounts(Motifs)
# HINT:   You need to use CountWithPseudocounts as a subroutine of ProfileWithPseudocounts
def CountWithPseudocounts(Motifs):
    count = {}
    aa = []
    cc = []
    gg = []
    tt = []
    z = 0 
    while z <= (len(Motifs[0])-1):
        somaa = 0
        for i in Motifs:
            if i[z] == 'A':
                somaa += 1
        somaa = somaa + 1
        aa.append(somaa)
    
        somac = 0
        for i in Motifs:
            if i[z] == 'C':
                somac += 1
        somac = somac + 1
        cc.append(somac)

        somag = 0
        for i in Motifs:
            if i[z] == 'G':
                somag += 1
        somag = somag + 1
        gg.append(somag)

        somat = 0
        for i in Motifs:
            if i[z] == 'T':
                somat += 1
        somat = somat + 1
        tt.append(somat)

    
        z += 1
    count['A'] = aa       
    count['C'] = cc
    count['G'] = gg
    count['T'] = tt
    return count


### GreedyMotifSearchWithPseudocounts 

def GreedyMotifSearchWithPseudocounts(Dna,k,t):
    BestMotifs = []
    for i in range(0, t):
        BestMotifs.append(Dna[i][0:k])
    n = len(Dna[0])
    for m in range(n-k+1):
        Motifs = []
        Motifs.append(Dna[0][m:m+k])
        for j in range(1, t):
            P = ProfileWithPseudocounts(Motifs[0:j])
            Motifs.append(ProfileMostProbablePattern(Dna[j], k, P))
        if Score(Motifs) < Score(BestMotifs):
            BestMotifs = Motifs
    return BestMotifs
       
                    
                    
                    
def ProfileWithPseudocounts(Motifs):
    profile = {}
    k = len(Motifs[0])
    t = len(Motifs)
    profile = CountWithPseudocounts(Motifs)
    
    for i in range(k):
        p = 0
        
        for symbol in "ACGT":
            p = p+profile[symbol][i]
        for symbol in "ACGT":
            profile[symbol][i] = profile[symbol][i]/p
    return profile

def Score(Motifs):
    consensus = Consensus(Motifs)
    count = CountWithPseudocounts(Motifs)
    
    k = len(Motifs[0])
    t = len(Motifs)  
    c = 0
    for Motif in Motifs:
        for i in range(k):
                if Motif[i] != consensus[i]:
                    c = c+1
    return c


def Consensus(Motifs):
    consensus = ""
    k = len(Motifs[0])
    t = len(Motifs)
    count = CountWithPseudocounts(Motifs)
    for j in range(k):
        m = 0
        frequentSymbol = ""
        for symbol in "ACGT":
            if count[symbol][j] > m:
                m = count[symbol][j]
                frequentSymbol = symbol
        consensus += frequentSymbol
    return consensus


def CountWithPseudocounts(Motifs):
    count = {}
    k = len(Motifs[0])
    for symbol in "ACGT":
        count[symbol] = []
        for j in range(k):
             count[symbol].append(1)
    t = len(Motifs)
    for i in range(t):
        for j in range(k):
            symbol = Motifs[i][j]
            count[symbol][j] = count[symbol][j]+1
    return count

def ProfileMostProbablePattern(Text, k, profile):
    n = len(Text)
    m = -1 #This is the bug - need to adjust it to -1 from 0
    x = Text[1:k]
    for i in range(n-k+1):
        Pattern =  Text[i:i+k]
        p = Pr(Pattern, profile)
        if p>m:
            m = p
            x = Pattern
            
    return x


def Pr(Text, Profile):
   
    p = 1
    for i in range(len(Text)):
        
        p1 = Profile[(Text[i])][i]
        p = p*p1
        
    return p


### GreedyMotifSearch with pseudocounts

def Count(Motifs):

    count = {}

    k = len(Motifs[0])

    for symbol in "ACGT":

        count[symbol] = []

        for j in range(k):

            count[symbol].append(1)

    t = len(Motifs)

    for i in range(t):

        for j in range(k):

            symbol = Motifs[i][j]

            count[symbol][j] += 1

    return count

def Consensus(Motifs):

    k = len(Motifs[0])

    count = Count(Motifs)

    consensus = ""

    for j in range(k):

        m = 0

        frequentSymbol = ""

        for symbol in "ACGT":

            if count[symbol][j] > m:

                m = count[symbol][j]

                frequentSymbol = symbol

        consensus += frequentSymbol

    return consensus

def Profile(Motifs):

    t = len(Motifs)

    k = len(Motifs[0])

    profile = Count(Motifs)

    for i in 'ACTG':

        for j in range(k):

            profile[i][j] = profile[i][j]/t  

    return profile

def Score(Motifs):

    # Insert code here

    score = 0

    k = len(Motifs[0])

    count = Count(Motifs)

    max_symbol = Consensus(Motifs)

    sum1 = 0

    for j in range(k):

        m = 0

        for symbol in "ATCG":

            if count[symbol][j] > m:

                sum1 += count[symbol][j]

    for j in range(k):

        m = 0

        for symbol in "AGTC":

            if count[symbol][j] > m:

                m = count[symbol][j]

        score += m  

    return sum1-score

def Pr(Text, Profile):

    p=1

    k = len(Profile["A"])

    for i in range(len(Text)):

        p=p*Profile[Text[i]][i]

    return p

def ProfileMostProbablePattern(text,k,profile):

    p=-1

    result=text[0:k]

    for i in range(len(text)-k+1):

        seq=text[i:i+k]

        pr=Pr(seq,profile)

        if pr>p:

            p=pr

            result=seq

    return result

def GreedyMotifSearch(Dna,k,t):

    BestMotifs = []

    for i in range(0, t):

        BestMotifs.append(Dna[i][0:k])

    n = len(Dna[0])

    for m in range(n-k+1):

        Motifs = []

        Motifs.append(Dna[0][m:m+k])

        for j in range(1, t):

            P = Profile(Motifs[0:j])

            Motifs.append(ProfileMostProbablePattern(Dna[j], k, P))

        if Score(Motifs) < Score(BestMotifs):

            BestMotifs = Motifs

    return BestMotifs

k = 15

t = 10

Motifs = GreedyMotifSearch(Dna, k, t)

## ADD Dna VARIABLE

print(Motifs)


print(Score(Motifs))

## RandomMotifs

import random

def RandomMotifs(Dna, k, t):

    t = len(Dna)
    l = len(Dna[0])
    RandomMotif =[]
    for i in range(t):
        r = random.randint(1,l-k) # 1 is not added as it is inclusive of last element also
        RandomMotif.append(Dna[i][r:r+k])
    return RandomMotif


### RandomMotifSearch

import random

def randomMotifs(dna,k,t):

    kmm = []
    sc = []
    k = 3
    D = {}
    for i in range(0,len(dna)):
        km = []
        for kk in range(len(dna[i])-k+1):
            km += [dna[i][kk:kk+k]]
        D[i] = km
    for m in range(0,t):
        ran = random.randint(0,len(D[0])-1)
        kmm += [D[m][ran]]

    return kmm


def ProfileWithPseudocounts(Motifs):


    t = len(Motifs)
    k = len(Motifs[0])
    profile = CountWithPseudocounts(Motifs) # output variable
    for symbol in profile:
        for kk in range(0,len(profile[symbol])):
            profile[symbol][kk] = profile[symbol][kk]/(len(Motifs) + 4)

    return profile


def CountWithPseudocounts(Motifs):

    count = {}
    for i in 'ACGT':
        count[i] = []
        for ii in range(len(Motifs[0])):
            count[i].append(1)
    for i in range(len(Motifs)):
        for j in range(len(Motifs[0])):
            symbol = Motifs[i][j]
            count[symbol][j] += 1

    return count



def Score(Motifs):


    count = 0
    L = Consensus(Motifs)
    for i in Motifs:
        for chr1, chr2 in zip(i,L):
            if chr1 != chr2:
                count += 1
    return count


def Consensus(Motifs):



    k = len(Motifs[0])
    count = Count(Motifs)
    consensus = ""
    for j in range(k):
        m = 0
        frequentSymbol = ""
        for symbol in "ACGT":
            if count[symbol][j] > m:
                m = count[symbol][j]
                frequentSymbol = symbol
        consensus += frequentSymbol

    return consensus


def Count(Motifs):

    count = {}
    for i in 'ACGT':
        count[i] = []
        for ii in range(len(Motifs[0])):
            count[i].append(0)
    for i in range(len(Motifs)):
        for j in range(len(Motifs[0])):
            symbol = Motifs[i][j]
            count[symbol][j] += 1

    return count


def RandomMotifs(dna,k,t):

    kmm = []
    sc = []
    D = {}
    for i in range(0,len(dna)):
        km = []
        for kk in range(len(dna[i])-k+1):
            km += [dna[i][kk:kk+k]]
        D[i] = km
    for m in range(0,t):
        ran = random.randint(0,len(D[0])-1)
        kmm += [D[m][ran]]

    return kmm

def Motifs(pf,dna):

    k = len(pf['A'])
    D = []
    for i in range(0,len(dna)):
        km = []
        sc = []
        for kk in range(len(dna[i])-k+1):
            km += [dna[i][kk:kk+k]]
        for i in km:
            sc += [Pr(i,pf)]
        D += [km[sc.index(max(sc))]]

    return D


def Pr(Text, Profile):

    p = 1
    for i in range(0,len(Text)):
        p *= Profile[Text[i]][i]

    return p


def RandomizedMotifSearch(Dna, k, t):

    M = RandomMotifs(Dna, k, t)
    BestMotifs = M

    while True:
        Profile = ProfileWithPseudocounts(M)
        M = Motifs(Profile, Dna)
        if Score(M) < Score(BestMotifs):
            BestMotifs = M
        else:
            return BestMotifs


### Normalize Values

def Normalize(P):

    d = {}
    for k,v in P.items():
        d[k] = P[k]/sum(P.values())
    return d

d = {'A': 0.1, 'C': 0.1, 'G': 0.1, 'T': 0.1}


### WeightedDie

def WeightedDie(Probabilities):
    n = random.uniform(0, 1)
    for p in Probabilities:
        n -= Probabilities[p]
        if n <= 0:
            return p


###

import random
from operator import itemgetter
# then, copy Pr, Normalize, and WeightedDie below this line
def WeightedDie(d):
    ran = random.uniform(0, 1)
    #print(ran,d)
    tot = 0
    for k, v in sorted(d.items(),key=itemgetter(1)):
        if tot <= ran < v + tot:
            return k
        tot += v
def Normalize(P):
    D = {}
    for k,v in P.items():
        D[k] = P[k]/sum(P.values())
    return D
def Pr(Text, Profile):
    p = 1
    for i in range(0,len(Text)):
        p *= Profile[Text[i]][i]
    return p
# Input:  A string Text, a profile matrix Profile, and an integer k
# Output: ProfileGeneratedString(Text, profile, k)
def ProfileGeneratedString(Text, profile, k):
    # your code here
    n = len(Text)
    probabilities = {}
    for i in range(0,n-k+1):
        probabilities[Text[i:i+k]] = Pr(Text[i:i+k], profile)

    probabilities = Normalize(probabilities)
    return WeightedDie(probabilities)


### GibbsSampler 

import random
# Input:  Integers k, t, and N, followed by a collection of strings Dna
# Output: GibbsSampler(Dna, k, t, N)
def GibbsSampler(Dna, k, t, N):
    BestMotifs = [] # output variable
    # your code here
    Motifs = RandomMotifs(Dna, k, t)
    BestMotifs = Motifs
    for j in range(1,N):
        i = random.randint(0,t-1)
        ReducedMotifs = []
        for j in range(0,t):
            if j != i:
                ReducedMotifs.append(Motifs[j])
        Profile = ProfileWithPseudocounts(ReducedMotifs)
        Motif_i = ProfileGeneratedString(Dna[i], Profile, k)
        Motifs[i] = Motif_i
        if Score(Motifs) < Score(BestMotifs):
                BestMotifs=Motifs
    return BestMotifs
# place all subroutines needed for GibbsSampler below this line
# Input:  A list of strings Dna, and integers k and t
# Output: RandomMotifs(Dna, k, t)
# HINT:   You might not actually need to use t since t = len(Dna), but you may find it convenient
def RandomMotifs(Dna, k, t):
    # place your code here.
    s = len(Dna[0])
    rm = []
    for i in range(0,t):
        init_index = random.randint(1,s-k)
        rm.append(Dna[i][init_index:init_index+k])    
    return rm
# Input:  A set of kmers Motifs
# Output: ProfileWithPseudocounts(Motifs)
def ProfileWithPseudocounts(Motifs):
    t = len(Motifs)
    k = len(Motifs[0])
    profile = {} # output variable
    # your code here
    c = CountWithPseudocounts(Motifs)
    for n in 'ACGT':
        p = []
        for i in range(0,k):
            p.append(c[n][i]/(t+4))
        profile[n] = p
    return profile
# Input:  A set of kmers Motifs
# Output: CountWithPseudocounts(Motifs)
def CountWithPseudocounts(Motifs):
    t = len(Motifs)
    k = len(Motifs[0])
    # insert your code here
    count = {} # initializing the count dictionary
    for symbol in "ACGT":
        count[symbol] = []
        for j in range(k):
             count[symbol].append(1)
    for i in range(t):
        for j in range(k):
             symbol = Motifs[i][j]
             count[symbol][j] += 1
    return count 
#tests in which of the intervals defined by list ar the number r lies
def testinterval(ar,r):
    ar.sort()
    if r<= ar[0]:
      return ar[0]
    for i in range(1,len(ar)-1):
      if ar[i-1]<r<=ar[i]:
        return ar[i]
    if ar[len(ar)-2]< r:
      return ar[len(ar)-1]
# Input:  A dictionary Probabilities whose keys are k-mers and whose values are the probabilities of these kmers
# Output: A randomly chosen k-mer with respect to the values in Probabilities
def WeightedDie(Probabilities):
    # your code here
    sumprob = {}
    s = 0
    for p in Probabilities:
        s += Probabilities[p]
        sumprob[p] = s
    revprob = {}
    for q in sumprob:
      revprob[sumprob[q]] = q
    w = list(sumprob.values())
    r = random.uniform(0,1)
    kmer = revprob[testinterval(w,r)]
    return kmer
# Input:  A string Text, a profile matrix Profile, and an integer k
# Output: ProfileGeneratedString(Text, profile, k)
def ProfileGeneratedString(Text, profile, k):
    # your code here
    n = len(Text)
    probabilities = {} 
    for i in range(0,n-k+1):
        probabilities[Text[i:i+k]] = Pr(Text[i:i+k], profile)
    probabilities = Normalize(probabilities)
    return WeightedDie(probabilities)
# Input:  String Text and profile matrix Profile
# Output: Pr(Text, Profile)
def Pr(Text, Profile):
    # insert your code here
    p = 1
    for i in range(0,len(Text)):
        p *= Primport random

def randomMotifs(dna,k,t):

    kmm = []
    sc = []
    k = 3
    D = {}
    for i in range(0,len(dna)):
        km = []
        for kk in range(len(dna[i])-k+1):
            km += [dna[i][kk:kk+k]]
        D[i] = km
    for m in range(0,t):
        ran = random.randint(0,len(D[0])-1)
        kmm += [D[m][ran]]

    return kmm


def ProfileWithPseudocounts(Motifs):


    t = len(Motifs)
    k = len(Motifs[0])
    profile = CountWithPseudocounts(Motifs) # output variable
    for symbol in profile:
        for kk in range(0,len(profile[symbol])):
            profile[symbol][kk] = profile[symbol][kk]/(len(Motifs) + 4)

    return profile


def CountWithPseudocounts(Motifs):

    count = {}
    for i in 'ACGT':
        count[i] = []
        for ii in range(len(Motifs[0])):
            count[i].append(1)
    for i in range(len(Motifs)):
        for j in range(len(Motifs[0])):
            symbol = Motifs[i][j]
            count[symbol][j] += 1

    return count



def Score(Motifs):


    count = 0
    L = Consensus(Motifs)
    for i in Motifs:
        for chr1, chr2 in zip(i,L):
            if chr1 != chr2:
                count += 1
    return count


def Consensus(Motifs):



    k = len(Motifs[0])
    count = Count(Motifs)
    consensus = ""
    for j in range(k):
        m = 0
        frequentSymbol = ""
        for symbol in "ACGT":
            if count[symbol][j] > m:
                m = count[symbol][j]
                frequentSymbol = symbol
        consensus += frequentSymbol

    return consensus


def Count(Motifs):

    count = {}
    for i in 'ACGT':
        count[i] = []
        for ii in range(len(Motifs[0])):
            count[i].append(0)
    for i in range(len(Motifs)):
        for j in range(len(Motifs[0])):
            symbol = Motifs[i][j]
            count[symbol][j] += 1

    return count


def RandomMotifs(dna,k,t):

    kmm = []
    sc = []
    D = {}
    for i in range(0,len(dna)):
        km = []
        for kk in range(len(dna[i])-k+1):
            km += [dna[i][kk:kk+k]]
        D[i] = km
    for m in range(0,t):
        ran = random.randint(0,len(D[0])-1)
        kmm += [D[m][ran]]

    return kmm

def Motifs(pf,dna):

    k = len(pf['A'])
    D = []
    for i in range(0,len(dna)):
        km = []
        sc = []
        for kk in range(len(dna[i])-k+1):
            km += [dna[i][kk:kk+k]]
        for i in km:
            sc += [Pr(i,pf)]
        D += [km[sc.index(max(sc))]]

    return D


def Pr(Text, Profile):

    p = 1
    for i in range(0,len(Text)):
        p *= Profile[Text[i]][i]

    return p


def RandomizedMotifSearch(Dna, k, t):

    M = RandomMotifs(Dna, k, t)
    BestMotifs = M

    while True:
        Profile = ProfileWithPseudocounts(M)
        M = Motifs(Profile, Dna)
        if Score(M) < Score(BestMotifs):
            BestMotifs = M
        else:
            return BestMotifsionary Probabilities, where keys are k-mers and values are the probabilities of these k-mers (which do not necessarily sum up to 1)
# Output: A normalized dictionary where the probability of each k-mer was divided by the sum of all k-mers' probabilities
def Normalize(Probabilities):
    # your code here
    result = {}
    sum = 0
    for m in Probabilities:
        sum += Probabilities[m]
    for n in Probabilities:
        result[n]= Probabilities[n]/sum
    return result  
# Input:  A set of k-mers Motifs
# Output: The score of these k-mers.
def Score(Motifs):
    # Insert code here
    k = len(Motifs[0])
    t = len(Motifs)
    cs = ConsensusWithPseudocounts(Motifs)
    score = 0
    for j in range(0,k):
        for i in range(0,t):
            if Motifs[i][j] != cs[j]:
                score += 1
    return score
# Input:  A set of kmers Motifs
# Output: A consensus string of Motifs.
def ConsensusWithPseudocounts(Motifs):
    # insert your code here
    k = len(Motifs[0])
    count = CountWithPseudocounts(Motifs)
    consensus = ""
    for j in range(k):
        m = 0
        frequentSymbol = ""
        for symbol in "ACGT":
            if count[symbol][j] > m:
                m = count[symbol][j]
                frequentSymbol = symbol
        consensus += frequentSymbol
    return consimport random

def randomMotifs(dna,k,t):

    kmm = []
    sc = []
    k = 3
    D = {}
    for i in range(0,len(dna)):
        km = []
        for kk in range(len(dna[i])-k+1):
            km += [dna[i][kk:kk+k]]
        D[i] = km
    for m in range(0,t):
        ran = random.randint(0,len(D[0])-1)
        kmm += [D[m][ran]]

    return kmm


def ProfileWithPseudocounts(Motifs):


    t = len(Motifs)
    k = len(Motifs[0])
    profile = CountWithPseudocounts(Motifs) # output variable
    for symbol in profile:
        for kk in range(0,len(profile[symbol])):
            profile[symbol][kk] = profile[symbol][kk]/(len(Motifs) + 4)

    return profile


def CountWithPseudocounts(Motifs):

    count = {}
    for i in 'ACGT':
        count[i] = []
        for ii in range(len(Motifs[0])):
            count[i].append(1)
    for i in range(len(Motifs)):
        for j in range(len(Motifs[0])):
            symbol = Motifs[i][j]
            count[symbol][j] += 1

    return count



def Score(Motifs):


    count = 0
    L = Consensus(Motifs)
    for i in Motifs:
        for chr1, chr2 in zip(i,L):
            if chr1 != chr2:
                count += 1
    return count


def Consensus(Motifs):



    k = len(Motifs[0])
    count = Count(Motifs)
    consensus = ""
    for j in range(k):
        m = 0
        frequentSymbol = ""
        for symbol in "ACGT":
            if count[symbol][j] > m:
                m = count[symbol][j]
                frequentSymbol = symbol
        consensus += frequentSymbol

    return consensus


def Count(Motifs):

    count = {}
    for i in 'ACGT':
        count[i] = []
        for ii in range(len(Motifs[0])):
            count[i].append(0)
    for i in range(len(Motifs)):
        for j in range(len(Motifs[0])):
            symbol = Motifs[i][j]
            count[symbol][j] += 1

    return count


def RandomMotifs(dna,k,t):

    kmm = []
    sc = []
    D = {}
    for i in range(0,len(dna)):
        km = []
        for kk in range(len(dna[i])-k+1):
            km += [dna[i][kk:kk+k]]
        D[i] = km
    for m in range(0,t):
        ran = random.randint(0,len(D[0])-1)
        kmm += [D[m][ran]]

    return kmm

def Motifs(pf,dna):

    k = len(pf['A'])
    D = []
    for i in range(0,len(dna)):
        km = []
        sc = []
        for kk in range(len(dna[i])-k+1):
            km += [dna[i][kk:kk+k]]
        for i in km:
            sc += [Pr(i,pf)]
        D += [km[sc.index(max(sc))]]

    return D


def Pr(Text, Profile):

    p = 1
    for i in range(0,len(Text)):
        p *= Profile[Text[i]][i]

    return p


def RandomizedMotifSearch(Dna, k, t):

    M = RandomMotifs(Dna, k, t)
    BestMotifs = M

    while True:
        Profile = ProfileWithPseudocounts(M)
        M = Motifs(Profile, Dna)
        if Score(M) < Score(BestMotifs):
            BestMotifs = M
        else:
            return BestMotifsensus    

k = 8
t = 5
N = 100
Dna = ["CGCCCCTCTCGGGGGTGTTCAGTAAACGGCCA","GGGCGAGGTATGTGTAAGTGCCAAGGTGCCAG","TAGTACCGAGACCGAAAGAAGTATACAGGCGT","TAGATCAAGTTTCAGGTGCACGTCGGTGAACC","AATCCACCAGCTCCACGTGCAATGTTGGCCTA"]

'''

# Run Code Here #

# End Code #