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

#### Find frequency of a certain nucleotoide in a window
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

#### Code to show a plot for a nucleotoide's occurence in the genome with a window of length n//2

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


'''

# Run Code Here #

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

print(MinimumSkew('GATACACTTCCCGAGTAGGTACTG'))

# End Code #