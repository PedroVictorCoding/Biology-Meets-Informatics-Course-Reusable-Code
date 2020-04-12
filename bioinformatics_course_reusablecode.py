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

def SymbolArray(Genome, symbol):
    array = []
    n = len(Genome)
    ExtendedGenome = Genome + Genome[0:n//2]
    for i in range(n):
        count = PatternCount(ExtendedGenome[i:i+(n//2)], symbol)
        array.append(count)
    return array


'''

# Run Code Here #

def PatternCount(Pattern, Text):
    count = 0
    for i in range(len(Text)-len(Pattern)+1):
        if Text[i:i+len(Pattern)] == Pattern:
            count = count+1
    return count 

def SymbolArray(Genome, symbol):
    array = []
    n = len(Genome)
    ExtendedGenome = Genome + Genome[0:n//2]
    for i in range(n):
        count = PatternCount(ExtendedGenome[i:i+(n//2)], symbol)
        array.append(count)
    return array

# End Code #