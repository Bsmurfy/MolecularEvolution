'''
Brendan Murphy

A program to simulate the molecular evolution of the huntingtin protein, with
the ultimate goal of determining how long it would take for half of a population
to display Huntington's disease.
'''

import csv
import random
from datetime import datetime

# define dictionaries for transitions and transversions
transitions = {'A':'G', 'G':'A', 'C':'T', 'T':'C'};
transversions = {'A':random.choice('CT'), 'G':random.choice('CT'),
                 'C':random.choice('AG'), 'T':random.choice('AG')};

def main():
    run = 0
    for run in range(100):
        print('') # makes the output more reader-friendly
        print('Run ' + str(run + 1) + ':')
        bulk()
        run += 1

def bulk():
    # file containing generation count over multiple runs of the program
    file = open('generationcount.csv', 'a', newline='\n')
    csvwriter = csv.writer(file)
    
    n = 1000 # population size

    # transitions are (A<->G) and (C<->T); rate of 10^-7
    alpha = 0.0000001
    # transversions are (A<->C), (A<->T), (G<->C), and (G<->T); rate of 10^-8
    beta = 0.00000001
    # CAG repeat extension; rate of 10^-4
    trirate = 0.0001

    startTime = datetime.now() # set start time

    # we are going to trim down the gene significantly to help with run time
    beforeCAG = 'AAGTCCTTC' # three codons before the repeat section
    afterCAG = 'CAACAGCCGCCACCGCCGCCGCCGCCGCCGCCGCCTCCTCAGCTTCCTCAGCCGCCGCCGCAGGCACAGCCGCTGCTGCCTCAGCCGCAGCCGCCCCCGCCGCCGCCC'

    # a normal CAG-repeat count is 16-21
    minNormalRepeat = 16
    maxNormalRepeat = 21
    # to display Huntington's, the CAG-repeat count must be 40
    displayRepeat = 40

    # prep to build population
    population = []

    # one individual will have Huntington's at the beginning
    initial = ''
    initial += beforeCAG
    initial += 'CAG' * 40 # 40 repeats required to display Huntington's
    initial += afterCAG
    population.append(initial)

    # the rest of the individuals in the starting population have normal values
    for individual in range(0, n-1):
        initial = ''
        initial += beforeCAG
        # add a number of CAG repeats within the accepted normal range
        initial += 'CAG' * random.randint(minNormalRepeat, maxNormalRepeat)
        initial += afterCAG
        population.append(initial)

    huntCount = 1 # number of individuals with Huntington's
    threshold = n/2 # threshold is half the population size
    genCount = 0 # number of generations simulated
    
    # simulate molecular evolution until threshold is met
    # testcount = 0
    while (huntCount < threshold):
        population = reproduce(population, n)
        population = mutate(population, alpha, beta, trirate)
        genCount += 1
        # check number of indviduals with Huntington's (aka 40+ CAG repeats)
        huntCount = countDisplays(population)
        # testcount += 1
        # print('test: ' + str(testcount) + '   hunt: ' + str(huntCount))

    # add generation when threshold is exceeded to the csv file
    csvwriter.writerow([str(genCount)])

    # Print generation number when threshold is exceeded and the program's run time
    print('By generation ' + str(genCount) + ', half the population will have Huntington\'s.')
    print('Done in: ' + str(datetime.now() - startTime))

    file.close() # close the opened file

'''
mutate - mutates the population according to the Kimura-2-parameter model
param pop: the population at a given generation
param a: alpha, the mutation rate of transitions
param b: beta, the mutation rate of transversions
param tri, the rate of CAG-repeat extensions
return newpop: the population after mutation has occurred
'''
def mutate(pop, a, b, tri):
    newpop = []
    for individual in pop:
        newseq = ''
        for n in range(0, len(individual), 3):
            codon = individual[n:n+3]
            if (codon != 'CAG'): # no back-mutating for CAGs, for sake of simplicity
                for nucleotide in codon:
                    tocompare = random.random()
                    if (tocompare < a): # transition
                        newseq += transitions[nucleotide]
                    elif (tocompare < a + b): # transversion
                        newseq += transversions[nucleotide]
                    else:
                        newseq += nucleotide
            else:
                newseq += 'CAG'
        # reflect that the probability of an extension slightly increases with the number of repeats
        probchange = 0.00002 * (countRepeats(individual) - 16)
        tri += probchange
        tocompare = random.random()
        if (tocompare < tri): # CAG extension
            extension = 'CAG' * random.randint(1, 5)
            # add extension after first 9 characters (start of repeat section)
            newseq = newseq[:9] + extension + newseq[9:]
        newpop.append(newseq)
    return newpop

'''
reproduce - produces a new population, accounting for natural selection
param pop: the population at a given generation
param n: the number of sequences in the population
return newpop: the population at the next generation
'''
def reproduce(pop, n):
    newpop = []
    # assumes haploid individuals, may be changed in future experiments
    for i in range(n):
        newpop.append(random.choice(pop))
    return newpop            

'''
countDisplays - counts the number of individuals displaying Huntington's in a given population
param pop: the population at a given generation
return count: the number of individuals displaying Huntington's
'''
def countDisplays(pop):
    count = 0
    for individual in pop:
        if (countRepeats(individual) >= 40):
            count += 1
    return count 

'''
countRepeats - counts the number of CAG repeats in a sequence
param sequence: the sequence to be searched
return count: the number of CAG repeats in the sequence
'''
def countRepeats(sequence):
    count = 0
    for n in range(0, len(sequence), 3):
        codon = sequence[n:n+3]
        if (codon == 'CAG'):
            count += 1
    return count
    
main()
