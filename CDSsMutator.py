import random, math, numpy as np, time, csv
import matplotlib.pyplot as plt
import subprocess
from subprocess import DEVNULL
from Bio import SeqIO
import itertools

PopulationSize = 50
Generations = 50
MRate = 0.05  # Mutation Rate
t = 5

AminoAcids = {}
AminoAcids['Phe'] = ['TTT', 'TTC']
AminoAcids['Leu'] = ['TTA', 'TTG', 'CTT', 'CTC', 'CTA', 'CTG']
AminoAcids['Ile'] = ['ATT', 'ATC', 'ATA']
AminoAcids['Met'] = ['ATG']
AminoAcids['Val'] = ['GTT', 'GTC', 'GTA', 'GTG']
AminoAcids['Ser'] = ['TCT', 'TCC', 'TCA', 'TCG', 'AGT', 'AGC']
AminoAcids['Pro'] = ['CCT', 'CCC', 'CCA', 'CCG']
AminoAcids['Thr'] = ['ACT', 'ACC', 'ACA', 'ACG']
AminoAcids['Ala'] = ['GCT', 'GCC', 'GCA', 'GCG']
AminoAcids['Tyr'] = ['TAT', 'TAC']
AminoAcids['TER'] = ['TAA', 'TAG', 'TGA']
AminoAcids['His'] = ['CAT', 'CAC']
AminoAcids['Gln'] = ['CAA', 'CAG']
AminoAcids['Asn'] = ['AAT', 'AAC']
AminoAcids['Lys'] = ['AAA', 'AAG']
AminoAcids['Asp'] = ['GAT', 'GAC']
AminoAcids['Glu'] = ['GAA', 'GAG']
AminoAcids['Cys'] = ['TGT', 'TGC']
AminoAcids['Trp'] = ['TGG']
AminoAcids['Arg'] = ['CGT', 'CGC', 'CGA', 'CGG', 'AGA', 'AGG']
AminoAcids['Gly'] = ['GGT', 'GGC', 'GGA', 'GGG']

Bases = ["A", "T", "C", "G"]

CDSPositions = [[265, 21555], [265, 13461], [21563, 25384], [25393, 26220], [26245, 26472], [26523, 27191],
                [27202, 27387], [27394, 27759], [27756, 27887], [27894, 28259], [28274, 29533], [29558, 29674]]

refarray = []
reference = "/home/armin/Desktop/COVID19/RSCU-Human.csv"

global RefGenome

fasta_reference = SeqIO.parse(open("/home/armin/Desktop/COVID19/COVID-19Genome.fa"), 'fasta')
for i in fasta_reference:
    RefGenome = str(i.seq)

CDS = []
for i in CDSPositions:
    CDS.append(RefGenome[i[0]:i[1] + 1])


def RefProcess(address):
    csvfile = open(address, newline='')
    ref = csv.reader(csvfile, delimiter=',')
    for i in ref:
        tmp = np.array(i)
        tmp = tmp.astype(float)
        refarray.append(tmp)


def Compatibility(chromosome):
    Input = open("/home/armin/Desktop/COVID19/Temp/Input.dat", 'w+')
    c = 0

    Input.write(">CDS" + str(c) + "\n")
    s = ""
    w = s.join(chromosome[0:])
    Input.write(w + "\n")
    Input.close()

    bashCommand = "codonw Input.dat output.out output.blk -nomenu -nowarn -silent -cutot -rscu -machine"
    process = subprocess.Popen(bashCommand.split(), stdout=DEVNULL, stderr=DEVNULL,
                               cwd='/home/armin/Desktop/COVID19/Temp')
    output, error = process.communicate()
    FileOut = open("/home/armin/Desktop/COVID19/Temp/output.blk", 'r')
    fr = FileOut.readlines()
    tarray = []  # RSCU of our chromosome
    for line in fr:
        ta = line.split()
        tta = []
        for i in range(0, 16):
            tta.append(float(ta[i]))
        for j in range(0, 16, 4):
            tarray.append([tta[j], tta[j + 1], tta[j + 2], tta[j + 3]])
    tarray = np.array(tarray)
    dif = np.subtract(refarray, tarray)
    difsquared = np.power(dif, 2)
    sumdif = np.sum(difsquared)
    FileOut.close()
    return sumdif


def GCCalc(chromosome, ChrSize):
    gc = 0
    for i in range(0, ChrSize):
        e = chromosome[i]
        if (e == "G" or e == "C"):
            gc += 1
    return gc / ChrSize


def CpGD(chromosome, ChrSize):
    cpg = 0
    c = 0
    g = 0
    for i in range(0, ChrSize):
        e = chromosome[i]
        if i < ChrSize - 1:
            ne = chromosome[i + 1]
            if (e == "C" and ne == "G"):
                cpg += 1
        if (e == "G"):
            g += 1
        elif (e == "C"):
            c += 1
    return ((cpg / ChrSize) / ((c / ChrSize) * (g / ChrSize)))


def Identity(chromosome, ChrSize, cdsSeq):
    notmatch = 0
    for i in range(0, ChrSize):
        if chromosome[i] != cdsSeq[i]:
            notmatch += 1
    return (notmatch/ChrSize)*100


def Fitness(chromosome):
    ChrSize = len(chromosome)
    # Check CDS criteria
    comp = Compatibility(chromosome) / 7  # 222 is max difference
    if comp == 0:
        return 0
    # Check for Identity with COVID-19 Sequence

    score = CpGWeight * CpGD(chromosome, ChrSize) + CompWeight * comp

    return (1 / score)


# Selection Function

def Selection(Population, Childs):  # Tournament method with t as an option (default=3)
    TotalPop = Population + Childs
    for j in range(0, PopulationSize):
        tmpFit = []
        ap = []
        for i in range(0, t):
            rand = random.randint(0, len(TotalPop) - 1)
            ap.append(TotalPop[rand])
            tmpFit.append(Fitness(ap[i]))
        Population[j] = ap[tmpFit.index(max(tmpFit))]
        TotalPop.remove(ap[tmpFit.index(max(tmpFit))])
    return Population


# Mutation Function

def Mutation(chro):  # Mutation
    global newcodon
    for j in range(0, len(chro), 3):
        if (MRate > random.random()):
            lcodon = chro[j:j + 3]
            codon = ""
            codon = codon.join(lcodon)
            for a in AminoAcids.keys():
                if codon in AminoAcids[a]:
                    newcodon = random.choice(AminoAcids[a])
            chro[j] = newcodon[0]
            chro[j + 1] = newcodon[1]
            chro[j + 2] = newcodon[2]
    # chro[27386]='A'
    # chro[27387] = 'T'
    # chro[27388] = 'G'
    chro[0] = 'A'
    chro[1] = 'T'
    chro[2] = 'G'
    return chro


def CreateChr(cdsSeq):
    global newcodon
    ChrSize = len(cdsSeq)
    a = ChrSize * ['E']
    for j in range(0, ChrSize, 3):
        codon = cdsSeq[j:j + 3]
        for p in AminoAcids.keys():
            if codon in AminoAcids[p]:
                newcodon = random.choice(AminoAcids[p])
        a[j] = newcodon[0]
        a[j + 1] = newcodon[1]
        a[j + 2] = newcodon[2]
    return a


# Prepare Human Reference Codon Usage
RefProcess(reference)

CDStoProcess=[]
CDStoProcess.append(CDS[2])


for CpGWeight,CompWeight in itertools.product(range(1,11),range(1,11)):
    BestResult=[]
    for cdsSeq in CDStoProcess:
        print('CDS #', CDS.index(cdsSeq), ":","Length =",len(cdsSeq))
        print("CompWeight:",CompWeight,"CpGWeight:",CpGWeight)
        ChrSize = len(cdsSeq)
        stime = time.time()

        # CpGWeight = 2  # Weight of CpG in Fitness
        # CompWeight = 1  # Weight of Compatibility in Fitness

        # Create Initial Population
        Population = []
        k = 0
        c = 0
        Fits = []
        while k < PopulationSize:
            c += 1
            a = CreateChr(cdsSeq)
            Population.append(a)
            Fits.append(Fitness(a))
            k += 1

        # Print Best Answer in First Generation
        best = []
        best.append(Population[Fits.index(max(Fits))])
        best.append(max(Fits))

        iteration = 1
        # print("Generation", iteration, "-> Max Fitness is: ", best[1])

        grf = []
        # ------------------------- #
        # -------- Section 5 ------- #

        # Repeat Algorithm
        for ite in range(0, Generations - 1):
            iteration = iteration + 1
            # Create Childs - Mutation
            Childs = []
            for i in Population:
                tempComp = 0
                child = []
                while tempComp == 0:
                    child = Mutation(i)
                    tempComp = Compatibility(child)
                Childs.append(child)

            # Calculate Childs Fitnesses and Print Best of Them
            CFits = []
            for i in range(0, PopulationSize):
                CFits.append(Fitness(Childs[i]))

            if (max(CFits) >= best[1]):
                best[0] = Childs[CFits.index(max(CFits))]
                best[1] = max(CFits)
            # print("Generation", iteration, "-> Max Fitness is: ", best[1])
            grf.append(best[1])
            Population = Selection(Population, Childs)

        # ------------------------- #
        # -------- Section 6 ------- #

        etime = time.time()
        # print("The Best Answer is:\n", best[0], "\nThe Identity of The Best Answer is:", (100-Identity(best[0], ChrSize, cdsSeq)),
        #       "\nCpG Score:",CpGD(best[0], ChrSize),"\nRuntime(min): ", (etime - stime) / 60,"\n\n")
        print("Identity:",(100-Identity(best[0], ChrSize, cdsSeq)),"\n")
        Summary=[CompWeight,CpGWeight,best[0],best[1],(100-Identity(best[0], ChrSize, cdsSeq))]
        BestResult.append(Summary)
        '''
        XB = np.arange(1, Generations, 1)
        plt.plot(XB, grf)
        plt.ylabel('Fitness')
        plt.xlabel('Generation')
        title = "MRate=" + str(MRate) + "   t=" + str(t)
        plt.title(title)
        name = "/home/armin/Desktop/COVID19/PlotCodonWise.png"
        plt.savefig(name)
        plt.clf()
        '''
        # ------------------------- #

print("\n\n\n---- Results ----")
for i in BestResult:
    print("CompWeight:",i[0],"CpGWeight:",i[1],"Identity:",i[4])