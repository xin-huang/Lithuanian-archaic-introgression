# This file is modified based on the script from
# https://github.com/xzhang-popgen/maladapt/blob/main/empirical/1empirical_compute-features.py
# Several parameters were originally hard-coded and have now been adjusted
# EOF errors were fixed in get_vcf_modern and get_vcf_archaic
# Permission to redistribute this code has been granted by the original developer, Prof. Xinjun Zhang


import msprime, pyslim, os, random, itertools,argparse,gzip
import numpy as np
import traceback


parser = argparse.ArgumentParser(description="A script for computing summary statistics in 50kb windows across given chromosome, between modern human and archaic human.")
parser.add_argument('-c', '--chr', action="store", dest="chr_id",
                        help="which chromosome, default: 6",
                        default=6, type=int)
parser.add_argument('-p', '--pop', action="store", dest="pop_id",
                        help="which non-YRI 1000G population: 'BEB','CHB' etc, default: 'BEB'",
                        default=1, type=int)
 
# Directory arguments
parser.add_argument('--vcfm', action="store", dest="DIR_vcfm",
                    help="Directory for modern VCF output files",
                    default="results/maladapt/vcf_modern_out/")
parser.add_argument('--vcfa', action="store", dest="DIR_vcfa",
                    help="Directory for archaic VCF output files",
                    default="results/maladapt/vcf_archaic_out/")
parser.add_argument('--dir-pop', action="store", dest="DIR_pop",
                    help="Directory for population files",
                    default="results/maladapt/dir_pop/")
parser.add_argument('--out', action="store", dest="DIR_out",
                    help="Directory for output files",
                    default="results/maladapt/dir_out/")
parser.add_argument('--stat', action="store", dest="DIR_stat",
                    help="Directory for statistical output files",
                    default="results/maladapt/dir_stat_out/")
parser.add_argument('--archaic', action="store", dest="DIR_archaic",
                    help="Directory for archaic input files",
                    default="results/maladapt/vcf_archaic/")
parser.add_argument('--1000g', action="store", dest="DIR_1000g",
                    help="Directory for 1000G input files",
                    default="results/maladapt/vcf_modern/")
parser.add_argument('--range', action="store", dest="file_range",
                    help="File containing window range custom data",
                    default="results/maladapt/window_range_custom/CHR6_windows_overlap.txt")

args = parser.parse_args()

allpop = ['CHB','KHV','GIH','PJL','CHS','CEU','GBR','IBS','PEL','TSI','JPT','CDX','FIN','MXL','PUR','CLM','BEB','STU','ITU', 'LIT']

chrom = args.chr_id
population = allpop[args.pop_id - 1]

#hardcoded population LIT
population = 'LIT'
#population = 'BEB'
print(chrom)
print(population)

DIR_vcfm = args.DIR_vcfm
DIR_vcfa = args.DIR_vcfa
DIR_pop = args.DIR_pop
DIR_out = args.DIR_out
DIR_stat = args.DIR_stat
DIR_archaic = args.DIR_archaic
DIR_1000g = args.DIR_1000g
file_range = args.file_range

print("Directories and file paths:")
print(f"Modern VCF Directory: {DIR_vcfm}")
print(f"Archaic VCF Directory: {DIR_vcfa}")
print(f"Population Directory: {DIR_pop}")
print(f"Output Directory: {DIR_out}")
print(f"Stat Directory: {DIR_stat}")
print(f"Archaic Directory: {DIR_archaic}")
print(f"1000G Directory: {DIR_1000g}")
print(f"File Range: {file_range}")

os.makedirs(DIR_vcfm, exist_ok=True)
os.makedirs(DIR_vcfa, exist_ok=True)
os.makedirs(DIR_pop, exist_ok=True)
os.makedirs(DIR_out, exist_ok=True)
os.makedirs(DIR_stat, exist_ok=True)
os.makedirs(DIR_archaic, exist_ok=True)
os.makedirs(DIR_1000g, exist_ok=True)

#DIR_vcfm = "results/maladapt/vcf_modern_out/"
#DIR_vcfa = "results/maladapt/vcf_archaic_out/"

#DIR_pop = "results/maladapt/dir_pop/"
#DIR_out = "results/maladapt/dir_out/"
#DIR_stat = "results/maladapt/dir_stat_out/"

#archaic infolder
#DIR_archaic = "results/maladapt/vcf_archaic/"

#file_range = "results/maladapt/window_range_custom/CHR"+str(chrom)+"_windows_overlap.txt"

#modern infolder
#DIR_1000g = "results/maladapt/vcf_modern/"

#original pre- and suffixes - not used
file_prefix_1000g = "ALL.chr"
file_suffix_1000g = ".phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz"
file_prefix_archaic = "AltaiNea_Den_combined."
file_suffix_archaic = ".vcf.gz"


file_prefix_1000g = "merged_"
file_suffix_1000g = "_output.vcf.gz"


file_prefix_archaic = "chr"
file_suffix_archaic = "_mq25_mapab100.vcf.gz"


#####################
#read range info
file = open(file_range,"rt")
line = file.readline()
field = line.split()
start = field
print("start")
print(start)
line = file.readline()
field = line.split()
end = field
print("end")
print(end)
line = file.readline()
field = line.split()
length = field
print("length")
print(length)
line = file.readline()
field = line.split()
snps = field
print("snps")
print(snps)

file.close()

#####################
    
def get_vcf_modern (chromosome,region_start,region_end,outname):
    vcf_infile = DIR_1000g+file_prefix_1000g+str(chromosome)+file_suffix_1000g
    outfile = outname+".recode.vcf"
    vf = gzip.open(vcf_infile, 'rt')
    out = open(outfile, 'w')
    vcf_line = vf.readline()
    while "#" in vcf_line:
        out.write(vcf_line)
        vcf_line = vf.readline()
    vcf_line_split = vcf_line.split()
    vcf_pos = vcf_line_split[1]
    vcf_pos = int(vcf_pos)
    vcf_ref = vcf_line_split[3]
    vcf_alt = vcf_line_split[4]
    try:
        vcf_qual = float(vcf_line_split[5])
    except:
        vcf_qual = 100
    while vcf_pos < int(region_start):
        vcf_line = vf.readline()
        if vcf_line == '': break
        vcf_line_split = vcf_line.split()
        vcf_pos = vcf_line_split[1]
        vcf_pos = int(vcf_pos)
        vcf_ref = vcf_line_split[3]
        vcf_alt = vcf_line_split[4]

        try:
            vcf_qual = float(vcf_line_split[5])
        except:
            vcf_qual = 100
    while vcf_pos <= int(region_end):

        if (len(vcf_ref)==1) & (len(vcf_alt)==1) & (vcf_qual>=40):
            out.write(vcf_line)
        vcf_line = vf.readline()
        if vcf_line == '': break
        vcf_line_split = vcf_line.split()
        vcf_ref = vcf_line_split[3]
        vcf_alt = vcf_line_split[4]
        try:
            vcf_qual = float(vcf_line_split[5])
        except:
            vcf_qual = 100
        if len(vcf_line_split)> 1:
            vcf_pos = int(vcf_line_split[1]) 
            vcf_ref = vcf_line_split[3]
            vcf_alt = vcf_line_split[4]
            try:
                vcf_qual = float(vcf_line_split[5])
            except:
                vcf_qual = 100
        else:
            vcf_pos += 1
    vf.close()
    out.close()

def get_vcf_archaic (chromosome,region_start,region_end,outname):
    vcf_infile = DIR_archaic+file_prefix_archaic+str(chromosome)+file_suffix_archaic
    outfile = outname+".recode.vcf"
    vf = gzip.open(vcf_infile, 'rt')
    out = open(outfile, 'w')
    vcf_line = vf.readline()
    while "#" in vcf_line:
        out.write(vcf_line)
        vcf_line = vf.readline()
    vcf_line_split = vcf_line.split()
    vcf_pos = vcf_line_split[1]
    vcf_pos = int(vcf_pos)
    vcf_ref = vcf_line_split[3]
    vcf_alt = vcf_line_split[4]
    #vcf_qual = float(vcf_line_split[5])
    try:
        vcf_qual = float(vcf_line_split[5])
    except:
        vcf_qual = 100

    vcf_filter = vcf_line_split[6]
    while vcf_pos < int(region_start):
        vcf_line = vf.readline()
        if vcf_line == '': break # EOF
        vcf_line_split = vcf_line.split()
        vcf_pos = vcf_line_split[1]
        vcf_pos = int(vcf_pos)
        vcf_ref = vcf_line_split[3]
        vcf_alt = vcf_line_split[4]
        try:
            vcf_qual = float(vcf_line_split[5])
        except:
            vcf_qual = 100
        vcf_filter = vcf_line_split[6]
    while vcf_pos <= int(region_end):
        if (len(vcf_ref)==1) & (len(vcf_alt)==1) & (vcf_qual>=100) & (vcf_filter == "."):
            out.write(vcf_line)
        vcf_line = vf.readline()
        if vcf_line == '': break # EOF
        vcf_line_split = vcf_line.split()
        vcf_ref = vcf_line_split[3]
        vcf_alt = vcf_line_split[4]
        #vcf_qual = float(vcf_line_split[5])
        try:
            vcf_qual = float(vcf_line_split[5])
        except:
            vcf_qual = 100

        vcf_filter = vcf_line_split[6]
        if len(vcf_line_split)> 1:
            vcf_pos = int(vcf_line_split[1])
            vcf_ref = vcf_line_split[3]
            vcf_alt = vcf_line_split[4]
            #vcf_qual = float(vcf_line_split[5])
            try:
                vcf_qual = float(vcf_line_split[5])
            except:
                vcf_qual = 100
            vcf_filter = vcf_line_split[6]
        else:
            vcf_pos += 1
    vf.close()
    out.close()

def get_pop (pop,pop2,chrom,start,end,vcf_kg): #pop = string = "YRI"
    vcf_infile=vcf_kg
    if pop == "YRI":
    	outvcfname = DIR_out+pop+pop2+"_"+"chr"+str(chrom)+"_"+str(start)+"-"+str(end)+".recode.vcf"
    else:
    	outvcfname = DIR_out+pop+"_"+"chr"+str(chrom)+"_"+str(start)+"-"+str(end)+".recode.vcf"
    outfile = outvcfname    
    popname = DIR_pop+pop+".txt"    
    indlist = list()
    with open(popname) as f:
        for line in f:
            line = line.split("\n")[0]
            indlist.append(line)
    
    vcf = open(vcf_infile, 'rt')
    out=open(outfile,"w")
    start = 0
    startg = 0    
    for count,line in enumerate(vcf):
        if line[0:6]=="#CHROM":
            start +=1
            startg = count
            fields = line.split()
            these = []
            for ind in indlist:
                if ind in fields:
                    these.append(fields.index(ind))

            these.extend(list(range(9)))
            these.sort()
            newfields = [fields[i] for i in these]
            out.writelines(str(i)+"\t" for i in newfields)
            out.writelines("\n")            
        if (start == 1) & (startg!=0) & (count>startg):
            fields = line.split()
            newfields = [fields[i] for i in these]
            out.writelines(str(i)+"\t" for i in newfields)
            out.writelines("\n")
    vcf.close()
    out.close()
    
    
def vcf2genos (vcfpath,archaic_id): 
    try:
        print("vcf2genos")
        print(archaic_id)
        print(vcfpath)
        if archaic_id !=0:
                try:
                    vcf = open(vcfpath)
                except:
                    print("open vcfpath ERROR!")
                start = 0
                startg = 0
                pos = []
                geno = []
                for count,line in enumerate(vcf):
                    if line[0:6]=="#CHROM":
                        start +=1
                        startg = count
                    if (start == 1) & (startg!=0) & (count>startg):
                        fields = line.split()

                        #if (len(fields[3])==1) & (len(fields[4])==1)  & (fields[3]!=".") & (fields[4]!=".")& (float(fields[5])>=100)& ("." not in fields[9][0:3])& ("." not in fields[10][0:3]):

                        if (len(fields[3])==1) & (len(fields[4])==1)  & (fields[3]!=".") & (fields[4]!=".")& (float(fields[5])>=100)& ("." not in fields[9][0:3]):

                            pos.append(int(fields[1]))
                            geno.append(fields[9:])
                vcf.close()
        elif archaic_id ==0:

            vcf = open(vcfpath)
            start = 0
            startg = 0
            pos = []
            geno = []
            for count,line in enumerate(vcf):
                #print("enum 0")
                if line[0:6]=="#CHROM":
                    start +=1
                    startg = count
                elif start == 1:
                    fields = line.split()

                    #
                    #if (len(fields[3])==1) & (len(fields[4])==1)  & (fields[3]!=".") & (fields[4]!=".")& (float(fields[5])>=100)& ("." not in fields[9][0:3])& ("." not in fields[10][0:3]):

                    pos.append(int(fields[1]))
                    geno.append(fields[9:])

            vcf.close()

        geno1 = []

        if (archaic_id==1):
            for n in range(len(geno)):
                geno1.append([geno[n][x][0]+"|"+geno[n][x][2] for x in [0]])
        elif (archaic_id ==2):
            for n in range(len(geno)):
                geno1.append([geno[n][x][0]+"|"+geno[n][x][2] for x in [1]])
        else:
            for n in range(len(geno)):
                geno1.append([geno[n][x][0]+"|"+geno[n][x][2] for x in range(len(geno[0]))])


        num_pos = len(geno1)
        num_ind = len(geno1[1])


        geno_ref = np.empty([num_pos,num_ind], dtype=object)
        geno_alt = np.empty([num_pos,num_ind], dtype=object)

        for n in range(num_pos):



            try:
                #geno_ref[n] = [int(i.split('|', 1)[0]) for i in geno1[n]]
                #geno_alt[n] = [int(i.split('|', 1)[1]) for i in geno1[n]]
                geno_ref[n] = [(i.split('|', 1)[0]) for i in geno1[n]]
                geno_alt[n] = [(i.split('|', 1)[1]) for i in geno1[n]]
            #except:
            except Exception as error:
                '''
                print("An error occurred:", error)

                print([(i.split('|', 1)[0]) for i in geno1[n]])
                print([(i.split('|', 1)[1]) for i in geno1[n]])
                '''
                pass



        geno_ref = np.transpose(geno_ref)
        geno_alt = np.transpose(geno_alt)

        hapmat = np.append(geno_ref,geno_alt,axis=0)

    except Exception:
        pass
        #traceback.print_exc()
        
    return hapmat,geno1,pos


def insert_anc_alleles (allpos,pos,hap): #insert non-variant alleles to an existing haplotype matrix
    # List of all locations in allpos where pos has a value
    pos_locs = np.isin(allpos, pos)
    # Transposed hap matrix to add and remove columns more easily
    hapT = hap.T
    # creating blank hap matrix, Transposed to more easilly assign colums
    new_hap = np.zeros((hap.shape[0], pos_locs.shape[0])).T
    #counter to keep track of current column being placed
    hap_counter = 0
    # finds all indexes where pos_loc is 1, and if so make that column in the new hap the next column in the original hap
    for i in range(0, pos_locs.shape[0]):
        if pos_locs[i] == 1:
            #print (hap[hap_counter])
            new_hap[i] = hapT[hap_counter]
            hap_counter = hap_counter + 1
    # since at the end of the for loop, pos would just be allpos, I return allpos, and I transpose the new haplotype matrix to get the actual matrix
    return allpos, new_hap.T

def insert_anc_alleles_geno (allpos,pos,geno): #insert non-variant alleles to an existing haplotype matrix
    # List of all locations in allpos where pos has a value
    newgeno = []
    
    for i in range(len(allpos)):
        if allpos[i] in pos:
            newgeno.append(geno[pos.index(allpos[i])])
        else:
            newgeno.append(['0|0']*len(geno[0]))

    
    # since at the end of the for loop, pos would just be allpos, I return allpos, and I transpose the new haplotype matrix to get the actual matrix
    return newgeno



def calc_freq (pop_hap):
    popfreq = np.sum(pop_hap, axis=0)
    popfreq = popfreq/ float(pop_hap.shape[0])
    return popfreq

def vSumFunc(source_hap,other_hap, currentArchi):
    current_hap = np.array([source_hap[currentArchi,]])
    div = np.zeros(other_hap.shape)
    ones = np.ones((other_hap.shape[0],1))
    current_hap = current_hap
    current_hap_extended = np.dot(ones, current_hap)
    div = np.logical_xor(current_hap_extended == 1, other_hap == 1)
    return np.add.reduce(div, 1)   


def calc_stats (vcf_adm,vcf_arc,vcf_yri,len_genome,archaic_id):
    try:
        hap_adm, geno_adm,pos_adm = vcf2genos(vcf_adm,0)
        hap_arc, geno_arc,pos_arc = vcf2genos(vcf_arc,archaic_id)
        hap_yri, geno_yri,pos_yri = vcf2genos(vcf_yri,0)

    except:
        pass


    p1_pos = pos_arc
    p2_pos = pos_yri
    p3_pos = pos_adm
    
    p1_hap = hap_arc
    p2_hap = hap_yri
    p3_hap = hap_adm
    p1_geno = geno_arc
    p2_geno = geno_yri
    p3_geno = geno_adm
            
    p1size = 2
    p2size = p2_hap.shape[0]
    p3size = p3_hap.shape[0]


    
    all_pos = np.unique(np.concatenate((p1_pos,p2_pos,p3_pos))) #len=244
    
    p1_geno = insert_anc_alleles_geno (all_pos,p1_pos,p1_geno)
    p2_geno = insert_anc_alleles_geno (all_pos,p2_pos,p2_geno)
    p3_geno = insert_anc_alleles_geno (all_pos,p3_pos,p3_geno)
    #3. insert non-segregating sites/ancestral alleles 0s to the hap matrices and pos lists  
    p1_pos,p1_hap = insert_anc_alleles(all_pos,p1_pos,p1_hap)
    p2_pos,p2_hap = insert_anc_alleles(all_pos,p2_pos,p2_hap)
    p3_pos,p3_hap = insert_anc_alleles(all_pos,p3_pos,p3_hap)


    #D
    p1_freq = calc_freq (p1_hap)
    p2_freq = calc_freq (p2_hap)
    p3_freq = calc_freq (p3_hap)

    abbavec = (1.0 - p2_freq)*p3_freq*p1_freq
    babavec = p2_freq*(1.0 - p3_freq)*p1_freq
    abba = np.sum(abbavec)
    baba = np.sum(babavec)

    if (abba + baba > 0):
        Dstat = (abba - baba) / (abba + baba)
    else:
        Dstat = float('nan')
    #fD
    checkfd1 = (p3_freq > p1_freq)
    abbafd1 = (1.0 - p2_freq)*p3_freq*p3_freq
    babafd1 = p2_freq*(1.0 - p3_freq)*p3_freq
    checkfd2 = (p3_freq < p1_freq)
    abbafd2 = (1.0 - p2_freq)*p1_freq*p1_freq
    babafd2 = p2_freq*(1.0 - p1_freq)*p1_freq
    abbafd = checkfd1 * abbafd1 + checkfd2 * abbafd2
    babafd = checkfd1 * babafd1 + checkfd2 * babafd2
    abbafd = np.sum(abbafd)
    babafd = np.sum(babafd)
    if (abbafd + babafd > 0):
        fD = (abba - baba) / (abbafd - babafd)
    else:
        fD = float('nan')

    #Heterozygosity
    hetvec = 2 * p3_freq * (1.0 - p3_freq)
    Het = np.sum(hetvec) / len(all_pos)

    #Divergence
    #div = divergence ratio between pairs of individual sequences (admixed/source), 
    #and average across all pairs of individuals
    divratio = []

    for archi in range(0, p1_hap.shape[0]): #iterate over 0-99 haps; 100 total
        divarchintro = vSumFunc(p1_hap,p3_hap, archi)
        divarchintro = divarchintro.astype("float")
        divarchnonintro = vSumFunc(p1_hap,p2_hap, archi)
        
        divarchnonintro = divarchnonintro.astype("float") #took the inversion here so that the multiplying below is really actually dividing; probably not necessary       
        #pairwise combos of divarchintro/divarchnonintro -> 100*100 comparisons
        for comb in itertools.product(divarchintro,divarchnonintro): #pairwise combos of divarchintro/divarchnonintro
            #print(comb)
            if comb[1] != 0:
                #print(comb)
                divratio.append(comb[0]/comb[1])
    divratioavg = float(sum(divratio)) / float(len(divratio)) #len(divratio) = 100* (100*100)
    
    #Q Dist
    Arc100 = (p1_freq == 1)
    NonAdm1 = (p2_freq < 0.01) 
    NonAdm10 = (p2_freq < 0.1)
    Arc100NonAdm1 = (Arc100 & NonAdm1)
    Arc100NonAdm10 = (Arc100 & NonAdm10)
    Freqs_Arc100NonAdm10 = p3_freq[np.where(Arc100NonAdm10 == True)]
    Freqs_Arc100NonAdm1 = p3_freq[np.where(Arc100NonAdm1 == True)]
    if Freqs_Arc100NonAdm10.size > 0:
        Q_10_100_q95 = np.percentile(Freqs_Arc100NonAdm10,95)
        Q_10_100_q90 = np.percentile(Freqs_Arc100NonAdm10,90)
        Q_10_100_max = np.max(Freqs_Arc100NonAdm10)
    else:
        Q_10_100_q95 = float('nan')
        Q_10_100_q90 = float('nan')
        Q_10_100_max = float('nan')
    if Freqs_Arc100NonAdm1.size > 0:
        Q_1_100_q95 = np.percentile(Freqs_Arc100NonAdm1,95)
        Q_1_100_q90 = np.percentile(Freqs_Arc100NonAdm1,90)
        Q_1_100_max = np.max(Freqs_Arc100NonAdm1)
    else:
        Q_1_100_q95 = float('nan')
        Q_1_100_q90 = float('nan')
        Q_1_100_max = float('nan')


    #Uniquely Shared Alleles
    U_10_0_100 = ( Arc100NonAdm10 & (p3_freq > 0) )
    U_10_20_100 = ( Arc100NonAdm10 & (p3_freq > 0.2) )
    U_10_50_100 = ( Arc100NonAdm10 & (p3_freq > 0.5) )
    U_10_80_100 = ( Arc100NonAdm10 & (p3_freq > 0.8) )
		
    U_1_0_100 = ( Arc100NonAdm1 & (p3_freq > 0) )
    U_1_20_100 = ( Arc100NonAdm1 & (p3_freq > 0.2) )
    U_1_50_100 = ( Arc100NonAdm1 & (p3_freq > 0.5) )
    U_1_80_100 = ( Arc100NonAdm1 & (p3_freq > 0.8) )
		
    U_10_0_100 = np.sum(U_10_0_100)
    U_10_20_100 = np.sum(U_10_20_100)
    U_10_50_100 = np.sum(U_10_50_100)
    U_10_80_100 = np.sum(U_10_80_100)
        
    U_1_0_100 = np.sum(U_1_0_100)
    U_1_20_100 = np.sum(U_1_20_100)
    U_1_50_100 = np.sum(U_1_50_100)
    U_1_80_100 = np.sum(U_1_80_100)
    
    
    p1_S, p2_S, p3_S, thetaW1, thetaW2, thetaW3,thetapi1, thetapi2, thetapi3,thetaH1, thetaH2, thetaH3,p1_H1, p1_H2,p1_H12, p1_H2H1,p2_H1, p2_H2,p2_H12, p2_H2H1,p3_H1, p3_H2,p3_H12, p3_H2H1 = stats_segthetasH12 (p1_hap,p2_hap,p3_hap, p1size,p2size,p3size)

    Zns1 = Zns (p1_hap,p1_geno)
    Zns2 = Zns (p2_hap,p2_geno)
    Zns3 = Zns (p3_hap,p3_geno)

    return Dstat, fD, Het, divratioavg,Q_10_100_q95,Q_10_100_q90,Q_10_100_max,Q_1_100_q95,Q_1_100_q90,Q_1_100_max,U_10_0_100,U_10_20_100,U_10_50_100,U_10_80_100, U_1_0_100,U_1_20_100,U_1_50_100,U_1_80_100,p1_S, p2_S, p3_S, thetaW1, thetaW2, thetaW3,thetapi1, thetapi2, thetapi3,thetaH1, thetaH2, thetaH3,p1_H1, p1_H2,p1_H12, p1_H2H1,p2_H1, p2_H2,p2_H12, p2_H2H1,p3_H1, p3_H2,p3_H12, p3_H2H1,Zns1,Zns2,Zns3 


def numseg (pop_hap):
    pop_freq = calc_freq (pop_hap)
    
    pop_S = np.sum((pop_freq>0)&(pop_freq<1))
    
    S_freq = pop_freq[(pop_freq>0)&(pop_freq<1)]
    
    return pop_S,S_freq 

def thetaW (numS, popsize):
    x = list(range(1, popsize*2))
    denom = [1/float(i) for i in x]
    
    thetaW = numS/sum(denom)
    return thetaW

def thetapi (S_freq, popsize):
    pi = sum(2*S_freq*(1.-S_freq))
    
    n=popsize *2
    
    thetapi = pi*n/(n-1)

    return thetapi

def thetaH (pophap, popsize):
    hap = np.transpose(pophap)
    n= popsize*2
    Si = [sum(i) for i in hap]
    Si = np.asarray([sum(i) for i in hap])
    num = [x**2*sum(Si == x) for x in list(range(1,n))]
    thetaH = sum(num)*2/(n*(n-1))
    return thetaH

    
def hap_homo (pophap,popsize):
    n=popsize*2
    
    unique,count = np.unique(pophap, return_counts=True,axis=0)
    if len(unique) > 1:
        
        freq = count/n
        freq[::-1].sort()
        p1 = freq[0]
        p2 = freq[1]
    
        H1 = sum(freq**2)
        H2 = H1-p1**2
        H12 = H1+2*p1*p2
        H2H1 = H2/H1
    
    else:H1,H2,H12,H2H1=1,0,1,0
        
    
    return H1,H2,H12,H2H1

def stats_segthetasH12 (p1_hap,p2_hap,p3_hap, p1size,p2size,p3size):
    p1_S, p1_Sfreq= numseg(p1_hap)
    p2_S, p2_Sfreq= numseg(p2_hap)
    p3_S, p3_Sfreq= numseg(p3_hap)

    thetaW1, thetaW2, thetaW3 = thetaW (p1_S, p1size),thetaW (p2_S, p2size),thetaW (p3_S, p3size)

    thetapi1, thetapi2, thetapi3 = thetapi(p1_Sfreq,p1size),thetapi(p2_Sfreq,p2size),thetapi(p3_Sfreq,p3size)

    thetaH1, thetaH2, thetaH3 = thetaH (p1_hap, p1size),thetaH (p2_hap, p2size),thetaH (p3_hap, p3size)

    p1_H1, p1_H2,p1_H12, p1_H2H1 =  hap_homo (p1_hap, p1size)     
    p2_H1, p2_H2,p2_H12, p2_H2H1 =  hap_homo (p2_hap, p2size)
    p3_H1, p3_H2,p3_H12, p3_H2H1 =  hap_homo (p3_hap, p3size)
         
    return p1_S, p2_S, p3_S, thetaW1, thetaW2, thetaW3,thetapi1, thetapi2, thetapi3,thetaH1, thetaH2, thetaH3,p1_H1, p1_H2,p1_H12, p1_H2H1,p2_H1, p2_H2,p2_H12, p2_H2H1,p3_H1, p3_H2,p3_H12, p3_H2H1
    #24 total

def r2(pAB,paB,pAb,pab):
    pA = (pAB+pAb)
    pB = (pAB+paB)
    pa = (paB+pab)
    pb = (pAb+pab)
    D = pAB-pA*pB
    Dab = pab-pa*pb
    if (pA*(1-pA)*pB*(1-pB) != 0 ):
        r2 = (D**2)/(pA*(1-pA)*pB*(1-pB))
    else:
        r2= (Dab**2)/(pa*(1-pa)*pb*(1-pb))
    return r2

def Zns (pophap,popgeno):
    pop_freq = calc_freq (pophap)
    
    geno = np.transpose(popgeno)
    S_geno = np.asarray([g[(pop_freq>0)&(pop_freq<1)] for g in geno])
        
    S = len(S_geno[0])
    if S>1:
        ls = range(1,S+1)
        pairs = list(itertools.combinations(ls,2))
    
        geno_list = np.asarray([[g[n].replace('|', '') for n in range(0,len(S_geno[0]))] for g in S_geno])
        genos_combos = [[list(zip(g[x-1],g[y-1])) for g in geno_list] for x,y in pairs]
    
        geno = [[''.join(i) for n in g for i in n] for g in genos_combos]
        AB = [sum([y == '11' for y in x]) for x in geno]
        aB = [sum([y == '01' for y in x]) for x in geno]
        Ab = [sum([y == '10' for y in x]) for x in geno]
        ab = [sum([y == '00' for y in x]) for x in geno]
        ps = [[w/float((w+x+y+z)),x/float((w+x+y+z)),y/float((w+x+y+z)),z/float((w+x+y+z))] for w,x,y,z in list(zip(AB,aB,Ab,ab))]

        r2s = [r2(w,x,y,z) for w,x,y,z in ps]
    
        Zns = sum(r2s)*2./(S*(S-1))
    else: 
        Zns = float("Nan")
    
    return Zns


#####################


for pop in [population]:
    outfile_nea = DIR_stat+str(chrom)+"output_nea_"+pop+".txt"
    windowfile_nea = DIR_stat+str(chrom)+"window_nea_"+pop+".txt"
    outfile_den = DIR_stat+str(chrom)+"output_den_"+pop+".txt"
    windowfile_den = DIR_stat+str(chrom)+"window_den_"+pop+".txt"
    
    windowfile1 = open(windowfile_nea,'w')
    windowfile2 = open(windowfile_den,'w')
    
    stats_name = ["chr","start","end","length","n_snps","Dstat","fD","Het","divratioavg","Q_10_100_q95","Q_10_100_q90","Q_10_100_max","Q_1_100_q95","Q_1_100_q90","Q_1_100_max","U_10_0_100","U_10_20_100","U_10_50_100","U_10_80_100","U_1_0_100","U_1_20_100","U_1_50_100","U_1_80_100","p1_S","p2_S","p3_S","thetaW1","thetaW2","thetaW3","thetapi1","thetapi2","thetapi3","thetaH1","thetaH2","thetaH3","p1_H1","p1_H2","p1_H12","p1_H2H1","p2_H1","p2_H2","p2_H12","p2_H2H1","p3_H1","p3_H2","p3_H12","p3_H2H1","Zns1","Zns2","Zns3"]

    windowfile1.writelines(i+"\t" for i in stats_name)
    windowfile1.write("\n")
    
    windowfile2.writelines(i+"\t" for i in stats_name)
    windowfile2.write("\n")



    for i in range(len(start)):
        start_pos = int(start[i])
        end_pos = int(end[i])
        len_genome = int(length[i])
        n_snps = int(snps[i])
        chrom = int(chrom)
    
        outname_kg = DIR_vcfm+"modern_chr"+str(chrom)+"_"+str(start_pos)+"-"+str(end_pos)+str(pop)
        outname_arc = DIR_vcfa+"archaic_chr"+str(chrom)+"_"+str(start_pos)+"-"+str(end_pos)+str(pop)
        get_vcf_modern (chrom,start_pos,end_pos,outname_kg)
        get_vcf_archaic (chrom,start_pos,end_pos,outname_arc)
        
    
        vcf_kg = outname_kg+".recode.vcf"
        get_pop (pop,pop,chrom,start_pos,end_pos,vcf_kg)       
        get_pop ("YRI",pop,chrom,start_pos,end_pos,vcf_kg)
        
        vcf_adm = DIR_out+pop+"_"+"chr"+str(chrom)+"_"+str(start_pos)+"-"+str(end_pos) + ".recode.vcf"
        vcf_yri = DIR_out+"YRI"+pop+"_"+"chr"+str(chrom)+"_"+str(start_pos)+"-"+str(end_pos) + ".recode.vcf"
        vcf_arc = outname_arc+".recode.vcf"
        
        archaic_id =1
        
        try:

        	Dstat, fD, Het, divratioavg,Q_10_100_q95,Q_10_100_q90,Q_10_100_max,Q_1_100_q95,Q_1_100_q90,Q_1_100_max,U_10_0_100,U_10_20_100,U_10_50_100,U_10_80_100, U_1_0_100,U_1_20_100,U_1_50_100,U_1_80_100,p1_S, p2_S, p3_S, thetaW1, thetaW2, thetaW3,thetapi1, thetapi2, thetapi3,thetaH1, thetaH2, thetaH3,p1_H1, p1_H2,p1_H12, p1_H2H1,p2_H1, p2_H2,p2_H12, p2_H2H1,p3_H1, p3_H2,p3_H12, p3_H2H1,Zns1,Zns2,Zns3 =calc_stats (vcf_adm,vcf_arc,vcf_yri,len_genome,archaic_id)


        
        	windowfile1.write("%d\t%d\t%d\t%d\t%d\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t" % (chrom,start_pos,end_pos,len_genome,n_snps,Dstat, fD, Het, divratioavg,Q_10_100_q95,Q_10_100_q90,Q_10_100_max,Q_1_100_q95,Q_1_100_q90,Q_1_100_max,U_10_0_100,U_10_20_100,U_10_50_100,U_10_80_100, U_1_0_100,U_1_20_100,U_1_50_100,U_1_80_100))
        	windowfile1.write("%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t" % (p1_S, p2_S, p3_S, thetaW1, thetaW2, thetaW3,thetapi1, thetapi2, thetapi3,thetaH1, thetaH2, thetaH3,p1_H1, p1_H2,p1_H12, p1_H2H1,p2_H1, p2_H2,p2_H12, p2_H2H1,p3_H1, p3_H2,p3_H12, p3_H2H1,Zns1,Zns2,Zns3))
        	windowfile1.write("\n")
        	windowfile1.flush()
        #except:
        except Exception as error:
            #print("An error occurred in 1:", error)

        	pass
        


        archaic_id =2
        
        try:
        	Dstat, fD, Het, divratioavg,Q_10_100_q95,Q_10_100_q90,Q_10_100_max,Q_1_100_q95,Q_1_100_q90,Q_1_100_max,U_10_0_100,U_10_20_100,U_10_50_100,U_10_80_100, U_1_0_100,U_1_20_100,U_1_50_100,U_1_80_100,p1_S, p2_S, p3_S, thetaW1, thetaW2, thetaW3,thetapi1, thetapi2, thetapi3,thetaH1, thetaH2, thetaH3,p1_H1, p1_H2,p1_H12, p1_H2H1,p2_H1, p2_H2,p2_H12, p2_H2H1,p3_H1, p3_H2,p3_H12, p3_H2H1,Zns1,Zns2,Zns3 =calc_stats (vcf_adm,vcf_arc,vcf_yri,len_genome,archaic_id)
        
        	windowfile2.write("%d\t%d\t%d\t%d\t%d\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t" % (chrom,start_pos,end_pos,len_genome,n_snps,Dstat, fD, Het, divratioavg,Q_10_100_q95,Q_10_100_q90,Q_10_100_max,Q_1_100_q95,Q_1_100_q90,Q_1_100_max,U_10_0_100,U_10_20_100,U_10_50_100,U_10_80_100, U_1_0_100,U_1_20_100,U_1_50_100,U_1_80_100))
        	windowfile2.write("%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t" % (p1_S, p2_S, p3_S, thetaW1, thetaW2, thetaW3,thetapi1, thetapi2, thetapi3,thetaH1, thetaH2, thetaH3,p1_H1, p1_H2,p1_H12, p1_H2H1,p2_H1, p2_H2,p2_H12, p2_H2H1,p3_H1, p3_H2,p3_H12, p3_H2H1,Zns1,Zns2,Zns3))
        	windowfile2.write("\n")
        	windowfile2.flush()
        #except:
        except Exception as error:
            #print("An error occurred in 2:", error)

        	pass
        
        os.system('rm '+vcf_kg)
        os.system('rm '+vcf_adm)
        os.system('rm '+vcf_arc)
        os.system('rm '+vcf_yri)
    
    windowfile1.close()
    windowfile2.close()
    

        

    

    
    
