from Generic_Methods import *
from Class_Methods import *
import PrimerIndex

#Can be removed later
import AmpliconGeneration
import time
import numpy as np
import multiprocessing
from math import ceil
from multiprocessing import shared_memory
import argparse
import random

def run_comparison(args):
    #initialize variables to store information
    runtimes = []
    cur_time = time.time()

    #Read sequences
    st = time.time()
    sequences = generate_sequences(args.metadata, args.sequences)
    with open(args.output + '/runtimes_' + str(args.seed) + '.txt', 'w+') as f:
        f.write('Time spent generating sequences: ' + str(time.time() - st) + '\n')
    
    #Randomize sequences
    random.seed(args.seed)
    random.shuffle(sequences)
    sequences = sequences[:args.n_sequences]
    #Note that we should change the numerical sequence ids
    for i in range(len(sequences)):
        sequences[i].id_num = i
    
    #Generate comparison matrix
    comparison_matrix = generate_opportunistic_matrix()

    #Sequence preprocessing
    st = time.time()
    #Check of variant locations is supplied and if list of variants to include is supplied
    if args.variants_location and args.variants:
        variants = []
        for variant in args.variants:
            if variant == 'a' and 'Alpha' not in variants:
                variants.append('Alpha')
            elif variant == 'b' and 'Beta' not in variants:
                variants.append('Beta')
            elif variant == 'c' and 'Gamma' not in variants:
                variants.append('Gamma')
            elif variant == 'd' and 'Delta' not in variants:
                variants.append('Delta')
            elif variant == 'e' and 'Epsilon' not in variants:
                variants.append('Epsilon')
            elif variant == 'z' and 'Zeta' not in variants:
                variants.append('Zeta')
            elif variant == 'n' and 'Eta' not in variants:
                variants.append('Eta')
            elif variant == 'k' and 'Kappa' not in variants:
                variants.append('Kappa')
            elif variant == 'm' and 'Mu' not in variants:
                variants.append('Mu')
            elif variant == 'o' and 'Omicron' not in variants:
                variants.append('Omicron')
        sequences, lb, ub, feasible_amplicons, relevant_nucleotides = preprocess_sequences(sequences, args.search_width, variants_location=args.variants_location, variants=variants, amplicon_width=args.amplicon_width, misalign_threshold=args.misalign_threshold)
    
        with open(args.output + '/sequences_included_' + str(args.seed) + '.txt', 'w+') as f:
            for sequence in sequences:
                f.write(sequence.id + '\n')
        
        with open(args.output + '/runtimes_' + str(args.seed) + '.txt', 'a') as f:
            f.write('Variants considered:\n')
            for variant in variants:
                f.write(variant + '\n')
            f.write('Total sequences = ' + str(len(sequences)) + '\n')
        
    else:
        sequences, lb, ub, feasible_amplicons, relevant_nucleotides = preprocess_sequences(sequences, args.search_width, amplicon_width=args.amplicon_width, misalign_threshold=args.misalign_threshold)
        
        with open(args.output + '/sequences_included_' + str(args.seed) + '.txt', 'w+') as f:
            for sequence in sequences:
                f.write(sequence.id + '\n')
        
        with open(args.output + '/runtimes_' + str(args.seed) + '.txt', 'a') as f:
            f.write('Variants considered:\n')
            for variant in ['Alpha','Beta','Gamma','Delta','Epsilon','Zeta','Eta','Kappa','Mu','Omicron']:
                f.write(variant + '\n')
            f.write('Total sequences = ' + str(len(sequences)) + '\n')

    with open(args.output + '/runtimes_' + str(args.seed) + '.txt', 'a') as f:
        f.write('Time spent pre-processing sequences and determining feasible amplicons: ' + str(time.time() - st) + ', number of feasible amplicons: ' + str(len(feasible_amplicons)) + '\n')
    
    #Generate primer index
    st = time.time()
    PI = PrimerIndex.PrimerIndex.generate_index_mp(sequences, args.primer_width, comparison_matrix, processors=args.cores)
    PI.remove_redundant()
    with open(args.output + '/runtimes_' + str(args.seed) + '.txt', 'a') as f:
        f.write('Time spent generating primer index and filtering for feasible primers: ' + str(time.time() - st) + '\n')

    #Generate amplicons
    st = time.time()
    amplicons, diffs_per_amplicon = generate_amplicons_mp_hybrid(sequences, args.amplicon_width, comparison_matrix, amplicon_threshold=args.amplicon_threshold, feasible_amplicons=feasible_amplicons, relevant_nucleotides=relevant_nucleotides, processors=args.cores)
    with open(args.output + '/runtimes_' + str(args.seed) + '.txt', 'a') as f:
        f.write('Time spent generating amplicon differentiation ' + str(time.time() - st) + '\n')
        f.write('Total feasible amplicons: ' + str(len(amplicons)) + '\n')

    st = time.time()
    logs, amplicons, result_amplicons = greedy(sequences, amplicons, diffs_per_amplicon, args.primer_width, args.search_width, PI, comparison_matrix, args.amplicons, args.coverage, 5, logging=True, mulitplex=args.multiplex)
    with open(args.output + '/runtimes_' + str(args.seed) + '.txt', 'a') as f:
        f.write('Time spent running greedy algorithm: ' + str(time.time() - st) + '\n')
    
    #Run final optimization
    if args.multiplex:
        st = time.time()
        cur_primers = check_primer_feasibility(sequences, result_amplicons, PI, optimize=1, coverage=args.coverage)
        with open(args.output + '/runtimes_' + str(args.seed) + '.txt', 'a') as f:
            f.write('Time spent doing final primer optimization: ' + str(time.time() - st) + '\n')
            f.write('Forward primers\n')
            for fwd in cur_primers['forward']:
                f.write(fwd + '\n')
            f.write('Reverse primers\n')
            for rev in cur_primers['reverse']:
                f.write(rev + '\n')
    else:
        st = time.time()
        for amplicon in result_amplicons:
            cur_primers = check_primer_feasibility(sequences, [amplicon], PI, optimize=1, coverage=args.coverage)
            with open(args.output + '/runtimes_' + str(args.seed) + '.txt', 'a') as f:
                f.write('Amplicon: ' + str(amplicon.id) + '\n')
                f.write('Forward primers\n')
                for fwd in cur_primers['forward']:
                    f.write(fwd + '\n')
                f.write('Reverse primers\n')
                for rev in cur_primers['reverse']:
                    f.write(rev + '\n')
                
    with open(args.output + '/logfile_' + str(args.seed) + '.txt', 'w') as f:
        for line in logs:
            f.write(line + '\n')

def sm_test(sequences, amplicons):
    OG_copy = np.zeros((len(amplicons), len(sequences), len(sequences)), dtype=np.int8)

    shm = shared_memory.SharedMemory(create=True, size=(len(amplicons), len(sequences), len(sequences)))
    dst = np.ndarray(OG_copy.shape, dtype=np.int8, buffer=shm.buf)
    dst[:] = OG_copy[:]

    print(dst.nbytes / 10**9)

    sequence_pairs = list(itertools.combinations(sequences,2))
    seqs_partitioned = [sequence_pairs[i:i+ceil(len(sequence_pairs)/4)] for i in range(0, len(sequence_pairs), ceil(len(sequence_pairs)/4))]

    with mp.Pool(4) as pool:
        pool.starmap(AmpliconGeneration.shm_test, zip(itertools.repeat(amplicons), itertools.repeat(sequences), seqs_partitioned, itertools.repeat(shm.name), [1,2,3,4]))

    shm.close()
    shm.unlink()

def translate_to_numeric(sequences, amplicons, relevant_nucleotides, comparison_matrix):
    
    chars = ['a','c','t','g','u','r','y','k','m','s','w','b','d','h','v','n','-']
    char_comp = np.zeros((len(chars), len(chars)), dtype=np.int8)
    chars2num = {}
    seqs_num = np.zeros((len(sequences), sequences[0].length), dtype=np.int32)
    amplicons_num = np.zeros((len(amplicons)), dtype=np.int32)
    amplicons_lb = np.zeros((len(amplicons)), dtype=np.int32)
    amplicons_ub = np.zeros((len(amplicons)), dtype=np.int32)
    
    for char_index in range(len(chars)):
        chars2num[chars[char_index]] = char_index
        
    for c1 in range(len(chars)):
        for c2 in range(len(chars)):
            if not comparison_matrix[(chars[c1], chars[c2])][0]:
                char_comp[c1][c2] = 1
                
    for s in range(len(sequences)):
        for i in range(sequences[s].length):
            seqs_num[s][i] = chars2num[sequences[s].sequence[i]]
            
    for a in range(len(amplicons)):
        amplicons_num[a] = amplicons[a][0]
        cur = np.where(relevant_nucleotides < amplicons[a][0])[0]
        if cur.shape[0] > 0:
            amplicons_lb[a] = cur[-1]
        cur = np.where(relevant_nucleotides < amplicons[a][1])[0]
        if cur.shape[0] > 0:
            amplicons_ub[a] = cur[-1]
                
    return chars2num, char_comp, seqs_num, amplicons_num, amplicons_lb, amplicons_ub

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Run the greedy amplicon and primer selection algorithm.')
    #Input data
    parser.add_argument('sequences', type=str, help='File location of the (aligned) sequences')
    parser.add_argument('metadata', type=str, help='File location of the metadata for the sequences')
    #Output data
    parser.add_argument('output', type=str, help='Folder where results will be stored')
    #Amplicon parameters
    parser.add_argument('-aw', '--amplicon_width', default=200, type=int, help='Amplicon size')
    parser.add_argument('-at', '--amplicon_threshold', default=0, type=int, help='Number of allowed mismatches in an amplicon')
    parser.add_argument('-mis', '--misalign_threshold', default=5, type=int, help='Number of allowed misalign characters in an amplicon')
    #Primer parameters
    parser.add_argument('-pw', '--primer_width', default=25, type=int, help='Primer size')
    parser.add_argument('-sw', '--search_width', default=50, type=int, help='Search window for finding primers')
    parser.add_argument('-cov', '--coverage', default=1.0, type=float, help='Fraction of sequences that should be covered by both a forward and reverse primer for every amplicon')
    #Greedy parameters
    parser.add_argument('-amps', '--amplicons', default=5, type=int, help='Number of amplicons to find')
    parser.add_argument('-mp', '--multiplex', action='store_true', help='If supplied, will find primers that work in a single batch (multiplex)')
    #General parameters
    parser.add_argument('-c', '--cores', default=1, type=int, help='Number of processor cores to use when using multiprocessing')
    parser.add_argument('-vl', '--variants_location', default=None, type=str, help='File location containing the lineages per variant')
    parser.add_argument('-v', '--variants', default=None, type=str, help='Letters representing the variants to consider')
    parser.add_argument('-n', '--n_sequences', default=10**9, type=int, help='Number of sequences to include')
    parser.add_argument('-sd', '--seed', default=0, type=int, help='Seed to use when randomizing sequences')
    parser.add_argument('-spl', '--sequences_per_lineage', default=0, type=int, help='Minimum number of sequences a lineage should include')
    parser.add_argument('-ll', '--lineages_location', default=None, type=str, help='File location containing the number of sequences per lineage')
    #Runtype parameter
    parser.add_argument('-run', '--run_type', default='greedy', type=str, help='Type of program to run')
    args = parser.parse_args()
    
    if args.run_type == 'runtime_comparison':
        run_comparison(args)
    
    elif args.run_type == 'greedy':
        print('test')
        #run_greedy(args)



    '''
    sequences = generate_sequences('/Users/jaspervanbemmelen/Documents/Wastewater/source_code/amplivar/testing', '/Users/jaspervanbemmelen/Documents/Wastewater/source_code/amplivar/testing')
    #sequences = generate_sequences('/tudelft.net/staff-umbrella/SARSCoV2Wastewater/jasmijn/ref_sets_gisaid_2022_08_18/global_all_time_N0_L29000', '/tudelft.net/staff-umbrella/SARSCoV2Wastewater/jasmijn/ref_sets_gisaid_2022_08_18/global_all_time_N0_L29000')

    sequences = sequences[:500]
    sequences, lb, ub, feasible_amplicons, relevant_nucleotides = preprocess_sequences(sequences, 50, amplicon_width=400, misalign_threshold=10)
    
    seqs = [s.sequence for s in sequences]
    lineages = [s.lineage_num for s in sequences]
    ids = [s.id_num for s in sequences]

    amps = list(feasible_amplicons)
    amps.sort(key = lambda x : x[0])
    amps = amps[:]

    amplicon_width = 400
    M = generate_opportunistic_matrix()
    c, C, S, A, A_lb, A_ub = translate_to_numeric(sequences, amps, relevant_nucleotides, M)
    
    rel = relevant_nucleotides[relevant_nucleotides < amps[-1][1]]

    runtimes = 0

    #st = time.time()
    #amps1 = generate_amplicons_mp_smartest(sequences, 400, M, amplicon_threshold=1, feasible_amplicons=amps, relevant_nucleotides=rel, processors=4)
    #print('New runtime: ' + str(time.time() - st))

    #st = time.time()
    #generate_amplicons_mp_exp_cy(sequences, 400, M, amplicon_threshold=1, feasible_amplicons=amps, processors=4)
    #print('Old runtime: ' + str(time.time() - st))

    st = time.time()
    amps2, X = generate_amplicons_mp_hybrid(sequences, 400, M, amplicon_threshold=1, feasible_amplicons=amps, relevant_nucleotides=rel, processors=4)
    print('Newest runtime: ' + str(time.time() - st))

    same = True
    if len(amps1) != len(amps2):
        same = False
        print('Unequal lengths: ' + str(len(amps1)) + ', ' + str(len(amps2)))
    for i in range(len(amps1)):
        if amps1[i].differences != amps2[i].differences:
            print(str(len(amps1[i].differences)) + ', ' + str(len(amps2[i].differences)))
            same = False
    print(same)

    st = time.time()
    amps1 = AmpliconGeneration.determine_differences_cy(amps, seqs, lineages, ids, 1, M)
    print('Initial method: ' + str(time.time() - st))
    print(len(amps1))
    

    st = time.time()
    amps1 = generate_amplicons_mp_exp_cy(sequences[:500], 400, M, lb=None, ub=None, amplicon_threshold=1, feasible_amplicons=amps, processors=2)
    print('Split on amplicons: ' + str(time.time() - st))
    
    st = time.time()
    res = []
    amps2 = AmpliconGeneration.generate_amplicons_smarter_cy(A, A_lb, A_ub, amplicon_width, A.shape[0], S, S.shape[0], lineages, C, rel, rel.shape[0], 1)
    for amp in amps2:
        res.append(Amplicon(amp[0], amp[1]))
        res[-1].differences = amps2[amp]
    print('New method: ' + str(time.time() - st))
    print(len(amps2))

    same = True
    print(len(amps1) == len(res))
    for i in range(len(amps1)):
        if amps1[i].differences != res[i].differences:
            same = False
    print(same)


    for a in amps1:
        if amps1[a] != amps2[a]:
            same = False
            print(a)
            print(amps1[a] - amps2[a])
            print('-'*200)
            print(amps2[a] - amps1[a])
            #print(str(len(amps1[a])) + ', ' + str(len(amps2[a])))
            break
    print(same)
    '''
    '''
    n_seqs = 100
    n_amps = 2000
    amplicon_threshold = 1
    comparison_matrix = generate_opportunistic_matrix()

    seqs = [s.sequence for s in sequences]
    amplicons = list(feasible_amplicons)
    amplicons.sort(key = lambda x : x[0])
    lineages = [s.lineage_num for s in sequences]

    ids = [s.id_num for s in sequences]

    chars = ['a','c','t','g','u','r','y','k','m','s','w','b','d','h','v','n','-']

    st = time.time()
    #amps1 = AmpliconGeneration.determine_differences_cy(amplicons[:n_amps], seqs[:n_seqs], lineages[:n_seqs], ids[:n_seqs], 1, comparison_matrix)
    amps1 = generate_amplicons_mp_exp_cy(sequences[:n_seqs], 200, comparison_matrix, lb=None, ub=None, amplicon_threshold=1, feasible_amplicons=amplicons[:n_amps], processors=4)
    print('Split on amplicons: ' + str(time.time() - st))

    st = time.time()
    amps2 = generate_amplicons_mp_smart(sequences[:n_seqs], 200, comparison_matrix, lb=None, ub=None, amplicon_threshold=1, feasible_amplicons=amplicons[:n_amps], relevant_nucleotides=relevant_nucleotides, processors=4)
    print('Split on sequences: ' + str(time.time() - st))
    """
    st = time.time()
    #amps2 = AmpliconGeneration.determine_differences_cy2(amplicons[:n_amps], seqs[:n_seqs], lineages[:n_seqs], ids[:n_seqs], 1, comparison_table, np.where(final_to_check == 1)[0])
    amps2 = generate_amplicons_mp_exp_cy2(sequences[:n_seqs], 200, comparison_matrix, lb=None, ub=None, amplicon_threshold=1, feasible_amplicons=amplicons[:n_amps], relevant_nucleotides=relevant_nucleotides, processors=2)
    print(time.time() - st)
    """
    print(len(amps1) == len(amps2))
    same = True
    for i in range(len(amps1)):
        if not amps1[i].differences == amps2[i].differences:
            print(amps1[i].differences - amps2[i].differences)
            print(amps2[i].differences - amps1[i].differences)
            print(len(amps1[i].differences))
            print(len(amps2[i].differences))
            same = False
            break
    print(same)
    '''