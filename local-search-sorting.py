import random as rd
import time

# swap 2 random elements in seq
def perturb(seq):
    swapped = seq.copy()
    a, b = rd.sample(range(len(seq)), 2)
    swapped[a], swapped[b] = swapped[b], swapped[a]
    return swapped

# fitness score of seq
def fitness(seq):
    return sum(1 for i in range(len(seq)) for j in range(i + 1, len(seq)) if seq[i] <= seq[j])

# find best neighbour seq from all swaps
def best_neighbor_all_pairs_swap(seq):
    curr_best_nbor = seq
    curr_best_val = fitness(seq)
    for i in range(len(seq)):
        for j in range(i + 1, len(seq)):
            swapped = seq.copy()
            swapped[i], swapped[j] = swapped[j], swapped[i]
            swapped_fitness = fitness(swapped)
            if swapped_fitness > curr_best_val:
                curr_best_nbor = swapped
                curr_best_val = swapped_fitness
    return curr_best_nbor

# find best neighbour seq by considering swaps for the first k elements
def best_neighbor_k_first_swaps(seq, k):
    curr_best_nbor = seq
    curr_best_val = fitness(seq)
    for i in range(k):
        for j in range(i + 1, len(seq)):
            swapped = seq.copy()
            swapped[i], swapped[j] = swapped[j], swapped[i]
            swapped_fitness = fitness(swapped)
            if swapped_fitness > curr_best_val:
                curr_best_nbor = swapped
                curr_best_val = swapped_fitness
    return curr_best_nbor

# sort seq using local search
def sorting_local_search(in_seq):
    current_sol = in_seq
    best_sol = current_sol

    while True:
        current_sol = perturb(current_sol)

        while True:
            candidate = best_neighbor_k_first_swaps(current_sol, 2)
            if fitness(candidate) > fitness(current_sol):
                current_sol = candidate
            else:
                break

        if fitness(current_sol) > fitness(best_sol):
            best_sol = current_sol
        else:
            break

    return best_sol

sequence = list(range(5, 0, -1))
print("Fitness score before sorting: ", fitness(sequence))
print("Sorting sequence...")

tic = time.perf_counter()
sorted_sequence = sorting_local_search(sequence)
toc = time.perf_counter()

print(f"Sorted sequence in {toc - tic:0.4f} seconds")
print("Fitness score after sorting: ", fitness(sorted_sequence))
print("Maximum possible score: ", 0.5 * len(sorted_sequence) * (len(sorted_sequence) - 1))
