import random
from math import *
import Add_func
import copy

def NSGAII_Calc(collective):
    if (collective.was_erased):
        Rank_Crowding_Distance_Assign(collective.indiv, collective.fit_index[0])
        collective.was_erased = False

    off_pop = collective.Crossover_NSGAII(Select(collective.indiv,collective.fit_index[0]))

    collective.Mutation_NSGAII(off_pop)


    mix_pop = off_pop + collective.indiv


    off_pop = []

    Nondominated_Sort_Fill(mix_pop, off_pop, collective.size, collective.fit_index[0])

    collective.indiv = copy.deepcopy(off_pop)

    nfes = collective.size

    #collective.save()

    collective.Fitness_Calc_Col()

    return nfes



def Select(old_pop, fit_indexes):
    pop_size = len(old_pop)
    a1 = []
    a2 = []
    for i in range(0,pop_size):
        a1.append(i)
        a2.append(i)

    mult_size = 0

    while (len(a1)%4 != 0):
        a1.append(random.randint(0,pop_size-1))
        a2.append(random.randint(0,pop_size-1))
        mult_size+=1

    a_size = len(a1)

    for i in range(0,a_size):
        rand = random.randint(0,a_size-1)
        temp = a1[rand]
        a1[rand] = a1[i];
        a1[i] = temp;
        rand = random.randint(0,a_size-1)
        temp = a2[rand]
        a2[rand] = a2[i];
        a2[i] = temp;
    off_pop = []
    for i in range(0,pop_size, 4):
        off_pop.append(copy.deepcopy(Tournament(old_pop[a1[i]], old_pop[a1[i + 1]], fit_indexes)))
        off_pop.append(copy.deepcopy(Tournament(old_pop[a1[i + 2]], old_pop[a1[i + 3]], fit_indexes)))
        off_pop.append(copy.deepcopy(Tournament(old_pop[a2[i]], old_pop[a2[i + 1]], fit_indexes)))
        off_pop.append(copy.deepcopy(Tournament(old_pop[a2[i + 2]], old_pop[a2[i + 3]], fit_indexes)))

    for i in range(0,mult_size):
        off_pop.pop(-1)

    if (len(off_pop) != pop_size):
        raise SystemError("error")

    return off_pop

def Tournament(ind1, ind2, fit_indexes):
    flag = Add_func.Dominance_check(ind1, ind2, fit_indexes)


    if (flag == 1):
        return ind1
    elif (flag == -1):
        return ind2
    if (ind1.crowd_dist > ind2.crowd_dist):
        return ind1
    elif (ind1.crowd_dist < ind2.crowd_dist):
        return ind2

    if (random.random() <= 0.5):
        return ind1
    else:
        return ind2
    

def Rank_Crowding_Distance_Assign(population, fit_index):
    
    pop_size = len(population)

    assigned_rank = [0] * pop_size
    rank = 1;

    while(True):
        dominated = [-1]* pop_size
        count = 0;
        
        for i in range(0,pop_size):
            if (dominated[i] == 0 or assigned_rank[i] != 0):
                continue
            for j in range (i+1, pop_size):
                if (dominated[j] == 0 or assigned_rank[j] != 0):
                    continue
                flag = Add_func.Dominance_check(population[i],population[j],fit_index)
                if (flag == 0):
                    continue
                elif (flag == 1):
                    dominated[j] = 0
                    count +=1
                else:
                    dominated[i] = 0
                    count +=1
                    break
        if (count == 0):
            break
        curr_rank_pop = []
        for i in range(0,pop_size):
            if (dominated[i] == -1 and assigned_rank[i] == 0):
                population[i].rank = rank
                assigned_rank[i] = rank
                curr_rank_pop.append(population[i]);


        Crowding_Distance_Assign(curr_rank_pop)
        rank +=1


def Nondominated_Sort_Fill(mixed_pop, off_pop, size, fit_index):

   
    rank = 1
    size_mixed = len(mixed_pop)

    assigned_rank = [0] * size_mixed 

    while (len(off_pop) < size):
        dominated = [-1]* size_mixed
        for i in range(0,size_mixed):
            if (dominated[i] == 0 or assigned_rank[i] != 0):
                continue
            for j in range (i+1, size_mixed):
                if (dominated[j] == 0 or assigned_rank[j] != 0):
                    continue
                flag = Add_func.Dominance_check(mixed_pop[i],mixed_pop[j], fit_index)
                if (flag == 0):
                    continue
                elif (flag == 1):
                    dominated[j] = 0
                else:
                    dominated[i] = 0
                    break
        curr_rank = []
        for i in range(0,size_mixed):
            if (dominated[i] == -1 and assigned_rank[i] == 0):
                mixed_pop[i].rank = rank
                assigned_rank[i] = rank
                curr_rank.append(mixed_pop[i]);

        Crowding_Distance_Assign(curr_rank)
        
        if (len(curr_rank) + len(off_pop) <= size):
           
           for i in range (0, len(curr_rank)):
               off_pop.append(curr_rank[i])
           rank += 1
        else:
            Crowding_Fill(off_pop,  curr_rank, size)


    if (len(off_pop) != size):
        raise SystemError("error")


        

def Crowding_Fill(off_pop,  curr_rank, size):

    amt = size - len(off_pop)
    dist_size = len(curr_rank)
    dist = []
    for i in range(0,dist_size):
        dist.append([i, curr_rank[i].crowd_dist])

    Dist_Quicksort(dist)

    for i in range(0, amt):
        off_pop.append(curr_rank[dist[dist_size - i - 1][0]])



def Crowding_Distance_Assign(curr_rank):
    front_size = len(curr_rank)
    nobj = len(curr_rank[0].fitness)

    if (front_size <=2):
        for i in range(0,front_size):
            curr_rank[i].crowd_dist = 1.0e14
        return

    dist = []
    for i in range(0, front_size):
        dist.append(i)



    obj_array = []
    for i in range(0, nobj):
        obj_array.append(dist)
        Front_Obj_Quicksort(curr_rank, i, obj_array[i], front_size)
    
    for j in range(0, front_size):
        curr_rank[dist[j]].crowd_dist = 0
    
    for i in range(0, nobj):
    
        curr_rank[obj_array[i][0]].crowd_dist = 1.0e14;
    
    for i in range(0, nobj):
    
        for j in range(1, front_size-1):
        
            if (curr_rank[obj_array[i][j]].crowd_dist != 1.0e14):
                if (curr_rank[obj_array[i][front_size - 1]].fitness[i] != curr_rank[obj_array[i][0]].fitness[i]):
                    temp = curr_rank[obj_array[i][j]].crowd_dist + (curr_rank[obj_array[i][j + 1]].fitness[i] - curr_rank[obj_array[i][j - 1]].fitness[i]) / (curr_rank[obj_array[i][front_size - 1]].fitness[i] - curr_rank[obj_array[i][0]].fitness[i]);
                    
                    curr_rank[obj_array[i][j]].crowd_dist = temp;
        
    for j in range(0, front_size):
    
        if (curr_rank[dist[j]].crowd_dist != 1.0e14):
        
            curr_rank[dist[j]].crowd_dist = (curr_rank[dist[j]].crowd_dist) / nobj;


def Front_Obj_Quicksort(curr_rank, i, obj_array, obj_array_size):
    Front_Obj_Quicksort_Actual(curr_rank, i, obj_array, 0, obj_array_size - 1);


def Front_Obj_Quicksort_Actual(pop, objcount, obj_array, left, right):
   
    if (left<right):
        index = random.randint(left, right);
        temp = obj_array[right];
        obj_array[right] = obj_array[index];
        obj_array[index] = temp;
        pivot = pop[obj_array[right]].fitness[objcount];
        i = left - 1;
        for j in range(left,right):
            if (pop[obj_array[j]].fitness[objcount] <= pivot):
                i += 1;
                temp = obj_array[j];
                obj_array[j] = obj_array[i];
                obj_array[i] = temp;
            
        index = i + 1;
        temp = obj_array[index];
        obj_array[index] = obj_array[right];
        obj_array[right] = temp;
        Front_Obj_Quicksort_Actual(pop, objcount, obj_array, left, index - 1);
        Front_Obj_Quicksort_Actual(pop, objcount, obj_array, index + 1, right);

def Dist_Quicksort(dist):
    Dist_Quicksort_Actual(dist, 0, len(dist)- 1);


def Dist_Quicksort_Actual(dist, left, right):
    if (left<right):
    
        index = random.randint(left, right);
        temp = dist[right];
        dist[right] = dist[index];
        dist[index] = temp;
        pivot = dist[right][1];
        i = left - 1;
        for j in range(left,right):
        
            if (dist[j][1] <= pivot):
                i += 1;
                temp = dist[j];
                dist[j] = dist[i];
                dist[i] = temp;
            
        
        index = i + 1;
        temp = dist[index];
        dist[index] = dist[right];
        dist[right] = temp;
        Dist_Quicksort_Actual(dist, left, index - 1);
        Dist_Quicksort_Actual(dist, index + 1, right);
    