import parameters as par
import random as rand
import copy
import NSGAII
import math
import Class
import MOEAD

clonepopulation = []
Archive = []
DE_rate = 0.5
sort_fit_ix = 0
col_ix = 0
def HEIA_Calc(col, PF_storage, i_fes):
    global col_ix
    col_ix = col.index - 1
    clonesize = int(col.size / 5)

    if (col.was_erased):
        HEIA_Pop_Init(col)
        col.was_erased = False

    SBXpop = []
    DEpop = []
    offspringpopulation = Clone_Operation(col)

    if (len(offspringpopulation) == 0):
        raise SystemError("error")

    for i in range(0, len(offspringpopulation)):
        if (par.HEIA_crossover_prob < rand.random()):
            SBXpop.append(offspringpopulation[i])
        else:
            DEpop.append(offspringpopulation[i])

    offspringpopulation = []

    nfes = DE_Update(DEpop,col, i_fes)
    nfes += SBX_Update(SBXpop,col, i_fes + nfes)

    offspringpopulation = DEpop + SBXpop

    Archive_Update(offspringpopulation,col.size,clonesize,col.fit_index[0])

    if (len(offspringpopulation) != col.size):
        raise SystemError("error")

    col.indiv = copy.deepcopy(offspringpopulation)

    col.Fitness_Calc_Col()
    PF_storage[:] = copy.deepcopy(Archive[col_ix])


    return nfes

def HEIA_Init(n_col, pop):
    global clonepopulation
    global Archive


    Archive = []
    clonepopulation = []
    for i in range(0,n_col):
        temp =[]
        Archive.append(copy.deepcopy(temp))
    clonepopulation = copy.deepcopy(Archive)
    pop.Fitness_Calc()

    return pop.size


def HEIA_Pop_Init(col):
    global col_ix
    col_ix = col.index - 1

    global Archive
    global clonepopulation
    global sort_fit_ix
    Archive[col_ix] = []
    clonepopulation[col_ix] = []

    clonesize = int(col.size / 5)

    col_pop = col.indiv

    NSGAII.Rank_Crowding_Distance_Assign(col_pop,col.fit_index[0])

    index = 0
    sort_fit_ix = 0
    col_pop.sort(key = Sort_Fit)

    for i in range(0,len(col_pop)):
        if (col_pop[i].rank == 1):
            col_pop[i].table = []
            col_pop[i].table.append(index)
            index += 1
            Archive[col_ix].append(copy.deepcopy(col_pop[i]))


    temp_arch = copy.deepcopy(Archive[col_ix])

    temp_arch.sort(key = Sort_Crowd, reverse = True)
    for i in range(0,len(Archive[col_ix])):
        if (i >= clonesize):
            break
        clonepopulation[col_ix].append(copy.deepcopy(temp_arch[i]))



def Clone_Operation(col):
    size = col.size
    offspring = []
    parents = copy.deepcopy(clonepopulation[col_ix])
    parents_size = len(parents)
    min_dist = 0.0
    max_dist = 1.0
    sum_dist = 0.0
   
    for k in range (0,parents_size):
        if (parents[k].crowd_dist != 1.0e14):
            max_dist = parents[k].crowd_dist
            min_dist = parents[-1].crowd_dist
            for l in range (0,k):
                parents[l].crowd_dist = parents[k].crowd_dist
            break

    if (parents[0].crowd_dist == 1.0e14):
        for l in range (0,parents_size):
                parents[l].crowd_dist = 1.0

    for k in range (0,parents_size):
        sum_dist += parents[k].crowd_dist

    clones = [0.0]*parents_size

    for k in range (0,parents_size):
        if (sum_dist == 0):
            clones[k] = math.ceil(float(size)/parents_size)
        else:
            clones[k] = math.ceil(size*parents[k].crowd_dist / sum_dist)

    remain = size
    i = 0

    for k in range (0,parents_size):
        for l in range(0,clones[k]):
            if (remain > 0):
                offspring.append(copy.deepcopy(parents[k]))
                remain -= 1
            i += 1
        if (remain == 0):
            break
    if (len(offspring) == 0):
        raise SystemError("error")

    return offspring

def DE_Update(DEpop, col, i_fes):

    clonepop_size = len(clonepopulation[col_ix])
    nfes = 0
    for i in range(0, len(DEpop)):
    
        if (nfes + i_fes >= par.max_iterations):
            break
        parents = []
        parents.append(DEpop[i])
        if (clonepop_size < 20):
            if (clonepop_size > 1):
                perm = []
                Random_Perm(perm, clonepop_size)
                parents.append(clonepopulation[col_ix][perm[0]])
                parents.append(clonepopulation[col_ix][perm[1]])
            else:
                parents.append(clonepopulation[col_ix][0])
                parents.append(clonepopulation[col_ix][0])
        else:
            if (0.1 < rand.random()):            
                neighbours = DEpop[i].table[0]
                perm = []
                Random_Perm(perm, 20)
                selected = perm[0]
                selected2 = perm[1]
                arch_size = len(Archive[col_ix])
                if (neighbours < 10):
                    parents.append(Archive[col_ix][selected2])
                    parents.append(Archive[col_ix][selected])
                elif (neighbours > (arch_size - 10)):
                    parents.append(Archive[col_ix][arch_size - 20 + selected2])
                    parents.append(Archive[col_ix][arch_size - 20 + selected])
                else:                
                    parents.append(Archive[col_ix][neighbours - 10 + selected2])
                    parents.append(Archive[col_ix][neighbours - 10 + selected])
            else:
                perm = []
                Random_Perm(perm, clonepop_size)
                parents.append(clonepopulation[col_ix][perm[0]])
                parents.append(clonepopulation[col_ix][perm[1]])
        offspring = MOEAD.Diff_Evo_Xover_B(parents[0], parents[1], parents[2], DE_rate, col.f_code, col.c_code)
        col.Mutation_Indi(offspring)
        offspring.Fitness_Calc(col.f_code)
       
        nfes+=1
        DEpop[i] = copy.deepcopy(offspring)
    return nfes
   
def SBX_Update(SBXpop, col, i_fes):
    nfes = 0
    for i in range(0, len(SBXpop)):
        if (nfes + i_fes >= par.max_iterations):
            break

        parents = []
        parents.append(SBXpop[i])
        random_index = rand.randint(0, len(clonepopulation[col_ix]) - 1)
        parents.append(clonepopulation[col_ix][random_index])
        offspring = col.Crossover_NSGAII(parents)[0]
        col.Mutation_Indi(offspring)
        offspring.Fitness_Calc(col.f_code)
       
        nfes+=1
        SBXpop[i] = copy.deepcopy(offspring)

    return nfes


def Archive_Update(merge_pop, popsize, clonesize, fit_indexes):
    global Archive
    global clonepopulation
    global sort_fit_ix
    solution_union = copy.deepcopy(Archive[col_ix])+copy.deepcopy(merge_pop)

    Suppress(solution_union, fit_indexes)

    NSGAII.Rank_Crowding_Distance_Assign(solution_union, fit_indexes)

    Archive[col_ix] = []
    clonepopulation[col_ix] = []

    front = []
    for i in range (0, len(solution_union)):
        if (solution_union[i].rank == 1):
            front.append(solution_union[i])

    Suppress(front, fit_indexes)
    front.sort(key = Sort_Crowd, reverse = True)

    while (len(front)> popsize):
        front.pop(-1)
        Crowding_Dist_Update(front, fit_indexes)
        front.sort(key = Sort_Crowd, reverse = True)
    
    sort_fit_ix = 0
    front.sort(key = Sort_Fit)
    for i in range(0,len(front)):    
        front[i].table = []
        front[i].table.append(i)
        Archive[col_ix].append(copy.deepcopy(front[i]))
    

    front.sort(key = Sort_Crowd, reverse = True)
    for i in range(0,len(front)):
        if (i >= clonesize):
            break
        clonepopulation[col_ix].append(copy.deepcopy(front[i]))
        

def Crowding_Dist_Update(front, fit_indexes):
    global sort_fit_ix
    nobj = len(fit_indexes)
    size = len(front)

    if (size <= 2):    
        for i in range(0,size):
            front[i].crowd_dist = 1.0e14
        return
    

    for i in range(0,size):
        front[i].crowd_dist = 0.0

    for i in range(0,nobj):
    
        fit_ix = fit_indexes[i] - 1  
        sort_fit_ix = fit_ix
        front.sort(key = Sort_Fit)
        objetiveMinn = front[0].fitness[fit_ix]
        objetiveMaxn = front[size - 1].fitness[fit_ix]
           
        front[0].crowd_dist = 1.0e14
        front[size - 1].crowd_dist = 1.0e14


        for j in range(1,size-1):        
            distance = front[j + 1].fitness[fit_ix] - front[j - 1].fitness[fit_ix] 
            if (objetiveMaxn != objetiveMinn):
                distance = distance / (objetiveMaxn - objetiveMinn)
            distance += front[j].crowd_dist
            front[j].crowd_dist = distance
         

def Random_Perm(perm, size):
    index = []
    flag = []

    for i in range(0,size):    
        index.append(i)
        flag.append(True)
        perm.append(0)
    
    num = 0
    while (num < size):    
        start = rand.randint(0, size - 1)
        while (True):        
            if (flag[start]):            
                perm[num] = index[start]
                flag[start] = False
                num+=1
                break
            
            if (start == (size - 1)):
                start = 0
            else:
                start+=1


def Suppress(pop, fit_indexes):

    nobj = len(fit_indexes)
    pop_size = len(pop)
    k = 0
    while (k < pop_size):
        l = k+1
        while (l < pop_size):        
            m = 0
            while (m < nobj):            
                fit_ix = fit_indexes[m] - 1
                diff = abs(pop[k].fitness[fit_ix] - pop[l].fitness[fit_ix])
                if (diff>0.000001):
                    break
                m+=1
                
            if (m == nobj):
                pop.pop(l)
                l-=1
                pop_size-=1
            l +=1
        k +=1
            
        
    


def Sort_Fit(c1):
    return c1.fitness[sort_fit_ix]


def Sort_Crowd(c1):
    return c1.crowd_dist



