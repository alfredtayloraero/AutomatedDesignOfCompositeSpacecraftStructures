import parameters as par
import math
import random as rand
import copy
import Class
import Add_func
import parameters

idealpoint = []
nadirpoint = []


def MOEAD_Calc(collective, pareto_front, i_fes):
    if (collective.was_erased):
        Neighbourhood_Init(collective)
        collective.was_erased = False

    nfes = MOEAD_Evolve(collective, pareto_front, i_fes)

    collective.Fitness_Calc_Col()

    return nfes

def MOEAD_Init(n_obj, population):
    n_var = len(population.indiv[0].code)

    population.Fitness_Calc()

    psize = population.size

    weight_vector = Add_func.Uniform_Weights_Generate(psize, n_obj, 2)

    for i in range(0, psize):       #psize
        population.indiv[i].namda = copy.deepcopy(weight_vector[i])         #weight vector is different length than the population
        population.indiv[i].saved_fitness = copy.deepcopy(population.indiv[i].fitness)

    return psize


def Neighbourhood_Init(collective):
    csize = collective.size
    dist = [0.0]*csize
    idx = [0]*csize
    nobj = collective.f_code.num_objs

    MOEAD_niche = int(csize * par.MOEAD_niche_multi)

    if (MOEAD_niche <= 3):
        MOEAD_niche = 4

    for i in range(0,csize):
        if (len(collective.indiv[i].saved_fitness) != nobj):
            collective.indiv[i].saved_fitness = copy.deepcopy(collective.indiv[i].fitness)

        collective.indiv[i].table = []
        for j in range(0,csize):
            dist[j] = Dist_Vector(collective.indiv[i].namda,collective.indiv[j].namda) # this sucks
            idx[j] = j;

        Min_Fast_Sort(dist,idx, csize, MOEAD_niche)

        for k in range(0,MOEAD_niche):
            collective.indiv[i].table.append(idx[k])

def Utility_Comp(collective, i_fes):

    csize = collective.size
    fit_ix = collective.fit_index[0]
    MODE = collective.mode

    for n in range(0,csize):
        f1 =  Fitness_Function(collective.indiv[n].fitness, collective.indiv[n].namda,fit_ix,i_fes, MODE)

    f2 =  Fitness_Function(collective.indiv[n].saved_fitness, collective.indiv[n].namda,fit_ix,i_fes, MODE)

    delta = (f2-f1)/f2
    if (delta>0.001):
        collective.indiv[n].utility = 1.0
    else:
        collective.indiv[n].utility *= (0.95 + 0.05*delta/0.001)

    collective.indiv[n].saved_fitness = copy.deepcopy(collective.indiv[n].fitness)

def Dist_Vector(vec1, vec2):
    dim = len(vec2)
    sum = 0
    for n in range(0,dim):
        sum += math.pow(vec1[n]-vec2[n],2)

    return math.sqrt(sum)


def Min_Fast_Sort(x, idx , n, m):
    for i in range (0,m):
        for j in range (i+1,n):
            if(x[i]>x[j]):
                temp = x[i]
                x[i] = x[j]
                x[j] = temp
                id = idx[i];
                idx[i] = idx[j];
                idx[j] = id;

def MOEAD_Evolve(col, pareto_front, i_fes):
    order = []
    nfes = 0
    Tournament(10,order,col.size,col)

    for sub in range(0,len(order)):
        c_sub = order[sub]

        rnd = rand.random()

        if (rnd < par.MOEAD_mating_chance):
            type = 1
        else:
            type = 2

        plist = []

        Mate_Selection(plist, c_sub, 2, type, col)

        rate2 = 0.5
        
        child = Diff_Evo_Xover_B(col.indiv[c_sub], col.indiv[plist[0]], col.indiv[plist[1]], rate2, col.f_code, col.c_code)


        plist = []

        col.Mutation_Indi(child)

        col.Fitness_Calc_Indi(child)

        Idealpoint_Update(child.fitness)

        if (col.mode == "MOEADMSF" or col.mode == "MOEADPSF"):
            Problem_Update_Global(child, col, c_sub, type, i_fes+nfes, pareto_front)
        else:
            Problem_Update(child, col, c_sub, type, i_fes+nfes, pareto_front)
        nfes+=1

        if (nfes + i_fes >= par.max_iterations):
            break

    return nfes

def Tournament(depth, selected, csize, col):
    candidate = []
    nobj = col.f_code.num_objs
    for k in range(0,nobj):
        selected.append(k)
    for k in range(nobj,csize):
        candidate.append(k)

    while (len(selected) < int(float(csize)/par.MOEAD_new_sol_multi)):
        best_idd = rand.randint(0,len(candidate)-1)
        best_sub = candidate[best_idd]

        for i in range(1,depth):
            i2 = rand.randint(0,len(candidate)-1)
            s2 = candidate[i2]
            if (col.indiv[s2].utility > col.indiv[best_sub].utility):
                best_idd = i2
                best_sub = s2
        selected.append(best_sub)
        candidate.pop(best_idd)

def Mate_Selection(list, cid, size, type, col):
    ss = len(col.indiv[cid].table)

    while (len(list) < size):
        if (type == 1):
            id = rand.randint(0,ss-1)
            parent = col.indiv[cid].table[id]
        else:
            parent = rand.randint(0,col.size-1)
        flag = True

        for i in range (0, len(list)):
            if (list[i] == parent):
                flag = False
                break

        if (flag):
            list.append(parent)


def Diff_Evo_Xover_B(ind0, ind1, ind2, rate,fcode,ccode):
    nvar = fcode.num_vars
    idx_rnd = rand.randint(0,nvar-1)

    child_code = [0.0] * nvar
    
    for n in range (0,nvar):
        upper_b = fcode.bound_max[n]
        lower_b = fcode.bound_min[n]

        rnd1 = rand.random()
        if (rnd1 < 1.0 or n == idx_rnd):
            child_code[n] = ind0.code[n] + rate*(ind2.code[n] - ind1.code[n])
        else:
            child_code[n] = ind0.code[n]

        if (child_code[n] < lower_b):
            rnd = rand.random()
            child_code[n] = lower_b + rnd*(ind0.code[n] - lower_b)
        elif (child_code[n] > upper_b):
            rnd = rand.random()
            child_code[n] = upper_b - rnd*(upper_b - ind0.code[n])

    out_ind = Class.Individual(fcode)
    out_ind.code = copy.deepcopy(child_code)
    return out_ind

def Problem_Update(indiv,col,id,type,i_fes, pareto_front):

    time = 0
    save_temp = 0

    if (type == 1):
        size = len(col.indiv[id].table)
    else:
        size = col.size

    perm = []

    for k in range(0,size):
        perm.append(k)

    Permutation(perm)

    MOEAD_limit = col.size * par.MOEAD_limit_multi
    if (MOEAD_limit == 1):
        MOEAD_limit = 2

    fit_ix = col.fit_index[0]
    MODE = col.mode
    for i in range(0,size):
        if (type == 1):
            k = col.indiv[id].table[perm[i]]
        else:
            k = perm[i]

        f1 =  Fitness_Function(col.indiv[k].fitness, col.indiv[k].namda, fit_ix, i_fes, MODE)
        f2 =  Fitness_Function(indiv.fitness, col.indiv[k].namda, fit_ix, i_fes, MODE)

        if (f2 < f1):
            indiv.namda = copy.deepcopy(col.indiv[k].namda)
            indiv.table = copy.deepcopy(col.indiv[k].table)
            indiv.saved_fitness = copy.deepcopy(col.indiv[k].saved_fitness)
            
            col.indiv[k] = copy.deepcopy(indiv)
            if (save_temp == 0 and (MODE == "MOEAD" or MODE == "MOEADMSF" or MODE == "MOEADPSF")):
            
                pareto_front.append(copy.deepcopy(indiv))
                save_temp+=1
            
            
            time+=1;

        if (time >= MOEAD_limit and type != 3):
            perm = []
            return
    perm = []




def Problem_Update_Global(indiv,col,id,type,i_fes, pareto_front):
    pops = col.size
    time = 0
    save_temp = 0

    x = [0.0]*pops
    idx = [0]*pops

    MOEAD_limit = pops * par.MOEAD_limit_multi

    niche = len(col.indiv[0].table)

    fit_ix = col.fit_index[0]
    MODE = col.mode
    
    for i in range(0,pops):
        idx[i] = i
        x[i] = Fitness_Function(indiv.fitness, col.indiv[i].namda, fit_ix, i_fes, MODE)

    Min_Fast_Sort(x, idx, pops, niche)

   
    for j in range(0,niche):
        
        fj =  Fitness_Function(col.indiv[idx[j]].fitness, col.indiv[idx[j]].namda, fit_ix, i_fes, MODE)
        

        if (fj > x[j]):
            indiv.namda = copy.deepcopy(col.indiv[idx[j]].namda)
            indiv.table = copy.deepcopy(col.indiv[idx[j]].table)
            indiv.saved_fitness = copy.deepcopy(col.indiv[idx[j]].saved_fitness)
            
            col.indiv[idx[j]] = copy.deepcopy(indiv)
            if (save_temp == 0 and (MODE == "MOEAD" or MODE == "MOEADMSF" or MODE == "MOEADPSF")):
            
                pareto_front.append(copy.deepcopy(indiv))
                save_temp+=1
            
            
            time+=1;

        if (time >= MOEAD_limit and type != 3):
            break

def Permutation(perm):
    size = len(perm)

    index = []
    for i in range(0,size):
        index.append(i)

    for i in range(0,size):
        j = rand.randint(0,size-1)
        while(True):
            if (index[j]>=0):
                perm[i] = index[j]
                index[j] = -1
                break
            elif(j == size - 1):
                j = 0
            else:
                j +=1


def Fitness_Function(fit,namda, fit_indexes, nfes, MODE):
    max_fun = 0

    if (MODE == "MOEAD" or MODE == "BCE" ):
    
        max_fun = Fitness_Function_TCH(fit, namda, fit_indexes)
    
    elif (MODE == "MOEADPSF"):
    
        max_fun = Fitness_Function_PSF(fit, namda, fit_indexes, nfes)
    
    elif (MODE == "MOEADMSF"):
    
        max_fun = Fitness_Function_MSF(fit, namda, fit_indexes, nfes)
    
    else:
        raise SystemError("error")

    return max_fun


def Fitness_Function_TCH(fit, namda, fit_indexes):

    max_fun = -1.0e+30

    nobj = len(fit_indexes)


    for n in range(0,nobj):
        fit_ix = fit_indexes[n]-1
        diff = abs(fit[fit_ix]-idealpoint[fit_ix])

        if(namda[fit_ix] == 0.0):
            feval = 0.0001*diff
        else:
            feval = diff*namda[fit_ix]
        if (feval > max_fun):
            max_fun = feval
    return max_fun

def Fitness_Function_PSF(fit, namda, fit_indexes, nfes):
    max_fun = -1.0e+30
    min_fun = 1.0e+30
    au = 1
    wmax = -1.0e+30
    wmin = 1.0e+30
    eps = 1.0e-6

    gen_val = float(nfes)/par.max_iterations

    nobj = len(fit_indexes)

    for n in range(0,nobj):
        fit_ix = fit_indexes[n]-1
        diff = abs(fit[fit_ix]-idealpoint[fit_ix])

        if(namda[fit_ix] == 0.0):
            feval = diff/eps
        else:
            feval = diff/namda[fit_ix]

        if (feval > max_fun):
            max_fun = feval
        if (feval < min_fun):
            min_fun = feval
        if (namda[fit_ix]>wmax):
            wmax = namda[fit_ix]
        if (namda[fit_ix] < wmin):
            wmin = namda[fit_ix]

    nd = norm_vector(namda)

    realA = [0.0]*nobj
    temp_namda = []

    for n in range(0,nobj):
        fit_ix = fit_indexes[n]-1
        realA[n] = (fit[fit_ix] - idealpoint[fit_ix])
        temp_namda.append(namda[fit_ix])

    A1 = norm_vector(realA)
    AB = abs(innerproduct(realA,temp_namda))
    d2 = pow(A1,2.0) - pow(AB,2.0)/pow(nd,2.0)
    d2 = pow(abs(d2),0.5)
    au = nobj*wmin

    penalty = 1.0-gen_val
    func_val = max_fun +10.0*au*penalty*d2
    return func_val


def Fitness_Function_MSF(fit, namda, fit_indexes, nfes):
    nobj = len(fit_indexes)
    max_fun = -1.0e+30
    min_fun = 1.0e+30
    gen_val = float(nfes)/par.max_iterations
    au = 1
    wmax = -1.0e+30
    wmin = 1.0e+30
    eps = 1.0e-6
    alpha = 1
    if (nobj > 2):
       eps = 1.0e-3
    for n in range(0,nobj):
       fit_ix = fit_indexes[n]-1
       diff = abs(fit[fit_ix])

       if(namda[fit_ix] == 0.0):
           feval = diff/eps
       else:
           feval = diff/namda[fit_ix]

       if (feval > max_fun):
           max_fun = feval
       if (feval < min_fun):
           min_fun = feval
       if (namda[fit_ix]>wmax):
           wmax = namda[fit_ix]
       if (namda[fit_ix] < wmin):
           wmin = namda[fit_ix]

    au = nobj *wmin
    penalty = 1.0 - gen_val
    func_val = pow(max_fun / min_fun, alpha*au*penalty)*max_fun

    return func_val
    



def norm_vector(x):
    sum = 0
    for i in range(0,len(x)):
        sum +=pow(x[i],2)
    return math.sqrt(sum)
def innerproduct(vec1, vec2):
    sum = 0
    for i in range(0,len(vec1)):
        sum += vec1[i]*vec2[i]

    return sum


def Clear_Ideal_and_Nadir(nobj):
    global nadirpoint
    global idealpoint
    nadirpoint = []
    idealpoint = []
    idealpoint = [1.0e+30] * nobj


def Create_Nadir(nobj):

    global nadirpoint
    nadirpoint = [-1.0e+30] * nobj


def Update_Nadirpoint(pop, flag):

    pops = len(pop)
    nobj = len(pop[0].fitness)

    global nadirpoint

    for n in range(0,nobj):
    
        if (flag == 1):
            nadirpoint[n] = -1.0e+30
        for i in range(0,pops):
        
            if (pop[i].fitness[n]>nadirpoint[n]):
            
                nadirpoint[n] = pop[i].fitness[n]
            
def Idealpoint_Update(ind):
    global idealpoint

    nobj = len(ind)

    for n in range(0,nobj):
    
        if (ind[n]<idealpoint[n]):        
            idealpoint[n] = ind[n];