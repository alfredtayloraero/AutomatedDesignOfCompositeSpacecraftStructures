import Fit_Function
import random
import copy
import parameters
import Add_func
import GA_Operators
import MOEAD

class Individual(object):
    """description of class"""
    def __init__(self, fcode):
        self.code= [];
        self.fitness = [];
        self.cons_val = [];
        self.cons_violation = False;

        self.crowd_dist = 0.;
        self.rank = 0;

        self.namda = [];
        self.table = [];
        self.saved_fitness = [];
        self.utility = 0.;

        for i in range (0, fcode.num_vars):
            x = random.random();
            val = fcode.bound_min[i] + x*(fcode.bound_max[i] - fcode.bound_min[i]);
            self.code.append(val)

    def Fitness_Calc(self, fcode):
        self.fitness = fcode.Fitness_C(self.code)
            
        

    def Cons_Viol_Calc(self,fcode):
        self.cons_val = fcode.Cons_Calc(self.code,self.fitness)
        self.cons_violation = False

        for i in range (0, fcode.num_cons):
            if (self.cons_val[i] < parameters.cons_check_param):
                cons_violation = True;  
        return

class population(object):
    
    def __init__(self, fcode, mcode, scode, ccode, amt):
        self.indiv = [] 
        self.f_code = fcode
        self.s_code = scode
        self.m_code = mcode
        self.c_code = ccode
        self.fit_index = [[6,3],[3,6]] #this is buggy and only works with uf1 and uf2
        self.size = amt
        self.sort_vect = []

        for i in range (0,amt):
            self.indiv.append(Individual(self.f_code))

    def Fitness_Calc(self):
        for i in range (0, self.size):
            self.indiv[i].Fitness_Calc(self.f_code);

    def Fitness_Calc_Indi(self,indi):
        indi.Fitness_Calc(self.f_code)

    def Selection(self):
        temp_fit_index = self.fit_index[0]
        return self.s_code.Select(self.indiv, self.size, self.f_code.num_cons, temp_fit_index)

    def Crossover(self):
        selected = self.Selection()

        new_pop = self.c_code.Crossover(selected,self.f_code)

        if (len(new_pop) != len(self.indiv)):
            raise SystemError("error")
        return new_pop

    def Crossover_NSGAII(self, selected):
        new_pop = self.c_code.Crossover(selected,self.f_code)


        return new_pop

    def Mutation(self):
        for i in range (0, self.size):
            self.m_code.Mutate(self.indiv[i].code, self.f_code);

    def Mutation_Indi(self,indi):
        self.m_code.Mutate(indi.code, self.f_code)

    def Mutation_NSGAII(self,pop):
        for i in range (0, len(pop)):
            self.m_code.Mutate(pop[i].code, self.f_code)
            pop[i].Fitness_Calc(self.f_code)

    def Sort_Individuals(self):
        self.sort_vect = self.fit_index[0];
        self.indiv.sort(key = self.sort2)

    def Sort_Individuals2(self):
        self.sort_vect = self.fit_index[1];
        self.indiv.sort(key = self.sort2)

    def Add(self,indi):
        self.indiv.append(copy.deepcopy(indi))
        self.size +=1
        
    def Add_Pop(self,indi):
        for i in range (0, len(indi)):
            self.Add(indi[i])

    def Remove(self,ix):
        del self.indiv[ix]
        self.size -=1

    #def Remove_Last(self):
        #del self.indiv[self.size - 1]
        #self.size -=1

    def sort2(self,x):
        fit_size = len(self.sort_vect)

        val = 0;

        for i in range (0,fit_size):
            val += x.fitness[self.sort_vect[i]-1]
        return val;


class collective(population):
    def __init__(self, pop, label, ix, m):
        population.__init__(self,pop.f_code, pop.m_code, pop.s_code, pop.c_code, 0)
        self.index = ix;
        self.index_vid = ix
        for i in range (0, pop.size):
            if (label[i] == ix):
                self.indiv.append(copy.deepcopy(pop.indiv[i]))
                self.size +=1
        self.elite = []
        self.was_erased = True
        self.mode = m
        self.fitness = 0

    def Elite_Create(self):
        esize = int(self.size*parameters.MLSGA_elite_size)

        if (esize == 0):
            esize = 1;

        self.Sort_Individuals();
        self.elite = [];
        for i in range (0,esize):
            self.elite.append(copy.deepcopy(self.indiv[i]))

    def Elite_Replace(self):
        for i in range(0,len(self.elite)):
            self.indiv.append(copy.deepcopy(self.elite[i]))

        self.Sort_Individuals()
        i_size = len(self.indiv)
        while (i_size > self.size):
            del self.indiv[i_size-1]
            i_size -= 1

        if (len(self.indiv) != self.size):
            raise SystemError("error")

        self.Fitness_Calc_Col()

    def Fitness_Calc_Col(self):
        
        pop_size = self.size
        self.fit_index=[[1,2,3,4,5,6],[1,2,3,4,5,6]] #big bodge
        temp_col_fit = 0.
        for i_ind in range(0,pop_size):
            fit_size = len(self.fit_index[1])
            fit_temp = 0.
            for j in range(0,fit_size):
                fit_temp += self.indiv[i_ind].fitness[self.fit_index[1][j]-1]/fit_size
            temp_col_fit += fit_temp;
        temp_col_fit /= pop_size
        self.fitness = temp_col_fit

    def Erase(self):
        self.indiv = []
        self.size = 0
        self.fitness = 0
        self.elite = []
        self.was_erased = True

class pareto_front(object):
    def __init__(self):
        self.size = 0
        self.indiv = []

    def Pareto_Search(self,indiv):      #take a list of individuals to be checked
        if (len(indiv) == 0):
            return;

        for i_ind in range (0, len(indiv)):
            self.indiv.append(copy.deepcopy(indiv[i_ind]))
        
        self.size = len(self.indiv)
        n_func = len(self.indiv[0].fitness)
        
        ncons = len(self.indiv[0].cons_val)


        if (ncons > 0):
            for i in range (0,par_size):
                if (self.indiv[i].cons_violation == True):
                    del self.indiv[i]
                    self.size -= 1
        i = 0
        while (i < self.size):
            j = i+1
            while (j < self.size):
                m_i = 0;
                m_j = 0;
                for k in range (0,n_func):
                    val = self.indiv[i].fitness[k] - self.indiv[j].fitness[k]
                    if (val < pow(10, -parameters.PF_res )):
                        m_i +=1
                    else:
                        m_j +=1
                if (m_i == n_func):
                    del self.indiv[j]
                    self.size -= 1
                    j -= 1
                elif (m_j == n_func):
                    del self.indiv[i]
                    self.size -=1
                    i -= 1
                    j = self.size
                j +=1
            i+=1
        if (self.size > parameters.PF_refine_size + 200):
            self.Pareto_Refine(parameters.PF_refine_size)
        

    def Pareto_Search_Indi(self,indi):
        if (indi.cons_violation == True):
            return -1;
        
        n_func = len(self.indiv[0].fitness)


        i = 0
        while (i < self.size):
            m_i = 0;
            m_j = 0;
            for k in range (0,n_func):
                val = self.indiv[i].fitness[k] -indi.fitness[k]
                if (val < pow(10, -parameters.PF_res )):
                    m_i +=1
                else:
                    m_j +=1
            if (m_i == n_func):
                return -1;
            elif (m_j == n_func):
                del self.indiv[i]
                self.size -=1
                i -= 1
            i+=1
        self.indiv.append(copy.deepcopy(indi))
        self.size +=1;

        if (self.size > parameters.MTS_PF_refine_size + 200):
           self.Pareto_Refine(parameters.MTS_PF_refine_size)

        return self.size - 1;

    def Pareto_Refine(self,size=parameters.PF_size):
        if (self.size < size):
            return

        n_obj = len(self.indiv[0].fitness)

        if (n_obj == 2):
            self.indiv.sort(key = lambda x1: x1.fitness)

            while (self.size > size):
                temp_min_storage = [parameters.INF,0]

                for i in range (1, self.size - 1):
                    temp_min = Add_func.Distance(self.indiv[i].fitness,self.indiv[i-1].fitness);
                    temp_min2 = Add_func.Distance(self.indiv[i].fitness,self.indiv[i+1].fitness);

                    if (temp_min == 0 or temp_min2 == 0):
                        raise SystemError("error")

                    if (temp_min2 > temp_min):
                        temp_min = temp_min2
                    if (temp_min < temp_min_storage[0]):
                        temp_min_storage[0] = temp_min
                        temp_min_storage[1] = i
                del self.indiv[temp_min_storage[1]]
                self.size -=1;
        return