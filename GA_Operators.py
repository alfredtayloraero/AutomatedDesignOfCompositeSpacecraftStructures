import parameters
import random
import copy

class Crossover_SBX(object):
    def __init__(self):
        self.name = "SBX"

    def Crossover(self,selected,fcode):
        cross_indi = []
        psize = len(selected)
        psize_odd = False
        vars = fcode.num_vars
        Di_c = parameters.Di_c

        if (psize % 2 == 1):
            psize_odd = True
            psize -=1

        for j in range (0,psize,2):
            parent1_code = selected[j].code
            parent2_code = selected[j+1].code

            if (random.random() <= parameters.cross_prob):
                code1 = []
                code2 = []

                for i in range(0,vars):
                    if (random.random() <= 0.5):
                        parent1 = parent1_code[i]
                        parent2 = parent2_code[i]
                        y2 = parent1
                        y1 = parent2
                        upper_b = fcode.bound_max[i]
                        lower_b = fcode.bound_min[i]

                        if(abs(parent1 - parent2)> pow(10,-14)):
                            if(parent2 > parent1):
                                y2 = parent2
                                y1 = parent1

                            rnd = random.random()

                            beta = 1.0 + 2.0 * ((y1 - lower_b) / (y2 - y1));
                            alpha = 2.0 - pow(beta, -(Di_c + 1.0));
                            if (rnd < 1./alpha):
                                betaq = pow((rnd*alpha), (1.0 / (Di_c + 1)));
                            else:
                                betaq = pow((1.0 / (2.0 - rnd*alpha)), (1.0 / (Di_c + 1)));
                            c1 = 0.5*((y1 + y2) - betaq*(y2 - y1));

                            beta = 1.0 + 2.0 * ((upper_b - y2) / (y2 - y1));
                            alpha = 2.0 - pow(beta, -(Di_c + 1.0));
                            if (rnd < 1./alpha):
                                betaq = pow((rnd*alpha), (1.0 / (Di_c + 1)));
                            else:
                                betaq = pow((1.0 / (2.0 - rnd*alpha)), (1.0 / (Di_c + 1)));
                            c2 = 0.5*((y1 + y2) + betaq*(y2 - y1));

                            if (c1 < lower_b):
                                c1 = lower_b;
                            elif (c1 > upper_b):
                                c1 = upper_b;
                            if (c2 < lower_b):
                                c2 = lower_b;
                            elif (c2 > upper_b):
                                c2 = upper_b;

                            if (random.random() <= 0.5):
                                code1.append(c2)
                                code2.append(c1)
                            else:
                                code1.append(c1)
                                code2.append(c2)
                        else:
                            code1.append(parent1)
                            code2.append(parent2)
                    else:
                        code1.append(parent1_code[i])
                        code2.append(parent2_code[i])
                indi1 = copy.deepcopy(selected[j])
                indi1.code = code1
                indi2 = copy.deepcopy(selected[j+1])
                indi2.code = code2
                cross_indi.append(indi1)
                cross_indi.append(indi2)
            else:
                indi1 = copy.deepcopy(selected[j])
                indi1.rank = 0
                indi1.crowd_dist = 0.
                indi2 = copy.deepcopy(selected[j+1])
                indi2.rank = 0
                indi2.crowd_dist = 0.
                cross_indi.append(indi1)
                cross_indi.append(indi2)
        if (psize_odd == True):
            cross_indi.append(indi1)

        if (len(cross_indi)!= len(selected)):
            raise SystemError("error")

        return cross_indi

class Roulette_Wheel(object):
    def __init__(self):
        self.name = "Roulette_Wheel_MLS"

    def Select(self,indiv_pop,psize,ncons,fit_indexes):
        indiv = []
        indiv_ix = []
        for i in range (0,len(indiv_pop)):
            if (indiv_pop[i].cons_violation == False):
                indiv.append(indiv_pop[i].fitness)
                indiv_ix.append(i)
        indiv_size = len(indiv)


        if (indiv_size == 0):
            for i in range (0,len(indiv_pop)):
                indiv.append(indiv_pop[i].fitness)
                indiv_ix.append(i)
                indiv_size = len(indiv)
        index = []
        invers_fitn_tot = 0.
        inverse_fitn =[]

        fit_size = len(fit_indexes)

        for i in range (0,indiv_size):
            temp = 0.;
            for j in range (0,fit_size):
                temp += indiv_pop[i].fitness[fit_indexes[j] - 1] / fit_size;

            inverse_fitn.append(1/temp)
            invers_fitn_tot+= 1/temp
        
        inverse_fitn[0] /= invers_fitn_tot

        for i in range (1,indiv_size):
            inverse_fitn[i] /= invers_fitn_tot
            inverse_fitn[i] += inverse_fitn[i-1]

        inverse_fitn[indiv_size - 1] = 1.

        for i in range (0, psize):
            index_h = 0
            wheel_r = random.random()

            for j in range(0,indiv_size):
                if (inverse_fitn[j] >= wheel_r):
                    index_h = j;
                    break

            index.append(copy.deepcopy(indiv_pop[indiv_ix[index_h]]))
        if (len(index) != psize):
            raise SystemError("error")

        return index


class Mutation_poly(object):
    def __init__(self):
        self.name = "Polynomial"

    def Mutate(self,code,fcode):
        
        Di_M = parameters.Di_m
        for i in range(0,len(code)):
            if (random.random() <= parameters.mut_prob):
                code_val = code[i]
                delta = 0.
                upper_b = fcode.bound_max[i]
                lower_b = fcode.bound_min[i]

                if (code_val > lower_b and code_val < upper_b):
                    if ((code_val - lower_b) < (upper_b - code_val)):
                        delta = (code_val - lower_b) / (upper_b - lower_b);
                    else:
                        delta = (upper_b - code_val) / (upper_b - lower_b);
                    
                    indi = 1. / (Di_M +1.)
                    rnd1 = random.random()
                    deltaq = 0.

                    if (rnd1 <= 0.5):
                        val = 2 * rnd1 + (1.0 - 2.0 * rnd1)*pow((1.0 - delta), (Di_M + 1.));
                        deltaq = pow(val, indi) - 1.0;
                    else:
                        val = 2 * (1.0 - rnd1) + 2.0 * (rnd1 - 0.5)*pow((1.0 - delta), (Di_M + 1.));
                        deltaq = 1.0 - pow(val, indi);

                    code[i] +=(deltaq*(upper_b - lower_b));

                    if (code[i] < lower_b):
                        code[i] = lower_b
                    elif (code[i] > upper_b):
                        code[i] = upper_b;

                else:
                    code[i] = lower_b + random.random()*(upper_b-lower_b)
        return


