import Add_func


def IGD_calc(Pareto_F,real_PF):
    Pareto = Pareto_F.indiv;

    size = len(Pareto);

    if (size == 0):
        return 1e10;

    IGD_val = 0;

    di_min_GA = [];

    for i in range(0,len(real_PF)):
        for j in range(0,size):
            temp_min_di = Add_func.Distance(Pareto[j].fitness,real_PF[i]);

            if (j == 0):
                di_min_GA.append(temp_min_di);
            else:
                if (di_min_GA[i] > temp_min_di):
                    di_min_GA[i] = temp_min_di;

    IGD_val = sum(di_min_GA)/len(real_PF);

    return IGD_val;