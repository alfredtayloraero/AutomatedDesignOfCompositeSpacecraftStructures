import math
import parameters
import copy

def Sort_Namda(x):
    return x.namda[0]

def Sort_MTS(x):
    return x.utility

def Distance(point1, point2):
    dist = 0;
    #dist += pow(float(point1[i])-float(point2[i]),2)
    if (len(point1) != len(point2)):
        raise SystemError("error")
        #systemError

    for i in range (0,len(point1)):
        dist += pow(float(point1[i])-float(point2[i]),2)
        #list index out of range
    return math.sqrt(dist);


def MLS_Col_Fit_Create(col,MLS):
    nobj = col[0].f_code.num_objs;
    ncol = len(col);


    if (MLS == 0):
        temp_fit = [];

        for i in range (0,nobj):
            temp_fit.append(i+1);

        temp_fit_index = [temp_fit,temp_fit];

        col[0].fit_index = temp_fit_index;
        return;

    fit_index_sel = parameters.fit_index_sel;
    fit_index_col_sel = parameters.fit_index_col_sel;

    if (nobj == 2):
        if (MLS == 7):
            temp_fit_index1 = [ [ fit_index_sel, fit_index_col_sel ],[ fit_index_sel, fit_index_col_sel ] ];
            temp_fit_index2 = [ [ fit_index_sel, fit_index_col_sel ],[fit_index_col_sel ] ];
            temp_fit_index2_n =[ [ fit_index_sel],[ fit_index_col_sel ] ];
            temp_fit_index3 = [ [ fit_index_sel, fit_index_col_sel ],[ fit_index_sel ] ];
            temp_fit_index3_n = [ [fit_index_col_sel ],[ fit_index_sel ] ];

            for i_col in range(0,ncol):
                col_index = col[i_col].index;

                if (col_index % 3 == 1):
                    col[i_col].fit_index = temp_fit_index1
                elif (col_index % 3 == 0):
                    if (col[i_col].mode == "Normal"):
                        col[i_col].fit_index = temp_fit_index2_n;
                    else:
                        col[i_col].fit_index = temp_fit_index2;
                else:
                    if (col[i_col].mode == "Normal"):
                        col[i_col].fit_index = temp_fit_index3_n;
                    else:
                        col[i_col].fit_index = temp_fit_index3;
        else:
            raise SystemError("not implemented")
    elif (nobj < 7):
        if (MLS == 7):
            temp_fit = [];
            for i_o in range(0,nobj):
                temp_fit.append(i_o+1);
            temp_fit_index1 = [ temp_fit,temp_fit ];

            for i_col in range(0,ncol):
                col_index = col[i_col].index - 1;

                switch = col_index % (nobj+1);


                if (switch == 0):
                    col[i_col].fit_index = temp_fit_index1
                else:
                    col[i_col].fit_index = [temp_fit,switch];



        else:
            raise SystemError("not implemented")
    else:
        raise SystemError("not implemented")


def Output_First_Lines(index_file,time_file,GA_data_file):
    #write index file
    index_file.write("Modes,");
    index_file.write("Termination,");
    index_file.write("T value,");
    index_file.write("Selection,");
    index_file.write("Crossover,");
    index_file.write("Crossover rate,");
    index_file.write("Mutation,");
    index_file.write("Mutation rate,");
    index_file.write("Function,");
    index_file.write("Population,");
    if (parameters.MLSGA_Hybrid == True):
    
        index_file.write("Collectives,");
        index_file.write("Fit for indi,");
        index_file.write("Fit for col,");
        index_file.write("Elitism rate (%),");
        index_file.write("MLS,");
        index_file.write("Col elim limit,");
    
    index_file.write("Number of constrains,");
    index_file.write("PF_Refine\n");

    #write second line
    index_file.write(str(parameters.GA_mode)  + ",");
    index_file.write("Iterations,");
    index_file.write(str(parameters.max_iterations) + ",");
    index_file.write(str(parameters.Selection.name) + ",");
    index_file.write(str(parameters.Crossover.name) + ",");
    index_file.write(str(parameters.cross_prob)  + ",");
    index_file.write(str(parameters.Mutation.name) + ",");
    index_file.write(str(parameters.mut_prob)  + ",");
    index_file.write(str(parameters.Function.name_func) + ",");
    index_file.write(str(parameters.Pop_size)  + ",");
    if (parameters.MLSGA_Hybrid == True):
    
        index_file.write(str(parameters.MLSGA_n_col)  + ",");
        index_file.write(str(parameters.fit_index_sel)  + ",");
        index_file.write(str(parameters.fit_index_col_sel)  + ",");
        index_file.write(str(parameters.MLSGA_elite_size)  + ",");
        index_file.write(str(parameters.MLS)  + ",");
        index_file.write(str(parameters.MLSGA_n_col_elim)  + ",");
    
    index_file.write(str(parameters.Function.num_cons) + ",");
    index_file.write(str(parameters.PF_REFINE) + "\n");

    #write time file
    time_file.write("Time,");
    if (parameters.IGD_on):
        time_file.write("IGD,");
    if (parameters.HV_on):
        time_file.write("HV,");
    time_file.write("End_Generation,");
    time_file.write("End_Iteration\n");

    #write data file
    GA_data_file.write("Min time (run),");
    GA_data_file.write("Max time (run),");
    GA_data_file.write("Average time (run),");
    GA_data_file.write("Std deviation (time run),");
    if (parameters.IGD_on):
        GA_data_file.write("Min IGD,");
        GA_data_file.write("Max IGD,");
        GA_data_file.write("Average IGD,");
        GA_data_file.write("Std deviation (IGD),");
            
    if (parameters.HV_on):
        GA_data_file.write("Min HV,");
        GA_data_file.write("Max HV,");
        GA_data_file.write("Average HV,");
        GA_data_file.write("Std deviation (HV),");
        
    GA_data_file.write("Min generation,");
    GA_data_file.write("Max generation,");
    GA_data_file.write("Average generation,");
    GA_data_file.write("Std deviation (generations)\n");


def Output_Add_Time(time_file, time, IGD, HV, nfes, index_r,end_generation):
    time_file.write(str(time) + ",");
    if (parameters.IGD_on):
        time_file.write(str(IGD) + ",");
    if (parameters.HV_on):
        time_file.write(str(HV) + ",");
    time_file.write(str(end_generation) + ",");
    time_file.write(str(nfes) + "\n");

def Output_Add_Summary(GA_data_file, run_data):

    GA_data_file.write(str(run_data.time_struct[0])+ ",");
    GA_data_file.write(str(run_data.time_struct[1])+ ",");
    GA_data_file.write(str(run_data.time_struct[2])+ ",");
    GA_data_file.write(str(run_data.time_struct[3])+ ",");
    if (parameters.IGD_on):
        GA_data_file.write(str(run_data.IGD_struct[0])+ ",");
        GA_data_file.write(str(run_data.IGD_struct[1])+ ",");
        GA_data_file.write(str(run_data.IGD_struct[2])+ ",");
        GA_data_file.write(str(run_data.IGD_struct[3])+ ",");
            
    if (parameters.HV_on):
        GA_data_file.write(str(run_data.HV_struct[0])+ ",");
        GA_data_file.write(str(run_data.HV_struct[1])+ ",");
        GA_data_file.write(str(run_data.HV_struct[2])+ ",");
        GA_data_file.write(str(run_data.HV_struct[3])+ ",");
        
    GA_data_file.write(str(run_data.generation_struct[0])+ ",");
    GA_data_file.write(str(run_data.generation_struct[1])+ ",");
    GA_data_file.write(str(run_data.generation_struct[2])+ ",");
    GA_data_file.write(str(run_data.generation_struct[3])+ "\n");






def Dominance_check(ind1, ind2, fit_indexes):
	flag1 = flag2 = 0;
	nobj = len(ind1.fitness)

	cons_size = len(ind1.cons_val)
	if (cons_size > 0):
	
		constr_val_ind1 = 0
		constr_val_ind2 = 0
		for i in range(0,cons_size):
		
			if (ind1.cons_val[i] < 0):
				constr_val_ind1 += ind1.cons_val[i]
			if (ind2.cons_val[i] < 0):
				constr_val_ind2 += ind2.cons_val[i]
		
		if (constr_val_ind1 < 0 and constr_val_ind2 < 0):
		
			if (constr_val_ind1 > constr_val_ind2):
				return 1
			elif (constr_val_ind1 < constr_val_ind2):
				return -1
			else:
				return 0
		
		elif (constr_val_ind1 < 0. and constr_val_ind2 == 0.):
			return -1;
		elif (constr_val_ind1 == 0. and constr_val_ind2 < 0.):
			return 1;
	


	for i in range(0,len(fit_indexes)):
	
		fit_ind1 = 0
		fit_ind2 = 0

		fit_ind_ix = fit_indexes[i]-1;

		
		fit_ind1 = ind1.fitness[fit_ind_ix]
		fit_ind2 = ind2.fitness[fit_ind_ix]
		
												
		if (fit_ind1 < fit_ind2):
			flag1 = 1
		elif (fit_ind1 > fit_ind2):
			flag2 = 1
	
	if (flag1 == 1 and flag2 == 0):
		return (1);
	elif (flag1 == 0 and flag2 == 1):
		return (-1)
	else:
		return (0);



def Dominance_check_NCons(ind1, ind2, fit_indexes):
	flag1 = flag2 = 0;
	nobj = len(ind1.fitness)


	for i in range(0,len(fit_indexes)):
	
		fit_ind1 = 0
		fit_ind2 = 0

		fit_ind_ix = fit_indexes[i]-1;

		
		fit_ind1 = ind1.fitness[fit_ind_ix]
		fit_ind2 = ind2.fitness[fit_ind_ix]
		
												
		if (fit_ind1 < fit_ind2):
			flag1 = 1
		elif (fit_ind1 > fit_ind2):
			flag2 = 1
	
	if (flag1 == 1 and flag2 == 0):
		return (1);
	elif (flag1 == 0 and flag2 == 1):
		return (-1)
	else:
		return (0);



def Uniform_Weights_Generate(pop_size, nobj, type):  #moead type 3
    val_w = recursion_point(pop_size+1,nobj) #pop_size+1


    if(type == 1):
        for i in range (0, len(val_w)):
        
            temp_val = 0
            for j in range (0, len(val_w)):
                temp_val += pow(val_w[i][j], 2);
            temp_val = (temp_val)*0.5;

            for j in range (0, len(val_w)):
                val_w[i][j] /= temp_val;
        
    elif(type == 2):
    
        for i in range (0, len(val_w)):
            temp_val = 0;
            for j in range (0, len(val_w)):
                temp_val += (val_w[i][j])**0.5;
            temp_val = temp_val/(pow(temp_val,3));

            for j in range (0, len(val_w)):
                val_w[i][j] *= temp_val;
      
    return val_w


def recursion_point(pop_size, nobj):
    temp_nobj = nobj - 1
    layer = 0
    wnum = 0

    while (wnum < pop_size):
    
        layer+=1;
        wnum = Get_Weight_Num(layer, temp_nobj, 0)
    
    layer-=1

    wnum = Get_Weight_Num(layer, temp_nobj, 0)
    out = []
    for i in range(0, wnum):
        out.append([0.0]*nobj)

    step = 1.0 / (layer - 1.0)

    Get_Weight(layer, temp_nobj, step, out)

    for i in range(0,len(out)):
    
        temp_sum = 0
        for j in range(0,temp_nobj):
            temp_sum += out[i][j]
        if (temp_sum >= 1):
            temp_sum = 1

        out[i][temp_nobj] = 1. - temp_sum
    
    return out;

def Get_Weight_Num(layer, od, wnum):

    out = wnum;
    if (od == 1):
        out = wnum + layer
    else:
        for i in range(1,layer+1):        
            out = Get_Weight_Num(i, od - 1, out)
    return out


def Get_Weight(layer, od, step, out):

    if (len(out) == 0):
        return

    for i in range(1,layer+1):     
    
        n1 = Get_Weight_Num(i - 1, od, 0)
        n2 = Get_Weight_Num(i, od, 0)

        for j in range(n1,n2):        
            out[j][od - 1] = step*(layer - i)
        
        if (od > 1):
            temp = copy.deepcopy(out[n1:n2])
            Get_Weight(i, od - 1, step, temp)

            k = 0
            for j in range(n1,n2):      
                out[j] = temp[k]
                k+=1