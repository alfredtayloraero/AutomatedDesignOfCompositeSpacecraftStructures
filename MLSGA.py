import time
import os
import copy


import parameters
import Class
import GA_data
import Add_func
import IGD
import HV
import kmeans
import NSGAII
import MOEAD
import HEIA



Debug = True
#**MSLGA***/
#**Original idea by Adam Sobey***/
#**Written by Przemyslaw Grudniewski***/
#**Updates and additional ideas by Przemyslaw Grudniewski***/


                                                                          #Elitism
##std::ofstream file
##std::ofstream graph, *graph_v								#File for saving graph data
##std::ofstream Debug										#File for saving debugging data


##std::vector<std::ofstream> graph_v_vector;
##std::ofstream PF, PF2, IGD_stream, HV_stream								#Output files

#if (Add_func.Error_Check() == False):
 #   raise SystemError("Wrong parameters setting")
    
#get the comment to files and folders names
comment	= input("Comment (to the folder and results names): ");								#The comment to the folder and results name

    
#Get the starting date and time of the program
date = time.localtime(time.time());
comment = str(date[0]) + "_" + str(date[1])+ "_" + str(date[2])+ "_"+str(date[3])+ "_"+ str(date[4]) + "_" + comment;

print(comment);
comment2 = "Results/" + comment;
os.mkdir(comment2);

## FILE * plotPipe = _popen("c:\\gnuplot\\bin\\gnuplot.exe", "w")								#Pipe for GNUplot

if (parameters.MLSGA_Hybrid == False):
    parameters.MLSGA_n_col = 1;
    parameters.MLSGA_n_col_elim = 0;
    parameters.MLS = 0;


#intialise run independent parameters
coevol_amt = len(parameters.GA_mode);
n_col = parameters.MLSGA_n_col
n_col_elim = parameters.MLSGA_n_col_elim;
pop_size = parameters.Pop_size;
skip_gen = parameters.MLSGA_col_elim_limit
Func = parameters.Function												
Crossover = parameters.Crossover 
Mutation = parameters.Mutation
Selection = parameters.Selection
MODE_vector = parameters.GA_mode
nobj = Func.num_objs
nvar = Func.num_vars

termination_crit = parameters.max_iterations

run_data = GA_data.GA_data(parameters.IGD_on,parameters.HV_on);
IGD_storage = [];
HV_storage = [];
real_PF = Func.Plot_PF(parameters.PF_real_size)


#create output files (result file)
results_index = open(comment2 + "/" + comment + "_results_index.csv","w");
results_time= open(comment2 + "/" + comment + "_results_time.csv","w");
results_GA_Data = open(comment2 + "/" + comment + "_results_GA_Data.csv","w");

#create the first two lines:
Add_func.Output_First_Lines(results_index,results_time,results_GA_Data)

for i_run in range(1, parameters.n_runs + 1):

    print ("Run: #",i_run," started\n")

    MOEAD.Clear_Ideal_and_Nadir(Func.num_objs)


    run_t = time.perf_counter();
    #initialise parameters
    run_time = 0
    nfes = 0; #number of function evaluations for MOEAD and MTS
    cons_viol_count = 0;	#How many times constrains have been violated in the current run
    
    flag_stored = False;
    
    IGD_val_final = 0;												#Final IGD value
    HV_val_final = 0;										#Final HV value
    min_fitness = 0.;											#minimum fitness for current run - only for MLS1
    end_generation = 0;										#generation at which solution was found
    end_fitness = 0;
    GA_success = 0;

    #Open the file for saving PF results
    PF_data = open(comment2 + "/#" + str(i_run) + "_PF.csv","w");
    MOEAD.Clear_Ideal_and_Nadir(nobj);
    ##Open the file for saving temprary graph data

    #Actual GA intialisation
    coevol_subpop_size = int(pop_size / coevol_amt);

    
    #Initialise the Pareto Front
    Pareto_F = Class.pareto_front()

    #initialise the vector of collectives
    col = []

    #initialise the vector of temporary collectives
    coevol_col_vect = []

    for i_coevol in range (0, coevol_amt):
        MODE = MODE_vector[i_coevol];
        col_temp_size = coevol_subpop_size

        if (i_coevol + 1 == coevol_amt):            #I don't understand why this is here
            col_temp_size =  pop_size #int(pop_size - ((coevol_amt - 1) * coevol_subpop_size)); #This is just equal to 1 as coevol amt is equal to 1

        p1 = Class.population(Func,Mutation,Selection,Crossover,col_temp_size)
        #Initialise the population
        if (MODE == "MOEAD" or MODE == "MOEADMSF" or MODE == "MOEADPSF"):
           #initialise the MOEAD
            nfes += MOEAD.MOEAD_Init(nobj,p1);
        elif (MODE == "MTS"):
            MTS.MTS_Init(p1,Pareto_F,n_col);        
        elif (MODE == "HEIA"):  
            nfes += HEIA.HEIA_Init(n_col,p1);
        elif (MODE == "NSGAII" or MODE == "Normal"):
            True
        else:
            raise SystemError("error")
        
        
        #Initialise the vector for storaging real labels
        real_label = []
        #Do SVM
        if (n_col > 1):
            real_label = kmeans.Clustering(p1,n_col);
        elif (n_col == 1):
            real_label = [1]*pop_size
        else:
            raise SystemError("error")

        print("SVM Finished")

        #**** Collevtive generation *****/

        temp_col = []

        biggest_col = [0,0,0,0];



        for iCol in range (0, n_col):
            #Assign the individuals to the collective
            temp_col.append(Class.collective(p1,real_label,iCol + 1,MODE))

            if (n_col != 1):
                if (temp_col[iCol].size > biggest_col[1]):
                    biggest_col[2] = biggest_col[0];
                    biggest_col[3] = biggest_col[1];
                    biggest_col[0] = iCol;
                    biggest_col[1] = temp_col[iCol].size;
        #Rearrange collectives to remove too small ones
        for iCol in range (0, n_col):
            temp_i = temp_col[iCol].size;
            if (n_col != 1):
                if (temp_i <= (col_temp_size / parameters.min_col_size_multi) or temp_i <= 2):
                    if (biggest_col[1] >= biggest_col[3]):
                        for i in range (temp_i, int(col_temp_size /  parameters.new_col_size_multi)):
                            #Get the individual to be cut
                            temp_indi = temp_col[biggest_col[0]].indiv[0]
                            #Add the individual to the new collective
                            temp_col[iCol].Add(temp_indi);
                            #Remove it from the old one
                            temp_col[biggest_col[0]].Remove(0)
                            #Correct the size of the biggest collective
                            biggest_col[1] -=1
                    else:
                        for i in range (temp_i, int(col_temp_size /  parameters.new_col_size_multi)):
                            #Get the individual to be cut
                            temp_indi = temp_col[biggest_col[2]].indiv[0]
                            #Add the individual to the new collective
                            temp_col[iCol].Add(temp_indi);
                            #Remove it from the old one
                            temp_col[biggest_col[2]].Remove(0)
                            #Correct the size of the biggest collective
                            biggest_col[3] -=1
        
        n_col_target = int(n_col / coevol_amt);
        if (n_col % coevol_amt != 0):
            if (n_col % coevol_amt >= i_coevol + 1):
                n_col_target +=1;
        while (len(temp_col) > n_col_target):
            ix_small = 0;
            ix_small_size = temp_col[0].size + temp_col[1].size;
            #find the smalles pair- they have to be near each other in order to maintain
            for iCol in range (1,len(temp_col)-1):
                if ( temp_col[iCol].size < temp_col[iCol + 1].size):
                    ix_small = iCol

            temp_col[ix_small].Add_Pop(temp_col[ix_small + 1].indiv)

            #remove the added collective
            del temp_col[ix_small + 1];
        #check the size of population  
        temp_size = 0
        for iCol in range(0,len(temp_col)):
            temp_col[iCol].index = iCol+1
            temp_size += temp_col[iCol].size

        if (temp_size != col_temp_size):
            raise SystemError("error")

        coevol_col_vect.append(temp_col)

    if (coevol_amt == 1):
        col = coevol_col_vect[0]

    else:
        ix_col = 0
        for i_coevol in range (0, coevol_amt):
            for i_col in range (0, len(coevol_col_vect[i_coevol])):
                col.append(coevol_col_vect[i_coevol][i_col])
                col[ix_col].index_vid = ix_col + 1;	
                ix_col +=1
    
    if (n_col != len(col)):
        raise SystemError("error")
    Add_func.MLS_Col_Fit_Create(col,parameters.MLS)

    for iCol in range (0, n_col):
        MODE = col[iCol].mode
        if (coevol_amt != 1):
            col[iCol].index = col[iCol].index_vid;

        if (MODE == "Normal" or MODE == "NSGAII"):
            col[iCol].Fitness_Calc()
            nfes += col[iCol].size
        elif (MODE == "MOEAD" or MODE == "MOEADMSF" or MODE == "MOEADPSF"):
            True
        elif (MODE == "MTS"):
            MTS.MTS_Init_Col(col[iCol]);
        elif (MODE == "HEIA"):
            HEIA.HEIA_Pop_Init(col[iCol]);

        if (MODE == "MOEADPSF" or MODE == "MOEADMSF" or parameters.MLSGA_norm_obj == True):
            if (iCol == 0):
                MOEAD.Create_Nadir(nobj);
                MOEAD.Update_Nadirpoint(col[iCol].indiv, 1);
            elif (col[iCol - 1].mode != "MOEADPSF" and col[iCol-1].mode != "MOEADMSF" and parameters.MLSGA_norm_obj != True):
                MOEAD.Create_Nadir(nobj);
                MOEAD.Update_Nadirpoint(col[iCol].indiv, 1);
            else:
                MOEAD.Update_Nadirpoint(col[iCol].indiv, 2);

        Pareto_F.Pareto_Search(col[iCol].indiv)
        #end of iCol loop

    if (n_col != len(col)):
        raise SystemError("error")
    elif (Debug == True):
        col_size = 0;
        for iCol in range(0,n_col):
            col_size += col[iCol].size
        if (col_size != pop_size):
            raise SystemError("error")

    ##reinit

    print("Collective generation - complete!")

    PF_storage = []

    termination_index = 0
    iGen = 1
    while (nfes < termination_crit):

        if(Debug == True):
            termination_index = nfes;
            print (termination_index)
        iGen+=1

        for iCol in range (0,n_col):
            MODE = col[iCol].mode


            if (MODE == "Normal"):
                #Create the elite population
                col[iCol].Elite_Create()
                #Create the offspring collective via selection and crossover
                col[iCol].indiv = col[iCol].Crossover()
                #Mutate the offspring collective
                col[iCol].Mutation()
                #Calculate the fitness of the offspring population
                col[iCol].Fitness_Calc()
                nfes += col[iCol].size
                #find pareto front
                Pareto_F.Pareto_Search(col[iCol].indiv)
                #Insert the elite population
                col[iCol].Elite_Replace()
            elif (MODE == "NSGAII"):
                #Calculate the offspring collective
                nfes += NSGAII.NSGAII_Calc(col[iCol]);
                #find pareto front
                Pareto_F.Pareto_Search(col[iCol].indiv)

            elif (MODE == "MOEAD" or MODE == "MOEADMSF" or MODE == "MOEADPSF"):
                if (n_col == 1):
                    PF_storage = []
                temp_PF_storage = []
                nfes +=MOEAD.MOEAD_Calc(col[iCol], temp_PF_storage, nfes)

                for i_PF in range (0, len(temp_PF_storage)):
                    PF_storage.append(copy.deepcopy(temp_PF_storage[i_PF]))

                if (iGen -1 % 50 == 0):
                    MOEAD.Utility_Comp(col[iCol], iGen,nfes)
            elif (MODE == "MTS"):
                nfes += MTS.MTS_Calc(col[iCol],nfes)
            elif (MODE == "HEIA"):
                if (n_col == 1):
                    PF_storage = []
                temp_PF_storage = []
                nfes += HEIA.HEIA_Calc(col[iCol],temp_PF_storage, nfes)
                
                for i_PF in range (0, len(temp_PF_storage)):
                    PF_storage.append(copy.deepcopy(temp_PF_storage[i_PF]))
            else:
                raise SystemError("error")


            if (n_col == 1 and MODE != "Normal"and MODE != "NSGAII"):
                PF_storage = []
            elif (len(PF_storage) >= 200 or nfes >= termination_crit ):
                Pareto_F.Pareto_Search(PF_storage)
                PF_storage = []
            
            if (nfes>100000 and not flag_stored):
                PF_data.write("100K Iterations");
                PF_data.write("\n");
                for iPF in range (0,Pareto_F.size):
                    #Save objectives
                    for iObj in range (0,nobj):
                        PF_data.write(str(Pareto_F.indiv[iPF].fitness[iObj]));
                        PF_data.write(",");
            
                    for iVar in range(0,nvar):
                        PF_data.write(",");
                        PF_data.write(str(Pareto_F.indiv[iPF].code[iVar]));
                    PF_data.write("\n");
                flag_stored=True
            
            if (nfes >= termination_crit):  #Breaks the loop
                break

            if (MODE == "MOEADPSF" or MODE == "MOEADMSF" or parameters.MLSGA_norm_obj == True):
                if (iCol == 0):
                    #Add_func.Create_Nadir(nobj);
                    MOEAD.Update_Nadirpoint(col[iCol].indiv, 2);
                elif (col[iCol - 1].mode != "MOEADPSF" and col[iCol-1].mode != "MOEADMSF" and parameters.MLSGA_norm_obj != True):
                    #Add_func.Create_Nadir(nobj);
                    MOEAD.Update_Nadirpoint(col[iCol].indiv, 2);
                else:
                    MOEAD.Update_Nadirpoint(col[iCol].indiv, 2);

            #end of loop for collectives evolutionn
        
        if ((iGen-1) % skip_gen == 0):

            col.sort(key = lambda x: x.fitness)
            for iCol in range (0, n_col):
                col[iCol].Sort_Individuals2()

            iIndivSource = -1
            iIndivSourceMax = int(pop_size / (n_col * 3))
            for iElim in range (0,parameters.MLSGA_n_col_elim):
                if (abs(col[n_col - 1].fitness - col[n_col - n_col_elim - 1].fitness) < pow(10, -parameters.pow_eq_zero)):
                    break

                col_elim_index = n_col - 1 - iElim
                col_size_temp = col[col_elim_index].size

                MODE = col[col_elim_index].mode
                temp_namda = []
                if (MODE == "MOEAD" or MODE == "MOEADMSF" or MODE == "MOEADPSF"):
                    for i in range (0, col_size_temp):
                        temp_namda.append(col[col_elim_index].indiv[i].namda)

                #Erase the worst collective
                col[col_elim_index].Erase();
                iIndivNew = 0
                k = 0;
                temp_const_limit = 0;
                while (iIndivNew < col_size_temp):
                    if (k == 0):
                        iIndivSource+=1
                    if (iIndivSource >= iIndivSourceMax):
                        iIndivSource = 0
                    if (iIndivSource < col[k].size):
                        if (col[k].indiv[iIndivSource].cons_violation == False or temp_const_limit > pop_size / 2.):
                            col[col_elim_index].Add(col[k].indiv[iIndivSource])
                            if (MODE == "MOEAD" or MODE == "MOEADMSF" or MODE == "MOEADPSF" ):
                                col[col_elim_index].indiv[iIndivNew].namda = temp_namda[iIndivNew]
                            iIndivNew+=1
                        else:
                            temp_const_limit+=1
                    elif (k != (n_col - n_col_elim - 1)):
                        k+=1
                        continue
                    else:
                        k = 0
                        continue
                    k = iIndivNew % (n_col - n_col_elim)

                if (col[col_elim_index].size <= 0 or col[col_elim_index].size != col_size_temp):
                    raise SystemError("error")
        #end of elimination
        ##if (parameters.IGD_on):
       # #calculate the IGD
           
            ##IGD_val = IGD.IGD_calc(Pareto_F,real_PF);
            ##print(" " ,IGD_val, " ");

        ##Calculate the IGD every generation and save PF
        ##if (parameters.PERF_VAL_GEN == True):
          ##   Pareto_F.Pareto_Search(PF_storage)
            ## PF_storage = []
            


    #end of loop for number of generations

    if (n_col == 1 and MODE_vector[0] != "Normal"):
        Pareto_F.Pareto_Search(col[0].indiv)
    if (parameters.PF_REFINE == True):
        Pareto_F.Pareto_Refine()
    PF_data.write("Final Pareto Front");
    PF_data.write("\n");
    #save the results
    for iPF in range (0,Pareto_F.size):
        #Save objectives
        for iObj in range (0,nobj):
            PF_data.write(str(Pareto_F.indiv[iPF].fitness[iObj]));
            PF_data.write(",");

        for iVar in range(0,nvar):
            PF_data.write(",");
            PF_data.write(str(Pareto_F.indiv[iPF].code[iVar]));
        PF_data.write("\n");

    if (parameters.IGD_on):
        #calculate the IGD
        IGD_val = IGD.IGD_calc(Pareto_F,real_PF);

        IGD_val_final = IGD_val;
        IGD_storage.append(IGD_val);

    if (parameters.HV_on):
        #calculate the IGD
        HV_val = HV.HV_calc(Pareto_F,real_PF);

        HV_val_final = HV_val;
        HV_storage.append(HV_val);

    #close data files
    PF_data.close();

    #save the time
    run_time = time.perf_counter() - run_t;

    print("Run: #",i_run," finished time:")
    print(run_time,"s\n")
    if (parameters.IGD_on):
        print("IGD:", IGD_val_final, "\n");

    if (parameters.HV_on):
        print("HV:", HV_val_final, "\n");
    #save the run data to files
    Add_func.Output_Add_Time(results_time,run_time,IGD_val_final,HV_val_final,nfes,i_run,iGen);

    #save the run data
    run_data.Add(run_time,iGen,IGD_val_final,HV_val_final);


#end of the loop for the number of runs

#Calculate the standard deviation and average values for the current run
run_data.Std_Dev_Calculation();
#save the data to the file
Add_func.Output_Add_Summary(results_GA_Data,run_data)


##Calculate the average IGD over the runs for each generation number and save it to file and make graph



results_index.close();
results_time.close();
results_GA_Data.close();

print("Optimisation Finished\n")
wait = input("PRESS ENTER TO CONTINUE.")    


 #x_go = [];
  #          y_go = [];
   #         x1_go = [];
    #        y1_go = [];
     #      
      #      for icol in range(0,n_col):
      #
             #   for ind in range(0,col[icol].size):
              #      x1_go.append(col[icol].indiv[ind].fitness[0]);
               #     y1_go.append(col[icol].indiv[ind].fitness[1]);
 #           for iPF in range (0,Pareto_F.size):
  #              
   #     #Save objectives
     #           x_go.append(Pareto_F.indiv[iPF].fitness[0]);
    #            y_go.append(Pareto_F.indiv[iPF].fitness[1]);
     #       trace = go.Scatter( x = x_go, y = y_go, mode = 'markers',);
      #      trace1 = go.Scatter( x = x1_go, y = y1_go, mode = 'markers',);
       #     layout = go.Layout(xaxis=dict(range=[0,2]),yaxis=dict(range=[0,2]));
        #    data = [trace1,trace];
         #   fig = go.Figure(data = data,layout=layout)
          # 
           # py.plot(fig,filename = 'basic-scatter');