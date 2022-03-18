import Fit_Function
import GA_Operators

Function = Fit_Function.CGRID()
GA_mode = ["MOEAD"];
MLSGA_Hybrid = True;
MLSGA_norm_obj = False;		#Default false			#Define if the objectives for fitness evaluations should be normalised (for MSF and PSF, but true makes PSF inoperative)

MLSGA_n_col = 6;	#Default 8					#Numbers of collectives  ***with a choice of 4, 6, 8 
MLSGA_n_col_elim = 1;	#Default 1				#Number of collectives elimiated ***with a choice of 1, 2  (is overriden to 0 if MLSGA_Htbrid == false)
MLSGA_col_elim_limit = 10; #Default 4			#Define how many genrations have to pass before col elimination occurs - if 1 elimiantion will occur every generation - MTS have override that it change to 1 - line 494 in MLSGA.cpp
MLSGA_elite_size = 0.1; #Default 0.1					#How many individuals become elite <0.0-1.0>

#Parameter setting*/
Pop_size = 1200;	#Default 1000	#Population size
MLS = 7; #Default 7					#Index of the MLS type used
n_runs = 30; #Default 30 CEC'09		#Number of runs for each algorithm (90 total)
mut_prob = 0.1; #Default 0.08		#mutation probability - for real coded it is the absolute value, for binary it is divided by number of variables
cross_prob = 1;	#Default 0.70		#crossover probability
Crossover = GA_Operators.Crossover_SBX()
Mutation = GA_Operators.Mutation_poly()
Selection = GA_Operators.Roulette_Wheel()


#Termination condition*/
max_iterations = 300000;	#Goal: 300000			#max iterations

#Performance evaluation*/
HV_on = False;	
IGD_on = False;

#Dynamic functions parameter*/
n_steps_dyn = 10;	#Default 5							#Number of distinct steps; represents the severity of change
T_dyn = 10;	#Default 5								#Window where the dynamic problem remains cosntant
Reinit_Mode = "VP"; #list of reinitialisation modes that have to be used "None", "BR", "CER", "VP"
gen_memory_size = 2;	#Default 2						#How many generations are saved and used for reinitialisation
JY__window = 0;	#Default 100					#Window where JY functions remain atant


#Additional GA parameters*/
Di_c = 20;	#Default 20					#Distribution index of crossover
Di_m = 20;	#Default 20					#Distribution index of mutation - suttipong mut
Di_m2 = 0.00001;									#Distribution index of mutation - PG mut	(1. - random mutation, 0. - no mutation)
fit_index_sel = 5;  #Default 2			#Which fitness is used for indiv selection
fit_index_col_sel = 2;  #Default 1		#Which fitness is used for col selection 
Binary_string_size = 1;	#Default					#Lenght of each variable in binary string
MPCross_size = 2;		#Default 2						#Defines how many separation pos are in multi-po crossover (1 is two po crossver)

#GA parameters - MOEAD*/
MOEAD_niche_multi = 0.1;	#Default 0.1			#multiplier of the neighbour size (* pop size)
MOEAD_limit_multi = 0.01;	#Default 0.01			#multiplier of the maximal number of solutions replaced(* pop size)
MOEAD_mating_chance = 0.9;	#Default 0.9			#probability chance for mating
MOEAD_new_sol_multi = 5;		#Default 5				#How many new individuals are generated each generation = pop_size/MOEAD_new_sol_multi

#GA parameters - MTS*/
MTS_LocalSearch_test_amount = 1;	#Default 5				#How many times each local search is tested each generation
MTS_LocalSearch_amount = 1;	#Default 45				#How many local searches are made each generation
MTS_foreground_multp = 0.125;		#Default 0.125		#Size of the foreground
MTS_MPL1 = 9;	#Default 9									#Bonus 1
MTS_MPL2 = 2;	#Default 2									#Bonus 2
MTS_PF_refine_size = 1000;	#Default 100000		#Size of the PF after normal refining, when PF is too big (not the final one) for MTS only

#GA parameters - HEIA*/
HEIA_crossover_prob = 0.5;	#Default 0.5f		#Crossover probabilty for HEIA (between SBX and DE)

#Precision*/
prec = 6;	#Default 6						#Precision of results showing
PF_prec = 16;	#Default 16				#Precision of PF
PF_real_size = 10000; #Default 10000					#Size of the real PF - for the IGD calculation
PF_size = 1000;  #Default 100 - CEC'09				#Output size of the calculated PF / When PF_Refine is true
PF_refine_size = 1000;	#Default 700					#Size of the PF after normal refining, when PF is too big (not the final one) - should be at least 2x PF_size - not for MTS
PF_res = 6; #Default 5					#Pareto front search resolution
c_plot_res = 500; #Default 500						#Resolution of contour plot generation (higher -> better, but slower). Define size of the net
c_plot_deg = 10; #Default 10				#Degree of contour plot - increase difference between best and worst pos on contour plot, but high values may cover some of results
min_col_size_multi = 30;#Default 30		#Min size multiplier for collective size; col_min_size = pop_size / min_pop_size_multi
new_col_size_multi = 15;#Default 15		#Min size multiplier for collective size; col_min_size = pop_size / min_pop_size_multi
pow_eq_zero = 14;	#Default 14			#which decimal will decide if two numbers are equal
cons_check_param = -0.000001;	#Default -0.000001 - CEC09	#When rain is violated for IGD check


#Random*/
rnd_seed = 2;	#Default 1								#Seed for the pseudo random generator

#Output*/
real_PF_out = "PF_Real";						#File name for real Pareto Front
temp_data_size = 100;	#Default 500					#Define amount of runs after which the temp storage will be cleaned


#Input*/
input_name = "Input/LHCGA.csv";				#File name and path for distribution input

#Video*/
frame_h = 720; #Default 720				#Frame height for video output
frame_w = 1280;  #Default 1280			#Frame width for video output
FPS = 10;	#Default 3						#Frame per second
frame_img_multi = 1;	#Default 1			#How many times one frame is added - increase if want to make video longer and find a suitable frame
frame_skip = 5;	

#constants
pi = 3.141592653589793238462643383279502884197;
e = 2.71828182845904523536028747135266249775724709369995;
INF = 1.0e14;

#define
PF_REFINE = True;           #If POF should be refined to predefined size
PERF_VAL_GEN = False;       #if IGD and HV should be calculated every generation