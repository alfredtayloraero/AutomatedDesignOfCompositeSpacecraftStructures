import math
import numpy as np
import parameters

#Fixed parameters
FibreThermalCoeff, MatrixThermalCoeff = -0.38e-6,71e-6 #googles
FibreDensity, MatrixDensity = 1.8e3,1.17e3 #kg/m3 guess from feeg1002 textbook
MaxTempGradient = 260 #Prelaunch conditions (ambient-cryogenic LOX)
FibreYield, MatrixYield = 4.9e9,80e6 #guess from feeg1002 textbook
FibreSpecCost, MatrixSpecCost = 500,25
FibreModulus, MatrixModulus = 230e9,3e9 #guess from feeg1002 textbook
BucklingFac, Length, Radius = 1,41,1.85 # ,F9 Diameter/2
InternalP, AxialP = 3.447e5,3.789e6 #3.447bar taken from an estimation 
                        #Axial load is 4g acceleration on 97570kg second stage


class Function(object): #Function class definition
    """description of class"""
    def __init__(self, fname, varsn, objsn, consn):
        self.name_func = fname; #Name of the function
        self.num_vars = varsn;  #Number of variables
        self.num_objs = objsn;  #Number of objectives
        self.num_cons = consn;  #Number of constraints
        self.bound_min = [];    #Minimum bounds of functions
        self.bound_max = [];    #Maximum bounds of functions

class CGRID(Function):
    def __init__(self):
        Function.__init__(self,"CGRID",5,6,0) #variables objectives constraints
        #Inputs : Helical rib angle, Distance between circular ribs, Distance between helical ribs, 
        #Circular rib thickness, Helical rib thickness, Fibre volume fraction (Matrix volume fraction)
        #Height of stiffeners
        self.Bound_Set();   #Sets variable bounds
        #self.num_objs=6;
        #self.num_vars=7;

    def Bound_Set(self):     #Automatically sets bounds of functions between zero and one
        self.bound_max = [1.571];
        self.bound_min = [0.001];
        for i in range(1,self.num_vars):
            self.bound_min.append(0.001);
        self.bound_max.append(1)
        self.bound_max.append(1)
        self.bound_max.append(1)
        self.bound_max.append(1)

    def Fitness_C(self, code):  #Fitness function for UF1
        HelAngle, HelRatio, CircRatio, FibreVolumeFrac, Height = code[0], code[1], code[2] , code[3], code[4]    
        CompositeModulus = FibreModulus*FibreVolumeFrac+MatrixModulus*(1-FibreVolumeFrac)
        CompositeYieldStress=FibreYield*FibreVolumeFrac+MatrixYield*(1-FibreVolumeFrac)
        ShellBucklingForce= 2*math.pi*Height**2*math.cos(HelAngle)**2*math.sqrt(CompositeModulus**2*(2/3*CircRatio*HelRatio))
        LocalBucklingForce= BucklingFac*math.pi**3/3*CompositeModulus*Height*(HelRatio)**3*(math.sin(2*HelAngle))**2*(math.cos(HelAngle))**2
        CriticalHelicalForce= 4*math.pi*CompositeYieldStress*Radius*Height*(HelRatio)*math.cos(HelAngle)**2
        PressureStress=(2*HelRatio+HelRatio)*InternalP/(math.pi*Radius*Height**2)
        ThermalStress=CompositeModulus*MaxTempGradient*(MatrixThermalCoeff-(MatrixThermalCoeff-FibreThermalCoeff)*((2-FibreVolumeFrac)*FibreVolumeFrac*FibreModulus-
                                                                       (1+FibreVolumeFrac)*(CompositeModulus-MatrixModulus*(1-FibreVolumeFrac)))/(1-2*FibreVolumeFrac))
        Mass=2*math.pi*Radius*Length*Height*(FibreVolumeFrac*FibreDensity+(1-FibreVolumeFrac)*MatrixDensity)*(2*HelRatio+HelRatio)
        Cost=Mass*(FibreSpecCost*FibreVolumeFrac+MatrixSpecCost*(1-FibreVolumeFrac))
        
        fitness1, fitness2, fitness3, fitness4, fitness5, fitness6= 0.0, 0.0, 0.0, 0.0, 0.0, 0.0;
        '''fitness1 = 1.
        if (ShellBucklingForce<AxialP):
            fitness1=0.
        fitness2 = 1.
        if (LocalBucklingForce<AxialP):
            fitness2 = 0.
        fitness3 = 1.
        if (CriticalHelicalForce<AxialP):
            fitness3 = 0.
        fitness4 = 1.
        if (PressureStress+ThermalStress>CompositeYieldStress):
            fitness4 = 0.
        '''
        safety=[]
        safety.append(ShellBucklingForce/AxialP)
        safety.append(LocalBucklingForce/AxialP)
        safety.append(CriticalHelicalForce/AxialP)
        safety.append(ShellBucklingForce/AxialP)
        fitness1 = np.exp(-max(AxialP/(ShellBucklingForce),0.2))/4 #minimising inverse of the safety factor divided by a fudge factor
        fitness2 = np.exp(-max(AxialP/(LocalBucklingForce),0.2))/4 #as the results were super skewed towards safety factor
        fitness3 = np.exp(-max(AxialP/(CriticalHelicalForce),0.2))/4
        fitness4 = np.exp(-max((PressureStress+ThermalStress)/(CompositeYieldStress),0.2))/4
        fitness5 = np.exp(-Mass/1e3)
        fitness6 = np.exp(-Cost/1e6)
        
        fitness = [fitness1,fitness2,fitness3,fitness4,fitness5,fitness6]; #returns the fitness values for the problem

        return fitness;

    #This is probably where it's going wrong
    def Plot_PF(self,size): #
        Real_PF = [];

        for i in range (0,size):
            temp_f1 = i/(size-1.);
            temp_f2 = 1- math.sqrt(temp_f1);

            PF_val = [temp_f1, temp_f2]
            Real_PF.append(list(PF_val));

        return Real_PF




class UF1(Function):
    def __init__(self):
        Function.__init__(self,"UF1",30,2,0) #variables objectives constants
        self.Bound_Set();   #Sets variable bounds

    def Bound_Set(self):     #Automatically sets bounds of functions between zero and one
        self.bound_max = [1];
        self.bound_min = [0];
        for i in range(1,self.num_vars):
            self.bound_max.append(1);
            self.bound_min.append(-1);

    def Fitness_C(self, code):  #Fitness function for UF1
        SumJ1 = 0.;
        SumJ2 = 0.;
        sizeJ1 = 0;
        sizeJ2 = 0;
        for i in range (2, self.num_vars + 1):
            y = code[i - 1] - math.sin((6. * parameters.pi * code[0]) + (i * parameters.pi) / self.num_vars)
            if (i % 2 == 1):
                sizeJ1 += 1;
                SumJ1 += pow(y, 2)
            else:
                sizeJ2 +=1;
                SumJ2 += pow(y, 2)
        
        fitness1 = code[0] + (2. / sizeJ1)*SumJ1; #Calculates the first fitness function using J1 params
        fitness2 = 1. - math.sqrt(code[0]) + (2. / sizeJ2)*SumJ2; #Calculates the 2nd fitness function using J2 params

        fitness = [fitness1,fitness2]; #returns the fitness values for the problem

        return fitness;

    def Plot_PF(self,size): #
        Real_PF = [];

        for i in range (0,size):
            temp_f1 = i/(size-1.);
            temp_f2 = 1- math.sqrt(temp_f1);

            PF_val = [temp_f1, temp_f2]
            Real_PF.append(list(PF_val));

        return Real_PF



class UF2(Function): #repeat of UF1 for a different problem
    def __init__(self):
        Function.__init__(self,"UF2",30,2,0)
        self.Bound_Set();

    def Bound_Set(self):
        self.bound_max = [1];
        self.bound_min = [0];
        for i in range(1,self.num_vars):
            self.bound_max.append(1);
            self.bound_min.append(-1);
            
    def Fitness_C(self, code):
        SumJ1 = 0.;
        SumJ2 = 0.;
        sizeJ1 = 0;
        sizeJ2 = 0;

        for i in range (2, self.num_vars + 1):
            if (i % 2 == 1):
                sizeJ1 += 1;
                SumJ1 += pow(code[i-1] - (0.3*pow(code[0], 2)*math.cos(24 * parameters.pi*code[0] + 4. * i*parameters.pi / self.num_vars) + 0.6*code[0])*math.cos(6 * parameters.pi*code[0] + i*parameters.pi / self.num_vars), 2);
            else:
                sizeJ2 +=1;
                SumJ2 += pow(code[i-1] - (0.3*pow(code[0], 2)*math.cos(24 * parameters.pi*code[0] + 4. * i*parameters.pi / self.num_vars) + 0.6*code[0])*math.sin(6 * parameters.pi*code[0] + i*parameters.pi / self.num_vars), 2);
        
        fitness1 = code[0] + (2. / sizeJ1)*SumJ1;
        fitness2 = 1. - code[0]**0.5 + (2. / sizeJ2)*SumJ2;

        fitness = [fitness1,fitness2];

        return fitness;

    def Plot_PF(self,size):
        Real_PF = [];

        for i in range (0,size):
            temp_f1 = i/(size-1.);
            temp_f2 = 1- math.sqrt(temp_f1);

            PF_val = [temp_f1, temp_f2]
            Real_PF.append(list(PF_val));

        return Real_PF