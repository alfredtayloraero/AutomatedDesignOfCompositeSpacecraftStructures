import math

class GA_data(object):
    def __init__(self, IGD_o, HV_o):
        
        self.HV_on = HV_o
        self.IGD_on = IGD_o
        
        self.time = []
        self.time_struct = [0.,0.,0.,0.] #min,max,avereage,std
        if (IGD_o):
            self.IGD = []
            self.IGD_struct = [0.,0.,0.,0.]
        if (HV_o):
            self.HV = []
            self.HV_struct = [0.,0.,0.,0.]
        self.generation = []
        self.generation_struct =[0,0,0,0]

        return

    def Add(self, time_val, generation_val, IGD_val, HV_val):
        
        self.time.append(time_val)
        self.generation.append(generation_val)

        if (self.IGD_on):
            self.IGD.append(IGD_val)
        if (self.HV_on):
            self.HV.append(HV_val)

    def Average_Calc(self):
        size = len(self.time)

        self.time_struct[2] = sum(self.time)/size
        self.generation_struct[2] = sum(self.generation)/size

        if (self.IGD_on):
            self.IGD_struct[2] = sum(self.IGD)/size
        if (self.HV_on):
            self.HV_struct[2] = sum(self.HV)/size
        return


    def Std_Dev_Calculation(self):
        self.Average_Calc()

        #get the min max

        self.time_struct[0] = min(self.time)
        self.time_struct[1] = max(self.time)
        self.generation_struct[0] = min(self.generation)
        self.generation_struct[1] = max(self.generation)

        if (self.IGD_on):
            self.IGD_struct[0] = min(self.IGD)
            self.IGD_struct[1] = max(self.IGD)
        if (self.HV_on):
            self.HV_struct[0] = min(self.HV)
            self.HV_struct[1] = max(self.HV)



        size = len(self.time)

        time_std_temp = 0
        generation_std_temp = 0
        IGD_std_temp = 0
        HV_std_temp = 0    
        

        for i in range (0,size):
            time_std_temp += pow(self.time[i] - self.time_struct[2], 2);
            generation_std_temp += pow(self.generation[i] - self.generation_struct[2], 2);

            if (self.IGD_on):
                IGD_std_temp += pow(self.IGD[i] - self.IGD_struct[2], 2);  

            if (self.HV_on):
                HV_std_temp += pow(self.HV[i] - self.HV_struct[2], 2);

        self.time_struct[3] = math.sqrt(time_std_temp / size)
        self.generation_struct[3] = math.sqrt(generation_std_temp / size)
        
        if (self.IGD_on):
            self.IGD_struct[3] = math.sqrt(IGD_std_temp / size)
        if (self.HV_on):
            self.HV_struct[3] = math.sqrt(HV_std_temp / size)
        return
             

