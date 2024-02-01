# Beam 3: Design optimization of RC beam 3 under simple optimization.

from MOSA.RCbeamSOB import rcb
import pandas as pd
import numpy as np
import plotly
import time

#First adjust the values of contraints and constructability function as per requirement
dff= []
rdw= []
rdc= []
rdco= []
sdw= []
sdc= []
sdco= []

min_wt=[]
min_cost=[]
min_co2=[]
t1= time.time()
for i in range (1):
    r1= rcb(width=300, depth=450, length=7, bending_moment=[303.93, 115.94, 298.9], shear_force=[179.37, 90.09, 177.56], cover=30, left_end_continous= True,right_end_continous= True)      #Creating beam object
    r1.beam_optimization()              #Performing multiobjective optimization

    

    df= r1.optimization_result                         #Storing Pareto front results
    dff.append(df) 
    min_weight_index= df["Weight (Kg)"].idxmin()

    min_cost_index= df["Cost (INR)"].idxmin()
    min_co2_index= df["CO2 (KgCO2e)"].idxmin()
    rd1= r1.rd_list[min_weight_index] 
    rd2= r1.rd_list[min_cost_index] 
    rd3= r1.rd_list[min_co2_index] 

    rdw.append(rd1)
    rdc.append(rd2)
    rdco.append(rd3)

    sd1= r1.sd_list[min_weight_index] 
    sd2= r1.sd_list[min_cost_index] 
    sd3= r1.sd_list[min_co2_index] 
    
    sdw.append(sd1)
    sdc.append(sd2)
    sdco.append(sd3)

    min_wt.append(df["Weight (Kg)"].min())
    min_cost.append(df["Cost (INR)"].min())
    min_co2.append(df["CO2 (KgCO2e)"].min())
    
t2= time.time()
comp_time= t2-t1
print (comp_time)

with pd.ExcelWriter('./Result_Beam3_WC_df.xlsx',) as writer:
     for i in range (1):
        dff[i].to_excel(writer, sheet_name=f'Opti{i+1}')


with pd.ExcelWriter('./Result_Beam3_WC_rd.xlsx',) as writer:
     for i in range (1):
        rdw[i].to_excel(writer, sheet_name=f'RDweight{i+1}')
        rdc[i].to_excel(writer, sheet_name=f'RDCost{i+1}')
        rdco[i].to_excel(writer, sheet_name=f'RDCO2{i+1}')

with pd.ExcelWriter('./Result_Beam3_WC_sd.xlsx',) as writer:
     for i in range (1):
        sdw[i].to_excel(writer, sheet_name=f'SDweight{i+1}')
        sdc[i].to_excel(writer, sheet_name=f'SDCost{i+1}')
        sdco[i].to_excel(writer, sheet_name=f'SDCO2{i+1}')
    