# Beam 1: Design optimization of RC beam 1 under market practise based constructability function.

from MOSA.RCbeamMOB import rcb
import pandas as pd
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

t0=time.time()
for i in range (3):
    r2= rcb(width=300, depth=450, length=4, bending_moment=[138.26, 64.07, 138.26], shear_force=[121.55, 76.94, 121.55], cover=30)     #Creating beam object
    r2.beam_optimization(nearest_value=5)              #Performing multiobjective optimization

    df= r2.optimization_result                         #Storing Pareto front results
    dff.append(df)
    min_weight_index= df["Weight (Kg)"].idxmin()
    min_cost_index= df["Cost (INR)"].idxmin()
    min_co2_index= df["CO2 (KgCO2e)"].idxmin()
    rd1= r2.rd_list[min_weight_index] 
    rd2= r2.rd_list[min_cost_index] 
    rd3= r2.rd_list[min_co2_index] 

    rdw.append(rd1) 
    rdc.append(rd2)
    rdco.append(rd3)

    sd1= r2.sd_list[min_weight_index] 
    sd2= r2.sd_list[min_cost_index] 
    sd3= r2.sd_list[min_co2_index] 
    
    sdw.append(sd1)
    sdc.append(sd2)
    sdco.append(sd3)

    min_wt.append(df["Weight (Kg)"].min())
    min_cost.append(df["Cost (INR)"].min())
    min_co2.append(df["CO2 (KgCO2e)"].min())
   

tf= time.time()
comp_time_f= tf- t0

print (comp_time_f/3)

with pd.ExcelWriter('./Result_Beam1_df.xlsx',) as writer:
     for i in range (3):
        dff[i].to_excel(writer, sheet_name=f'Opti{i+1}')
    
with pd.ExcelWriter('./Result_Beam1_rd.xlsx',) as writer:
     for i in range (3):
        rdw[i].to_excel(writer, sheet_name=f'RDweight{i+1}')
        rdc[i].to_excel(writer, sheet_name=f'RDCost{i+1}')
        rdco[i].to_excel(writer, sheet_name=f'RDCO2{i+1}')

with pd.ExcelWriter('./Result_Beam1_sd.xlsx',) as writer:
     for i in range (3):
        sdw[i].to_excel(writer, sheet_name=f'SDweight{i+1}')
        sdc[i].to_excel(writer, sheet_name=f'SDCost{i+1}')
        sdco[i].to_excel(writer, sheet_name=f'SDCO2{i+1}')