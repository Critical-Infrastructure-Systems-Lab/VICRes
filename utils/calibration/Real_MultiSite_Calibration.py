# LIBRARIES
import os
import time
import shutil
import datetime
import numpy as np
import pandas as pd
import multiprocessing
from indices3 import NSE3, TRMSE3, MSDE, ROCE
from platypus import Real, Problem, EpsNSGAII, ProcessPoolEvaluator, Hypervolume, nondominated


def viccall(vars):
    print(f"Vars values: {vars}")

    global VIC_count, sim_period, hotstart, obj_fn, var_no, zones
    dtdr = '/mnt/d/VIC_and_VICRes/Mekong6_multisite'  #'/mnt/d/VIC_and_VICRes/Mekong1'    ,  'D:/VIC_and_VICRes/Mekong1'
    
    # Save the current solution for tracking purposes
    vars_file_path = os.path.join(dtdr, 'Calibration', 'performance_tracking', f'Vars_Round{VIC_count}.txt')
    np.savetxt(vars_file_path, vars, fmt='%s', delimiter='\n')
    site_names = os.listdir(os.path.join(dtdr, 'Routing', 'Parameters', 'OBS_Discharge'))  # List of the observed discharge sites  
    site_count = len(site_names)                                                           # Number of observed discharge sites 
    
    # Get the rank and process ID
    rank = multiprocessing.current_process()._identity[0]
    print("Thread no: ", rank, " Round: ", VIC_count) 
    
    # # =================== Modify soil parameters    
    os.chdir(os.path.join(dtdr, 'Calibration'))  
    zone_lines = pd.read_csv('zone.txt', delimiter='\t', header=None)
    zone_lines.columns = [f'Col{i+1}' for i in range(zone_lines.shape[1])]
        
    os.chdir(os.path.join(dtdr, 'Rainfall_Runoff', 'Parameters'))  
    soil_lines = pd.read_csv('soil.txt', delimiter='\t', header=None)
    soil_lines.columns = [f'Col{i+1}' for i in range(soil_lines.shape[1])]
    
    var_lines = pd.read_csv(vars_file_path, delimiter='\t', header=None)
    var_lines.columns = [f'Col{i+1}' for i in range(var_lines.shape[1])]
    
    group_size = var_no//zones
    all_selected_values = []
    zon = 0
    for i in range(zones):
        # Get a group of 9 values from the DataFrame
        selected_values = var_lines.iloc[i * group_size:(i + 1) * group_size].values.flatten()
        selected_values = [float(value) for value in selected_values]
        all_selected_values.append(selected_values)
        zon+=1
        filtered_rows = zone_lines[zone_lines['Col5'] == zon]
        
        # Find rows in soil_lines that match the latitude (Col3) and longitude (Col4) of filtered_rows
        matched_indices = soil_lines[
            (soil_lines['Col3'].isin(filtered_rows['Col3'])) &
            (soil_lines['Col4'].isin(filtered_rows['Col4']))
        ].index
        
        # Columns to be updated with selected_values
        columns_to_update = ['Col5', 'Col6', 'Col7', 'Col8', 'Col9', 'Col19', 'Col20']

        # Update matched rows in soil_lines with selected_values
        for col, value in zip(columns_to_update, selected_values[:7]):
            soil_lines.loc[matched_indices, col] = round(float(value), 6)
            
    # Save the updated soil back to a file
    soil_lines.to_csv(os.path.join(dtdr, 'Calibration', 'soil_tracking', 'Soilparam_Round' + str(VIC_count) + '.txt'), sep='\t', index=False, header=False)
    soil_lines.to_csv(os.path.join(dtdr, 'Rainfall_Runoff', 'Parameters', 'soil.txt'), sep='\t', index=False, header=False)
    
    all_selected_values = np.array(all_selected_values)
    zone_mean = np.mean(all_selected_values, axis=0)
    
    # # =================== Modify routing parameters
    os.chdir(os.path.join(dtdr, 'Routing', 'Parameters'))
    text_file = open('configurationOP6.txt','r')			    # Modify flow routing file
    lines = text_file.read().splitlines()						# 1st running, lines in the configuration file does not contain \n
    text_file.close()
    if len(lines) < 37:
        raise ValueError("configuration file does not have enough lines.")
    year_line = lines[27].strip()
    parts = year_line.split()
    start_year = int(parts[0])
    end_year = int(parts[2])
    # with open('configurationOP6.txt','w') as my_csv:			# Modify when needed
        # for i in range(37):
            # if ((i==5) and (dec_var[8]>0)):
                # my_csv.write("%f\n"%(zone_mean[7]))				# velocity
            # elif ((i==8) and (dec_var[9]>0)):
                # my_csv.write("%f\n"%(zone_mean[8]))				# diffusivity
            # else:
                # my_csv.write("%s\n"%(lines[i]))

                                         
    #================Rainfall-Runoff and Routing ============           
    # Run Rainfall_Runoff model 
    os.chdir(os.path.join(dtdr, 'Rainfall_Runoff', 'Output'))         
    os.system('rm -f *')
    sourcecode_dir = os.path.join(dtdr, 'Rainfall_Runoff', 'SourceCode')
    os.chdir(sourcecode_dir)    
    os.system('./vicNl -g ../Parameters/globalparam.txt')     
    
                
    # Run Routing model 
    sourcecode_dir = os.path.join(dtdr, 'Routing', 'SourceCode')
    os.chdir(sourcecode_dir)                              
    os.system('rm -f *.uh_s')
    os.system('./rout ../Parameters/configurationOP6.txt')
    #========================================================

    
    os.chdir(os.path.join(dtdr, 'Routing', 'Results', 'ResOP6'))
    all_stations = 'OUTPUT_Round'+str(VIC_count)+'.day'
    shutil.copy('OUTPUT.day', dtdr + f'/Calibration/station_tracking/{all_stations}') 
    
    sum_nash = 0
    sum_trmse = 0
    sum_msde = 0
    sum_roce = 0
    nash_values, trmse_values, msde_values, roce_values = [], [], [], []
    for site in site_names:
        
        station = site[:-4]
        station_name = f'{station}_Round'+str(VIC_count)+'.day'
        os.chdir(os.path.join(dtdr, 'Routing', 'Results', 'ResOP6'))
        shutil.copy(f'{station}.day', station_name)
        shutil.copy(f'{station}.day', dtdr + f'/Calibration/station_tracking/{station_name}')
    
        # Read observed discharge (values)
        os.chdir(os.path.join(dtdr, 'Routing', 'Parameters', 'OBS_Discharge'))
        rs_q = pd.read_csv(f'{station}.csv')
        subset = rs_q[(rs_q['Year'] >= start_year) & (rs_q['Year'] <= end_year)]
        RS_data = subset.to_numpy()
        gaudata = [float(RS_data[i][3]) for i in range(len(RS_data))]
    
        # Read modelled discharge (values)
        os.chdir(os.path.join(dtdr, 'Routing', 'Results', 'ResOP6'))
        sim_q = open(station_name,'r') 
        lines = sim_q.read().split('\n')
        Sim_data = [[0 for x in range(4)] for y in range(len(lines))]
        count_no = 0
        for line in lines:
            try:
               year,month,day,flows = filter(None,line.split(' '))
               Sim_data[count_no][0] = int(year)
               Sim_data[count_no][1] = int(month)
               Sim_data[count_no][2] = int(day)
               Sim_data[count_no][3] = float(flows)
               count_no+=1
            except:
               print("...")                                        
        sim_q.close()
        VIC_data = [[0 for x in range(4)] for y in range(len(RS_data))]
        for i in range(len(RS_data)):
            VIC_data[i][0] = int(RS_data[i][0])
            VIC_data[i][1] = int(RS_data[i][1])           
            VIC_data[i][2] = int(RS_data[i][2])
            for j in range(len(lines)):
               if Sim_data[j][0] == VIC_data[i][0] and Sim_data[j][1] == VIC_data[i][1] and Sim_data[j][2] == VIC_data[i][2]:
                VIC_data[i][3] = float(Sim_data[j][3])	
    
        # Calculate objective functions (NSE, TRMSE, MSDE, ROSE)
        nortrmse = 0
        normsde = 0
        meangautrmse = 0
        for i in range(len(RS_data)):
            meangautrmse+=gaudata[i]
            normsde+=pow((gaudata[i]-gaudata[i-1])*2,2)
        meangautrmse = meangautrmse /(len(RS_data))
        for i in range(len(RS_data)):
            zobs = (pow((abs(gaudata[i])+1),0.3)-1)/0.3
            nortrmse+=pow(zobs-meangautrmse,2)
        nash  = NSE3(gaudata,VIC_data,3)   
        trmse = TRMSE3(gaudata,VIC_data,3)
        msde  = MSDE(gaudata,VIC_data,3)
        roce  = ROCE(gaudata,VIC_data,3)
        
        # Store actual objective function values for each site
        nash_values.append(nash)
        trmse_values.append(trmse)
        msde_values.append(msde)
        roce_values.append(roce)
        
        sum_nash += nash
        sum_trmse += trmse
        sum_msde += msde
        sum_roce += roce
    
    # Calculate average of each objective function
    avg_nash = sum_nash / site_count
    avg_trmse = sum_trmse / site_count
    avg_msde = sum_msde / site_count
    avg_roce = sum_roce / site_count 
        
    # Normalized NSE
    if avg_nash<=0:
        standnash = 1
    else:
        standnash = 1 - avg_nash
        
    # Normalized TRMSE
    standtrmse = avg_trmse/pow(nortrmse/(sim_period-hotstart),0.5)
    if (standtrmse>1):
        standtrmse = 1
    elif (standtrmse<0):
        standtrmse = 0
        
    # Normalized MSDE
    standmsde = avg_msde / normsde * (sim_period-hotstart)
    if (standmsde>1):
        standmsde = 1
    elif (standmsde<0):
        standmsde = 0
        
    # Normalized ROCE
    standroce = 0
    if (avg_roce > 1):
        standroce = 1
    elif (avg_roce<0):
        standroce = 0   
   
    calibrationresults = []
    if (obj_fn[0]>0):
        calibrationresults.append(standnash)
        print("NSE = ",avg_nash)
    if (obj_fn[1]>0):
        calibrationresults.append(standtrmse)
        print("TRMSE = ",avg_trmse)
    if (obj_fn[2]>0):
        calibrationresults.append(standmsde)
        print("MSDE = ",avg_msde)
    if (obj_fn[3]>0):
        calibrationresults.append(standroce)
        print("ROCE = ",avg_roce)
                
    # Tracking
    track = [nash_values, trmse_values, msde_values, roce_values, vars]
    os.chdir(os.path.join(dtdr, 'Calibration', 'performance_tracking'))
    np.savetxt('Performance_' + 'Round' + str(VIC_count) + '.txt', track, fmt='%s')
        
    VIC_count += 1
    return calibrationresults


# ======  MAIN  =====  MAIN  ======  MAIN  =========  MAIN  =======
# Start simulations
if __name__ == '__main__':

    dtdr = '/mnt/d/VIC_and_VICRes/Mekong6_multisite'
    start = datetime.datetime.now()
    VIC_count = 1

    # Create required directories
    for dir_name in ['performance_tracking', 'station_tracking', 'soil_tracking']:
        directory = os.path.join(dtdr, 'Calibration', dir_name)
        os.makedirs(directory, exist_ok=True)
            
    # Define optimization problem (objective functions, decision variables)
    # Read operation file (calibration.txt)
    os.chdir(dtdr+'/Calibration')
    cali_setup = open('multicalibration.txt','r')
    lines = cali_setup.read().split('\n')
    sim_period = int(lines[2].split('\t')[0])   # simulation period (days)
    hotstart   = int(lines[3].split('\t')[0])   # hotstart period (days)
    iteration  = int(lines[4].split('\t')[0])   # number of iterations/evaluations
    pop_size   = int(lines[5].split('\t')[0])   # population size
    cores      = int(lines[6].split('\t')[0])   # number of cores
    stations   = int(lines[24].split('\t')[0])
    zones      = int(lines[26].split('\t')[0])

    
    # Define decision variables
    dec_var = [0 for x in range(10)]
    var_no  = 0     # number of selected decision variables (up to 10)
    sel_var = []    # selected variabes     
    for zo in range(zones):
        for i in range(10):
            dec_var[i] = int(lines[i+8].split('\t')[0])
            if (dec_var[i]==1):
                var_no += 1
                if (i==0):
                    sel_var.append(Real(0.05, 0.9))		# b
                elif (i==1):
                    sel_var.append(Real(0.000001, 0.95))	# Ds
                elif (i==2):
                    sel_var.append(Real(2, 30))			# Dmax
                elif (i==3):
                    sel_var.append(Real(0.05, 0.999999))	# Ws
                elif (i==4):
                    sel_var.append(Real(1, 3))			# c
                elif (i==5):
                    sel_var.append(Real(0.05, 0.25))		# d1
                elif (i==6):
                    sel_var.append(Real(0.3, 1.5))		# d2
                elif (i==7):
                    sel_var.append(Real(0.04, 0.052))		# n
                elif (i==8):
                    sel_var.append(Real(2, 5))			# v
                elif (i==9):
                    sel_var.append(Real(250, 4000))		# df									

                
    # Define objective functions          
    obj_fn = [0 for x in range(4)]            
    obj_no = 0      # number of selected objective functions (up to 4)        
    min_obj = []
    max_obj = []
    eps = []
    for i in range(4):
        obj_fn[i] = int(lines[i+19].split('\t')[0])
        if (obj_fn[i]==1):
            obj_no += 1
            min_obj.append(0)   # normalize objective function values 
            max_obj.append(1)   # with range from 0 to 1
            eps.append(0.001)
    
    # Setup optimization parameters
    problem = Problem(var_no, obj_no)
    problem.types[:] = sel_var
    problem.function = viccall
    hyp = Hypervolume(minimum = min_obj, maximum = max_obj)
    x = []
    generation = 1
    
    with ProcessPoolEvaluator(cores) as evaluator:
        algorithm = EpsNSGAII(problem, eps, population_size=pop_size, evaluator=evaluator)
        print(f"Starting optimization. Population size: {pop_size}, Max evaluations: {iteration}")
        while algorithm.nfe < iteration:  
            print(f"Generation = {generation}") 
            algorithm.step() 
            y = hyp.calculate(algorithm.result)
            x.append(y)
            generation +=1
        print(f"Completed {algorithm.nfe} evaluations")  
    nondominated_solutions = nondominated(algorithm.result)

    # Finish and write out results
    os.chdir(dtdr+'/Calibration')
    print("Finish running simulations! See vicres_calibration.obj.txt and vicres_calibration.var.txt for results.")
    np.savetxt("vicres_calibration.obj.txt",[s.objectives[:] for s in nondominated_solutions],fmt="%s")
    np.savetxt("vicres_calibration.var.txt",[s.variables[:] for s in nondominated_solutions],fmt="%s")
    print("Hypervolume indicator:")
    print(x)
    end = datetime.datetime.now()
    print("Start",start)
    print("End",end)
    print(f"Total time: {end - start}")

# End of file