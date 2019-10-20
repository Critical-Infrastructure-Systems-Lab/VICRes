import os, sys
import numpy
# Run rainfall-runoff model
os.chdir('Rainfall-runoffSetup/Results')
#os.system('rm fluxes*.*')
#os.system('rm snow*.*')
os.chdir('../../RoutingSetup/input')
#os.system('rm fluxes*.*')
os.chdir('../../Rainfall-runoff')
#os.system('./vicNl -g ../Rainfall-runoffSetup/globalparam.txt')    
#files = [f for f in os.listdir("../Rainfall-runoffSetup/Results") if os.path.isfile(os.path.join("../Rainfall-runoffSetup/Results",f))]
#for file in files:             
#        if (len(file)>=22):
#            if (len(file)==22):
#                toadox = file[7:14]
#                toadoy = file[15:22]
#           else:
#               toadox = file[7:15]
#                toadoy = file[16:23]
#            newfile = "../RoutingSetup/input/fluxes_"+toadoy+"_"+toadox            
#            file = "../Rainfall-runoffSetup/Results/" + file            
#            os.rename(file,newfile)    
# Run routing model
os.chdir('../Routing/SourceCode')
os.system('rm *.uh_s')
os.system('./rout ../../RoutingSetup/configuration.txt')