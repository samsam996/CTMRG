
import numpy as np
import ncon as nc
import matplotlib.pyplot as plt
import argparse
import os 
import sys

from src import CTMRG_C4v
from src import CTMRG_nosymm
from src import Ising2D
from src import qStatePotts
from src import AshkinTeller

# from src import one_site_obs
# from src import partition_function


if __name__ == '__main__':

    parser = argparse.ArgumentParser(description="All variables are set in the data.txt file")
    parser.add_argument('--input', type=str, required=True, help='Input file path.')
    # parser.add_argument('--output', type=str, required=True, help='Output file path.')
    # parser.add_argument('--verbose', action='store_true', help='Enable verbose output.')
    # parser.add_argument('--chi', type=int, required=True, help='Value of chi')
    # parser.add_argument('--temp_low', type=float, required=True, help='Lower temperature')
    # parser.add_argument('--temp_high', type=float, required=True, help='Value of chi')
    # parser.add_argument('--nbstep', type=int, required=True, help='Lower temperature')

    args = parser.parse_args()


    Tc = 2/(np.log(np.sqrt(2)+1))


    with open(args.input, 'r') as file:
        data = file.readlines()

    # Parse the file content and set variables
    variables = {}
    with open(args.input, 'r') as file:
        for line in file:
            line = line.strip()  # Remove leading/trailing whitespace
            if line.startswith('#') or not line:  # Ignore comments or empty lines
                continue
            key, value = line.split('=')  # Split by '=' to get key-value pairs
            key = key.strip()
            value = value.strip()
            
            # Try to determine the type of the value
            try:
                # Attempt to convert to an integer
                if '.' in value:  # Check if the value contains a decimal point
                    variables[key] = float(value)  # Store as float
                else:
                    variables[key] = int(value)  # Store as int
            except ValueError:
                # If conversion fails, treat it as a string
                variables[key] = value

    chi = variables['chi']
    temp_low = variables['temp_low']
    temp_high = variables['temp_high']
    nbstep = variables['nbstep']
    bc = variables['boundary_condition']
    model = variables['model']
    q = variables['q']
    J = variables['J']
    lmbda = variables['lambda']

    temp = np.linspace(temp_low,temp_high,nbstep)
    magne = []
    free_energy = []

    for i in range(len(temp)):
  
        if model == "Ising":
            a,b = Ising2D(temp[i],J)
        elif model=="Potts":
            a,b = qStatePotts(temp[i],J,q)
        elif model=="AshkinTeller":
            a,b = AshkinTeller(temp[i],J,lmbda)
        else: 
            print("Wrong choice of model, please chose Ising, Potts or AshkinTeller") 
            sys.exit("Exiting the program due to invalid model.")



        d = np.shape(a)[0]
        my_ctm = CTMRG_C4v(chi,d)
        my_ctm2 = CTMRG_nosymm(chi,d)
        my_ctm2.fixed_boundary_condition(a)
        my_ctm2.evolution(a)
  
       

        if bc == "fixed":
            my_ctm.fixed_boundary_condition(a)
        elif bc == "open":
            my_ctm.open_boundary_condition(a)
        else:
            print("Wrong choice of boundary conditions, please chose set or fixed") 
            sys.exit("Exiting the program due to invalid boundary condition.")

        free_ener = 0
        err = 1
        it = 0


        while abs(err) > 1e-8 and it < 100:
            it += 1
            # my_ctm.evolution(a)
            # my_ctm2.evolution(a)
            my_ctm2.evolution(a)

            tmp = free_ener
            free_ener = -temp[i]*np.log(my_ctm2.partition_function(a))
            free_ener = my_ctm2.partition_function(a)

            # free_ener = -temp[i]*np.log(my_ctm2.partition_function(a))
            err = free_ener - tmp
            print('it : ',it, ', error : ',err)

        magne.append(my_ctm.one_site_obs(a,b))
        free_energy.append(free_ener)



    with open('output.txt', 'w') as file:
        for i in range(len(temp)):
            file.write(f"{temp[i]} {magne[i]} {free_energy[i]}\n")

    plt.plot(temp, magne, 'o')
    plt.xlabel('temperature')
    plt.ylabel('magnetisation')

    # Specify the folder where the figure is saved and create it if it doesn't exist
    folder = 'figures'  
    if not os.path.exists(folder):
        os.makedirs(folder)

    # Save the figure in the specified folder
    plt.savefig(os.path.join(folder, 'magnetisation.png'))  
