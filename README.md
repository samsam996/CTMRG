# CTMRG for the 2D Ising model

Simulate a series of 2D classical models with the Corner Transfer Matrix Renormalisation Group (CTMRG) algorithm [1]. The parameters are fixed in the parameters.txt file and the magnetisation is saved in the figure folder. The temperature, magnetisation and free energy per site are stored in the output.txt file. 

We can simulate three different models, the Ising model, the q-state Potts model and the Ashkin-Teller model, whose Hamiltonian are given by:

$$
H_{Ising} = -J\sum_{i,j}
$$

To run the code, use the following command:

```
python main.py --input parameters.txt
```

Required packages are stored in the requirements.txt file

[1] https://arxiv.org/abs/cond-mat/9507087
