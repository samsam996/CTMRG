# CTMRG for the 2D Ising model

Simulate a series of 2D classical models with the Corner Transfer Matrix Renormalisation Group (CTMRG) algorithm [1]. The parameters are fixed in the parameters.txt file and the magnetisation is saved in the figure folder. The temperature, magnetisation and free energy per site are stored in the output.txt file. 

We can simulate three different models, the Ising model, the q-state Potts model and the Ashkin-Teller model, whose Hamiltonian are given by:

$$
H_{Ising} = -J\sum_{\langle i,j\rangle } \sigma_i \sigma_j \qquad \sigma \in {\pm 1}  \\
H_{qPotts} = -J\sum_{ \langle i,j \rangle } \delta_{\sigma_i \sigma_j} \qquad \sigma \in {1,.., q}  \\
H_{AT} = - \sum_{ \langle i,j \rangle } (\sigma_i \sigma_j + \tau_i \tau_j + \lambda (\sigma_i \sigma_j \tau_i \tau_j ) \qquad \sigma \in {\pm 1}, \tau \in {\pm 1}  
$$

To run the code, use the following command:

```
python main.py --input parameters.txt
```

Required packages are stored in the requirements.txt file

[1] https://arxiv.org/abs/cond-mat/9507087
