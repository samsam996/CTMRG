# CTMRG for the 2D Ising model

Simulate the 2D ising model with the Corner Transfer Matrix Renormalisation Group (CTMRG) algorithm [1]. The parameters are fixed in the parameters.txt file and the magnetisation is saved in the figure folder. The temperature, magnetisation and free energy per site are stored in the output.txt file. 

To run the code, use the following command:

```
python main.py --input parameters.txt
```

Required packages are stored in the requirements.txt file

[1] https://arxiv.org/abs/cond-mat/9507087
