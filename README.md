# Diffusion Capacity of Single and Interconnected Networks

Tiago A. Schieber, Laura C. Carpi, Panos M. Pardalos, Cristina Masoller, Albert Díaz-Guilera and Martín G. Ravetti

This repository represents an R framework of functions and examples for computing the Diffusion Capacity of Single and Interconnected Networks proposed by Schieber et al. in 2023.

## Running and testing

Run functions.r file.

load("g_example.RData")

For the Single Layer Diffusion Capacity:

sldc(g_example)

returns the diffusion capacity of each network's vertice.

Multilayered Case:

load("layers.RData")

mldc(g_example,layers)

returns the Multilayer Diffusion Capacity or each network vertice in each layer. The "layers" variable is a list of indices in each layer. In this case, 2 layers the first containing 10 vertices and the second 90 vertices.

## The Heat and Kuramoto Models

The heat and Kuramoto models can be easily implemented in R. 

See heat_model.r

See kuramoto.r

## Climate Network Data

The Climate Network Dataset can be freely downloaded from https://psl.noaa.gov/data/gridded/data.ncep.reanalysis.html 



