using SparseArrays
using LinearAlgebra
using ITensors, ITensorMPS
using ProgressBars
using KrylovKit
using Plots
using LaTeXStrings
using ClassicalOrthogonalPolynomials
using Random

let 
    N = 100
    sites = siteinds("Qudit", N)
    U = 1
    J = 1 

    # Trying to build the Hamiltonian
    os = OpSum()
    for j=1:N-1
        os += -J,"A",j,"Adag",j+1
        os += - J,"Adag",j,"A",j+1
        os += U/2,"n",j,"n",j 
        os += -U/2,"n",j
    end
    H = MPO(os,sites)

    "psi0 = random_mps(sites;linkdims=10)

    nsweeps = 5
    maxdim = [10,20,100,100,200]
    cutoff = [1E-10]

    energy,psi = dmrg(H,psi0;nsweeps,maxdim,cutoff)"

    return
end