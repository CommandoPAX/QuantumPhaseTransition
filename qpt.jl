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
    sites = siteinds("Arg",N)
    U = 1
    J = 1 

    # Trying to build the Hamiltonian
    os = OpSum()
    for j=1:N-1
        os += - J,"b",j,"b dag",j+1
        os += - J,"b dag",j,"b",j+1
        os += U/2,"b dag", "b", j,("b dag", "b", j, -1)
    end
    H = MPO(os,sites)

    psi0 = random_mps(sites;linkdims=10)

    nsweeps = 5
    maxdim = [10,20,100,100,200]
    cutoff = [1E-10]

    energy,psi = dmrg(H,psi0;nsweeps,maxdim,cutoff)

    return

end