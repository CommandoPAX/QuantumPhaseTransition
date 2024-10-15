using SparseArrays
using LinearAlgebra
using ITensors, ITensorMPS
using ProgressBars
using KrylovKit
using Plots
using LaTeXStrings
using ClassicalOrthogonalPolynomials
using Random

function initstate(N,max_per_site)
    state = zeros(N)
    j = 1
    while j < N
        r = rand(1:max_per_site)
        state[j] = r
        if j+r < N
            j += r
        else
            state[N] += N-sum(state)
            j += r
        end
    end
    return state
end

let 
    # Initializes N bosons sites
    N = 10
    sites = siteinds("Qudit", N, dim=N+1;conserve_number=true)
    U = 100
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
    TestState = ["1","1","1","1","1","1","1","1","1","1"]
    psi_test = MPS(sites, TestState)
    psi0 = random_mps(sites, TestState;linkdims=10)
    @show(expect(psi0, "n"; sites=1:N))

    nsweeps = 5 # number of sweeps : 5
    maxdim = [10,20,100,100,200] # bonds dimension
    cutoff = [1E-10] # truncation error

    #energy,psi = dmrg(H,psi0;nsweeps,maxdim,cutoff)
    #average = expect(psi, "n"; sites=1:N)
    #@show(average)
    "for i=1:N
        for j=1:N"
            #average = expect(psi, "Adag",i,"A",j; sites=j)
            "@show(average)
        end
    end"
    return
end
