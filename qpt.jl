using SparseArrays
using LinearAlgebra
using ITensors, ITensorMPS
using ProgressBars
using KrylovKit
using Plots
using LaTeXStrings
using ClassicalOrthogonalPolynomials
using Random

#https://www.overleaf.com/2663516136vbjstqbfdvgk#c9d482

function true_rng(N,max_per_site)
    state = [rand(0:max_per_site) for j in 1:N]
    tot = sum(state)
    while tot > N
        j = rand(1:N)
        if state[j] != 0
            sub = rand(1:state[j])
            if tot-sub < N
                sub += N-tot
            end
            state[j] -= sub
            tot -= sub
        end
    end
    tot = sum(state)
    while tot < N
        r = rand(1:N)
        if state[r] == 0
            state[r] = N-tot
        end
    end
    return state
end

function Googoogaga()
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

function plot(N,proba)
    surface(1:N,1:N,proba,xlabel="i",ylabel="j",zlabel="proba i to j")
end

N=100
plot(N,true_rng(N,5))
