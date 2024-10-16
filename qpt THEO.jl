using ITensors
using Plots
using Random
using ITensorMPS
using LaTeXStrings

print("START OF SIMULATION \n")
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
    result = map(string,state)
    return result
end

function SPDM(sites, psi)
    sing_density=zeros(N,N)
    for j in 1:N
        for k in 1:N
            temp=deepcopy(psi) # deepcopy of psi to avoid accidentaly messing with the files
            temp=apply(op("A",sites,k),temp) # applies the annihilation operator on site k
            temp=apply(op("Adag",sites,j),temp) # applies the creation operator on site j 
            sing_density[j,k]= inner(psi,temp) # computes the inner product of psi and the single particle density operator
        end
    end
    Partic_Numb = sum([sing_density[j, j] for j in 1:N])
    return sing_density, Partic_Numb
end

function Plot_3D(N,DATA)
    col_grad = cgrad([:orange, :blue], [0.1, 0.3, 0.8])
    Plots.surface(1:N,1:N,DATA,xlabel="i",ylabel="j",zlabel="Proba",color=col_grad)
end

function HITMAP(N,DATA)
    xs = 1:N
    ys = 1:N
    plt = Plots.heatmap(xs,ys,DATA)
    display(plt)
end

function Run_Simulation(N, U, J)
    # Initializes N bosons sites
    sites = siteinds("Qudit", N, dim=N+1;conserve_number=true, conserve_qns = true)
    
    # Variables needed for the dmrg algorithm
    nsweeps = 50 # number of sweeps : 5
    maxdim = [10,20,100,100,200] # bonds dimension
    cutoff = [1E-10] # truncation error

    # Creates the Bose-Hubbard model Hamiltonian
    os = OpSum()
    for j=1:N-1
        os += -J,"A",j,"Adag",j+1
        os += - J,"Adag",j,"A",j+1
    end
    for j=1:N 
        os += U/2,"n",j,"n",j 
        os += -U/2,"n",j
    end
    H = MPO(os,sites)

    # Intialises the random state from a given distribution of states
    Init_State = ["2", "0", "1", "1", "1", "1", "1", "1", "1", "1"]
    psi0 = random_mps(sites, Init_State;linkdims=10)

    # Executes the DRMG algorithm
    energy,psi = dmrg(H,psi0;nsweeps,maxdim,cutoff)

    # Gets the single particle density matrix from the DMRG results
    Single_Particle_Density, Particle_Number = SPDM(sites, psi)
end
    
let 
    plot(10, Run_Simulation(10, 7, 1))
end
