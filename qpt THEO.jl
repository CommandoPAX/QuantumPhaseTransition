using ITensors
using Plots
using Random
using ITensorMPS
using LaTeXStrings
using ForwardDiff
using EasyFit

#https://www.overleaf.com/2663516136vbjstqbfdvgk#c9d482

function true_rng(N,max_per_site)
    state = [rand(0:max_per_site) for j in 1:N] #creation of random boson patches
    tot = sum(state)
    while tot > N #Suppression of bosons if we have too many
        j = rand(1:N)
        if state[j] != 0
            tot = sum(state)
            sub = rand(1:state[j]) #subtract a random value
            if tot-sub < N #check if said value doesn't bring the total under the N total of boson required
                sub = tot-N
            end
            state[j] -= sub
        end
    end
    while tot < N #Addition of bosons if we have not enough
        j = rand(1:N)
        if state[j] < max_per_site
            add = rand(state[j]:max_per_site)
            tot = sum(state)
            if tot+add > N #check if said value doesn't bring the total above the N total of boson required
                add = N-tot
            end
            state[j] += add
        end
    end
    result = map(string,state)
    return result
end

function SPDM(sites, psi, N)
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

function isLinear(x,y,thresh)
    fit = fitlinear(x,y)
    if(fit.R < thresh)
        return false
    else 
        return true
    end
end

function Plot_3D(N,DATA)
    col_grad = cgrad([:orange, :blue], [0.1, 0.3, 0.8])
    Plots.surface(1:N,1:N,DATA,xlabel="i",ylabel="j",zlabel="Proba",color=col_grad)
end

#Plot the density for the site j in log-log scale
function Plot_one_site_density(single_density,j)
    one_site_density = log.(single_density[:,j])
    N=length(one_site_density)
    plt=Plots.plot(1:N,one_site_density,xlabel="site #",ylabel="one-site density",title="One site density for site " * string(j),
    legend=false, linewidth=2,linecolor=[:black])
    display(plt)
end 

function Hitmap(N,DATA)
    xs = 1:N
    ys = 1:N
    plt = Plots.heatmap(xs,ys,DATA)
    display(plt)
end

function Run_Simulation(N, U, J)

    # Initializes N bosons sites
    print("Initializing...\n")
    sites = siteinds("Qudit", N, dim=N+1;conserve_number=true, conserve_qns = true)
    
    # Variables needed for the dmrg algorithm
    nsweeps = 20 # number of sweeps
    maxdim = [10,20,100,100,200] # bonds dimension
    cutoff = [1E-10] # truncation error

    # Creates the Bose-Hubbard model Hamiltonian
    print("Computing Hamiltonian...\n")
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
    print("Computing initial state...\n")
    Init_State = true_rng(N, 5)
    psi0 = random_mps(sites, Init_State;linkdims=10)

    # Executes the DRMG algorithm
    print("Applying DRMG...\n")
    energy,psi = dmrg(H,psi0;nsweeps,maxdim,cutoff)

    # Gets the single particle density matrix from the DMRG results
    print("Getting single particle densities...\n")
    Single_Particle_Density, Particle_Number = SPDM(sites, psi, N)

    # Creates the plots and display them (with the particle number at the end to verify)
    print("Final particle number : ", Particle_Number)
    #Plot_3D(N, Single_Particle_Density)
    #Hitmap(N, Single_Particle_Density)
    Plot_one_site_density(Single_Particle_Density, Int(N/2))

    #Check for the phase we're in 
    return isLinear(log.(Single_Particle_Density[Int(N/2):N,Int(N/2)]) , [1.0*i for i in Int(N/2):N],0.995)

end
    
let 
    Run_Simulation(28, 2, 1)
end
