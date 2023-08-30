function expOnce(dx,dy,Gamma0,T0,time,dt,trunc,nsim,PBC=true)
    Gamma = Gamma0*dt
    T = T0*dt
    numt = convert(Int64,time/dt)

    lattice_corr = zeros(ComplexF64,numt-trunc,dx,dy)
    vortice_corr = zeros(ComplexF64,numt-trunc,dx,dy)

    for sim_i = 1:nsim
        latticeHist = zeros(Float64,numt,dx,dy,PBC)

        for it = 1:numt-1
            updateStep(it,dx,dy,latticeHist,Gamma,T,PBC)
        end

        vorticeHist = calcVortices(numt,dx,dy,latticeHist,PBC)

        lattice_corr = lattice_corr + findCorrelation(numt,dx,dy,latticeHist,trunc,PBC)
        vortice_corr = vortice_corr + findCorrelation(numt,dx,dy,vorticeHist,trunc,PBC)
    end
    lattice_corr = lattice_corr/nsim
    vortice_corr = vortice_corr/nsim

    return lattice_corr, vortice_corr
end

function expOnce(dx,dy,Gamma0,T0,time,dt,trunc,nsim,PBC=true,showProg=true)
    Gamma = Gamma0*dt
    T = T0*dt
    numt = convert(Int64,time/dt)

    lattice_corr = zeros(ComplexF64,numt-trunc,dx,dy)
    vortice_corr = zeros(ComplexF64,numt-trunc,dx,dy)

    @showprogress for sim_i = 1:nsim
        latticeHist = zeros(Float64,numt,dx,dy,PBC)

        for it = 1:numt-1
            updateStep(it,dx,dy,latticeHist,Gamma,T,PBC)
        end

        vorticeHist = calcVortices(numt,dx,dy,latticeHist,PBC)

        lattice_corr = lattice_corr + findCorrelation(numt,dx,dy,latticeHist,trunc,PBC)
        vortice_corr = vortice_corr + findCorrelation(numt,dx,dy,vorticeHist,trunc,PBC)
    end
    lattice_corr = lattice_corr/nsim
    vortice_corr = vortice_corr/nsim

    return lattice_corr, vortice_corr
end


function scanT0(T0_start,T0_end,T0_step,dx,dy,Gamma0,time,dt,trunc,nsim,PBC=true)
end