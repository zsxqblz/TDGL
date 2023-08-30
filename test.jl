include("dependencies.jl")
include("sim.jl")
include("exp.jl")
include("pp.jl")

# run one simulation
let 
    run(`clear`)
    dx = 20
    dy = 20
    Gamma0 = 1
    T0 = 5
    time = 10
    dt = 0.1
    trunc = 10
    
    Gamma = Gamma0*dt
    T = T0*dt
    numt = convert(Int64,time/dt)

    latticeHist = zeros(Float64,numt,dx,dy)

    for it = 1:numt-1
        updateStep(it,dx,dy,latticeHist,Gamma,T)
    end

    vorticeHist = calcVortices(numt,dx,dy,latticeHist)

    corr = findCorrelation(numt,dx,dy,latticeHist,trunc)
    vortice_corr = findCorrelation(numt,dx,dy,vorticeHist,trunc)

    x_l = collect(1:dx)
    y_l = collect(1:dy)
    t_l = collect(range(dt,time,step=dt))
    save3DData(t_l,x_l,y_l,latticeHist,"data/230829/test_1")
    save3DData(t_l,x_l,y_l,vorticeHist,"data/230829/vortice_1")
    saveComplex3DData(t_l[1:numt-trunc],x_l,y_l,corr,"data/230829/corr_1")
    saveComplex3DData(t_l[1:numt-trunc],x_l,y_l,vortice_corr,"data/230829/vortice_corr_1")
end

# run one exp
let 
    run(`clear`)
    dx = 20
    dy = 20
    Gamma0 = 1
    T0 = 5
    time = 10
    dt = 0.1
    trunc = 10
    nsim = 1

    numt = convert(Int64,time/dt)

    lattice_corr, vortice_corr = expOnce(dx,dy,Gamma0,T0,time,dt,trunc,nsim,true,true)

    x_l = collect(1:dx)
    y_l = collect(1:dy)
    t_l = collect(range(dt,time,step=dt))
    saveComplex3DData(t_l[1:numt-trunc],x_l,y_l,lattice_corr,"data/230829/lattice_corr_1")
    saveComplex3DData(t_l[1:numt-trunc],x_l,y_l,vortice_corr,"data/230829/vortice_corr_1")
end