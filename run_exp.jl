include("dependencies.jl")
include("sim.jl")
include("exp.jl")
include("pp.jl")

dx = parse(Int64,ARGS[1])
dy = parse(Int64,ARGS[2])
Gamma0 = parse(Float64,ARGS[3])
T0 = parse(Float64,ARGS[4])
time = parse(Float64,ARGS[5])
dt = parse(Float64,ARGS[6])
trunc = parse(Int64,ARGS[7])
nsim = parse(Int64,ARGS[8])
filename = ARGS[9]

numt = convert(Int64,time/dt)

lattice_corr, vortice_corr = expOnce(dx,dy,Gamma0,T0,time,dt,trunc,nsim,true,true)

x_l = collect(1:dx)
y_l = collect(1:dy)
t_l = collect(range(dt,time,step=dt))
saveComplex3DData(t_l[1:numt-trunc],x_l,y_l,lattice_corr,filename*"lattice_corr")
saveComplex3DData(t_l[1:numt-trunc],x_l,y_l,vortice_corr,filename*"vortice_corr")