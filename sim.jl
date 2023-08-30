function genLattice(dx,dy,randomPhase=false)
    if randomPhase
        return zeros(Float64,dx,dy)
    else
        return pi*rand(Float64, (dx, dy))
    end
end

function updateStep(it,dx,dy,latticeHist,Gamma,T,PBC=true)
    lattice = latticeHist[it,:,:]
    if PBC
        for y = 1:dy, x = 1:dx
            del = (
                sin(lattice[mod1(x+1,end),y]-lattice[mod1(x,end),y])
                -sin(lattice[mod1(x,end),y]-lattice[mod1(x-1,end),y])
                +sin(lattice[x,mod1(y,end)]-lattice[x,mod1(y,end)])
                -sin(lattice[x,mod1(y,end)]-lattice[x,mod1(y-1,end)])
            )
            latticeHist[it+1,x,y] = mod2pi(lattice[x,y]+Gamma*del+rand(Truncated(Normal(0,T),-3*T,3*T)))
        end
    else
        for xy = 1:dy
            for x = 1:dx
                #TODO
            end
        end
    end
end

function delxy(dx,dy,lattice,pbc=true)
    delx = zeros(Float64,dx,dy)
    dely = zeros(Float64,dx,dy)
    if pbc
        for y = 1:dy
            for x = 1:dx
                delx[x,y] = (
                    rem2pi(lattice[mod1(x+1,end),y]-lattice[mod1(x,end),y], RoundNearest)
                    +rem2pi(lattice[mod1(x,end),y]-lattice[mod1(x-1,end),y], RoundNearest)
                )
                dely[x,y] = (
                    rem2pi(lattice[x,mod1(y,end)]-lattice[x,mod1(y,end)], RoundNearest)
                    +rem2pi(lattice[x,mod1(y,end)]-lattice[x,mod1(y-1,end)], RoundNearest)
                )
            end
        end
    else
        for xy = 1:dy
            for x = 1:dx
                delx[x,y] = (
                    rem2pi(lattice[mod1(x+1,end),y]-lattice[mod1(x,end),y], RoundNearest)
                    +rem2pi(lattice[mod1(x,end),y]-lattice[mod1(x-1,end),y], RoundNearest)
                )
                dely[x,y] = (
                    rem2pi(lattice[x,mod1(y,end)]-lattice[x,mod1(y,end)], RoundNearest)
                    +rem2pi(lattice[x,mod1(y,end)]-lattice[x,mod1(y-1,end)], RoundNearest)
                )
            end
        end
    end
    return delx, dely
end