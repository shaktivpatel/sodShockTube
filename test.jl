using Plots

function shockTube(; n=201, tMax=0.285, CFL=0.95)
    # --- Setup computational domain ---
    xmin, xmax = -0.5, 0.5
    dx = (xmax - xmin) / (n - 1)
    x = range(xmin, xmax; length=n)

    gamma = 1.4

    # --- Pre-allocate arrays ---
    rNum  = zeros(n)
    uNum  = zeros(n)
    eNum  = zeros(n)
    pNum  = zeros(n)
    cNum  = zeros(n)
    
    U      = zeros(n,3)
    F      = zeros(n,3)
    Ucalc  = zeros(n,3)
    Ubar   = zeros(n,3)
    Fbar   = zeros(n,3)
    Ubar2  = zeros(n,3)

    uStar  = zeros(n)
    pStar  = zeros(n)
    eStar  = zeros(n)
    
    # --- Initial conditions ---
    @inbounds for i in 1:n
        if x[i] <= 0
            rNum[i] = 1.0
            pNum[i] = 1.0
        else
            rNum[i] = 0.125
            pNum[i] = 0.1
        end
        uNum[i] = 0.0
        eNum[i] = pNum[i]/(gamma-1) + 0.5*rNum[i]*uNum[i]^2
        cNum[i] = sqrt(gamma*pNum[i]/rNum[i])
    end

    Ucalc[:,1] .= rNum
    Ucalc[:,2] .= rNum .* uNum
    Ucalc[:,3] .= eNum

    t = 0.0
    tStar = 0.75 * tMax

    # --- Time loop ---
    while t < tStar
        # CFL time step
        dt = CFL*dx / maximum(cNum .+ abs.(uNum))

        # --- Compute flux vector ---
        @inbounds begin
            U[:,1] .= rNum
            U[:,2] .= rNum .* uNum
            U[:,3] .= eNum

            F[:,1] .= U[:,2]
            F[:,2] .= U[:,2].*uNum .+ pNum
            F[:,3] .= uNum .* U[:,3] .+ pNum .* uNum
        end

        # --- Predictor step ---
        @inbounds @simd for i in 2:n-1
            Ubar[i,:] .= U[i,:] .- (dt/dx) .* (F[i+1,:] .- F[i,:])
        end

        # Compute primitive variables at predicted step
        @inbounds @simd for i in 2:n-1
            uStar[i] = Ubar[i,2]/Ubar[i,1]
            pStar[i] = (gamma-1)*(Ubar[i,3] - 0.5*Ubar[i,2]^2/Ubar[i,1])
            eStar[i] = (Ubar[i,3] - 0.5*Ubar[i,1]*uStar[i]^2)/Ubar[i,1]

            Fbar[i,1] = Ubar[i,1]*uStar[i]
            Fbar[i,2] = Ubar[i,1]*uStar[i]^2 + pStar[i]
            Fbar[i,3] = uStar[i]*Ubar[i,3] + pStar[i]*uStar[i]
        end

        # --- Corrector step ---
        @inbounds @simd for i in 2:n-1
            Ubar2[i,:] .= U[i,:] .- (dt/dx) .* (Fbar[i,:] .- Fbar[i-1,:])
        end

        # --- Update conserved variables ---
        @inbounds begin
            Ucalc[2:end-1,:] .= 0.5 .* (Ubar[2:end-1,:] .+ Ubar2[2:end-1,:])
            # Boundary conditions
            Ucalc[1,:] .= Ucalc[2,:]
            Ucalc[end,:] .= Ucalc[end-1,:]
        end

        # --- Update primitive variables ---
        rNum .= Ucalc[:,1]
        uNum .= Ucalc[:,2] ./ rNum
        eNum .= Ucalc[:,3]
        @. pNum = max((gamma-1)*(eNum - 0.5*rNum*uNum^2), 1e-12)
        @. cNum = sqrt(gamma*pNum/rNum)

        t += dt
    end

    return x, rNum, uNum, eNum, pNum
end

# --- Call function and plot ---
x, r, u, e, p = shockTube()

plot(x, r, label="Density")
plot!(x, u, label="Velocity")
plot!(x, p, label="Pressure")
