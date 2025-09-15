using Plots

function shockTube(; n=201, tMax=0.285, CFL=0.5)
    # Domain
    xmin, xmax = -0.5, 0.5
    dx = (xmax - xmin)/(n-1)
    x = range(xmin, xmax, length=n)

    # Initialize variables
    gamma = 1.4
    t = 0.0
    tStar = 0.75 * tMax

    rNum = zeros(Float64, n)   # density
    uNum = zeros(Float64, n)   # velocity
    eNum = zeros(Float64, n)   # total energy
    pNum = zeros(Float64, n)   # pressure
    cNum = zeros(Float64, n)   # speed of sound

    # Conserved variables
    U = zeros(Float64, n, 3)      # [rho, rho*u, e]
    Ucalc = zeros(Float64, n, 3)
    Ubar = zeros(Float64, n, 3)
    Ubar2 = zeros(Float64, n, 3)
    F = zeros(Float64, n, 3)
    Fbar = zeros(Float64, n, 3)

    # Auxiliary arrays
    uStar = zeros(Float64, n)
    pStar = zeros(Float64, n)

    # Initial conditions
    for i in 1:n
        if x[i] <= 0
            rNum[i] = 1.0
            pNum[i] = 1.0
            uNum[i] = 0.0
        else
            rNum[i] = 0.125
            pNum[i] = 0.1
            uNum[i] = 0.0
        end
        eNum[i] = pNum[i]/(gamma - 1) + 0.5*rNum[i]*uNum[i]^2
        cNum[i] = sqrt(gamma * pNum[i] / rNum[i])
    end

    # Initialize conserved variables
    Ucalc[:,1] .= rNum
    Ucalc[:,2] .= rNum .* uNum
    Ucalc[:,3] .= eNum

    # Main time loop
    while t < tStar
        # CFL condition
        dt = CFL * dx / maximum(abs.(uNum) .+ cNum)

        # Compute flux
        U .= Ucalc  # copy current state
        F[:,1] .= U[:,2]                   # rho*u
        F[:,2] .= U[:,2].^2 ./ U[:,1] .+ (gamma-1)*(U[:,3] - 0.5*U[:,2].^2 ./ U[:,1])  # rho*u^2 + p
        F[:,3] .= (U[:,3] .+ (gamma-1)*(U[:,3] - 0.5*U[:,2].^2 ./ U[:,1])) .* (U[:,2]./U[:,1])  # energy flux

        # Predictor step
        @inbounds for i in 2:n-1
            Ubar[i,:] .= U[i,:] .- dt/dx .* (F[i+1,:] .- F[i,:])
        end

        # Compute primitive variables for predictor
        uStar[2:n-1] .= Ubar[2:n-1,2] ./ Ubar[2:n-1,1]
        pStar[2:n-1] .= (gamma-1) .* (Ubar[2:n-1,3] .- 0.5 .* Ubar[2:n-1,2].^2 ./ Ubar[2:n-1,1])
        pStar .= max.(pStar, 1e-12)  # pressure floor

        # Predictor flux
        Fbar[:,1] .= Ubar[:,1] .* uStar
        Fbar[:,2] .= Ubar[:,1] .* uStar.^2 .+ pStar
        Fbar[:,3] .= uStar .* Ubar[:,3] .+ uStar .* pStar

        # Corrector step
        @inbounds for i in 2:n-1
            Ubar2[i,:] .= U[i,:] .- dt/dx .* (Fbar[i,:] .- Fbar[i-1,:])
        end

        # Update
        Ucalc[2:n-1,:] .= 0.5 .* (Ubar[2:n-1,:] .+ Ubar2[2:n-1,:])

        # Apply simple boundary conditions (copy neighbor)
        Ucalc[1,:] .= Ucalc[2,:]
        Ucalc[end,:] .= Ucalc[end-1,:]

        # Update primitive variables
        rNum .= Ucalc[:,1]
        uNum .= Ucalc[:,2] ./ rNum
        eNum .= Ucalc[:,3]
        pNum .= (gamma-1) .* (eNum .- 0.5 .* rNum .* uNum.^2)
        pNum .= max.(pNum, 1e-12)   # pressure floor
        cNum .= sqrt.(gamma .* pNum ./ rNum)

        # Advance time
        t += dt
        print(t)
    end

    return x, rNum, uNum, pNum
end

# Example usage
x, rho, u, p = shockTube()

plot(x, rho, label="Density", lw=2)
plot!(x, u, label="Velocity", lw=2)
plot!(x, p, label="Pressure", lw=2)
