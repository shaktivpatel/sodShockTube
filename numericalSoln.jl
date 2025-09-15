using Plots


#sod shock tube numerical solution using MacCormack (2 step Lax-Wendroff Scheme) scheme
function shockTube(; n=201, tMax=0.04, CFL=0.95)
    #set up computational domain
    xmin = -.5
    xmax = .5
    dx = (xmax-xmin)/(n-1)
    x = xmin:dx:xmax

    t = 0.0
    tStar = .75*tMax

    #create empty variables and specificy values

    rNum = zeros(n) #density
    uNum = zeros(n) #velocity
    eNum = zeros(n) #specific energy
    pNum = zeros(n) #pressure
    cNum = zeros(n) #local speed of sound
    U = zeros(n,3) #U vector [r, ru, re]
    F = zeros(n,3) #flux vector [ru, ru^2 + p, rue + pu]
    Ucalc = zeros(n,3) #updated U vector
    Ubar = zeros(n,3) #Ubar vector

    uStar = zeros(n) 
    pStar = zeros(n) 
    eStar = zeros(n) 
    Fbar = zeros(n,3) 
    Ubar2 = zeros(n,3) 


    gamma = 1.4

    #set initial conditions
    for i = 1:n
        if x[i] <= 0
            rNum[i] = 1
            pNum[i] = 1
            uNum[i] = 0
        elseif x[i] > 0
            rNum[i] = .125
            pNum[i] = .1
            uNum[i] = 0
        end

        eNum[i] = pNum[i]/(gamma - 1) + 0.5*rNum[i]*(uNum[i]^2)
        cNum[i] = sqrt(gamma*pNum[i]/rNum[i])
    end

    Ucalc[:,1] = rNum
    Ucalc[:,2] = rNum.*uNum
    Ucalc[:,3] = eNum

    while(t < tStar)

        dt = CFL*dx./(maximum(cNum.+(uNum.*uNum)))

        for i = 1:n

            #calculate flux

            U[i,1] = rNum[i]
            U[i,2] = rNum[i].*uNum[i]
            U[i,3] = eNum[i]

            F[i,1] = rNum[i].*uNum[i]
            F[i,2] = rNum[i].*uNum[i].^2 + pNum[i]
            F[i,3] = uNum[i].*(eNum[i] + pNum[i])
        end

        for i = 2:(n-1)

            #predictor
            Ubar[i,1] = U[i,1] - (dt/dx)*(F[i+1,1] - F[i,1])
            Ubar[i,2] = U[i,2] - (dt/dx)*(F[i+1,2] - F[i,2])
            Ubar[i,3] = U[i,3] - (dt/dx)*(F[i+1,3] - F[i,3])

            uStar[i] = Ubar[i,2]/Ubar[i,1]
            pStar[i] = (gamma-1)*(Ubar[i,3]-(.5*Ubar[i,2]^2))/Ubar[i,1]
            eStar[i] = (Ubar[i,3]-Ubar[i,1]*(.5*uStar[i]^2))/Ubar[i,1]

            Fbar[i,1] = Ubar[i,1] .* uStar[i]
            Fbar[i,2] = Ubar[i,1] .* uStar[i].^2 + pStar[i]
            Fbar[i,3] = uStar[i] .* Ubar[i,3] + pStar[i]

            #corrector
            Ubar2[i,1] = U[i,1] - (dt/dx)*(Fbar[i,1] - Fbar[i-1,1])
            Ubar2[i,2] = U[i,2] - (dt/dx)*(Fbar[i,2] - Fbar[i-1,2])
            Ubar2[i,3] = U[i,3] - (dt/dx)*(Fbar[i,3] - Fbar[i-1,3])

            #update
            Ucalc[i,1] = .5 * (Ubar[i,1] + Ubar2[i,1])
            Ucalc[i,2] = .5 * (Ubar[i,2] + Ubar2[i,2])
            Ucalc[i,3] = .5 * (Ubar[i,3] + Ubar2[i,3])
            Ucalc[1,:] = Ucalc[2,:]
            Ucalc[n,:] = Ucalc[n-1,:]


        end

        Fbar[1,:] = Fbar[2,:]
        Fbar[end,:] = Fbar[end-1,:]

        rNum = Ucalc[:,1]
        uNum = Ucalc[:,2]./rNum
        eNum = Ucalc[:,3]
        pNum = (gamma-1) * (eNum-0.5.*rNum.*(uNum.*uNum))
        pNum = max.(pNum, 1e-12) 
        cNum = sqrt.(gamma.*pNum./rNum)

        t += dt

        println(t)
        end

    return x, rNum, uNum, eNum, pNum
end

#call function
x, r, u, e, p = shockTube()

#print(x)
#plot properties vs. position

scatter(x, r, label="Density", ms=2, ma=0.5)
scatter!(x, u, label="Velocity", ms=2, ma=0.5)
scatter!(x, p, label="Pressure", ms=2, ma=0.5)
