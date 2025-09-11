#sod shock tube numerical solution using MacCormack (2 step Lax-Wendroff Scheme) scheme
function shocktube()
    #set up computational domain
    n = 201
    xmin = -.5
    xmax = .5
    dx = (xmax-xmin)/(n-1)
    x = xmin:dx:xmax

    t = 0.0
    tMax = .285
    tStar = .75*tMax

    #create empty variables and specificy values

    rNum = zeros(n) #density
    uNum = zeros(n) #velocity
    eNum = zeros(n) #specific energy
    pNum = zeros(n) #pressure
    cNum = zeros(n) #local speed of sound
    U = zeros(n,3) #U vector
    F = zeros(n,3) #flux vector
    Ucalc = zeros(n,3) #Ubar vector

    gamma = 1.4
    CFL = .95

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

        eNum[i] = pNum[i]/(gamma - 1) + .5*rNum[i]*uNum[i]^2
        cNum[i] = sqrt(gamma*pNum[i]/rNum[i])
    end

    Ucalc[:,1] = rNum
    Ucalc[:,2] = rNum.*uNum
    Ucalc[:,3] = eNum

    while(t < tStar)

        dt = CFL*dx/(maximum(cNum.+abs.(uNum)))

        for i = 1:n

            #calculate flux

            U[i,1] = rNum[i]
            U[i,2] = rNum[i].*uNum[i]
            U[i,3] = eNum[i]

            F[i,1] = rNum[i].*uNum[i]
            F[i,2] = rNum[i].*uNum[i].^2 + pNum[i]
            F[i,3] = uNum[i].*(eNum[i] + pNum[i])
        end

        t += dt
    end

end