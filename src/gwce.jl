using DrWatson
@quickactivate "Bsnq.jl"

using Gridap
using Printf
using LineSearches: BackTracking
using CSV
using Tables
using Revise
using WaveSpec.Constants


function gwce_setup(domX, domY, dx, dy, ha, inletη, inletP,
    simT, simΔt, outΔt, probesxy, probname)


    # Generate Cartesian Domain 2DH 
    domain = (domX[1], domX[2], domY[1], domY[2])
    partition = ( Int((domX[2]-domX[1])/dx), 
        Int((domY[2]-domY[1])/dy))
    model = CartesianDiscreteModel(domain, partition)


    # Label sides of the domain
    labels = get_face_labeling(model)
    add_tag_from_tags!(labels, "inlet", ["tag_7","tag_1","tag_3"])
    add_tag_from_tags!(labels, "outlet", ["tag_8","tag_2","tag_4"])
    add_tag_from_tags!(labels, "sideWall", ["tag_5","tag_6"])
    writevtk(model, probname*"_model")


    # Define Test Fnc
    # ---------------------Start---------------------
    order = 2
    reffeη = ReferenceFE(lagrangian, Float64, order)
    Ψ1 = TestFESpace(model, reffeη, 
        conformity=:H1, dirichlet_tags=["inlet"])

    reffeJ = ReferenceFE(lagrangian, VectorValue{2,Float64}, order)
    Ψ2 = TestFESpace(model, reffeJ, 
        conformity=:H1)

    reffeP = ReferenceFE(lagrangian, VectorValue{2,Float64}, order)
    Ψ3 = TestFESpace(model, reffeP, 
        conformity=:H1, 
        dirichlet_tags=["inlet", "outlet","sideWall"],
        dirichlet_masks=[(true,false), (true,false), (false,true)])
    # ----------------------End----------------------


    # Define Trial Space
    # ---------------------Start---------------------    
    Eta = TransientTrialFESpace(Ψ1, inletη)

    J = TrialFESpace(Ψ2) #NOTE: Check if J should be transient
    
    g3b(x, t::Real) = VectorValue(0.0, 0.0)
    g3b(t::Real) = x -> g3b(x,t)
    g3c(x, t::Real) = VectorValue(0.0, 0.0)
    g3c(t::Real) = x -> g3c(x,t)
    P = TransientTrialFESpace(Ψ3, [inletP,g3b,g3c])

    Y = MultiFieldFESpace([Ψ1, Ψ2, Ψ3])
    X = TransientMultiFieldFESpace([Eta, J, P])
    #NOTE: J is not transient
    # ----------------------End----------------------

    
    # Water depth        
    h0 = interpolate_everywhere(ha, FESpace(model,reffeη))


    # Tau field 
    τₘ = 0.0
    function fτ(x) 
        lh0 = h0(x)
        if(lh0 > 200.0)
            lτ = 0.005          
        elseif(lh0 < 1.0)
            lτ = 1.0
        else
            lτ = 1.0/lh0
        end              
        return lτ*τₘ
    end
    τ = interpolate_everywhere(fτ, FESpace(model,reffeη))
    

    # Define integration space
    Ω = Triangulation(model)
    dΩ = Measure(Ω, 2*order)


    # Intermediate Functions
    sc2a(p,h) = 1.0/h*(∇⋅p) - 1.0/(h*h)*(p⋅∇(h))


    # Weak form    
    # ---------------------Start---------------------
    res(t, (η, j, p), (ψ1, ψ2, ψ3)) =
        ∫( ∂tt(η)*ψ1 )dΩ + 
        ∫( ∂t(η)*τ*ψ1 + g*(h0+η)*∇(ψ1)⋅∇(η) 
            - ψ1*∇(τ)⋅p - ∇(ψ1)⋅j )dΩ + 
        ∫( ψ2⋅j + 0.5*g*(ψ2⋅∇(η*η)) - (ψ2⋅p)*τ )dΩ +
        ∫( ψ2⋅(∇(p)'⋅p)/(h0+η) )dΩ +
        ∫( ψ2⋅(sc2a(p,(h0+η))*p) )dΩ +
        ∫( ∂t(p)⋅ψ3 - ψ3⋅j + (ψ3⋅p)*τ + g*(h0+η)*(ψ3⋅∇(η)) )dΩ    
         
    op_AD = TransientFEOperator(res, X, Y; order=2)    
    # ----------------------End----------------------


    # Initial soln
    t0 = 0.0
    x0 = interpolate_everywhere(
        [0.0, VectorValue(0.0,0.0), VectorValue(0.0,0.0)], X(t0))
    x0t = interpolate_everywhere(
        [0.0, VectorValue(0.0,0.0), VectorValue(0.0,0.0)], X(t0))
    x0tt = interpolate_everywhere(
        [0.0, VectorValue(0.0,0.0), VectorValue(0.0,0.0)], X(t0))
    
    
    # NL Solver    
    nls = NLSolver(show_trace=true, 
        method=:newton, linesearch=BackTracking(), iterations=10)
    #ode_solver = Newmark(nls, simΔt, 0.5, 0.25)        
    ode_solver = GeneralizedAlpha(nls, simΔt, 0.0)    
    solnht = solve(ode_solver, op_AD, (x0,x0t,x0tt), t0, simT)    

    
    # Probe setup
    numP = length(probesxy)    
    probeDa = zeros(Float64, 1, numP*5+1)    
    lDa = zeros(Float64, 1, numP*5+1)
    probeDa[1,1] = t0;    


    # for (uh_tn, tn) in solnht        
    #     println(tn)        
    # end

    
    createpvd(probname) do pvd
        ηh, jh, ph = x0
        tval = @sprintf("%5.3f",t0)                
        println("Time : $tval")
        tval = @sprintf("%d",floor(Int64,t0*1000))                
        pvd[t0] = createvtk(Ω,probname*"_$tval"*".vtu",
            cellfields=["eta"=>ηh, "P"=>ph, "J"=>jh, "h0"=>-h0])
    end


    # Execute
    outMod = floor(Int64,outΔt/simΔt);
    createpvd(probname, append=true) do pvd    
        cnt=0
        for (solh, t) in solnht                       
            cnt = cnt+1
            ηh, jh, ph = solh
            tval = @sprintf("%5.3f",t)                
            println("Time : $tval")
            tval = @sprintf("%d",floor(Int64,t*1000))                

            lDa[1] = t
            for ip in 1:numP
                ipxy = probesxy[ip];
                ik = (ip-1)*5 + 1
                lDa[ik+1] = ipxy[1]
                lDa[ik+2] = ipxy[2]
                lDa[ik+3] = ηh(ipxy)
                ipph = ph(ipxy)
                lDa[ik+4] = ipph[1]
                lDa[ik+5] = ipph[2]
            end
            probeDa = vcat(probeDa, lDa);

            if(cnt%outMod != 0) 
                continue
            end
            pvd[t] = createvtk(Ω,probname*"_$tval"*".vtu",
                cellfields=["eta"=>ηh, "P"=>ph, "J"=>jh, "h0"=>-h0])
        end
    end    

    CSV.write(probname*"_probes.csv", Tables.table(probeDa))
end


#
# Run simulation
g = 9.81
T = 10.0; #s
h = 1.0; #m
H = 0.05; #m
L = dispersionRel(h, T)
k = 2*π/L;
ω = 2*π/T;
probname = datadir("sims_202401","gwce_run","gwce")

probesxy = [Point(12.0, 5.0)
            Point(24.0, 5.0)
            Point(36.0, 5.0)
            Point(48.0, 5.0)];

# Water Depth
ha(x) = h
# function ha(x)
#     if(x[1]<50.0)
#         rha = h
#     elseif(x[1]>100.0)
#         rha = h/2.0
#     else
#         rha = h - h/2.0*(x[1]-50.0)/50.0
#     end

#     return rha
# end

g1(x, t::Real) = H/2.0*sin(-ω*t)
g1(t::Real) = x -> g1(x,t)
g2a(x, t::Real) = VectorValue(H/2.0*sin(-ω*t)*ω/k, 0.0)
g2a(t::Real) = x -> g2a(x,t)
gwce_setup((0,300), (0,10), 1.25, 1.25, ha, g1, g2a,
    80, 0.2, 2, probesxy, probname)