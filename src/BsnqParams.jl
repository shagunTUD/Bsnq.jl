module Bsnq2D

using Revise
using DrWatson
@quickactivate "Bsnq.jl"

using Parameters
using Gridap
using Printf
using LineSearches: BackTracking
using CSV
using TickTock
using Bsnq.Constants
using Bsnq.WaveMaking

export case_setup, default_params

const B = BsnqC


function case_setup(params, ha)    

  @unpack domX, domY, dx, dy = params    
  @unpack simT, simΔt, outΔt = params 
  @unpack probname = params
  @unpack T, H, h0 = params
  @unpack absx = params


  ## Generate Cartesian Domain 2DH 
  domain = (domX[1], domX[2], domY[1], domY[2])
  partition = ( Int((domX[2]-domX[1])/dx), 
    Int((domY[2]-domY[1])/dy))
  model = CartesianDiscreteModel(domain, partition)


  ## Label sides of the domain
  labels = get_face_labeling(model)
  add_tag_from_tags!(labels, "Inlet", ["tag_7","tag_1","tag_3"])
  add_tag_from_tags!(labels, "Outlet", ["tag_8","tag_2","tag_4"])
  add_tag_from_tags!(labels, "SideWall", ["tag_5","tag_6"])
  writevtk(model, probname*"_model")


  ## Domains
  Ω = Interior(model)
  Γ = Boundary(model)
  nΓ = get_normal_vector(Γ)


  ## Sponge domain
  # ---------------------Start---------------------
  absc = VectorValue( absx[3], absx[4] )
  absn = VectorValue( absx[1], absx[2] )  
  absCoef = 30.0 / absx[6] / (exp(1.0) - 1.0)
  function mapSpongeLayer(xs)
    n = length(xs)
    x = 1/n*sum(xs)    

    (x - absc)⋅absn > 0.0
  end
  
  xΩ = get_cell_coordinates(Ω)
  Ωd_to_Ω = lazy_map(mapSpongeLayer, xΩ)    
  Ωd = Triangulation( Ω, Ωd_to_Ω )
  writevtk(Ωd,probname*"_model"*"_Od")
  
  function spngFnc(x)
    dx = ( (x - absc)⋅absn ) / absx[5]
    dx = min( max( 0.0, dx ), 1.0 )

    return absCoef * ( exp(dx*dx) - 1.0 )
  end
  # ----------------------End----------------------



  ## Inlet boundary conditions
  wv = WaveMaker(T, H, h0; 
    theory = Fourier3(), ramp = (0, 2*T))
  println("Wave-making theory is : ", wv.theory)  
  inletη(t::Real) = x -> inletη(x, t::Real)
  inletη(x, t::Real) = WaveInletη(wv, x, t)  
  inletP(t::Real) = x -> inletP(x, t::Real)
  inletP(x, t::Real) = WaveInletP(wv, x, t)  


  ## Define Test Fnc
  # ---------------------Start---------------------
  ordη = 1
  reffeη = ReferenceFE(lagrangian, Float64, ordη)

  Ψw = TestFESpace(Ω, reffeη, 
    conformity=:H1)

  Ψη = TestFESpace(Ω, reffeη, 
    conformity=:H1, dirichlet_tags=["Inlet"])
  
  ordP = 2
  reffeP = ReferenceFE(lagrangian, VectorValue{2,Float64}, ordP)
  
  ΨP = TestFESpace(Ω, reffeP, 
    conformity=:H1, 
    dirichlet_tags=["Inlet", "Outlet","SideWall"],
    dirichlet_masks=[(true,false), (true,false), (false,true)])
  # ----------------------End----------------------


  ## Define Trial Space
  # ---------------------Start---------------------    
  Uw = TrialFESpace(Ψw)    
  
  Uη = TransientTrialFESpace(Ψη, inletη)    
  
  g3b(x, t::Real) = VectorValue(0.0, 0.0)
  g3b(t::Real) = x -> g3b(x,t)
  g3c(x, t::Real) = VectorValue(0.0, 0.0)
  g3c(t::Real) = x -> g3c(x,t)
  UP = TransientTrialFESpace(ΨP, [inletP,g3b,g3c])

  Y = MultiFieldFESpace([Ψw, Ψη, ΨP])
  X = TransientMultiFieldFESpace([Uw, Uη, UP])
  #NOTE: Uw is not transient
  # ----------------------End----------------------

  
  ## Water depth        
  h = interpolate_everywhere(ha, FESpace(Ω,reffeη))
  

  ## Define measures
  dΩ = Measure(Ω, 2*ordP)
  dΩd = Measure(Ωd, 2*ordP)
  dΓ = Measure(Γ, 2*ordP)

  ## Weak form    
  # ---------------------Start---------------------

  # Intermediate fncs
  # Check log_bsnq_v001_math.ipynb for the split of conv term
  conv2(p, td) = 1.0/td*(∇⋅p) - 1.0/(td*td)*(p⋅∇(td))


  # Residual
  res(t, (w, η, p), (ψw, ψη, ψp)) =
    ∫( w*ψw + ∇(ψw) ⋅ (h*∇(η)) )dΩ +
    ∫( -ψw * h * ∇(η)⋅nΓ )dΓ +
    ∫( ∂t(η)*ψη + (∇⋅p)*ψη )dΩ + 
    ∫( ∂t(p)⋅ψp + 
      ψp ⋅ (∇(p)'⋅p) / (h+η) + #convec-part1
      ψp ⋅ (conv2(p,(h+η)) * p) + #convec-part2
      g * (h+η) * (ψp⋅∇(η)) )dΩ + 
    ∫( (B + 1/3) * ( (∇⋅ψp) * h*h * (∇⋅∂t(p)) ) )dΩ + 
    ∫( -(B + 1/3) * ( ψp * h*h * (∇⋅∂t(p)) ) ⋅ nΓ )dΓ + 
    ∫( (2*B + 1/2) * ( ψp ⋅ ( h * (∇⋅∂t(p)) * ∇(h) ) ) )dΩ + 
    ∫( (-1/6) * ( ψp ⋅ (∇(∂t(p))⋅∇(h)) ) )dΩ +
    ∫( -B*g * ψp ⋅ (h*h*∇(w)) )dΩ + 
    ∫( ψη * η * spngFnc + (ψp ⋅ p) * spngFnc )dΩd

        
  op_AD = TransientFEOperator(res, X, Y)    
  # ----------------------End----------------------


  ## Initial soln
  t0 = 0.0
  x0 = interpolate_everywhere(
    [0.0, 0.0, VectorValue(0.0,0.0)], X(t0))
  
  
  ## Solver setup
  # ---------------------Start---------------------
  # # Linear Solver
  # lin_solver = LUSolver()
  # ode_solver = ThetaMethod(lin_solver, simΔt, 0.5)

  # NL Solver    
  nls = NLSolver(show_trace=true, 
    method=:newton, linesearch=BackTracking(), 
    iterations=10, ftol = 1e-5)  
  ode_solver = ThetaMethod(nls, simΔt, 0.5)
  
  solnht = solve(ode_solver, op_AD, x0, t0, simT)    
  # ----------------------End----------------------

  
  createpvd(probname) do pvd
    wh, ηh, ph = x0
    tval = @sprintf("%5.3f",t0)                
    tprint = @sprintf("%08d",t0*1000)
    println("Time : $tval")
    pvd[t0] = createvtk(Ω,probname*"_$tprint"*".vtu",
      cellfields=["eta"=>ηh, "P"=>ph, "h"=>-h])
  end


  ## Execute
  outMod = floor(Int64,outΔt/simΔt);
  tick()
  createpvd(probname, append=true) do pvd    
    cnt=0
    for (solh, t) in solnht                            
      cnt = cnt+1
      wh, ηh, ph = solh
      tval = @sprintf("%5.3f",t)                
      tprint = @sprintf("%08d",t*1000)
      println("Time : $tval \t Counter : $cnt")                  

      println("-x-x-")
      tock()
      println()
      tick()

      if(cnt%outMod != 0) 
        continue
      end
      pvd[t] = createvtk(Ω,probname*"_$tprint"*".vtu",
        cellfields=["eta"=>ηh, "P"=>ph, "h"=>-h])
    end
  end    
  tock()
  
end



"""
Memb_params

Parameters for the VIV.jl module.
"""

@with_kw struct default_params

  domX = (0,50)  
  domY = (0,10)
  dx = 2.5
  dy = 2.5

  simT = 5
  simΔt = 0.2
  outΔt = 1

  # Wave parameters
  T = 2 #s
  H = 0.05 # wave height
  h0 = 1.0 #water-depth  

  # Sponge layer
  # (nx, ny, cx, cy, len, T)
  absx = (1.0, 0.0, 30.0, 0.0, 20, 2)

  probname = datadir("sims_202401","gwce_run","gwce")

end

end 