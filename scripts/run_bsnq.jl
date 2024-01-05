using Revise
using DrWatson
@quickactivate "Bsnq"

include(srcdir("BsnqParams.jl"))

using Parameters
using TimerOutputs
#using .Bsnq2D

const to = TimerOutput()

resDir = datadir("sims_202401","bsnq_run","bsnq")

params = Bsnq2D.default_params( 
  probname = resDir,
  simT = 0.8
)

ha(x) = params.h0

@timeit to "Setup Run" Bsnq2D.case_setup(params, ha)

println("Warmup Done")
println()

params = Bsnq2D.default_params( 

  domX = (0,36.576),
  domY = (0,6.096),
  dx = 0.2032,
  dy = 0.2032,

  simT = 40,
  simΔt = 0.08,
  outΔt = 0.8,

  # Wave parameters
  T = 2, #s
  H = 0.015, # wave height
  h0 = 0.4572, #water-depth

  # Sponge layer
  # (nx, ny, cx, cy, len, T)
  absx = (1.0, 0.0, 30.0, 0.0, 6.576, 2.0),
  
  probname = resDir
)

function haWhalin(x)
  G = sqrt( x[2] * (6.096-x[2]) )

  if( x[1] ≤ (10.67-G) )
    return 0.4572
  elseif( x[1] ≤ (18.29-G) )
    return 0.4572 + 1/25*(10.67-G-x[1])
  else
    return 0.1524
  end
end

@timeit to "Total Run" Bsnq2D.case_setup(params, haWhalin)

show(to)

