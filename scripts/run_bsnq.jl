using DrWatson
@quickactivate "Bsnq"

include(srcdir("bsnq_params.jl"))

using Revise
using Parameters
using TimerOutputs
using WaveSpec.Constants
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

  simT = 20,
  simΔt = 0.08,
  outΔt = 0.8,

  # Wave parameters
  h0 = 0.4572, #water-depth
  H = 0.015, # wave height
  ω = 2*π / 2, #T = 2s
  k = dispersionRelAng(0.4572, 2*π / 2; msg=true),
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

