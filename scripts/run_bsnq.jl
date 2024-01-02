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

# params = Bsnq2D.default_params( 
#   probname = resDir,
#   simT = 20 
# )

params = Bsnq2D.default_params( 

  domX = (0,36.576),
  domY = (0,6.096),
  dx = 0.2032,
  dy = 0.2032,

  simT = 8,
  simΔt = 0.08,
  outΔt = 0.8,

  # Wave parameters
  h0 = 0.4572, #water-depth
  H = 0.015, # wave height
  ω = 2*π / 2, #T = 2s
  k = dispersionRelAng(0.4572, 2*π / 2; msg=true),
  probname = resDir
)

@timeit to "Total Run" Bsnq2D.case_setup(params)

show(to)

