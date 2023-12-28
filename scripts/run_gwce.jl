using DrWatson
@quickactivate "Bsnq"

include(srcdir("gwce.jl"))

using Revise
using Parameters
using TimerOutputs
#using .GWCE

const to = TimerOutput()

resDir = datadir("sims_202401","gwce_run","gwce")

params = GWCE.GWCE_params( 
  probname = resDir 
)

@timeit to "Total Run" GWCE.gwce_setup(params)

show(to)

