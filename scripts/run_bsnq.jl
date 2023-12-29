using DrWatson
@quickactivate "Bsnq"

include(srcdir("bsnq_params.jl"))

using Revise
using Parameters
using TimerOutputs
#using .Bsnq2D

const to = TimerOutput()

resDir = datadir("sims_202401","bsnq_run","bsnq")

params = Bsnq2D.default_params( 
  probname = resDir 
)

@timeit to "Total Run" Bsnq2D.case_setup(params)

show(to)

