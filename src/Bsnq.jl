module Bsnq

using Revise 
using DrWatson
@quickactivate "Bsnq.jl"

include(srcdir("subroutines","Constants.jl"))
include(srcdir("subroutines","WaveMaking.jl"))

export Constants, WaveMaking

function main()
  return 0
end

end