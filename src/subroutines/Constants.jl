module Constants

using Revise
using Roots: find_zero
using Random
using Distributions

export g, BsnqC 
export dispersionRel, dispersionRelAng
export randomPhase
export simpsonInteg1D, trapzInteg1D, gaussQuad1D

g = 9.81 #accel due to gravity
BsnqC = 1 / 15


function randomPhase(ω; seed = 1234)

  # Seeding is needed to generate the same 
  # array of random numbers every time
  Random.seed!(seed) 
  ϕ = rand(Uniform(-π,π), length(ω))
  return ϕ
end


#----------- Dispersion Rel ----------#
function dispersionRel(h::Real, T::Real; msg::Bool=true)

  f(L) = L-(g/2/π*T*T*tanh(2*π/L*h))
  L = find_zero(f, g/2/π*T*T)

  if(msg)
    println("Time Period \t T \t ",T)
    println("Water Depth \t h \t ",h)
    println("Wave Length \t L \t ",L)
    println("Wave Celerity \t C \t ",L/T)
    println("Disp. Regime \t h/L \t ",round.(h/L; digits=3))
  end
  
  return L
end


function dispersionRelAng(h::Real, ω::Real; msg::Bool=true)

  f(k) = ω^2-(g*k*tanh(k*h))
  k = find_zero(f, ω^2/g)

  if(msg)
    println("Angular Freq \t ω \t ",ω)
    println("Water Depth \t h \t ",h)
    println("Wave number \t k \t ",k)
    println("Wave Celerity \t C \t ",ω/k)
    println("Disp. Regime \t h/L \t ",round.(h*k/2π; digits=3))
  end
  
  return k
end
#--------- End Dispersion Rel --------#


#------------ Integration ------------#
function simpsonInteg1D(y, dx)
  ind = 1:length(y)  
  w = ifelse.(iseven.(ind), 4, 2)  
  w[1] = 1
  w[end] = 1

  return sum(dx/3.0 * w .* y)
end


function trapzInteg1D(y, dx)
  ind = 1:length(y)  
  w = ifelse.(iseven.(ind), 2, 2)  
  w[1] = 1
  w[end] = 1

  return sum(dx/2.0 * w .* y)
end


function gaussQuad1D(y, dx)
  # Gauss-Legendre Quadrature 2 point
  q1dx = 0.5*(1 - 1/√(3))
  q2dx = 0.5*(1 + 1/√(3))

  ly = y[1:end-1]
  ry = y[2:end]

  yq1 = (ry - ly)*q1dx + ly
  yq2 = (ry - ly)*q2dx + ly

  return sum(yq1 + yq2)*dx/2
end
#---------- End Integration ----------#

end