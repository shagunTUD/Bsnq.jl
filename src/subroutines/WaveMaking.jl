module WaveMaking

using Revise
using DrWatson
@quickactivate "Bsnq.jl"

include(srcdir("subroutines","Constants.jl"))

using Gridap
using Roots: find_zero
using .Constants

export WaveInletη, WaveInletP
export WaveMaker, Airy, Fourier3


abstract type WaveTheory end




"""
Common
======

"""
# ---------------------Start---------------------
struct WaveMaker  

  theory::WaveTheory
  # Options are Airy(), Fourier3()

  T::Float64
  H::Float64
  h0::Float64
  ϕ0::Float64  
  θDeg::Float64  
  x0::NTuple{2,Float64}

  L::Float64
  linL::Float64
  ω::Float64
  k::Float64  
  kh::Float64
  θRad::Float64
  csθ::Float64
  snθ::Float64
  kx::Float64
  ky::Float64
  c::Float64
  
end

function WaveMaker( T::Real, H::Real, h0::Real, 
  ϕ0::Real = 0.0, θDeg::Real = 0.0, x0 = (0.0,0.0); 
  theory = Airy() )    

  WaveMaker( theory, T, H, h0, 
    ϕ0, θDeg, x0 )
end  


function getPhase(wv::WaveMaker, x, t::Real)
  ϕ = 
    wv.kx * (x[1] - wv.x0[1]) +
    wv.ky * (x[2] - wv.x0[2]) +
    - wv.ω * t + 
    wv.ϕ0

  return ϕ
end


function WaveInletη( wv::WaveMaker, x, t::Real)
  return WaveInletη( wv, wv.theory, x, t )
end


function WaveInletP( wv::WaveMaker, x, t::Real)
  return WaveInletP( wv, wv.theory, x, t )  
end
# ----------------------End----------------------




"""
Airy (Linear) WaveMaking
========================

"""
# ---------------------Start---------------------
struct Airy <: WaveTheory  end


function WaveMaker( theory::Airy, T, H, h0, 
  ϕ0, θDeg, x0 )
  
  θRad = θDeg*π/180.0
  csθ = cos(θRad)
  snθ = sin(θRad)
  ω = 2*π/T
  k = dispersionRelAng(h0, ω, msg=true)
  kh = k*h0
  L = 2*π/k
  linL = L
  kx = k * cos(θRad)
  ky = k * sin(θRad)
  c = L/T

  WaveMaker( theory, T, H, h0, ϕ0, θDeg, x0,
    L, linL, ω, k, kh, θRad, csθ, snθ, kx, ky, c)
end


function WaveInletη( wv::WaveMaker, theory::Airy, x, t::Real)
  ϕ = getPhase(wv, x, t)    

  return wv.H/2.0*cos(ϕ)  
end


function WaveInletP( wv::WaveMaker, theory::Airy, x, t::Real)
  ϕ = getPhase(wv, x, t)  
  
  # Integrating from -h to 0
  pn = wv.H/2.0 * wv.ω / wv.k * cos(ϕ)

  # Using Wheeler stretching
  η = WaveInletη( wv, theory, x, t )
  pn = pn * (wv.h0 + η)/wv.h0

  p = pn * wv.csθ
  q = pn * wv.snθ

  return VectorValue(p, q)
end

# ----------------------End----------------------




"""
Fourier 3 WaveMaking
====================

Reference
  Madsen, P. A., & Sørensen, O. R. (1993). 
  Bound waves and triad interactions in shallow water. 
  Ocean Engineering, 20(4), 359–388. 
  https://doi.org/10.1016/0029-8018(93)90002-Y
"""
# ---------------------Start---------------------
struct Fourier3 <: WaveTheory 
  G2::Float64
  G3::Float64

  α1::Float64
  α2::Float64
  α3::Float64

  Fourier3(G2::Real, G3::Real, α1::Real, α2::Real, α3::Real) = 
    new(G2, G3, α1, α2, α3)
  Fourier3() = new(0.0, 0.0, 0.0, 0.0, 0.0)
end


function WaveMaker( theory::Fourier3, T, H, h0, 
  ϕ0, θDeg, x0 )
  
  θRad = θDeg*π/180.0
  csθ = cos(θRad)
  snθ = sin(θRad)
  ω = 2*π/T

  linL = 2*pi / dispersionRelAng(h0, ω, msg=false)

  # Solve for kh^2
  σ3(y) = 9.0/16.0/y * ( 1 + 2 * (BsnqC + 1/9) * y ) /
    ( 1 + (2*BsnqC + 1/3) * y )
  f(y) = sqrt( g/h0 * ( y * (1 + BsnqC*y) ) / ( 1 + (BsnqC + 1/3)*y ) ) *
    ( 1 + H*H/4.0/h0/h0 * σ3(y) ) - ω
  y = find_zero(f, (2*π/linL*h0)^2 )
  
  kh = sqrt(y)
  k = kh/h0
  L = 2*π/k  
  kx = k * cos(θRad)
  ky = k * sin(θRad)
  c = L/T

  khsq = kh*kh
  coefG2 = 3.0/4.0/khsq * ( 1.0 + (BsnqC+1/9)*khsq )

  coefG3 = 27.0/64.0/khsq/khsq * ( 1.0 + 2*(BsnqC+1/9)*khsq )  

  # Constants in eta expression
  α1 = H/2.0 
  α2 = H^2 / 4.0 / h0 * coefG2
  α3 = H^3 / 8.0 / (h0^2) * coefG3

  if(true)
    println("Time Period \t T \t ",T)
    println("Wave Height \t H \t ",H)
    println("Water Depth \t h0 \t ",h0)
    println("Wave Length \t L \t ",L)
    println("Wave Len-Lin \t linL \t ",linL)
    println("Wave Celerity \t C \t ",c)
    println("Disp. Regime \t h/L \t ",round.(h0/L; digits=3))
  end

  WaveMaker( Fourier3(coefG2, coefG3, α1, α2, α3), 
    T, H, h0, ϕ0, θDeg, x0,
    L, linL, ω, k, kh, θRad, csθ, snθ, kx, ky, c)
end


function WaveInletη( wv::WaveMaker, theory::Fourier3, x, t::Real)
  ϕ = getPhase(wv, x, t)  
  
  return wv.theory.α1*cos(ϕ) + wv.theory.α2*cos(2.0*ϕ) + 
    wv.theory.α3*cos(3.0*ϕ)
end


function WaveInletP( wv::WaveMaker, theory::Fourier3, x, t::Real)
  
  η = WaveInletη( wv, theory, x, t )
  pn = η * wv.c

  p = pn * wv.csθ
  q = pn * wv.snθ

  return VectorValue(p, q)
end

# ----------------------End----------------------


end