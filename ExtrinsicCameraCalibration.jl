###############################################################################
"""

# General

This module contains Julia code for performing extrinsic camera calibration, as
described in Schweitzer & Cowen, 2021 (references below). It follows the
methodology develocped by Holland et al. (1997).

# References
K. T. Holland, R. A. Holman, T. C. Lippmann, J. Stanley and N. Plant, "Practical
use of video imagery in nearshore oceanographic field studies," in IEEE Journal
of Oceanic Engineering, vol. 22, no. 1, pp. 81-92, Jan. 1997, doi:10.1109/48.557542. 

Schweitzer, S. A., & Cowen, E. A. (2021). "Instantaneous river-wide water surface
velocity field measurements at centimeter scales using infrared quantitative
image velocimetry". Water Resources Research, 57,
e2020WR029279. https://doi.org/10.1029/2020WR029279

# License
This software is released under the MIT license, however, the authors request to
be credited if you use the software, preferably by citing the above-mentioned
Schweitzer & Cowen (2021) paper.

MIT License

Copyright (c) 2022 Seth A. Schweitzer

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.

""" 
module ExtrinsicCameraCalibration
###############################################################################

###############################################################################
using Statistics
using Printf # This modules is only required to enable the optional text output of HollandEtAlExtrinsicCalibration()
###############################################################################

###############################################################################
export IntrinsicCalibrationParameters, ExtrinsicCalibrationParameters, # structs
    HollandEtAlExtrinsicCalibration, CalculateDLTcoefficients # functions
###############################################################################


###############################################################################
"""
    IntrinsicCalibrationParameters

Data type for intrinsic camera calibration parameters.
"""
struct IntrinsicCalibrationParameters <: AbstractFloat
    u₀::Float64                 # u-coordinate (horizontal) of the image center (pixels)
    v₀::Float64                 # v-coordinate (vertical) of the image center (pixels)
    λᵤ::Float64                 # horizontal magnification factor
    λᵥ::Float64                 # vertical magnification factor
    PixelPitch::Float64
    radial_distortion_poly_degrees::Array{Float64,1}
    radial_distortion_poly_coeffs::Array{Float64,1}
end
###############################################################################

###############################################################################
"""
    ExtrinsicCalibrationParameters

Data type for extrinsic calibration parameters.
"""
struct ExtrinsicCalibrationParameters <: AbstractFloat
    φ::Float64                  # azimuth. clockwise from the north (y-axis), radians.
    τ::Float64                  # tilt. angle away from the vertical of the camera centeral ray, radians.
    σ::Float64                  # roll. rotation of the camera around the centeral ray, radians.
    f::Float64                  # distance b/w FPA and optical center of camera, pixels. Note that this will often be a negative number (since the FPA is behind the focal point)
    x_c::Float64                # x coordinate of camera's optical center, meters
    y_c::Float64                # y coordinate of camera's optical center, meters
    z_c::Float64                # z coordinate of camera's optical center, meters
end
### Convenience constructor:
# Input as NamedTuple: 
ExtrinsicCalibrationParameters(x::NamedTuple) = ExtrinsicCalibrationParameters(x.φ, x.τ, x.σ, x.f, x.x_c, x.y_c, x.z_c)  
###############################################################################


##############################################################################
struct RotationMatrix <: AbstractFloat
    m11::Float64
    m12::Float64
    m13::Float64
    m21::Float64
    m22::Float64
    m23::Float64
    m31::Float64
    m32::Float64
    m33::Float64
end
""" Initialize a RotationMatrix from a set of ExtrinsicCalibrationParameters """
function RotationMatrix(χ::ExtrinsicCalibrationParameters) 
    return RotationMatrix(f_m11(χ), f_m12(χ), f_m13(χ),
                          f_m21(χ), f_m22(χ), f_m23(χ),
                          f_m31(χ), f_m32(χ), f_m33(χ))
end
##############################################################################

###############################################################################
"""
    CalculateDLTcoefficients(χ, ι)

Calculates the DLT coefficients from the camera calibration coefficients

Uses the equations in the appendix of Holland et al. 1997 to calculate the
Direct Linear Transformation coefficients directly from the physical parameters
of the camera's installation - location, orientation, magnification.
"""
function CalculateDLTcoefficients(χ::ExtrinsicCalibrationParameters, ι::IntrinsicCalibrationParameters)::NamedTuple
    Rot = RotationMatrix(χ)
    L = fill(NaN, 11)
    LL = -( χ.x_c * Rot.m31 + χ.y_c * Rot.m32 + χ.z_c * Rot.m33 )
    abs(LL) < 1.0 ? @warn("LL has a small value! - results might be unstable.") : nothing
    L₁ = ( ι.u₀ * Rot.m31 -χ.f * Rot.m11 ) / (ι.λᵤ * LL)
    L₂ = ( ι.u₀ * Rot.m32 -χ.f * Rot.m12 ) / (ι.λᵤ * LL)
    L₃ = ( ι.u₀ * Rot.m33 -χ.f * Rot.m13 ) / (ι.λᵤ * LL)
    L₄ = -(L₁ * χ.x_c + L₂ * χ.y_c + L₃ * χ.z_c)
    L₅ = ( ι.v₀ * Rot.m31 -χ.f * Rot.m21 ) / (ι.λᵥ * LL)
    L₆ = ( ι.v₀ * Rot.m32 -χ.f * Rot.m22 ) / (ι.λᵥ * LL)
    L₇ = ( ι.v₀ * Rot.m33 -χ.f * Rot.m23 ) / (ι.λᵥ * LL)
    L₈ = -(L₅ * χ.x_c + L₆ * χ.y_c + L₇ * χ.z_c)
    L₉ = Rot.m31 / LL
    L₁₀ = Rot.m32 / LL
    L₁₁ = Rot.m33 / LL
    return (; L₁, L₂, L₃, L₄, L₅, L₆, L₇, L₈, L₉, L₁₀, L₁₁)
end
###############################################################################


###############################################################################
"""
    HollandEtAlExtrinsicCalibration(GCPuv, GCPxyz, χ, ι)

Returns the extrinsic camera calibration coefficients, calculated from a set of
GCPs in real-world and pixel-space coordinates, using the methodology of Holland
et al. 1997.
"""
function HollandEtAlExtrinsicCalibration(GCPuv::Array{Float64,2}, GCPxyz::Array{Float64,2},
                                         χinitial::ExtrinsicCalibrationParameters,
                                         ι::IntrinsicCalibrationParameters;
                                         MaxIterations = 50,
                                         PrintoutStats=true)
    # in practical experience if the method doesn't converge after about 30
    # iterations it probably will not converge within 101 either.
    χ = deepcopy(χinitial)      # this prevents modification of the original input
    nGCPs = minimum([size(GCPuv,1), size(GCPxyz,1)])
    A = fill(NaN, nGCPs*2,7)
    b = fill(NaN, nGCPs*2)
    CorrectionTerm = fill(Inf, 7)
    PrevCorrectionTerm = copy(CorrectionTerm)
    iteration,ContinueFlag = (1,true)
    while (ContinueFlag==true)
        for gcp in 1:nGCPs
            x,y,z = GCPxyz[gcp,:]
            u,v = GCPuv[gcp,:]
            F₀, G₀, Matrix_row1, Matrix_row2 = UpdateConstants(x, y, z, u, v, χ, ι)
            A[gcp,:]       = Matrix_row1
            A[gcp+nGCPs,:] = Matrix_row2
            b[gcp]         = -F₀
            b[gcp+nGCPs]   = -G₀
        end
        CorrectionTerm = A\b
        χ = ExtrinsicCalibrationParameters(
            χ.φ + CorrectionTerm[1],
            χ.τ + CorrectionTerm[2],
            χ.σ + CorrectionTerm[3],
            χ.f + CorrectionTerm[4],
            χ.x_c + CorrectionTerm[5],
            χ.y_c + CorrectionTerm[6],
            χ.z_c + CorrectionTerm[7],
        )
        if (iteration > MaxIterations)
            @warn join(["exceeded maximum  iterations ($iteration)",
                        "max(|ε|)=$(map(y -> @sprintf("%.0e",y),maximum(abs.(CorrectionTerm))))",
                        "rms(ε₁-ε₀)=$(map(y -> @sprintf("%.0e",y),sqrt(mean(CorrectionTerm-PrevCorrectionTerm).^2)))"],
                       ", ")
            ContinueFlag = false
        end
        if ( (maximum(abs.(CorrectionTerm)) < 1e-10) || # convergenece conditions
             (sqrt(mean(CorrectionTerm-PrevCorrectionTerm).^2) <= 1e-10) )
            ContinueFlag = false
        end
        iteration <= 3 ? (ContinueFlag = true) : nothing # but makes sure there have been at least 3 iterations
        ContinueFlag == true ? PrevCorrectionTerm = copy(CorrectionTerm) : nothing # store this iteration's results unless we're done - in that case keep the last iteration's results for reporting.
        iteration += 1
    end
    if PrintoutStats == true
        println(join(["Finished extrinsic calibration:",
                      "$iteration iterations",
                      "max(|ε|)=$(map(y -> @sprintf("%.0e",y),maximum(abs.(CorrectionTerm))))",
                      "rms(ε₁-ε₀)=$(map(y -> @sprintf("%.0e",y),sqrt(mean(CorrectionTerm-PrevCorrectionTerm).^2)))"],
                     ", "))
    end
    return χ
end
###############################################################################

##############################################################################
""" helper function for HollandEtAlExtrinsicCalibration() """
function UpdateConstants(x, y, z, u, v, χ, ι)
    # ---------- distances from the GCP to the camera:
    Δx  = f_Δx(x,χ)
    Δy  = f_Δy(y,χ)
    Δz  = f_Δz(z,χ)
    # ---------- undistorted, scale corrected image coordinates of the GCP:
    uˢ  = f_uˢ(u,ι)
    vˢ  = f_vˢ(v,ι)
    # ---------- elements of the 3x3 orthogonal rotation matrix (aka direction cosines):
    m11  = f_m11(χ)
    m12  = f_m12(χ)
    m13  = f_m13(χ)
    # --   
    m21  = f_m21(χ)
    m22  = f_m22(χ)
    m23  = f_m23(χ)
    # --   
    m31  = f_m31(χ)
    m32  = f_m32(χ)
    m33  = f_m33(χ)
    # ---------- results using rotation matrix:
    o = m11 * Δx + m12 * Δy + m13 * Δz
    p = m21 * Δx + m22 * Δy + m23 * Δz
    q = m31 * Δx + m32 * Δy + m33 * Δz 
    # ---------- bᵢⱼ are a set of constants used in LHS of the iterative calibration calculation:
    # Note that as calculated here they are multiplied by q.
    # b₁ⱼ are in the u coordinate direction
    qb11 = uˢ * (Δx*m32 -Δy*m31 ) + χ.f  * (Δx * m12 -Δy * m11 )
    qb12 = uˢ * (Δx * cos(χ.τ) * sin(χ.φ) + Δy * cos(χ.φ) * cos(χ.τ) + Δz * sin(χ.τ)) + χ.f  * (-Δx * sin(χ.φ) * sin(χ.τ) * sin(χ.σ) - Δy * cos(χ.φ) * sin(χ.τ) * sin(χ.σ) + Δz * cos(χ.τ) * sin(χ.σ))
    qb13 = χ.f * (Δx * m21 + Δy * m22 + Δz * m23 )
    qb14 = o                   
    qb15 = uˢ * m31 + χ.f * m11  
    qb16 = uˢ * m32 + χ.f * m12 
    qb17 = uˢ * m33 + χ.f * m13   
    # ----------
    # b₂ⱼ are in the v coordinate direction
    qb21 = vˢ * (Δx * m32 -Δy * m31 ) + χ.f  * (Δx * m22 -Δy * m21 ) 
    qb22 = vˢ * (Δx * cos(χ.τ) * sin(χ.φ) + Δy * cos(χ.φ) * cos(χ.τ) + Δz * sin(χ.τ)) + χ.f  * (-Δx * sin(χ.φ) * sin(χ.τ) * cos(χ.σ)-Δy * cos(χ.φ) * sin(χ.τ) * cos(χ.σ) + Δz * cos(χ.τ) * cos(χ.σ))
    qb23 = -χ.f * (Δx * m11 + Δy * m12 + Δz * m13 ) # 
    qb24 = p                   # 
    qb25 = vˢ * m31 + χ.f * m21
    qb26 = vˢ * m32 + χ.f * m22
    qb27 = vˢ * m33 + χ.f * m23
    ###############################################################################
    # F₀, G₀ form the RHS of the system of equations:
    F₀ = uˢ + o/q * χ.f 
    G₀ = vˢ + p/q * χ.f 
    ###############################################################################
    Matrix_row1 = [qb11, qb12, qb13, qb14, -qb15, -qb16, -qb17] ./ q
    Matrix_row2 = [qb21, qb22, qb23, qb24, -qb25, -qb26, -qb27] ./ q
    return F₀, G₀, Matrix_row1, Matrix_row2
end
###############################################################################



###############################################################################
# 
# Functions to calculate parameters defined in Holland et al. 1997:
# 
# ---------- distances from the GCP to the camera:
f_Δx(x,χ) = x - χ.x_c                  # [meters]
f_Δy(y,χ) = y - χ.y_c                  # [meters]
f_Δz(z,χ) = z - χ.z_c                  # [meters]
# ---------- undistorted, scale corrected image coordinates of the GCP:
f_uˢ(u,ι) = (u - ι.u₀)ι.λᵤ               # [pixels]
f_vˢ(v,ι) = (v - ι.v₀)ι.λᵥ               # [pixels]
# ---------- elements of the 3x3 orthogonal rotation matrix (aka direction cosines):
f_m11(χ) =  cos(χ.φ)cos(χ.σ) + sin(χ.φ)cos(χ.τ)sin(χ.σ)  # [dimensionless]
f_m12(χ) = -sin(χ.φ)cos(χ.σ) + cos(χ.φ)cos(χ.τ)sin(χ.σ) # [dimensionless]
f_m13(χ) =  sin(χ.τ)sin(χ.σ)                       # [dimensionless]
# --
f_m21(χ) = -cos(χ.φ)sin(χ.σ) + sin(χ.φ)cos(χ.τ)cos(χ.σ) # [dimensionless]
f_m22(χ) =  sin(χ.φ)sin(χ.σ) + cos(χ.φ)cos(χ.τ)cos(χ.σ)  # [dimensionless]
f_m23(χ) =  sin(χ.τ)cos(χ.σ)                       # [dimensionless]
# --
f_m31(χ) =  sin(χ.φ)sin(χ.τ)                       # [dimensionless]
f_m32(χ) =  cos(χ.φ)sin(χ.τ)                       # [dimensionless]
f_m33(χ) = -cos(χ.τ)                            # [dimensionless]
# ---------- scaling coefficients from direction cosines and distances:
f_o(x,y,z,χ) = f_m11(χ)*f_Δx(x,χ) + f_m12(χ)*f_Δy(y,χ) + f_m13(χ)*f_Δz(z,χ) # [meters]
f_p(x,y,z,χ) = f_m21(χ)*f_Δx(x,χ) + f_m22(χ)*f_Δy(y,χ) + f_m23(χ)*f_Δz(z,χ) # [meters]
f_q(x,y,z,χ) = f_m31(χ)*f_Δx(x,χ) + f_m32(χ)*f_Δy(y,χ) + f_m33(χ)*f_Δz(z,χ) # [meters]
###############################################################################


###############################################################################
###############################################################################
###############################################################################
end                             # end of module
###############################################################################
###############################################################################
###############################################################################

