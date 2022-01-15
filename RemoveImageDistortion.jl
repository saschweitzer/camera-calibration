###############################################################################
"""
Provides utilities for removal of radial distortion, by predicting the true u,v
coordinates of a pixel from the distorted coordinates, and the intrinsic
calibration parameters of the camera.

Also provides utilities to simulate radial distoriton.
"""
module RemoveImageDistortion
###############################################################################

###############################################################################
import ExtrinsicCameraCalibration: IntrinsicCalibrationParameters
# using Statistics
###############################################################################

###############################################################################
export RemoveRadialDistortion, SimulateRadialDistortion
###############################################################################


###############################################################################
function RemoveRadialDistortion(u_d::Array, v_d::Array, ι::IntrinsicCalibrationParameters)
    return RemoveRadialDistortion(u_d::Array, v_d::Array, ι.radial_distortion_poly_degrees,
                                  ι.radial_distortion_poly_coeffs, ι.u₀, ι.v₀)
end
function RemoveRadialDistortion(u_d::Array, v_d::Array, Degrees, Coeffs, u₀, v₀)
    # Convert the *observed* pixel to polar coordinates:
    ρ_d = abs.((u_d .- u₀) .+ (v_d .- v₀)im)
    Θ_d = angle.((u_d .- u₀) .+ (v_d .- v₀)im)
    # Calculate the Δρ correction needed at each point:
    Δρ = Polyvalc(ρ_d, Coeffs, Degrees)
    Δu = cos.(Θ_d) .* Δρ
    Δv = sin.(Θ_d) .* Δρ
    # Apply the correction to ρ and convert back to pixel coordinates:
    u_p = u_d - Δu
    v_p = v_d - Δv
    return u_p,v_p
end
###############################################################################


###############################################################################
function SimulateRadialDistortion(uv_p::Array, ι::IntrinsicCalibrationParameters)
    return SimulateRadialDistortion(u_p, v_p, ι)
end
function SimulateRadialDistortion(u_p::Array, v_p::Array, ι::IntrinsicCalibrationParameters)
    return SimulateRadialDistortion(u_p, v_p,
                                    ι.radial_distortion_poly_degrees,
                                    ι.radial_distortion_poly_coeffs,
                                    ι.u₀, ι.v₀)
end
function SimulateRadialDistortion(u_p::Array, v_p::Array, Degrees, Coeffs, u₀, v₀)
    ρ_p = abs.((u_p .- u₀) .+ (v_p .- v₀)im)
    Θ_p = angle.((u_p .- u₀) .+ (v_p .- v₀)im)
    # Calculate the Δρ correction needed at each point:
    Δρ = Polyvalc(ρ_p, Coeffs, Degrees)
    Δu = cos.(Θ_p) .* Δρ
    Δv = sin.(Θ_p) .* Δρ
    # Apply the correction to ρ and convert back to pixel coordinates:
    u_d = u_p + Δu
    v_d = v_p + Δv
    return u_d,v_d
end
###############################################################################


###############################################################################
# Polyfitc is similar to polyfit, but allows choosing powers (5,3,1, etc.)
# (based on Matlab file exchange Polyfitc calculations):
function Polyfitc(X,Y,Degrees)
    lsqM = ones(length(X), length(Degrees))
    for n = 1:length(Degrees)
        lsqM[:, n] = X.^Degrees[n]
    end
    Coeffs = (lsqM \ Y)'
    return Coeffs[:]
end
# ----------
function Polyvalc(X, Coeffs, Degrees)
    val(x,Coeffs,Degrees) = sum( [ Coeffs[i]*x^(Degrees[i]) for i in 1:length(Coeffs) ] )
    Y = [ val(x,Coeffs,Degrees) for x in X ]
end
###############################################################################

###############################################################################
###############################################################################
###############################################################################
end                             # end of module
###############################################################################
###############################################################################
###############################################################################

