function orbit_DCM(i::Float64, Ω::Float64, ω::Float64)::Matrix{Float64}
    return [
        cos(Ω) (-sin(Ω) * cos(i)) (sin(Ω)*sin(i));
        sin(Ω) (cos(Ω)*cos(i)) -(cos(Ω)*sin(i));
        0.0 sin(i) cos(i)
    ] * [
        cos(ω) -sin(ω) 0.0;
        sin(ω) cos(ω) 0.0;
        0.0 0.0 1.0
    ]
end


function Kepler_invsolve_ellipse!(E::Vector{Float64}, M::Vector{Float64}, ecc::Float64, n_samples::Integer; ε::Float64 = 10.0^-13, max_iter::Integer = 128)
    local n_active::Integer = 0
    for i ∈ 1:max_iter
        @tturbo for n ∈ 1:n_samples
            ΔEn = -((E[n] - ecc * sin(E[n])) - M[n]) / (1.0 - ecc * cos(E[n]))
            should_updaten = ε ≤ abs(ΔEn)
            E[n] += should_updaten * ΔEn
            n_active += should_updaten
        end
        if n_active == false
            break
        end
        n_active = 0
    end
    return nothing
end


function Kepler_invsolve_hyperbola!(H::Vector{Float64}, M::Vector{Float64}, ecc::Float64, n_samples::Integer; ε::Float64 = 10.0^-13, max_iter::Integer = 128)
    local n_active::Integer = 0
    for i ∈ 1:max_iter
        @tturbo for n ∈ 1:n_samples
            ΔHn = -((ecc * sinh(H[n]) - H[n]) - M[n]) / (ecc * cosh(H[n]) - 1.0)
            should_updaten = ε ≤ abs(ΔHn)
            H[n] += should_updaten * ΔHn
            n_active += should_updaten
        end
        if n_active == false
            break
        end
        n_active = 0
    end
    return nothing
end


function sample_ellipse!(
    r::Vector{Float64}, x::Vector{Float64}, y::Vector{Float64}, z::Vector{Float64}, x̂::Vector{Float64}, ŷ::Vector{Float64}, ẑ::Vector{Float64},
    v_x::Vector{Float64}, v_y::Vector{Float64}, v_z::Vector{Float64}, t_from_epoch::Vector{Float64},
    M_i::Float64, μ::Float64, a::Float64, ł::Float64, e::Float64, DCM::Matrix{Float64}, n_samples::Integer
)
    local M::Vector{Float64} = deepcopy(t_from_epoch)
    @tturbo for n ∈ 1:n_samples
        M[n] = muladd(M[n], √(μ / a^3), M_i)
    end
    local E::Vector{Float64} = deepcopy(M)
    Kepler_invsolve_ellipse!(E, M, e, n_samples)
    @tturbo for n ∈ 1:n_samples
        νn = 2.0atan(√((1.0 + e) / (1.0 - e)) * tan(E[n] / 2.0))
        r̂_xn = cos(νn)
        r̂_yn = sin(νn)
        r_magn = ł / (1.0 + e * r̂_xn)
        x̂n = (DCM[1, 1] * r̂_xn) + (DCM[1, 2] * r̂_yn)
        ŷn = (DCM[2, 1] * r̂_xn) + (DCM[2, 2] * r̂_yn)
        ẑn = (DCM[3, 1] * r̂_xn) + (DCM[3, 2] * r̂_yn)
        ϕn = (e * r̂_yn) / (1.0 + e * r̂_xn)
        v̂_xn = -r̂_yn * cos(ϕn) + r̂_xn * sin(ϕn)
        v̂_yn = r̂_yn * sin(ϕn) + r̂_xn * cos(ϕn)
        v_magn = √(μ * ((2.0 / r_magn) - (1.0 / a)))
        r[n] = r_magn
        x[n] = r_magn * x̂n
        y[n] = r_magn * ŷn
        z[n] = r_magn * ẑn
        x̂[n] = x̂n
        ŷ[n] = ŷn
        ẑ[n] = ẑn
        v_x[n] = v_magn * ((DCM[1, 1] * v̂_xn) + (DCM[1, 2] * v̂_yn))
        v_y[n] = v_magn * ((DCM[2, 1] * v̂_xn) + (DCM[2, 2] * v̂_yn))
        v_z[n] = v_magn * ((DCM[3, 1] * v̂_xn) + (DCM[3, 2] * v̂_yn))
    end  
    return nothing
end


function sample_hyperbola!(
    r::Vector{Float64}, x::Vector{Float64}, y::Vector{Float64}, z::Vector{Float64}, x̂::Vector{Float64}, ŷ::Vector{Float64}, ẑ::Vector{Float64},
    v_x::Vector{Float64}, v_y::Vector{Float64}, v_z::Vector{Float64}, t_from_epoch::Vector{Float64},
    M_i::Float64, μ::Float64, a::Float64, ł::Float64, e::Float64, DCM::Matrix{Float64}, n_samples::Integer
)
    local M::Vector{Float64} = deepcopy(t_from_epoch)
    for n ∈ 1:n_samples
        M[n] = muladd(M[n], √(μ / -a^3), M_i)
    end
    local E::Vector{Float64} = deepcopy(M)
    Kepler_invsolve_hyperbola!(E, M, e, n_samples)
    @tturbo for n ∈ 1:n_samples
        νn = 2.0atan(√((e + 1.0) / (e - 1.0)) * tanh(E[n] / 2.0))
        r̂_xn = cos(νn)
        r̂_yn = sin(νn)
        r_magn = ł / (1.0 + e * r̂_xn)
        x̂n = (DCM[1, 1] * r̂_xn) + (DCM[1, 2] * r̂_yn)
        ŷn = (DCM[2, 1] * r̂_xn) + (DCM[2, 2] * r̂_yn)
        ẑn = (DCM[3, 1] * r̂_xn) + (DCM[3, 2] * r̂_yn)
        ϕn = (e * r̂_yn) / (1.0 + e * r̂_xn)
        v̂_xn = -r̂_yn * cos(ϕn) + r̂_xn * sin(ϕn)
        v̂_yn = r̂_yn * sin(ϕn) + r̂_xn * cos(ϕn)
        v_magn = √(μ * ((2.0 / r_magn) - (1.0 / a)))
        r[n] = r_magn
        x[n] = r_magn * x̂n
        y[n] = r_magn * ŷn
        z[n] = r_magn * ẑn
        x̂[n] = x̂n
        ŷ[n] = ŷn
        ẑ[n] = ẑn
        v_x[n] = v_magn * ((DCM[1, 1] * v̂_xn) + (DCM[1, 2] * v̂_yn))
        v_y[n] = v_magn * ((DCM[2, 1] * v̂_xn) + (DCM[2, 2] * v̂_yn))
        v_z[n] = v_magn * ((DCM[3, 1] * v̂_xn) + (DCM[3, 2] * v̂_yn))
    end  
    return nothing
end


function sample_parabola!(
    r::Vector{Float64}, x::Vector{Float64}, y::Vector{Float64}, z::Vector{Float64}, x̂::Vector{Float64}, ŷ::Vector{Float64}, ẑ::Vector{Float64},
    v_x::Vector{Float64}, v_y::Vector{Float64}, v_z::Vector{Float64}, t_from_epoch::Vector{Float64},
    μ::Float64, a::Float64, ł::Float64, DCM::Matrix{Float64}, n_samples::Integer
)
    local M::Vector{Float64} = deepcopy(t_from_epoch)
    for n ∈ 1:n_samples
        M[n] *= √(μ / (2.0 * a^3))
    end
    @tturbo for n ∈ 1:n_samples
        νn = 2.0atan(2.0sinh(asinh((3.0 / 2.0) * M[n]) / 3.0))
        r̂_xn = cos(νn)
        r̂_yn = sin(νn)
        r_magn = ł / (1.0 + r̂_xn)
        x̂n = (DCM[1, 1] * r̂_xn) + (DCM[1, 2] * r̂_yn)
        ŷn = (DCM[2, 1] * r̂_xn) + (DCM[2, 2] * r̂_yn)
        ẑn = (DCM[3, 1] * r̂_xn) + (DCM[3, 2] * r̂_yn)
        ϕn = (e * r̂_yn) / (1.0 + e * r̂_xn)
        v̂_xn = -r̂_yn * cos(ϕn) + r̂_xn * sin(ϕn)
        v̂_yn = r̂_yn * sin(ϕn) + r̂_xn * cos(ϕn)
        v_magn = √(μ * ((2.0 / r_magn) - (1.0 / a)))
        r[n] = r_magn
        x[n] = r_magn * x̂n
        y[n] = r_magn * ŷn
        z[n] = r_magn * ẑn
        x̂[n] = x̂n
        ŷ[n] = ŷn
        ẑ[n] = ẑn
        v_x[n] = v_magn * ((DCM[1, 1] * v̂_xn) + (DCM[1, 2] * v̂_yn))
        v_y[n] = v_magn * ((DCM[2, 1] * v̂_xn) + (DCM[2, 2] * v̂_yn))
        v_z[n] = v_magn * ((DCM[3, 1] * v̂_xn) + (DCM[3, 2] * v̂_yn))
    end  
    return nothing
end


function sample_conic!(
    r::Vector{Float64}, x::Vector{Float64}, y::Vector{Float64}, z::Vector{Float64}, x̂::Vector{Float64}, ŷ::Vector{Float64}, ẑ::Vector{Float64},
    v_x::Vector{Float64}, v_y::Vector{Float64}, v_z::Vector{Float64}, t_from_epoch::Vector{Float64},
    M_i::Float64, μ::Float64, a::Float64, e::Float64, DCM::Matrix{Float64}, n_samples::Integer
)
    local ł::Float64
    if e == 1
        ł = 2.0a
        sample_parabola!(r, x, y, z, x̂, ŷ, ẑ, v_x, v_y, v_z, t_from_epoch, μ, a, ł, DCM, n_samples) # no M_i needed since we just take periapsis as M_i
    else
        ł = a * (1.0 - e^2)
        if e > 1
            sample_hyperbola!(r, x, y, z, x̂, ŷ, ẑ, v_x, v_y, v_z, t_from_epoch, M_i, μ, a, ł, e, DCM, n_samples)
        else
            sample_ellipse!(r, x, y, z, x̂, ŷ, ẑ, v_x, v_y, v_z, t_from_epoch, M_i, μ, a, ł, e, DCM, n_samples)
        end
    end
    return nothing
end


