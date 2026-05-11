function Householder_solve!(x::Vector{Float64}, λ::Vector{Float64}, T_true::Vector{Float64}, M::Integer, tilesize::Integer; ε_1::Float64 = 10.0^-10, ε_2::Float64 = 10.0^-13, max_iter::Integer = 128)
    local ε::Float64 = ε_1
    if M > 0
        ε = ε_2
    end
    n_active::Integer = 0
    for i ∈ 1:max_iter
        @tturbo for n ∈ 1:tilesize
            λn = λ[n]
            xn = x[n]
            yn = √(1.0 - (λn^2 * (1.0 - xn^2)))
            is_ellipsen = xn < 1.0
            ψn = is_ellipsen * acos((xn * yn) + (λn * (1.0 - xn^2))) +
                !(is_ellipsen) * acosh((xn * yn) - (λn * (xn^2 - 1.0)))
            denomn = (1.0 - xn^2)^-1
            Tn = denomn * (((ψn + π * M) / √(abs(1.0 - xn^2))) - (xn + λn * yn))
            T′n = denomn * ((3.0Tn * xn) - 2.0 + (2.0λn^3 * (xn / yn)))
            T′′n = denomn * (3.0Tn + (5.0xn * T′n) + (2.0(1.0 - λn^2) * (λn^3 / yn^3)))
            T′′′n = denomn * ((7.0xn * T′′n) + 8.0T′n - (6.0(1.0 - λn^2) * λn^5 * (xn / yn^5)))
            Tn -= T_true[n]
            Δxn = -Tn * (
                (T′n^2 - ((Tn * T′′n) / 2.0)) / 
                ((T′n * (T′n^2 - (Tn * T′′n))) + ((T′′′n * Tn^2) / 6.0))
            )
            should_updaten = ε ≤ abs(Δxn)
            x[n] += should_updaten * Δxn
            n_active += should_updaten
        end
        if n_active == false
            break
        end
        n_active = 0
    end
    return nothing    
end


function Lambert_solve!(Δv_1::Matrix{Float64}, Δv_2::Matrix{Float64},
    μ::Float64, M_ref_o::Float64, M_ref_d::Float64, a_o::Float64, a_d::Float64, e_o::Float64, e_d::Float64, DCM_o::Matrix{Float64}, DCM_d::Matrix{Float64},
    t_departure::Float64, Δt_departure::Float64, tof_min::Float64, Δt_flight::Float64, tof_samples::Integer; tilesize::Integer = 4096
)
    departure_times::Vector{Float64} = collect(range(t_departure, step = Δt_departure, length = tilesize))
    arrival_times::Vector{Float64} = deepcopy(departure_times)
    tof::Float64 = tof_min
    @tturbo for n ∈ 1:tilesize
        arrival_times[n] += tof
    end
    
    
    # Memory for origin orbit parameters:
    r_o::Vector{Float64} = zeros(Float64, tilesize)
    x_o::Vector{Float64} = zeros(Float64, tilesize)
    y_o::Vector{Float64} = zeros(Float64, tilesize)
    z_o::Vector{Float64} = zeros(Float64, tilesize)
    x̂_o::Vector{Float64} = zeros(Float64, tilesize)
    ŷ_o::Vector{Float64} = zeros(Float64, tilesize)
    ẑ_o::Vector{Float64} = zeros(Float64, tilesize)
    vx_o::Vector{Float64} = zeros(Float64, tilesize)
    vy_o::Vector{Float64} = zeros(Float64, tilesize)
    vz_o::Vector{Float64} = zeros(Float64, tilesize)
    # Initialize the origin orbit:
    sample_conic!(
        r_o, x_o, y_o, z_o, x̂_o, ŷ_o, ẑ_o,
        vx_o, vy_o, vz_o, departure_times,
        M_ref_o, μ, a_o, e_o, DCM_o, tilesize
    )
    

    # Memory for destination orbit parameters:
    r_d::Vector{Float64} = zeros(Float64, tilesize)
    x_d::Vector{Float64} = zeros(Float64, tilesize)
    y_d::Vector{Float64} = zeros(Float64, tilesize)
    z_d::Vector{Float64} = zeros(Float64, tilesize)
    x̂_d::Vector{Float64} = zeros(Float64, tilesize)
    ŷ_d::Vector{Float64} = zeros(Float64, tilesize)
    ẑ_d::Vector{Float64} = zeros(Float64, tilesize)
    vx_d::Vector{Float64} = zeros(Float64, tilesize)
    vy_d::Vector{Float64} = zeros(Float64, tilesize)
    vz_d::Vector{Float64} = zeros(Float64, tilesize)
    # Memory for Householder solver parameters:
    λ::Vector{Float64} = zeros(Float64, tilesize)
    T::Vector{Float64} = zeros(Float64, tilesize)
    x::Vector{Float64} = zeros(Float64, tilesize)
    γ::Vector{Float64} = zeros(Float64, tilesize)
    ρ::Vector{Float64} = zeros(Float64, tilesize)
    σ::Vector{Float64} = zeros(Float64, tilesize)
    is_backwards::BitVector = falses(tilesize)


    for i ∈ 1:tof_samples
        sample_conic!(
            r_d, x_d, y_d, z_d, x̂_d, ŷ_d, ẑ_d,
            vx_d, vy_d, vz_d, arrival_times,
            M_ref_d, μ, a_d, e_d, DCM_d, tilesize
        )
        # Initialize the Householder solver:
        @tturbo for n ∈ 1:tilesize
            cn = √((x_d[n] - x_o[n])^2 + (y_d[n] - y_o[n])^2 + (z_d[n] - z_o[n])^2)
            sn = (cn + r_d[n] + r_o[n]) / 2.0
            γ[n] = √((μ * sn) / 2.0)
            ρn = (r_o[n] - r_d[n]) / cn
            ρ[n] = ρn
            σ[n] = √(1.0 - ρn^2)
            is_backwardsn = ((x_o[n] * y_d[n]) - (y_o[n] * x_d[n])) < 0.0 
            is_backwards[n] = is_backwardsn
            λn = √(1.0 - (cn / sn)) * (1.0 * !(is_backwardsn) - 1.0 * is_backwardsn)
            Tn = √(2.0μ / sn^3) * arrival_times[n]
            T_00n = acos(λn) + (λn * √(1.0 - λn^2))
            T_1n = (2.0 / 3.0) * (1.0 - λn^3)
            bool_1n = Tn ≥ T_00n
            bool_2n = Tn < T_1n
            λ[n] = λn
            T[n] = Tn
            x[n] = bool_1n * ((T_00n / Tn)^(2 / 3) - 1.0) + 
                !(bool_1n) * (
                    bool_2n * (((5.0 / 2.0) * (T_1n * (T_1n - Tn)) / (Tn * (1.0 - λn^5))) + 1.0)
                    + !(bool_2n) * (((T_00n / Tn)^(log2(T_1n / T_00n))) - 1.0)
                )
        end
        Householder_solve!(x, λ, T, 0, tilesize)
        @tturbo for n ∈ 1:tilesize
            xn = x[n]
            λn = λ[n]
            yn = √(1.0 - (λn^2 * (1.0 - xn^2)))
            γn = γ[n]
            ρn = ρ[n]
            σn = σ[n]
            c_minusn = muladd(λn, yn, -xn)
            c_plusn = muladd(λn, yn, xn)
            r_on = r_o[n]
            r_dn = r_d[n]
            vr_on = (γn * muladd(-ρn, c_plusn, c_minusn)) / r_on
            vr_dn = (-γn * muladd(ρn, c_plusn, c_minusn)) / r_dn
            vt_n = γn * (σn * muladd(λn, xn, yn))
            vt_on = vt_n / r_on
            vt_dn = vt_n / r_dn
            x̂_on = x̂_o[n]
            ŷ_on = ŷ_o[n]
            ẑ_on = ẑ_o[n]
            x̂_dn = x̂_d[n]
            ŷ_dn = ŷ_d[n]
            ẑ_dn = ẑ_d[n]
            is_backwardsn = is_backwards[n]
            fsignn = (1.0 * !(is_backwardsn) - 1.0 * is_backwardsn)
            n̂xn = fsignn * ((ŷ_on * ẑ_dn) - (ẑ_on * ŷ_dn))
            n̂yn = fsignn * ((ẑ_on * x̂_dn) - (x̂_on * ẑ_dn))
            n̂zn = fsignn * ((x̂_on * ŷ_dn) - (ŷ_on * x̂_dn))


            vx′_on = (vr_on * x̂_on) + (vt_on * ((n̂yn * ẑ_on) - (n̂zn * ŷ_on)))
            vy′_on = (vr_on * ŷ_on) + (vt_on * ((n̂zn * x̂_on) - (n̂xn * ẑ_on)))
            vz′_on = (vr_on * ẑ_on) + (vt_on * ((n̂xn * ŷ_on) - (n̂yn * x̂_on)))
            vx′_dn = (vr_dn * x̂_dn) + (vt_dn * ((n̂yn * ẑ_dn) - (n̂zn * ŷ_dn)))
            vy′_dn = (vr_dn * ŷ_dn) + (vt_dn * ((n̂zn * x̂_dn) - (n̂xn * ẑ_dn)))
            vz′_dn = (vr_dn * ẑ_dn) + (vt_dn * ((n̂xn * ŷ_dn) - (n̂yn * x̂_dn)))
            vx_on = vx_o[n]
            vy_on = vy_o[n]
            vz_on = vz_o[n]
            vx_dn = vx_d[n]
            vy_dn = vy_d[n]
            vz_dn = vz_d[n]
            Δv_1[n, i] = (vx′_on - vx_on)^2 + (vy′_on - vy_on)^2 + (vz′_on - vz_on)^2
            Δv_2[n, i] = (vx′_dn - vx_dn)^2 + (vy′_dn - vy_dn)^2 + (vz′_dn - vz_dn)^2
        end
        @tturbo for n ∈ 1:tilesize
            arrival_times[n] += Δt_flight
        end
    end
end






