using LoopVectorization
using GLMakie
include("KeplerSolvers.jl")
include("LambertSolvers.jl")


μ::Float64 = 1.32712440042 * 10.0^20
M_ref_o::Float64 = deg2rad(100.46435)
a_o::Float64 = 1.00000011 * 149597870700.0
e_o::Float64 = 0.01671022
M_ref_d::Float64 = deg2rad(355.45332)
a_d::Float64 = 1.52366231 * 149597870700.0
e_d::Float64 = 0.09341233


t_departure::Float64 = 0.0 * 86400.0
Δt_departure::Float64 = 0.2 * 86400.0
tof_min::Float64 = 50.0 * 86400.0
Δt_flight::Float64 = 0.2 * 86400.0
tof_samples::Integer = 4096
my_tilesize::Integer = 4096


DCM_o::Matrix{Float64} = orbit_DCM(deg2rad(0.00005), deg2rad(-11.26064), deg2rad(102.94719))
DCM_d::Matrix{Float64} = orbit_DCM(deg2rad(1.85061), deg2rad(49.57854), deg2rad(336.04084))


Δv_1::Matrix{Float64} = zeros(Float64, my_tilesize, tof_samples)
Δv_2::Matrix{Float64} = zeros(Float64, my_tilesize, tof_samples)
Lambert_solve!(Δv_1, Δv_2, μ, M_ref_o, M_ref_d, a_o, a_d, e_o, e_d, DCM_o, DCM_d, 
    t_departure, Δt_departure, tof_min, Δt_flight, tof_samples; tilesize = my_tilesize
)
cost::Matrix{Float64} = zeros(Float64, my_tilesize, tof_samples)
arrival_time::Matrix{Float64} = zeros(Float64, my_tilesize, tof_samples)
@tturbo for j ∈ 1:tof_samples
    for i ∈ 1:my_tilesize
        cost[i, j] = √(Δv_1[i, j])# + √(Δv_2[i, j])
        arrival_time[i, j] = (t_departure + Δt_departure * i + Δt_flight * j) / 86400.0
    end
end


departure_time_vec::Vector{Float64} = collect(range(t_departure, step = Δt_departure, length = my_tilesize)) / 86400.0
time_of_flight_vec::Vector{Float64} = collect(range(tof_min, step = Δt_flight, length = tof_samples)) / 86400.0
day_step::Integer = 50
arrival_time_vec::Vector{Float64} = collect(range(floor(Integer, t_departure / 86400.0), step = day_step, length = floor(Integer, maximum(arrival_time) / day_step)))


fig = Figure(size = (1600, 900))
ax = Axis(fig[1, 1], aspect = DataAspect(), title = "Earth to Mars Transfer Porkchop Plot (w/. aerobraking)", ylabel = "Time of Flight (Julian Days)", xlabel = "Departure Time Relative to J2000 Epoch (Julian Days)")
hm = heatmap!(ax, departure_time_vec, time_of_flight_vec, cost, colorscale = log10, colormap = :thermal)
Colorbar(fig[:, end + 1], hm, label = "Delta-V Cost")
contour!(ax, departure_time_vec, time_of_flight_vec, arrival_time, levels = arrival_time_vec, labels = true, color = :black, labelformatter = x -> "$x days")
contour!(ax, departure_time_vec, time_of_flight_vec, cost, levels = 10.0.^range(3.0, 5.0; length = 15), labels = true, color = :white, labelformatter = x -> "$x m / s")
save(normpath(joinpath((@__FILE__), raw"..\marswithaerobrake.png")), fig)
fig






