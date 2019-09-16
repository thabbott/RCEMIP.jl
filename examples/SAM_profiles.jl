# Example: generate input profiles for a simulation with SAM

import RCEMIP
using Printf

# Generate grid
z = RCEMIP.vertical_grid(66.0)

# Generate grd file
fid = open("examples/grd", "w+")
for zz in z
	write(fid, @sprintf(" %11.4f\n", zz))
end
close(fid)

# Generate snd file
T_s = 300.0
q_v = RCEMIP.q_v(z/1e3, T_s; use_extended = true)
T_v = RCEMIP.T_v(z/1e3, T_s)
T = RCEMIP.T(T_v, q_v)
p = RCEMIP.p(z/1e3, T_s)
θ = T .* (1000.0 ./ p).^(287.0./1004.0)
u = zeros(size(z))
v = zeros(size(z))
fid = open("examples/snd", "w+")
write(fid, "  z[m]  p[mb]  tp[K]  q[g/kg]  u[m/s]  v[m/s]\n")
write(fid, @sprintf("  %d   %d   %d\n", 0, length(z), RCEMIP.p_0))
for ii = 1:length(z)
	write(fid, 
		@sprintf("%11.4f %11.4f %11.4f %11.4f %11.4f %11.4f\n",
		z[ii], p[ii], θ[ii], 1e3*q_v[ii], u[ii], v[ii]))
end
write(fid, @sprintf("  %d   %d   %d\n", 1000, length(z), RCEMIP.p_0))
for ii = 1:length(z)
	write(fid, 
		@sprintf("%11.4f %11.4f %11.4f %11.4f %11.4f %11.4f\n",
		z[ii], p[ii], θ[ii], 1e3*q_v[ii], u[ii], v[ii]))
end
close(fid)

# Generate trc file
# Calculate ozone profile
mw_O3 = 3.0 * 16.0
mw_air = 29.0
z = collect(0.00:0.25:70.00)
p = RCEMIP.p(z, T_s)
η_O3 = RCEMIP.η_O3(p)
r_O3 = η_O3 .* mw_O3 ./ mw_air

# Calculate N2O profile
mw_N2O = 2.0 * 14.0 + 16.0
η_N2O = RCEMIP.η_N2O .* ones(size(p))
r_N2O = η_N2O .* mw_N2O ./ mw_air

# Calculate CH4 profile
mw_CH4 = 12.0 + 4.0 * 1.0
η_CH4 = RCEMIP.η_CH4 .* ones(size(p))
r_CH4 = η_CH4 * mw_CH4 / mw_air

# Calculate CFC profiles
r_CFC11 = zeros(size(p))
r_CFC12 = zeros(size(p))

# Create file
fid = open("examples/trc", "w+")
write(fid, 
	@sprintf("   %d ! Number of levels (RCEMIP sounding)\n",
		length(p)
	)
)
write(fid,
	" z(km)   p(mb)   o3(g/g)   n2o(g/g)   ch4(g/g)   cfc11(g/g)   cfc12(g/g)\n"
)
for ii = length(p):-1:1
	write(fid,
		@sprintf("%11.3f %11.3f %11.4e %11.4e %11.4e %11.4e %11.4e\n",
			z[ii], p[ii], r_O3[ii], r_N2O[ii], r_CH4[ii], r_CFC11[ii], r_CFC12[ii]
		)
	)
end
close(fid)