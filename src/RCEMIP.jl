module RCEMIP

# Parameters

# Table 1
Ω = 0.0
f = 0.0
R_E = 6371.0
g = 9.79764
Rd = 287.04
cpd = 1004.64
Rv = 461.50
cpv = 1846.0
Lv0 = 2.501e6
Lf0 = 3.337e5
Ls0 = 2.834e6

# Table 2 above the line
η_CO2 = 348e-6
η_CH4 = 1650e-9
η_N2O = 306e-9
η_CFC11 = 0.0
η_CFC12 = 0.0
η_CFC22 = 0.0
η_CCL4 = 0.0
g_1 = 3.6478e-6
g_2 = 0.83209
g_3 = 11.3515

# Table 2 below the line
z_t = 15.0
q_0295 = 12e-3
q_0300 = 18.65e-3
q_0305 = 24e-3
q_t = 1e-14
z_q1 = 4000e-3
z_q2 = 7500e-3
Γ = 0.0067
p_0 = 1014.8


# Profiles
"""
Calculate ozone profile

Parameters:
	p: pressure (hPa)

Returns:
	ozone molar concentration (nondim)
"""
function η_O3(p)
	return g_1 .* p.^g_2 .* exp(-p./g3)
end

"""
Calculate water vapor profile 

Parameters:
	z: height (km)
	T: surface temperature (K)

Keyword arguments:
	use_extended: calculate surface specific humidities
	appropriate for a wider range
	of surface temperature than used by the RCEMIP protocol.
	Default is false; must be set to true is T is not equal
	to 295, 300, or 305.

Returns:
	Water vapor mass concentration (nondim)
"""
function q_v(z, T; use_extended = false)

	if use_extended
		q0 = q_0_extended(T)
	elseif T == 295
		q_0 = q_0295
	elseif T == 300
		q_0 = q_0300
	elseif T == 305
		q_0 = q_0305
	else
		throw(ArgumentError(
			"must set use_extended = true if T is not one of 295, 300, and 305"))
	end

	q = q_0 .* exp(-z./z_q1) .* exp(-(z./z_q2).^2.0)
	strat = findall(z .> z_t)
	q[strat] = q_t
	return q

end

"""
Calculate surface specific humidity for non-standard surface temperature

Parameters:
	T: surface temperature (K)

Returns:
	Surface specific humidity (nondim) that gives a surface RH  of ≈0.8
"""
function q_0_extended(T)
	esat0 = 611.0
	T0 = 273.15
	RH_s = 0.8
	e = RH_s .* esat0 .* exp(-Lv0./Rd .* (1.0./T .- 1.0./T0))
	return (Rd / Rv) .* e ./ (1e2 * p_s)
end

"""
Calculate virtual temperature profile

Parameters:
	z: height (km)
	T: surface temperature (K)

Returns:
	Virtual temperature (K)
"""
function T_v(z, T)

	T_v0 = T
	T_v = T_v0 .- Γ.*z
	strat = findall(z .> z_t)
	T_t = T_v[strat[1]]
	T_v[strat] = T_t
	return T_v

end

"""
Calculate temperature

Parameters:
	T_v: virtual temperature (K)
	q: specific humidity (nondim)

Returns:
	Temperature (K)
"""
function T(T_v, q)
	return T_v ./ (1.0 .+ 0.608.*q)
end

"""
Calculate pressure

Parameters:
	z: height (km)
	T: surface temperature (K)

Returns:
	Pressure (hPa)
"""
function p(z, T)

	T_v0 = T
	T_v = T_v0 .- Γ.*z
	p = p_0 * (T_v./T_v0).^(g./(Rd.*Γ))
	strat = findall(z .> z_t)
	itp = strat[1]
	p_tp = p[itp]
	z_tp = z[itp]
	T_vtp = T_v[itp]
	p[itp:end] = p_tp .* exp(-g./(Rd.*T_vtp).*(z[itp:end] .- z_tp))
	return p 

end

end # module
