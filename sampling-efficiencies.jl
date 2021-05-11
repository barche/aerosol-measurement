### A Pluto.jl notebook ###
# v0.14.5

using Markdown
using InteractiveUtils

# ╔═╡ 1631d84e-adbf-11eb-1c2a-5540466cca90
using Plots

# ╔═╡ 5d985f72-22cb-43e3-a62c-75011a250797
using Unitful: cm, m, μm, L, minute, s, kg, NoUnits, °, K, N

# ╔═╡ b02afae3-a407-47cd-8870-ea8ead79fa22
using UnitfulRecipes

# ╔═╡ a4b8d878-01aa-4326-94fc-1911a716aee3
md"""
# Sampling efficiencies

Case study: [Palas Frog](https://www.palas.de/en/product/fidasfrog)

* Size range: 0.18 μm - 18 μm (or to 100 μm)
* Sampling tube length: 12 cm
* Sampling tube diameter: 5 mm.
* Flow rate: 1.4 l/min

"""

# ╔═╡ 25908de6-1500-413f-a34e-787bca795406
g = 9.81m/s^2

# ╔═╡ 57e6cab1-af6a-4e09-a38b-c308103d97fa
md"### Calculate sampling velocity u"

# ╔═╡ 69a117a3-5303-49c7-b4b9-17e378f4f90a
Q = 1.4L/minute

# ╔═╡ aab89b38-ca99-4cf2-ba1e-ac1813047b7d
d_tube = 0.5cm

# ╔═╡ 2480e551-a6a1-4851-b176-dfc208f1f128
L_tube = 12cm

# ╔═╡ 39d6cc4f-6892-4c13-8304-e43bc6329f97
A = π*d_tube^2/4

# ╔═╡ ec805ecc-05f4-4dc5-9011-79c26c5647f5
u = Q/A |> m/s

# ╔═╡ 152d407f-45b1-4a71-aeaf-8d898218ae67
md"### Particle properties"

# ╔═╡ 614c27b4-a9ad-4bc0-92e9-2cfc4fbe5b9b
μ = 1.82e-5kg/(m*s) # Air viscosity at 20°C

# ╔═╡ a0e9795c-288a-4465-8b67-e3817559c477
ρp = 1000kg/m^3

# ╔═╡ 09e7d6e4-d2ce-42ff-8bcf-f9c639a16d7e
md"Particle relaxation time:"

# ╔═╡ 0017c34e-6314-481c-9e6e-65c379986ed8
τp(d) = ρp*d^2/(18*μ) |> s

# ╔═╡ 6c288551-7a8c-41ed-a66c-435ee3c7d8bd
md"### Aspiration efficiencies"

# ╔═╡ 1ca9c0d3-a689-4af2-8b5b-5d090871f593
u0 = [1m/s u 1.5m/s]

# ╔═╡ d9431135-b9e3-4b9f-a2b6-95d0ca239f32
drange = 0.18μm:0.1μm:18μm

# ╔═╡ 136fa812-49ef-478a-ac92-1ee4df4d93eb
St = 0.001:0.001:10

# ╔═╡ 738fde5b-22f2-490d-ad00-d3bfe336419f
md"The formula of [Rader & Marple, 1987](https://doi.org/10.1080/02786828808959190):"

# ╔═╡ bd73e435-0df8-4347-bd3a-3eb72616042a
function η_asp_rader_marple(dp,u0)
	St = τp(dp)*u0/d_tube
	return 1+(u0/u-1)*(1-1/(1+3.77*St^0.883))
end

# ╔═╡ 43e65cfd-2c13-44a9-9849-c110710a10dc
plot(drange, η_asp_rader_marple.(drange, u0),xlabel="Particle size", ylabel="Aspiration efficiency")

# ╔═╡ 04d59cef-db45-4e70-b9f8-8270b3af647b
md"The correlation from [Durham & Lundgren, 1980](https://doi.org/10.1016/0021-8502(80)90033-6):"

# ╔═╡ f25169e9-cc4d-4e1d-bf56-4baed5a4636e
function η_asp_durgam_lundgren(St,u0_u,θ)
	stprime = St*exp(0.022*θ)
	return 1+(u0_u*cosd(θ)-1) * (1-(1+(2+0.617*u0_u)*stprime)^(-1))/(1-(1+2.617*stprime)^(-1)) * (1-(1+0.55*stprime*exp(0.25*stprime))^(-1))
end

# ╔═╡ 36f09930-46c7-4fd0-9c67-2a088f59d2a3
md"Sampling in calm air, following [Grinshpun, Willeke & Kalatoor, 1993](https://doi.org/10.1016/0960-1686(93)90132-I):"

# ╔═╡ 1738a078-c946-4cd8-a594-f34d1067e47f
function η_asp_grishpun(dp, φ)
	St = τp(dp)*u/d_tube
	vts = τp(dp)*g
	return vts/u * cos(φ) + exp(-(4*St^(1+sqrt(vts/u)))/(1+2*St))
end

# ╔═╡ 450deebf-cffc-4a5d-b879-955999b75de0
plot(drange, η_asp_grishpun.(drange, 0°),xlabel="Particle size", ylabel="Aspiration efficiency")

# ╔═╡ d76f6774-af38-4637-bd76-a1dc57c2630c
τp.(drange)*u/d_tube .|> NoUnits

# ╔═╡ ba6257cc-9b47-4d2c-8463-c42ce193e9f4
md"## Transmission efficiencies"

# ╔═╡ fc175159-ff89-42ec-8bee-a3603a48183a
ρ = 1.2kg/m^3

# ╔═╡ 29cb24a0-0376-405d-9c03-0acfb658d8dd
ν = μ/ρ

# ╔═╡ 6419d4a3-eb62-420f-91f1-aa13c54b1176
Re = u*d_tube/ν |> NoUnits

# ╔═╡ 70fcc13f-e442-48d8-9ea6-24b283666a5a
md"Horizontal pipe due to gravity from [Okazaki & Willeke, 1987](https://doi.org/10.1080/02786828708959164):"

# ╔═╡ ac28610b-c22f-41e7-9e59-9ea4d905c4e1
function η_trans_grav_okazaki(dp,u0)
	vts = τp(dp)*g
	Z = L_tube/d_tube*vts/u
	St = τp(dp)*u0/d_tube
	K = sqrt(Z*St)*Re^(-1//4)
	return exp(-4.7*K^(3//4))
end

# ╔═╡ de6f92ea-c5b1-425c-82f1-d72b908a7259
plot(drange, η_trans_grav_okazaki.(drange, [u/2 u 2*u]),xlabel="Particle size", ylabel="Gravitational efficiency")

# ╔═╡ b7d09280-4de7-425a-ad3d-cc083ea0e008
md"#### Diffusion"

# ╔═╡ b200606f-4b37-409d-8682-fff3d6523571
function cunningham_slip(dp)
	λ = 0.0665μm
	Kn = 2*λ/dp
	return 1+Kn*(1.142+0.558*exp(-0.999/Kn))
end

# ╔═╡ f7eb0297-697c-415b-8e3a-294c662bdc2a
function Dcoef(dp,T=293K)
	k = 1.38e-23N*m/K
	return k*T*cunningham_slip(dp)/(3*π*μ*dp) |> m^2/s
end

# ╔═╡ 3f412da9-5142-42d3-b4fe-ddbf896a32cb
md"Holman correlation for diffusion"

# ╔═╡ e0e54dfc-2185-4fc6-8571-d3feeef0ca6d
function η_diffusion_holman(dp)
	ξ = π*Dcoef(dp)*L_tube/Q
	Sh = 3.66 + 0.2672/(ξ+0.10079^(1//3))
	return exp(-ξ*Sh)
end

# ╔═╡ f5d14894-4097-48fc-9e0d-711c618436f8
η_diffusion_holman(0.1μm)

# ╔═╡ Cell order:
# ╠═1631d84e-adbf-11eb-1c2a-5540466cca90
# ╠═5d985f72-22cb-43e3-a62c-75011a250797
# ╠═b02afae3-a407-47cd-8870-ea8ead79fa22
# ╟─a4b8d878-01aa-4326-94fc-1911a716aee3
# ╠═25908de6-1500-413f-a34e-787bca795406
# ╟─57e6cab1-af6a-4e09-a38b-c308103d97fa
# ╠═69a117a3-5303-49c7-b4b9-17e378f4f90a
# ╠═aab89b38-ca99-4cf2-ba1e-ac1813047b7d
# ╠═2480e551-a6a1-4851-b176-dfc208f1f128
# ╠═39d6cc4f-6892-4c13-8304-e43bc6329f97
# ╠═ec805ecc-05f4-4dc5-9011-79c26c5647f5
# ╟─152d407f-45b1-4a71-aeaf-8d898218ae67
# ╠═614c27b4-a9ad-4bc0-92e9-2cfc4fbe5b9b
# ╠═a0e9795c-288a-4465-8b67-e3817559c477
# ╟─09e7d6e4-d2ce-42ff-8bcf-f9c639a16d7e
# ╠═0017c34e-6314-481c-9e6e-65c379986ed8
# ╠═6c288551-7a8c-41ed-a66c-435ee3c7d8bd
# ╠═1ca9c0d3-a689-4af2-8b5b-5d090871f593
# ╠═d9431135-b9e3-4b9f-a2b6-95d0ca239f32
# ╠═136fa812-49ef-478a-ac92-1ee4df4d93eb
# ╟─738fde5b-22f2-490d-ad00-d3bfe336419f
# ╠═bd73e435-0df8-4347-bd3a-3eb72616042a
# ╠═43e65cfd-2c13-44a9-9849-c110710a10dc
# ╟─04d59cef-db45-4e70-b9f8-8270b3af647b
# ╠═f25169e9-cc4d-4e1d-bf56-4baed5a4636e
# ╟─36f09930-46c7-4fd0-9c67-2a088f59d2a3
# ╠═1738a078-c946-4cd8-a594-f34d1067e47f
# ╠═450deebf-cffc-4a5d-b879-955999b75de0
# ╠═d76f6774-af38-4637-bd76-a1dc57c2630c
# ╟─ba6257cc-9b47-4d2c-8463-c42ce193e9f4
# ╠═fc175159-ff89-42ec-8bee-a3603a48183a
# ╠═29cb24a0-0376-405d-9c03-0acfb658d8dd
# ╠═6419d4a3-eb62-420f-91f1-aa13c54b1176
# ╟─70fcc13f-e442-48d8-9ea6-24b283666a5a
# ╠═ac28610b-c22f-41e7-9e59-9ea4d905c4e1
# ╠═de6f92ea-c5b1-425c-82f1-d72b908a7259
# ╠═b7d09280-4de7-425a-ad3d-cc083ea0e008
# ╠═b200606f-4b37-409d-8682-fff3d6523571
# ╠═f7eb0297-697c-415b-8e3a-294c662bdc2a
# ╠═3f412da9-5142-42d3-b4fe-ddbf896a32cb
# ╠═e0e54dfc-2185-4fc6-8571-d3feeef0ca6d
# ╠═f5d14894-4097-48fc-9e0d-711c618436f8
