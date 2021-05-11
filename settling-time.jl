### A Pluto.jl notebook ###
# v0.14.5

using Markdown
using InteractiveUtils

# ╔═╡ bd8544fa-1ef4-11eb-03cc-67a7b8dbb756
using Unitful: μm, m, s, °C, Pa, kg, mm, minute, hr

# ╔═╡ 508d3bb8-1ef5-11eb-0901-e398cefde302
md"""
# Particle settling basics

## Particle settling velocity and time

The time required for a particle to settle can be derived from the general particle equation of motion, as derived by e.g. Maxey and Riley (<https://doi.org/10.1063/1.864230>). Simplifying the equation to only account for the Stokes drag and gravitational force, the equilibrium free-fall velocity $v$ can be derived from:

$\frac{\rho_p - \rho_f}{\rho_p} g = \frac{18\mu}{\rho_p d_p^2} v$

Here, $\rho_p$ and $\rho_f$ are the particle and fluid densities, $d_p$ is the particle diameter and $\mu$ is the dynamic viscosity of the fluid. Solving for $v$ and assuming that the particle density is much greater than the fluid density yields: 

$v = \frac{\rho_p d_p^2}{18\mu}g$

Using this velocity, the time required for a particle to fall to the ground from a given height can be calculated. All this is under the assumption that the air is still and that other than gravity and Stokes drag no other forces act on the particle. It is also only valid (in air) for particles greater than about 1 $\mu m$.

### Numerical example

We consider a water droplet in air at 20°C:
"""

# ╔═╡ 871e3060-1ef5-11eb-2d30-dd398a7e85f9
μ = 1.82e-5kg/(m*s) # Air viscosity at 20°C

# ╔═╡ 8acf72be-1ef5-11eb-2ae6-a55199dca7d2
ρ_p = 1000kg/m^3 # Particle density

# ╔═╡ 8ad8cea4-1ef5-11eb-0cc0-cf4a73b9904f
g = 9.81m/s^2 # Gravitational acceleration;

# ╔═╡ 9a46f728-1ef5-11eb-06ba-635241932ddd
md"This function yields the settling velocity in mm/s, for a given diameter:"

# ╔═╡ 9e38dfb8-1ef5-11eb-0ae5-d16d6b49a53e
v(d) = g*ρ_p*d^2/(18*μ) |> mm/s

# ╔═╡ e3b259e6-1ef5-11eb-1435-29315e8a24e7
md"From this, we can compute the time it takes for a particle to fall 1.5 m:"

# ╔═╡ 03ff4fa4-1ef6-11eb-2fe2-7103efae1a87
tfall(d) = 1.5m/v(d) |> hr

# ╔═╡ be80bc26-1ef5-11eb-1bb3-a37d809b922f
md"Let's set a diameter:"

# ╔═╡ b7bb626a-1ef5-11eb-35a4-7bafab98e433
d = 3μm

# ╔═╡ d85dc2ba-1ef5-11eb-2324-8bf5eed6afce
md"""
For a diameter of $d, the particle:

* Reaches terminal velocity **$(v(d))**
* Takes **$(tfall(d))** to drop 1.5m
"""

# ╔═╡ a8e0a5a4-1ef8-11eb-066d-5147448c2c79
md"""
## Future work

The settling equations described here are combined with an evaporation equation in an article by Netz and Eaton: [https://www.pnas.org/content/117/41/25209](https://www.pnas.org/content/117/41/25209)

To get an idea of which particles evaporate and which particles fall to the floor, these equations could also be added to this notebook.
"""

# ╔═╡ Cell order:
# ╠═bd8544fa-1ef4-11eb-03cc-67a7b8dbb756
# ╟─508d3bb8-1ef5-11eb-0901-e398cefde302
# ╠═871e3060-1ef5-11eb-2d30-dd398a7e85f9
# ╠═8acf72be-1ef5-11eb-2ae6-a55199dca7d2
# ╠═8ad8cea4-1ef5-11eb-0cc0-cf4a73b9904f
# ╟─9a46f728-1ef5-11eb-06ba-635241932ddd
# ╠═9e38dfb8-1ef5-11eb-0ae5-d16d6b49a53e
# ╟─e3b259e6-1ef5-11eb-1435-29315e8a24e7
# ╠═03ff4fa4-1ef6-11eb-2fe2-7103efae1a87
# ╟─be80bc26-1ef5-11eb-1bb3-a37d809b922f
# ╠═b7bb626a-1ef5-11eb-35a4-7bafab98e433
# ╟─d85dc2ba-1ef5-11eb-2324-8bf5eed6afce
# ╟─a8e0a5a4-1ef8-11eb-066d-5147448c2c79
