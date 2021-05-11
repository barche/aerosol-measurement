### A Pluto.jl notebook ###
# v0.14.5

using Markdown
using InteractiveUtils

# ╔═╡ cba7e42f-9a13-4058-a055-51b190be061f
using jlmie

# ╔═╡ f13d1143-8a8b-42ea-8962-d2a6303d4e79
using Plots

# ╔═╡ 952fdfe3-c304-476e-9e6e-71d6654c631f
md"Implementation of Mie scattering from: [https://github.com/Hinamoooon/jlmie](https://github.com/Hinamoooon/jlmie)"

# ╔═╡ 0394e814-c5c9-405a-8912-539ddb9ae9dd
nmat = 1.5  # refractive index of a sphere

# ╔═╡ bd412ab5-d63d-4a24-9578-9728d5175c77
# refractive index of the environment
nenv = 1.00

# ╔═╡ 06400c78-e9f8-4c5e-9ad1-7f763e96b508
# wavelenth (nm)
λ = 683*1e-9

# ╔═╡ 99796b48-d714-4097-b4eb-470d8dbb9fef
# angular range
θ = range(0,π,length=180)

# ╔═╡ 2075b5e5-0b84-4f60-9984-c4135323f5d1
φ = 0

# ╔═╡ 3d25b5e6-4ac7-462f-aed6-347d482dd8bb
function mie(θ,ka)
	k = 2π/λ
	radius = ka / k
	nmax = -1
	Isff,_,_ = jlmie_Isff(nmat,radius,λ,nenv,θ,φ,nmax)
	return Isff
end

# ╔═╡ 597fc000-c15f-42e7-855e-1a7e78104061
function normalize(m)
	result = zeros(size(m))
	for j in 1:size(m,2)
		colmax = maximum(m[:,j])
		result[:,j] .= m[:,j] ./ colmax
	end
	return result
end

# ╔═╡ 0fa32b02-80bf-4ff2-9baa-64e26480396f
ka = [3 6 12]

# ╔═╡ 67cf5661-8ce0-491d-9f6f-02fc3297a8e5
mieresult = normalize(mie.(θ,ka))

# ╔═╡ c929c36f-c96a-4298-b952-78f0a17aeb7a
mieangles = plot(θ/π*180,
	mieresult,
	yaxis=:log,
	ylims=(1e-5,1),
	xticks=0:10:180,
	label="ka = " .* string.(ka),
	linewidth=2,
	xlabel="Scattering angle (°)",
	ylabel="Normalized scattering intensity"
)

# ╔═╡ c7ada24c-af26-47b6-a165-bcc323f8c6e0
#savefig(mieangles, "mieangles.svg")

# ╔═╡ 7d8b3fa2-2c0e-4b58-8fe0-b99b01b99629
function mieforward(dp,λ,n)
	θ = 0
	radius = dp/2
	nmax = -1
	Isff,_,_ = jlmie_Isff(n,radius,λ,nenv,θ,φ,nmax)
	return Isff
end

# ╔═╡ a5544416-ac12-43aa-bf09-87140d0fbe95
dplong = range(0.2e-6,3.0e-6,length=150)

# ╔═╡ abb16a3c-1fa3-40ce-8f51-45115877ed2e
miesizes = mieforward.(dplong, λ, nmat)

# ╔═╡ b4278f7d-c6e1-4d83-9295-f6c18cfc35e3
sizeresponse = plot(dplong.*1e6,
	miesizes,
	xaxis=:log,
	yaxis=:log,
	legend=nothing,
	linewidth=2,
	xticks=([0.2,0.5,1,2,3],[0.2,0.5,1,2,3]),
	xlabel="Particle diameter (μm)",
	ylabel="Scattering intensity")

# ╔═╡ c9fef4c4-4128-4142-9e9a-045a5acb0c47
#savefig(sizeresponse, "mie-sizeresponse.svg")

# ╔═╡ 0a829220-73d1-4377-8b20-3f4b228329e2
dp = range(0.5e-6,1.5e-6,length=150)

# ╔═╡ 41238f61-8725-4882-813d-09d6ce0ceff3
λrange = [380 500 750]

# ╔═╡ 174b91eb-7aca-4604-ba7d-477ba9d02e51
miewavelengths = plot(dp.*1e6,
	mieforward.(dp, λrange.*1e-9, nmat),
	xaxis=:log,
	yaxis=:log,
	legend=:topleft,
	linewidth=2,
	xticks=(0.5:0.1:1.5,0.5:0.1:1.5),
	xlabel="Particle diameter (μm)",
	ylabel="Scattering intensity",
	label="λ = " .* string.(λrange) .* " nm")

# ╔═╡ 25d017c6-e67b-463c-9a03-6b252a4b8406
#savefig(miewavelengths, "mie-wavelengths.svg")

# ╔═╡ 5da8c6c7-09ec-486a-878b-1cd4fb60f6b3
nrange = [1.3 1.4 1.5 1.6]

# ╔═╡ 8e49e929-a630-4d8c-b74f-02dc3fbf7af3
mierefidx = plot(dp.*1e6,
	mieforward.(dp, λ, nrange),
	xaxis=:log,
	yaxis=:log,
	legend=:topleft,
	linewidth=2,
	xticks=(0.5:0.1:1.5,0.5:0.1:1.5),
	xlabel="Particle diameter (μm)",
	ylabel="Scattering intensity",
	label="n = " .* string.(nrange)
)

# ╔═╡ d059b47c-7483-4a0f-9104-9825fe92ee19
#savefig(mierefidx, "mie-refidx.svg")

# ╔═╡ Cell order:
# ╠═cba7e42f-9a13-4058-a055-51b190be061f
# ╠═f13d1143-8a8b-42ea-8962-d2a6303d4e79
# ╠═952fdfe3-c304-476e-9e6e-71d6654c631f
# ╠═0394e814-c5c9-405a-8912-539ddb9ae9dd
# ╠═bd412ab5-d63d-4a24-9578-9728d5175c77
# ╠═06400c78-e9f8-4c5e-9ad1-7f763e96b508
# ╠═99796b48-d714-4097-b4eb-470d8dbb9fef
# ╠═2075b5e5-0b84-4f60-9984-c4135323f5d1
# ╠═3d25b5e6-4ac7-462f-aed6-347d482dd8bb
# ╠═597fc000-c15f-42e7-855e-1a7e78104061
# ╠═0fa32b02-80bf-4ff2-9baa-64e26480396f
# ╠═67cf5661-8ce0-491d-9f6f-02fc3297a8e5
# ╠═c929c36f-c96a-4298-b952-78f0a17aeb7a
# ╠═c7ada24c-af26-47b6-a165-bcc323f8c6e0
# ╠═7d8b3fa2-2c0e-4b58-8fe0-b99b01b99629
# ╠═a5544416-ac12-43aa-bf09-87140d0fbe95
# ╠═abb16a3c-1fa3-40ce-8f51-45115877ed2e
# ╠═b4278f7d-c6e1-4d83-9295-f6c18cfc35e3
# ╠═c9fef4c4-4128-4142-9e9a-045a5acb0c47
# ╠═0a829220-73d1-4377-8b20-3f4b228329e2
# ╠═41238f61-8725-4882-813d-09d6ce0ceff3
# ╠═174b91eb-7aca-4604-ba7d-477ba9d02e51
# ╠═25d017c6-e67b-463c-9a03-6b252a4b8406
# ╠═5da8c6c7-09ec-486a-878b-1cd4fb60f6b3
# ╠═8e49e929-a630-4d8c-b74f-02dc3fbf7af3
# ╠═d059b47c-7483-4a0f-9104-9825fe92ee19
