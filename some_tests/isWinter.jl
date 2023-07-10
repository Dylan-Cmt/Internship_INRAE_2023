### A Pluto.jl notebook ###
# v0.19.26

using Markdown
using InteractiveUtils

# ╔═╡ 73a0c096-1e1c-4f62-a4b5-587c5f18cdb1
#imports
using Parameters, StaticArrays

# ╔═╡ 42509056-371b-4d50-b52d-caeb086969cc
md"""
> À priori, je ne peux pas réutiliser `isWinter` de lmaillere car il ne s'agit pas d'une matrice mais d'un vecteur.
>
> Je suis donc obligé de procéder de faire une boucle.
"""

# ╔═╡ ba921fbc-48dc-4777-a2a0-af46f68c16c3
# TimeParam
struct TimeParam
	T::Int64
	τ::Int64
end

# ╔═╡ f4567a59-620f-4410-b1fe-12def48eaf48
begin
	# few objects for example
	tp = TimeParam(3, 1)
	svect = @SVector [1, 2, 3, 4, 5]
	mat = Matrix{SVector{5,Float64}}(undef, 2, 3)
	mat[1,1] = svect ; mat[2,1] = svect
end

# ╔═╡ a0e58f52-b504-4bd6-97d3-63cdd10dc1ae


# ╔═╡ ffcbe747-5375-40d0-8304-18e9a9a4bf51


# ╔═╡ e57becfe-1c8b-11ee-0032-79ae0303f048
# for vectors
isWinter(t,tp) =[mod(x, tp.T) < tp.τ ? 0 : 1 for x in t]

# ╔═╡ 1a27e8d1-f8f3-41f3-83aa-721f15d654b7
# for matrices
isWinter2(t,tp) =[isWinter(x,tp) for x in t[:,1]]

# ╔═╡ 8a154b29-3b69-47dd-b2b0-e0c6c8429e8d
isWinter3(t,tp) = [mod(x, tp.T) < tp.τ ? 0 : 1 for x in t]

# ╔═╡ e23b899b-fbfe-4369-b9a5-dd2c5cf1ccb3


# ╔═╡ 275cdb9e-600c-4564-93ce-16511b2f0aa6


# ╔═╡ a07e1d15-2973-46dc-9940-debeff38b6ba
# example with vectors

# ╔═╡ 53cb6476-5be3-4278-8377-11d3e5756425
isWinter(svect,tp)

# ╔═╡ 6646cf78-6069-4274-93be-3650b243db54
isWinter(mat[1,1],tp) # equivalent to isWinter(svect,tp)

# ╔═╡ bc636614-4e26-4bc2-a5e6-8d042b6469a6
# example with matrice

# ╔═╡ c557230a-1d2e-4e27-b095-1e50a6a4731f
isWinter2(mat[:,1],tp) # pareil mais sur chaque ligne

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
Parameters = "d96e819e-fc66-5662-9728-84c9c7592b0a"
StaticArrays = "90137ffa-7385-5640-81b9-e52037218182"

[compat]
Parameters = "~0.12.3"
StaticArrays = "~1.6.0"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.9.1"
manifest_format = "2.0"
project_hash = "bbfd7e164f63b3a36d0f5bbd02d685a214604917"

[[deps.Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"

[[deps.CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"
version = "1.0.2+0"

[[deps.Libdl]]
uuid = "8f399da3-3557-5675-b5ff-fb832c97cbdb"

[[deps.LinearAlgebra]]
deps = ["Libdl", "OpenBLAS_jll", "libblastrampoline_jll"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"

[[deps.OpenBLAS_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "4536629a-c528-5b80-bd46-f80d51c5b363"
version = "0.3.21+4"

[[deps.OrderedCollections]]
git-tree-sha1 = "d321bf2de576bf25ec4d3e4360faca399afca282"
uuid = "bac558e1-5e72-5ebc-8fee-abe8a469f55d"
version = "1.6.0"

[[deps.Parameters]]
deps = ["OrderedCollections", "UnPack"]
git-tree-sha1 = "34c0e9ad262e5f7fc75b10a9952ca7692cfc5fbe"
uuid = "d96e819e-fc66-5662-9728-84c9c7592b0a"
version = "0.12.3"

[[deps.Random]]
deps = ["SHA", "Serialization"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"

[[deps.SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"
version = "0.7.0"

[[deps.Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"

[[deps.StaticArrays]]
deps = ["LinearAlgebra", "Random", "StaticArraysCore"]
git-tree-sha1 = "0da7e6b70d1bb40b1ace3b576da9ea2992f76318"
uuid = "90137ffa-7385-5640-81b9-e52037218182"
version = "1.6.0"

    [deps.StaticArrays.extensions]
    StaticArraysStatisticsExt = "Statistics"

    [deps.StaticArrays.weakdeps]
    Statistics = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"

[[deps.StaticArraysCore]]
git-tree-sha1 = "6b7ba252635a5eff6a0b0664a41ee140a1c9e72a"
uuid = "1e83bf80-4336-4d27-bf5d-d5a4f845583c"
version = "1.4.0"

[[deps.UnPack]]
git-tree-sha1 = "387c1f73762231e86e0c9c5443ce3b4a0a9a0c2b"
uuid = "3a884ed6-31ef-47d7-9d2a-63182c4928ed"
version = "1.0.2"

[[deps.libblastrampoline_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850b90-86db-534c-a0d3-1478176c7d93"
version = "5.8.0+0"
"""

# ╔═╡ Cell order:
# ╟─42509056-371b-4d50-b52d-caeb086969cc
# ╠═73a0c096-1e1c-4f62-a4b5-587c5f18cdb1
# ╠═ba921fbc-48dc-4777-a2a0-af46f68c16c3
# ╠═f4567a59-620f-4410-b1fe-12def48eaf48
# ╟─a0e58f52-b504-4bd6-97d3-63cdd10dc1ae
# ╟─ffcbe747-5375-40d0-8304-18e9a9a4bf51
# ╠═e57becfe-1c8b-11ee-0032-79ae0303f048
# ╠═1a27e8d1-f8f3-41f3-83aa-721f15d654b7
# ╠═8a154b29-3b69-47dd-b2b0-e0c6c8429e8d
# ╟─e23b899b-fbfe-4369-b9a5-dd2c5cf1ccb3
# ╟─275cdb9e-600c-4564-93ce-16511b2f0aa6
# ╠═a07e1d15-2973-46dc-9940-debeff38b6ba
# ╠═53cb6476-5be3-4278-8377-11d3e5756425
# ╠═6646cf78-6069-4274-93be-3650b243db54
# ╠═bc636614-4e26-4bc2-a5e6-8d042b6469a6
# ╠═c557230a-1d2e-4e27-b095-1e50a6a4731f
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
