### A Pluto.jl notebook ###
# v0.19.26

using Markdown
using InteractiveUtils

# ╔═╡ 73a0c096-1e1c-4f62-a4b5-587c5f18cdb1
#imports
using Parameters, StaticArrays

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
	vect = [1, 2, 3, 4, 5]
	svect = @SVector [1, 2, 3, 4, 5]
	mat_SV = Matrix{SVector{5,Float64}}(undef, 2, 3)
	mat_V = Matrix{Vector{Float64}}(undef, 2, 2)
	mat_SV[1,1] = svect ; mat_SV[2,1] = svect
	mat_V[1,1] = vect ; mat_V[2,1] = vect ; mat_V[1,2] = vect ; mat_V[2,2] = vect
end

# ╔═╡ 2d5d7444-5eda-44db-be25-eccffc444fa4
mat_V

# ╔═╡ a0e58f52-b504-4bd6-97d3-63cdd10dc1ae


# ╔═╡ ffcbe747-5375-40d0-8304-18e9a9a4bf51


# ╔═╡ df69e50f-9329-41b3-b1f8-0b3ac5286d6c
md"""
> Le broadcast consiste à ajouter un "." après une fonction lors de son appel afin d'appliquer celle-ci sur un (S)Vecteur terme à terme.
"""

# ╔═╡ a8f7d66a-7dfc-4f3c-a190-8ca18ff0af0a
isWinterr(t) = mod(t, 15) < 10 ? 0 : 1

# ╔═╡ 22c8aafc-c131-4d73-91f4-ae6be4abe3ee
isWinterr.([0,3,6,9,12,15,18,21,24,27])

# ╔═╡ 9bcc28ad-0569-4486-ab40-e37e3d583a9f


# ╔═╡ 275cdb9e-600c-4564-93ce-16511b2f0aa6


# ╔═╡ 85fcd43e-d46f-4e4a-8b5d-ff729ca2bbd7
md"""
> Le broadcast fonctionne  avec les Vectors mais aussi les SVectors.
"""

# ╔═╡ 9ac45b01-b429-4fd9-8bc7-a0f874f0d23f
f(x) = exp(x*3)

# ╔═╡ 4051f5f9-2b29-415d-be7b-f5b47a3e15ad
exp.([0,3,6,9,12,15,18,21,24,27])

# ╔═╡ a9bce063-00b4-432f-8a93-41c6f5b98a67
f.(@SVector [0,3,6,9,12,15,18,21,24,27])

# ╔═╡ de104a2c-9606-471e-8d84-d93449a80123
md"""
> Cela fonctionne aussi lorsqu'on ajoute des arguments à broadcaster ou non.
"""

# ╔═╡ 938a4471-612d-43f6-8850-0a259df5f242
g(x,y) = exp(x*y)

# ╔═╡ 90c20b23-7446-43d8-af57-da13d1aef043
g.([0,1,2,3], [0,1,2,3])

# ╔═╡ 16e6e1de-86b9-4017-b9ad-8bff341c1198
g.([0,1,2,3], 3) # fait x*3 puis passe à l'exp

# ╔═╡ b9670e20-24de-4bd9-bcbb-15b7f5abe25b
md"""
> Mais c'est un peu plus délicat dans notre cas...
"""

# ╔═╡ ec1e935b-5313-4fa1-a481-4bbaa25e9425
h(x,tp) = exp(x*tp.T)

# ╔═╡ 4df34fda-1e76-47a1-a8e9-90bb7fd150d4
h.([0,1,2,3], TimeParam(T=2, τ=1))

# ╔═╡ 1f905ee5-deba-48bf-a1d0-b47add9a28d5


# ╔═╡ 0db36654-b594-4b4c-b192-667822870b35


# ╔═╡ 9913f3e2-514c-43e8-a18f-08b680d904cc
md"""
> Pour le moment, nous ne broadcastons pas.
>
> Sinon, une solution serait envisageable: définir isWinter dans la fonction `plot` de sorte à ce qu'il n'y ait qu'un seul argument `t`.
>
> Cependant, comme nous écrivons plusieurs fonctions `plots` nous n'allons pas procéder ainsi.
"""

# ╔═╡ 99a7fd3c-c0c4-4620-b40b-02f0cbeb6b23
# parcourt tous les elts du vecteur pour appliquer notre fct
isWinter_vect(t,tp) =[mod(x, 1) < tp.τ/tp.T ? 0 : 1 for x in t]

# ╔═╡ 1fef6e43-490e-45c0-8e4a-98554200595e
# parcourt la première colonne pour appliquer notre fct
isWinter(t,tp) =[isWinter_vect(x,tp) for x in t[:,1]]

# ╔═╡ 7ca12a1d-2490-4f6b-80c2-c066abf8581f


# ╔═╡ 51cd44cb-4e69-43af-a06d-0112ab2e133c
# example with vectors
isWinter_vect(svect,tp) # équivaut à isWinter_vect(mat_SV[1,1],tp)

# ╔═╡ 13aab6f6-8d79-4a87-a09f-8bba00196ee6
# example with matrice
isWinter(mat_SV[:,1],tp) # pareil mais sur chaque ligne

# ╔═╡ f994071e-465d-45e6-bb8a-e3f0aabc6d55
isWinter(mat_SV,tp) # mettre la matrice en entier ou que sa première colonne revient au même

# ╔═╡ 2ceb0cb7-2c47-482f-a1d8-8e0da28d7758
@time isWinter(Matrix{SVector{10,Float64}}(undef, 5, 5),tp);

# ╔═╡ 35dda99b-7011-438e-879c-cf7eba8b00a0


# ╔═╡ de42206f-2309-419b-9314-20f7696d7afb


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

julia_version = "1.9.0"
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
version = "5.7.0+0"
"""

# ╔═╡ Cell order:
# ╠═73a0c096-1e1c-4f62-a4b5-587c5f18cdb1
# ╠═ba921fbc-48dc-4777-a2a0-af46f68c16c3
# ╠═f4567a59-620f-4410-b1fe-12def48eaf48
# ╠═2d5d7444-5eda-44db-be25-eccffc444fa4
# ╟─a0e58f52-b504-4bd6-97d3-63cdd10dc1ae
# ╟─ffcbe747-5375-40d0-8304-18e9a9a4bf51
# ╟─df69e50f-9329-41b3-b1f8-0b3ac5286d6c
# ╠═a8f7d66a-7dfc-4f3c-a190-8ca18ff0af0a
# ╠═22c8aafc-c131-4d73-91f4-ae6be4abe3ee
# ╟─9bcc28ad-0569-4486-ab40-e37e3d583a9f
# ╟─275cdb9e-600c-4564-93ce-16511b2f0aa6
# ╟─85fcd43e-d46f-4e4a-8b5d-ff729ca2bbd7
# ╠═9ac45b01-b429-4fd9-8bc7-a0f874f0d23f
# ╠═4051f5f9-2b29-415d-be7b-f5b47a3e15ad
# ╠═a9bce063-00b4-432f-8a93-41c6f5b98a67
# ╟─de104a2c-9606-471e-8d84-d93449a80123
# ╠═938a4471-612d-43f6-8850-0a259df5f242
# ╠═90c20b23-7446-43d8-af57-da13d1aef043
# ╠═16e6e1de-86b9-4017-b9ad-8bff341c1198
# ╟─b9670e20-24de-4bd9-bcbb-15b7f5abe25b
# ╠═ec1e935b-5313-4fa1-a481-4bbaa25e9425
# ╠═4df34fda-1e76-47a1-a8e9-90bb7fd150d4
# ╟─1f905ee5-deba-48bf-a1d0-b47add9a28d5
# ╟─0db36654-b594-4b4c-b192-667822870b35
# ╟─9913f3e2-514c-43e8-a18f-08b680d904cc
# ╠═99a7fd3c-c0c4-4620-b40b-02f0cbeb6b23
# ╠═1fef6e43-490e-45c0-8e4a-98554200595e
# ╟─7ca12a1d-2490-4f6b-80c2-c066abf8581f
# ╠═51cd44cb-4e69-43af-a06d-0112ab2e133c
# ╠═13aab6f6-8d79-4a87-a09f-8bba00196ee6
# ╠═f994071e-465d-45e6-bb8a-e3f0aabc6d55
# ╠═2ceb0cb7-2c47-482f-a1d8-8e0da28d7758
# ╟─35dda99b-7011-438e-879c-cf7eba8b00a0
# ╟─de42206f-2309-419b-9314-20f7696d7afb
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
