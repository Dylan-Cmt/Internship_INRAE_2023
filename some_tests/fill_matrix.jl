### A Pluto.jl notebook ###
# v0.19.26

using Markdown
using InteractiveUtils

# ╔═╡ 1d78060c-74d6-4425-9cc3-fec449ac1790
using Parameters

# ╔═╡ 06319210-1b11-11ee-0abe-3b9d74e47ab9
t = [[1,2,3],[4,5,6]]

# ╔═╡ bcc22983-bcf7-4ed1-9479-d67069df96e6
typeof(t)

# ╔═╡ 9e9e8b77-79e9-4762-9830-596a29bf02b2
p = [[0,0,0],[1,1,1]]

# ╔═╡ 8177992d-9808-4761-b3b0-e4f7de057e28
m = [t p]

# ╔═╡ 58c9e2a8-504a-4657-8a04-3c0986c662af
m[:,1]

# ╔═╡ 1ed44598-64ce-4afc-aa82-f2cf4a91c0fa
vcat(m,[t p])

# ╔═╡ 7f13c77b-26b6-4d1d-8c32-b71bb3c1b723
m[1,2] = [ 2, 2, 2]

# ╔═╡ bbc8dd7d-18f5-4814-b430-a6e116ec6237
c = Array{Vector{Int64}}

# ╔═╡ e27db139-7b68-43a8-9b18-f0776aaf14b5
mat = Matrix{Vector{Int64}}(undef, 2, 2)

# ╔═╡ 8956d609-58d0-4104-a27c-2b3c554ce6a3
mat[1,2] = [ 2, 2, 2]

# ╔═╡ 8264110f-9af8-4f03-b742-dfd4a51b3bf6
mat

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
Parameters = "d96e819e-fc66-5662-9728-84c9c7592b0a"

[compat]
Parameters = "~0.12.3"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.9.1"
manifest_format = "2.0"
project_hash = "b4cda6050e9bc2a76ac8a5516fc984b386006961"

[[deps.OrderedCollections]]
git-tree-sha1 = "d321bf2de576bf25ec4d3e4360faca399afca282"
uuid = "bac558e1-5e72-5ebc-8fee-abe8a469f55d"
version = "1.6.0"

[[deps.Parameters]]
deps = ["OrderedCollections", "UnPack"]
git-tree-sha1 = "34c0e9ad262e5f7fc75b10a9952ca7692cfc5fbe"
uuid = "d96e819e-fc66-5662-9728-84c9c7592b0a"
version = "0.12.3"

[[deps.UnPack]]
git-tree-sha1 = "387c1f73762231e86e0c9c5443ce3b4a0a9a0c2b"
uuid = "3a884ed6-31ef-47d7-9d2a-63182c4928ed"
version = "1.0.2"
"""

# ╔═╡ Cell order:
# ╠═1d78060c-74d6-4425-9cc3-fec449ac1790
# ╠═06319210-1b11-11ee-0abe-3b9d74e47ab9
# ╠═bcc22983-bcf7-4ed1-9479-d67069df96e6
# ╠═9e9e8b77-79e9-4762-9830-596a29bf02b2
# ╠═8177992d-9808-4761-b3b0-e4f7de057e28
# ╠═58c9e2a8-504a-4657-8a04-3c0986c662af
# ╠═1ed44598-64ce-4afc-aa82-f2cf4a91c0fa
# ╠═7f13c77b-26b6-4d1d-8c32-b71bb3c1b723
# ╠═bbc8dd7d-18f5-4814-b430-a6e116ec6237
# ╠═e27db139-7b68-43a8-9b18-f0776aaf14b5
# ╠═8956d609-58d0-4104-a27c-2b3c554ce6a3
# ╠═8264110f-9af8-4f03-b742-dfd4a51b3bf6
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
