### A Pluto.jl notebook ###
# v0.19.31

using Markdown
using InteractiveUtils

# ╔═╡ fc6f163f-a7c3-4482-86d7-a2181f9492c0
md"### Data Input 
$\begin{gather}
\text{Only change the values of W, X, Y, Z to solve the task}
\end{gather}$"

# ╔═╡ 0f4138b0-a004-11ee-2141-75c9938539e1
begin
	# Please fill the following data
	W = 1.0 #Change with your number
	X = 5.0 #Change with you number
	Y = 0.0 #Change with your number
	Z = 0.0 #Change with your number
end

# ╔═╡ cf631dce-a4b1-4c1d-a1aa-ada9298990a7
md"### Data Given Calculations"

# ╔═╡ 02f1bff8-26ac-4a6c-af75-c83f527134ae
begin
	# Geometry
	a = 4 + ((X + Z)/10)
	b = 0.5 + ((W+Z)/20)
	c = 1 + ((X + Y)/20)
	h = 1 

	# Material
	E = (45 + W) * 10^3 
	ν = 0.15 + 0.01 * Y
	λ = (1.33 + 0.06 * Z) 
	αθ = (1 + (X/5)) * 10^-6

	# Loading
	Pₘₐₓ = a * 0.997 * 9.81
	θ₁=35
	θ₂ = (40 + Y)

	# Other
	θref = (25 - W)

	(
		a =   a , b = b , c = c , h = h,
		E = E, ν = ν , λ = λ , αθ = αθ,
		Pₘₐₓ = Pₘₐₓ , θ₁ = θ₁ , θ₂ = θ₂,
		θref = θref
	)
end

# ╔═╡ 51cc2bee-8a91-497d-92b8-7b175db4eb56
md"### Constitutive Laws"

# ╔═╡ 25330f13-c226-4b77-b03a-d6e0090cae29
begin
	ℂ = (E/((1+ν)*(1-2ν))) * [1-ν ν 0; ν 1-ν 0; 0 0 (1-2ν)/2]
	Λ = [λ 0; 0 λ]
	I = [1;1;0]
	(ℂ = ℂ , Λ = Λ , I = I)
end

# ╔═╡ 5e8132ac-0900-42fd-9762-64ff5063a8e4
md"### Nodal Coordinates 
$\begin{gather}
\begin{matrix} 
node & x_1 & x_2  \\
1 & b + c & 0 \\
2 & b & a \\
3 & 0 & 0
\end{matrix}
\end{gather}$"

# ╔═╡ f640b6cd-cfe2-4955-a060-d6bd789a086a
xᵉ = [b+c 0;b a; 0 0]

# ╔═╡ 43ecc6bb-72c9-4abe-b90c-beedf19a54c3
begin
	x1 = b+c 
	y1 = 0
	x2 = b
	y2 = a
	x3 = 0
	y3 = 0
end

# ╔═╡ 8dc5323a-37a3-4403-9275-841cddf4d981
A = 0.5 * abs(x1 * (y2 - y3) + x2 * (y3 - y1) + x3 * (y1 - y2))

# ╔═╡ 317583b2-d024-47f6-a02c-fb54cbd3dcd4
md"### For $\underline{K}^{global}_{{\theta}{\theta}}$
$\begin{gather}
\underline{K}^{global}_{{\theta}{\theta}} = \underline{B}_{θ}^T.\underline{λ}.\underline{B}_{θ}.A.h
\end{gather}$"

# ╔═╡ 0bca8691-12d7-45a3-86d4-822196808cc0
Bθ = (1/(2A)) *[y2 - y3 y3 - y1 y1 - y2; x3 - x2 x1 - x3 x2 - x1]

# ╔═╡ 28b32f7d-396e-418c-a11e-e4891fe2db39
Kθθ = transpose(Bθ) * Λ * Bθ * A * h

# ╔═╡ e52ed72c-9d10-4f40-be95-67e022508ce5
md"### For $\underline{K}^{global}_{uu}$
$\begin{gather}
\underline{K}^{global}_{uu} = \underline{B}_{ε}^T.\underline{ℂ}.\underline{B}_{ε}.A.h
\end{gather}$"

# ╔═╡ 89d66f50-1856-4a48-8e3d-6f09e191d9b1
# reduced B - Matrix
Bε = (1/(2*A))*[y2 - y3 y3-y1 0;0 0 x1 - x3; x3 - x2 x1 - x3 y3 - y1]

# ╔═╡ 8d677543-a20b-4ae9-bd0e-7b0575f09db1
Kuu = transpose(Bε) * ℂ * Bε * A * h

# ╔═╡ 0fee7b4e-aef3-4bec-9dd6-89815851bcc7
md"### For $\underline{K}^{global}_{u{\theta}}$
$\begin{gather}
\underline{K}^{global}_{u{\theta}} = \frac{Ah}{3} . \underline{B}_{ε}^T.\underline{ℂ}. α_{θ} . \underline{I} . \begin{matrix} [1 & 1 & 1] \end{matrix}
\end{gather}$"

# ╔═╡ c0cf7d2b-7ffc-477d-bb30-527a3d088910
Kuθ = ((A*h)/3) * transpose(Bε) * ℂ * αθ * I * [1 1 1]

# ╔═╡ 4fd00a4b-baf5-4d4f-b964-b5d1235830f5
md"### For $\underline{r}^{global}_{red}$
$\begin{gather}
\underline{r}^{global}_{red} = 
\begin{bmatrix} 
0 \\
\int_{-1}^{+1} \frac{1}{2} . (1 + \xi_2) . \frac{P_{max}}{2} . (1 - ξ_2).h.|\underline{J} |  \,d\xi_2  \\
0
\end{bmatrix}
\begin{matrix}
\scriptsize{δu_{11}} \\
\scriptsize{δu_{21}} \\
\scriptsize{δu_{22}}
\end{matrix}
\end{gather}$"

# ╔═╡ 18081f08-5c77-484b-974b-bf56a3ea13e9
begin
	#Numerical Integeration for r₂₁
	Lᵉ = sqrt(a^2+b^2)
	detJ = Lᵉ/ 2
	gp1 = (-1/sqrt(3))
	gp2 = (1/sqrt(3))
	α₁ = 1
	α₂ = 1
	r₂₁ = α₁ * 0.5 * (1+gp1) * (Pₘₐₓ/2) * (1-gp1) * h * detJ + α₂ * 0.5 * (1+gp2) * (Pₘₐₓ/2) * (1- gp2) * h * detJ 
	r = [0;r₂₁;0]
end

# ╔═╡ 851093ae-eaf1-4229-bdd8-e295712a5ef0
md"### Solve the system
$\begin{gather}
\begin{bmatrix}
\underline{K}^{global}_{uu} & -\underline{K}^{global}_{u\theta} \\
\underline{0} & \underline{K}^{global}_{θθ}
\end{bmatrix} .
\begin{bmatrix}
u^e \\ θ^e
\end{bmatrix} = 
\begin{bmatrix}
\underline{R}^u - \underline{K}^{global}_{uθ} . θ^{ref,e} \\
\underline{R}^T
\end{bmatrix}
\end{gather}$"

# ╔═╡ 311e4fa5-98f9-42e2-b039-db7082ef2154
begin
	rupper = r - Kuθ * [θref; θref; θref]
	zero  = zeros(3,3)
	Ktotal = [Kuu -Kuθ;zero Kθθ]
	rtotal = [rupper;0;0;0]
	# prescriped temperature
	θᵣ = [θ₁; θ₂]
	# load vector corresponding to unknown DOFs
	rᵤ  = [rtotal[1:3];rtotal[6]]
	# stiffnes matrix of unknown DOFs
	Kᵤ = [Kuu Ktotal[1:3,6];  transpose(Ktotal[6,1:3])  Ktotal[6,6]]
	# stiffness matrix of prescribed DOFs
	Kᵣ = [Ktotal[1:3,4:5] ; transpose(Ktotal[6,4:5])]
	u = inv(Kᵤ) * (rᵤ - Kᵣ * θᵣ)
end

# ╔═╡ Cell order:
# ╟─fc6f163f-a7c3-4482-86d7-a2181f9492c0
# ╠═0f4138b0-a004-11ee-2141-75c9938539e1
# ╟─cf631dce-a4b1-4c1d-a1aa-ada9298990a7
# ╠═02f1bff8-26ac-4a6c-af75-c83f527134ae
# ╟─51cc2bee-8a91-497d-92b8-7b175db4eb56
# ╠═25330f13-c226-4b77-b03a-d6e0090cae29
# ╟─5e8132ac-0900-42fd-9762-64ff5063a8e4
# ╠═f640b6cd-cfe2-4955-a060-d6bd789a086a
# ╠═43ecc6bb-72c9-4abe-b90c-beedf19a54c3
# ╠═8dc5323a-37a3-4403-9275-841cddf4d981
# ╟─317583b2-d024-47f6-a02c-fb54cbd3dcd4
# ╠═0bca8691-12d7-45a3-86d4-822196808cc0
# ╠═28b32f7d-396e-418c-a11e-e4891fe2db39
# ╟─e52ed72c-9d10-4f40-be95-67e022508ce5
# ╠═89d66f50-1856-4a48-8e3d-6f09e191d9b1
# ╠═8d677543-a20b-4ae9-bd0e-7b0575f09db1
# ╟─0fee7b4e-aef3-4bec-9dd6-89815851bcc7
# ╠═c0cf7d2b-7ffc-477d-bb30-527a3d088910
# ╟─4fd00a4b-baf5-4d4f-b964-b5d1235830f5
# ╠═18081f08-5c77-484b-974b-bf56a3ea13e9
# ╟─851093ae-eaf1-4229-bdd8-e295712a5ef0
# ╠═311e4fa5-98f9-42e2-b039-db7082ef2154
