### A Pluto.jl notebook ###
# v0.19.31

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
end

# ╔═╡ d81f70c0-51bb-4f94-bd72-0ca5863f5b22
using PlutoUI

# ╔═╡ fcf2c92e-bc96-11ee-2232-e9d269857ecd
md"# Homework 03: Beam Task"

# ╔═╡ ef205520-e255-468b-ab5e-c8b06480ee00
# W -> Change with your number
@bind W NumberField(0:9, default=1.0)

# ╔═╡ 060330c3-fef2-4d99-9f85-7e6a802c5dd0
# X -> Change with your number
@bind X NumberField(0:9, default=5.0)

# ╔═╡ d9e97680-9d32-419d-b08c-78824cdfb886
# Y -> Change with your number
@bind Y NumberField(0:9, default=0.0)

# ╔═╡ 85601ffa-0f4e-4ec4-9e4f-d26303905264
# Z -> Change with your number
@bind Z NumberField(0:9, default=0.0)

# ╔═╡ 4187d88a-a16b-4dd4-ba28-7ca0d52374d4
md"## Geometry Data 
$\begin{gather}
a = 1 + 0.05 ⋅ X \; (m) \\
b = 4 + 0.1 ⋅ Z \; (m) \\
c = 0.1 + 0.2 ⋅ Y \; (m) \\
d = 3 + 0.1 ⋅ W \; (m) \\
t_1 = 0.75 + 0.05 ⋅ X - 0.04 ⋅ Z \; (m) \\
t_2 = \frac{t_1}{1.5}
\end{gather}$"

# ╔═╡ b3c5407b-ed14-40a2-8ae5-eedf034436bd
begin
	a = 1 + 0.05 * X 
	b = 4 + 0.1 * Z
	c = 0.1 + 0.2 * Y
	d = 3 + 0.1 * W
	t1 = 0.75 + 0.05 * X - 0.04 * Z
	t2 = t1 / 1.5 
	(
		a = a, b = b, c = c, d = d,  t₁ = t1 , t₂ = t2
	)
end

# ╔═╡ 27602c4a-a57b-47ad-846e-7d09530fdf56
md"## Material Data 
$\begin{gather}
E_b = (4 + 0.25⋅X) ⋅ 10^6 \; (KN/m^2) \\
EA_{2/4} = E_b ⋅ t_2^2 \; (KN) \\
EI_{2/4} = E_b ⋅ \frac{t_2^4}{12} \; (KN.m^2) \\ 
EA_t = (3 + 0.15 ⋅ W) ⋅ 10^5 \; (KN)
\end{gather}$"

# ╔═╡ 2be6033c-9e6a-4379-9743-3a6b8c14fbcc
begin
	Eb = (4 + 0.25 * X) * 10^6 
	EA2_4 = Eb * t2^2 
	EI2_4 = Eb * (t2^4/12)
	EA_t = (3+0.15 * W) * 10^5 
	(Eb = Eb , EA2_4 = EA2_4 , EI2_4 = EI2_4 , EA_t )
end

# ╔═╡ 9e034970-98a4-4b11-8343-5ab0e87e75ef
md"## Loading Data 
$\begin{gather}
p = 80 + 2 ⋅ Y \; (KN/m) \\
Q^* = 210 + 5 ⋅ Z \; (KN)
\end{gather}$"

# ╔═╡ 725c95f4-a30a-4acc-afe7-1f267eea89c2
begin
	p = 80 + 2 * Y 
	Q⁺ = 210 + 5 * Z 
	( p = p , Q⁺ = Q⁺)
end

# ╔═╡ 25a88f18-81c2-4299-9fa8-f39ad46ecf2e
md"## Spring Data 
$\begin{gather}
K_s = (5 + 0.5 ⋅ Y) ⋅ 10^3 \; (KN/m)
\end{gather}$"

# ╔═╡ 69051927-03da-4205-b9c0-50fcdba34c92
begin
	Ks = (5 + 0.5 * Y) * 10^3 
	(Kₛ = Ks)
end

# ╔═╡ 2f8c24a4-1d56-4733-8ab5-50a75185c203
md"## Degree of Freedoms
$\begin{gather}
\text{unkowns: } u_2, w_2, ϕ_2, w_3, \phi_3^l \rightarrow (5 \; DOKI) \\
\text{from symmetry: } u_3 = 0, \; \phi_3^r = -\phi_3^l , \; u_2 = -u_5, \; \phi_2 = -\phi_5, \; w_2 = w_5
\end{gather}$"

# ╔═╡ 2db25eec-7461-4af9-8dfa-d731eb82fdc5
md"## Element Stiffness Matrix
### **For Beam 1** :
$\begin{gather}
from: Node \; (1) \rightarrow Node \; (2) \\
DOFs: u_2,w_2,\phi_2
\end{gather}$
#### Linear Approximation:
$\begin{gather}
N_1 = \frac{1}{2} ⋅ (1 - \xi) \\
N_2 = \frac{1}{2} ⋅ (1+ \xi)
\end{gather}$
#### Area Approximation:
$\begin{gather}
A_1 = t_1^2 \\
A_2 = t_2^2 \\
A(\xi) = N_1 A_1 + N_2 A_2
\end{gather}$
#### Inertia Approximation:
$\begin{gather}
I_1 = \frac{t_1^4}{12} \\ 
I_2 = \frac{t_2^4}{12} \\
I(\xi) = N_1 I_1 + N_2 I_2 
\end{gather}$

#### Local Stiffness Matrix:
$\begin{gather}
\underline{K}^{1,l} = 
\begin{bmatrix}
 \underline{K}^u & 0 \\
	0 & \underline{K}^w 
\end{bmatrix} \\
\underline{B}^u = \frac{1}{L} \; \; \& \;\; \underline{ℂ}^u = EA(\xi) \\
\underline{K}^u = \int_{-1}^{1} \underline{B}^{u,T} ⋅ \underline{ℂ}^u ⋅ \underline{B}^u ⋅ \frac{L}{2} d\xi \\
\underline{B}^w = \begin{bmatrix} -(\frac{2}{L})^2 ⋅ \frac{d^2N_2^w}{d\xi^2}  && -(\frac{2}{L})^2 ⋅ \frac{d^2N_2^\phi}{d\xi^2}\end{bmatrix} \\
\underline{ℂ}^w = EI(\xi) \\
\underline{K}^w = \int_{-1}^{1} \underline{B}^{w,T} ⋅ \underline{ℂ}^w ⋅ \underline{B}^w ⋅ \frac{L}{2} d\xi \\
\end{gather}$

#### Global Stiffness Matrix:
$\begin{gather}
C = cos(\alpha) \\
S = sin(\alpha) \\ 
\underline{T} = 
\begin{bmatrix}
C & -S & 0 \\
S & C & 0 \\
0 & 0 & 1
\end{bmatrix} \\
\underline{K}^{1,g} = \underline{T}^T ⋅ \underline{K}^{1,l} ⋅ \underline{T}
\end{gather}$
"

# ╔═╡ fde0a499-b19b-4fb8-9bfc-0e7894c4ae25
begin
	A1 = t1^2 
	A2 = t2^2 
	I1 = (t1^4 )/12
	I2 = (t2^4)/12
	(A₁ = A1, A₂ = A2, I₁ = I1, I₂ = I2 )
end

# ╔═╡ e5e6e4ec-13db-4bad-be6f-144225e5c0aa
function A(ξ)
	N1 = 0.5 * (1 - ξ)
	N2 = 0.5 * (1+ξ)
	A = N1 * A1 +N2 * A2 
end

# ╔═╡ 83327af2-b0ad-489d-8ab9-433c665de2ec
function I(ξ)
	N1 = 0.5 * (1 - ξ)
	N2 = 0.5 * (1+ξ)
	I = N1 * I1 +N2 * I2 
end

# ╔═╡ 5569c2e2-55d6-44ef-be60-e943dd43bea4
begin
	L = b
	Bᵘ = 1/L 
	E = Eb 
	# Numerical Integeration
	α = 2 
	ξ = 0
	Kᵘ = α * (Bᵘ * E * A(ξ) * Bᵘ) * (L/2)
	(Kᵘ = Kᵘ, L)
end

# ╔═╡ 3886c3ac-87f7-4604-aac5-b3c1f8191f79
function Bʷ(ξ)
	B = [(2/L)^2 * (3/2) * ξ  (2/L)^2 * (L/4) * (1+ 3* ξ) ]
end

# ╔═╡ d8368947-46e4-48cc-b286-c94451263877
begin
	# Numerical Integration
	α1 = 1 
	ξ1 = -1/sqrt(3)
	α2 = 1
	ξ2 = 1 / sqrt(3)
	Cʷ = (ξ) -> E * I(ξ)
	Kʷ = α1 * transpose(Bʷ(ξ1)) * Cʷ(ξ1) * Bʷ(ξ1) * (L/2) + α2 * transpose(Bʷ(ξ2)) * Cʷ(ξ2) * Bʷ(ξ2) * (L/2)
	(Kʷ = Kʷ , L)
end

# ╔═╡ 7b996bfc-3f79-4994-b39e-fb298835d4bb
begin
	K_1_local = [Kᵘ 0 0;0 Kʷ[1,1] Kʷ[1,2];0 Kʷ[2,1] Kʷ[2,2]]
	node1 = (-d,a+b)
	node2 = (-d,a)
	C1 = (node2[1] - node1[1])/b
	C2 = (node2[2] - node1[2])/b
	T1 = [C1 C2 0;
		-C2 C1 0;
		 0 0 1]
	K¹ = transpose(T1) * K_1_local * T1
	(K¹ = K¹ , T¹ = T1)
end

# ╔═╡ f8c25c03-8190-4b90-b108-124b937070cb
md"### **For Beam 2** :
$\begin{gather}
from: Node \; (2) \rightarrow Node \; (3) \\
DOFs: u_2,w_2,\phi_2, w_3, \phi_3^l
\end{gather}$
#### Nodal Coordinates:
$\begin{gather}
Node (2): (-d,a) \\
Node (3): (0,0) \\
α = tan^{-1}(a/d) \\
c = cos(α) \\
s = sin(α) \\ 
L = \sqrt{(d)^2 + (a)^2}
\end{gather}$

#### Global Stiffness Matrix:

$\begin{gather}
\underline{K}^{2,g} = 
\begin{bmatrix}
\frac{EAc^2}{L} + \frac{12EIs^2}{L^3} & \frac{-EAcs}{L} + \frac{12EIcs}{L^3} & \frac{-6EIs}{L^2} & \frac{EAcs}{L} - \frac{12EIcs}{L^3} & \frac{-6 EI s}{L^2} \\
\frac{-EAcs}{L} + \frac{12EIcs}{L^3} & \frac{EAs^2}{L} + \frac{12EIc^2}{L^3} & \frac{-6EIc}{L^2} & \frac{-EAs^2}{L} - \frac{12EIc^2}{L^3} & \frac{-6 EI c}{L^2} \\
\frac{-6EIs}{L^2} & \frac{-6EIc}{L^2} & \frac{4EI}{L} & \frac{6EIc}{L^2} & \frac{2EI}{L} \\
\frac{EAcs}{L} - \frac{12EIcs}{L^3} & \frac{-EAs^2}{L} - \frac{12EIc^2}{L^3} & \frac{6EIc}{L^2} & \frac{EAs^2}{L} + \frac{12EIc^2}{L^3} & \frac{6 EI c}{L^2}\\
\frac{-6 EI s}{L^2} & \frac{-6 EI c}{L^2} & \frac{2EI}{L} & \frac{6 EI c}{L^2} & \frac{4EI}{L}
\end{bmatrix}
\end{gather}$

"

# ╔═╡ d90f4d01-fb77-445a-97e1-e747db9613c2
begin
	node3 = (0,0)
	L2 = sqrt((node3[1] - node2[1])^2 + (node3[2] - node2[2])^2)
	C1_2 = (node3[1] - node2[1]) / L2
	C2_2 = (node3[2] - node2[2]) / L2
	T2 = [C1_2 C2_2 0 0 0 0;
		 -C2_2 C1_2 0 0 0 0;
			0    0  1 0 0 0;
			0    0 0  C1_2 C2_2 0;
			0    0  0  -C2_2 C1_2 0;
			0     0  0   0  0  1]
	
	EI = EI2_4
	EA = EA2_4
	K2_local = [EA/L2 0 0 -EA/L2 0 0;
				0 12 * EI/L2^3 -6 * EI/L2^2 0 -12 * EI/L2^3 -6*EI/L2^2;
				0 -6 * EI/L2^2 4 * EI / L2 0 6 * EI / L2^2 2*EI/L2;
				-EA/L2 0 0 EA/L2 0 0;
				0 -12 * EI/L2^3 6 * EI/L2^2 0 12 * EI / L2^3 6 * EI / L2^2;
				0 -6 * EI/L2^2 2*EI/L2 0 6 * EI/L2^2 4 * EI/L2]
	K² = transpose(T2) * K2_local * T2
	K² = K²[1:end .!= 4  ,1:end .!= 4]
	(K² = K² , L2)
	
end

# ╔═╡ 891eb96b-b563-46f7-9738-16b64caea5bd
md"### **For Truss 3** :
$\begin{gather}
from: Node \; (3) \rightarrow Node \; (4) \\
DOFs: w_3
\end{gather}$
#### Nodal Coordinates:
$\begin{gather}
Node (3): (0,0) \\
Node (4): (0,a + b + c) \\
\end{gather}$"

# ╔═╡ d8bd602a-b6a5-41ad-82c5-5de48469be19
begin
	node4 = (0,a+b+c)
	L3 = a+b+c
	C1_3 = (node4[1] - node3[1])/L3
	C2_3 = (node4[2] - node3[2])/L3

	T3 = [C1_3 C2_3;-C2_3 C1_3]
	K³ = transpose(T3) * (EA_t / (2 * L3)) * T3
	K³ = K³[2,2]
	(K³ = K³ , L3)
end

# ╔═╡ e5f071b8-11cb-4620-b69b-33e0a3cc61f3
md"### **For Spring** :
$\begin{gather}
from: Node \; (1) \rightarrow Node \; (2) \\
DOFs: u_2 , w_2 , \phi_2
\end{gather}$
#### Local Stiffness Matrix:
$\begin{gather}
\underline{K}^s = \underline{B}^T * K_s * \underline{B} \\
w_s = N_1^w w_1 + N_1^ϕ \phi_1 +N_2^w w_2 + N_2^ϕ \phi_2 
\end{gather}$"

# ╔═╡ 4b9bac2d-ca94-4789-8007-55a837a4052b
function xToxi(x)
	x2 = b
	ξ = (x/(0.5 * b)) - 1 
end

# ╔═╡ 9254cc5d-9ea5-4e39-a779-ebed9ad97d47
begin
	ξs = xToxi(b/3)
	L1 = L
	Bs = [0 0.25 * (1-ξs)^2 * (2+ξs) -(L1/8) * (1-ξs)^2 * (1+ξs)]
	Kˢ = transpose(T1) * (transpose(Bs) * Ks * Bs) * T1
	Kˢ = Kˢ[1:end .!= 2, 1:end .!= 2]
end

# ╔═╡ 3d4278b4-6b65-405a-9d16-a3154e3ec3b9
md"## Overall (Assembly) Stiffness Matrix
```math
\begin{align}
u_2 \; \; \;w_2 \; \; \; \phi_2 \;\;  w_3 \; \; \; \phi_3^l\; \; \; \; \; \; \;\\
\underline{K} = 
\begin{bmatrix}
- & - & - & - & - \\
- & - & - & - & -\\
- & - & - & - & - \\
- & - & - & - & -\\
- & - & - & - & -
\end{bmatrix}
\begin{matrix}
\delta u_2 \\ \delta w_2 \\ \delta \phi_2 \\ \delta w_3 \\ \delta \phi_3^l
\end{matrix}
\end{align}
```
"

# ╔═╡ 5456067e-b957-4037-8c09-d96e119624ae
begin

	# Beam 1: u2 , w2 , ϕ2
		K¹_all = [K¹ zeros(3) zeros(3);0 0 0 0 0 ; 0 0 0 0 0] 	
	# Beam2: u2, w2, ϕ2 , w3
	K²_all = K²
	# truss: w3 
	K³_all = zeros(5,5)
	K³_all[4,4] = K³
	K³_all
	# spring; u2, ϕ2
	Kˢ_all =[Kˢ[1,1] 0 Kˢ[1,2] 0 0;
			 0 0 0 0 0;
			 Kˢ[2,1] 0 Kˢ[2,2] 0 0;
			0 0 0 0 0;
	0 0 0 0 0 ]
	K = K¹_all + K²_all + K³_all + Kˢ_all
end

# ╔═╡ b0bcc6e0-82fd-4da4-93fe-37ac82e49bd5
md"## Load Vector
#### Distributed Load: 
```math
\begin{gather}
 \underline{r}^2_u = \int_{-1}^{+1} \underline{N}^T ⋅ \underline{b} ⋅ \frac{L}{2} d\xi \\
 \underline{N}^T = 
\begin{bmatrix}
N_2^u & 0 \\
0 & N_2^w \\
0 & N_2^ϕ \\
N_3^u & 0 \\
0 & N_3^w \\
0 & N_3^ϕ \\
\end{bmatrix}\\
\underline{b} = 
\begin{bmatrix}
 0 \\ q
\end{bmatrix}

\end{gather}
```
#### Point Load: 

```math
\begin{gather}
δW_{ext} = δw(\xi^*) ⋅ Q^* = \underline{δu}^{e,T} ⋅ \underline{B}^T ⋅ Q^* \\
\underline{B}^T(\xi^*) = 
\begin{bmatrix}
N_2^w & N_2^ϕ & N_3^w & N_3^ϕ
\end{bmatrix} \\
\underline{r}_p = \underline{B}^T ⋅ Q^*
\end{gather}
```
#### Total Load & Coordinate Transformation:

```math
\begin{gather}
\underline{r}^l = \underline{r}_u + \underline{r}_p \\
\underline{r}^g = \underline{T}^T ⋅ \underline{r}^l
\end{gather}
```
"

# ╔═╡ 980398f5-06c8-43c3-b310-84585271f01a
function loadNᵀ(ξ)
	N2_u = 0.5 *(1- ξ)
	N2_w = 0.25 * (1-ξ)^2 * (2+ξ)
	N2_ϕ = (-L2/8) * (1-ξ)^2 * (1+ξ)

	N3_u = 0.5 *(1+ ξ)
	N3_w = 0.25 * (1+ξ)^2 * (2-ξ)
	N3_ϕ = (L2/8) * (1+ξ)^2 * (1-ξ)

	Nᵀ = [N2_u 0 ;
		   0 N2_w;
		   0  N2_ϕ;
	       N3_u  0;
	        0   N3_w;
	        0    N3_ϕ]
end

# ╔═╡ 61c51495-e812-4984-b279-39c6efea112c
begin
	# Numerical Integration
	α11 = 1
	ξ11 = (-1/sqrt(3))
	α22 = 1
	ξ22 = (1/sqrt(3))
	load_b = [0; p]
	ru = α11 * loadNᵀ(ξ11) * load_b * (L2/2) +α22 *  loadNᵀ(ξ22)  * load_b * (L2/2)
	Bᵀ = (ξ) -> loadNᵀ(ξ)[:,2]
	rp =Bᵀ(0) * Q⁺
	(rᵤ = ru , rₚ = rp)
end

# ╔═╡ 45a20abc-5001-4b57-9ba1-39e2d2cd7f30
begin
	r_total = ru + rp
	rg = transpose(T2) * r_total
	r = [rg[1:3]; rg[4:5]]
	(r = r, T2)
end

# ╔═╡ 28e6021b-c936-43d6-a2a5-1a509fa12464
md"## Solution
```math
\begin{gather}
 \underline{K} ⋅ \underline{u} = \underline{r} \rightarrow \underline{u} = \underline{K}^{-1} ⋅ \underline{r}
\end{gather}
```
"

# ╔═╡ 41797e54-8956-45d0-a9a2-f9b399e813b8
u=inv(K) * r

# ╔═╡ 5e881f53-2d12-436b-b1a3-c80466beadf4
md"## Post Processing
```math
\begin{gather}
 ξ^m = \frac{-1}{2+Z} \\
κ = \frac{-d^2w}{dx^2} = \underline{B} ⋅ \underline{u}^e \\
κ = - (\frac{2}{L}) ⋅ 
\begin{bmatrix}
\frac{d^2N_2^w}{dξ^2} & \frac{d^2N_2^ϕ}{dξ^2} & \frac{d^2N_3^w}{dξ^2} & \frac{d^2N_3^ϕ}{dξ^2} 
\end{bmatrix} \\
M = EIκ
\end{gather}
```
"

# ╔═╡ 91d9ab2a-ad4a-44c7-be93-eea84da6d276
begin
	ξm =-1/(2+Z)
	twoLsq = - (2/L2)^2 
	M =EI * twoLsq * [1.5 * ξm (L2/4) * (1-3*ξm) -1.5 * ξm (-L2/4) * (1+3*ξm)] * u[2:5]
end

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
PlutoUI = "7f904dfe-b85e-4ff6-b463-dae2292396a8"

[compat]
PlutoUI = "~0.7.55"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.9.3"
manifest_format = "2.0"
project_hash = "f64cdffc70331b0a2f407efefd54fd84eb680773"

[[deps.AbstractPlutoDingetjes]]
deps = ["Pkg"]
git-tree-sha1 = "c278dfab760520b8bb7e9511b968bf4ba38b7acc"
uuid = "6e696c72-6542-2067-7265-42206c756150"
version = "1.2.3"

[[deps.ArgTools]]
uuid = "0dad84c5-d112-42e6-8d28-ef12dabb789f"
version = "1.1.1"

[[deps.Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"

[[deps.Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"

[[deps.ColorTypes]]
deps = ["FixedPointNumbers", "Random"]
git-tree-sha1 = "eb7f0f8307f71fac7c606984ea5fb2817275d6e4"
uuid = "3da002f7-5984-5a60-b8a6-cbb66c0b333f"
version = "0.11.4"

[[deps.CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"
version = "1.0.5+0"

[[deps.Dates]]
deps = ["Printf"]
uuid = "ade2ca70-3891-5945-98fb-dc099432e06a"

[[deps.Downloads]]
deps = ["ArgTools", "FileWatching", "LibCURL", "NetworkOptions"]
uuid = "f43a241f-c20a-4ad4-852c-f6b1247861c6"
version = "1.6.0"

[[deps.FileWatching]]
uuid = "7b1f6079-737a-58dc-b8bc-7a2ca5c1b5ee"

[[deps.FixedPointNumbers]]
deps = ["Statistics"]
git-tree-sha1 = "335bfdceacc84c5cdf16aadc768aa5ddfc5383cc"
uuid = "53c48c17-4a7d-5ca2-90c5-79b7896eea93"
version = "0.8.4"

[[deps.Hyperscript]]
deps = ["Test"]
git-tree-sha1 = "179267cfa5e712760cd43dcae385d7ea90cc25a4"
uuid = "47d2ed2b-36de-50cf-bf87-49c2cf4b8b91"
version = "0.0.5"

[[deps.HypertextLiteral]]
deps = ["Tricks"]
git-tree-sha1 = "7134810b1afce04bbc1045ca1985fbe81ce17653"
uuid = "ac1192a8-f4b3-4bfe-ba22-af5b92cd3ab2"
version = "0.9.5"

[[deps.IOCapture]]
deps = ["Logging", "Random"]
git-tree-sha1 = "8b72179abc660bfab5e28472e019392b97d0985c"
uuid = "b5f81e59-6552-4d32-b1f0-c071b021bf89"
version = "0.2.4"

[[deps.InteractiveUtils]]
deps = ["Markdown"]
uuid = "b77e0a4c-d291-57a0-90e8-8db25a27a240"

[[deps.JSON]]
deps = ["Dates", "Mmap", "Parsers", "Unicode"]
git-tree-sha1 = "31e996f0a15c7b280ba9f76636b3ff9e2ae58c9a"
uuid = "682c06a0-de6a-54ab-a142-c8b1cf79cde6"
version = "0.21.4"

[[deps.LibCURL]]
deps = ["LibCURL_jll", "MozillaCACerts_jll"]
uuid = "b27032c2-a3e7-50c8-80cd-2d36dbcbfd21"
version = "0.6.3"

[[deps.LibCURL_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "MbedTLS_jll", "Zlib_jll", "nghttp2_jll"]
uuid = "deac9b47-8bc7-5906-a0fe-35ac56dc84c0"
version = "7.84.0+0"

[[deps.LibGit2]]
deps = ["Base64", "NetworkOptions", "Printf", "SHA"]
uuid = "76f85450-5226-5b5a-8eaa-529ad045b433"

[[deps.LibSSH2_jll]]
deps = ["Artifacts", "Libdl", "MbedTLS_jll"]
uuid = "29816b5a-b9ab-546f-933c-edad1886dfa8"
version = "1.10.2+0"

[[deps.Libdl]]
uuid = "8f399da3-3557-5675-b5ff-fb832c97cbdb"

[[deps.LinearAlgebra]]
deps = ["Libdl", "OpenBLAS_jll", "libblastrampoline_jll"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"

[[deps.Logging]]
uuid = "56ddb016-857b-54e1-b83d-db4d58db5568"

[[deps.MIMEs]]
git-tree-sha1 = "65f28ad4b594aebe22157d6fac869786a255b7eb"
uuid = "6c6e2e6c-3030-632d-7369-2d6c69616d65"
version = "0.1.4"

[[deps.Markdown]]
deps = ["Base64"]
uuid = "d6f4376e-aef5-505a-96c1-9c027394607a"

[[deps.MbedTLS_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "c8ffd9c3-330d-5841-b78e-0817d7145fa1"
version = "2.28.2+0"

[[deps.Mmap]]
uuid = "a63ad114-7e13-5084-954f-fe012c677804"

[[deps.MozillaCACerts_jll]]
uuid = "14a3606d-f60d-562e-9121-12d972cd8159"
version = "2022.10.11"

[[deps.NetworkOptions]]
uuid = "ca575930-c2e3-43a9-ace4-1e988b2c1908"
version = "1.2.0"

[[deps.OpenBLAS_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "4536629a-c528-5b80-bd46-f80d51c5b363"
version = "0.3.21+4"

[[deps.Parsers]]
deps = ["Dates", "PrecompileTools", "UUIDs"]
git-tree-sha1 = "8489905bcdbcfac64d1daa51ca07c0d8f0283821"
uuid = "69de0a69-1ddd-5017-9359-2bf0b02dc9f0"
version = "2.8.1"

[[deps.Pkg]]
deps = ["Artifacts", "Dates", "Downloads", "FileWatching", "LibGit2", "Libdl", "Logging", "Markdown", "Printf", "REPL", "Random", "SHA", "Serialization", "TOML", "Tar", "UUIDs", "p7zip_jll"]
uuid = "44cfe95a-1eb2-52ea-b672-e2afdf69b78f"
version = "1.9.2"

[[deps.PlutoUI]]
deps = ["AbstractPlutoDingetjes", "Base64", "ColorTypes", "Dates", "FixedPointNumbers", "Hyperscript", "HypertextLiteral", "IOCapture", "InteractiveUtils", "JSON", "Logging", "MIMEs", "Markdown", "Random", "Reexport", "URIs", "UUIDs"]
git-tree-sha1 = "68723afdb616445c6caaef6255067a8339f91325"
uuid = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
version = "0.7.55"

[[deps.PrecompileTools]]
deps = ["Preferences"]
git-tree-sha1 = "03b4c25b43cb84cee5c90aa9b5ea0a78fd848d2f"
uuid = "aea7be01-6a6a-4083-8856-8a6e6704d82a"
version = "1.2.0"

[[deps.Preferences]]
deps = ["TOML"]
git-tree-sha1 = "00805cd429dcb4870060ff49ef443486c262e38e"
uuid = "21216c6a-2e73-6563-6e65-726566657250"
version = "1.4.1"

[[deps.Printf]]
deps = ["Unicode"]
uuid = "de0858da-6303-5e67-8744-51eddeeeb8d7"

[[deps.REPL]]
deps = ["InteractiveUtils", "Markdown", "Sockets", "Unicode"]
uuid = "3fa0cd96-eef1-5676-8a61-b3b8758bbffb"

[[deps.Random]]
deps = ["SHA", "Serialization"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"

[[deps.Reexport]]
git-tree-sha1 = "45e428421666073eab6f2da5c9d310d99bb12f9b"
uuid = "189a3867-3050-52da-a836-e630ba90ab69"
version = "1.2.2"

[[deps.SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"
version = "0.7.0"

[[deps.Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"

[[deps.Sockets]]
uuid = "6462fe0b-24de-5631-8697-dd941f90decc"

[[deps.SparseArrays]]
deps = ["Libdl", "LinearAlgebra", "Random", "Serialization", "SuiteSparse_jll"]
uuid = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"

[[deps.Statistics]]
deps = ["LinearAlgebra", "SparseArrays"]
uuid = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"
version = "1.9.0"

[[deps.SuiteSparse_jll]]
deps = ["Artifacts", "Libdl", "Pkg", "libblastrampoline_jll"]
uuid = "bea87d4a-7f5b-5778-9afe-8cc45184846c"
version = "5.10.1+6"

[[deps.TOML]]
deps = ["Dates"]
uuid = "fa267f1f-6049-4f14-aa54-33bafae1ed76"
version = "1.0.3"

[[deps.Tar]]
deps = ["ArgTools", "SHA"]
uuid = "a4e569a6-e804-4fa4-b0f3-eef7a1d5b13e"
version = "1.10.0"

[[deps.Test]]
deps = ["InteractiveUtils", "Logging", "Random", "Serialization"]
uuid = "8dfed614-e22c-5e08-85e1-65c5234f0b40"

[[deps.Tricks]]
git-tree-sha1 = "eae1bb484cd63b36999ee58be2de6c178105112f"
uuid = "410a4b4d-49e4-4fbc-ab6d-cb71b17b3775"
version = "0.1.8"

[[deps.URIs]]
git-tree-sha1 = "67db6cc7b3821e19ebe75791a9dd19c9b1188f2b"
uuid = "5c2747f8-b7ea-4ff2-ba2e-563bfd36b1d4"
version = "1.5.1"

[[deps.UUIDs]]
deps = ["Random", "SHA"]
uuid = "cf7118a7-6976-5b1a-9a39-7adc72f591a4"

[[deps.Unicode]]
uuid = "4ec0a83e-493e-50e2-b9ac-8f72acf5a8f5"

[[deps.Zlib_jll]]
deps = ["Libdl"]
uuid = "83775a58-1f1d-513f-b197-d71354ab007a"
version = "1.2.13+0"

[[deps.libblastrampoline_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850b90-86db-534c-a0d3-1478176c7d93"
version = "5.8.0+0"

[[deps.nghttp2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850ede-7688-5339-a07c-302acd2aaf8d"
version = "1.48.0+0"

[[deps.p7zip_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "3f19e933-33d8-53b3-aaab-bd5110c3b7a0"
version = "17.4.0+0"
"""

# ╔═╡ Cell order:
# ╠═d81f70c0-51bb-4f94-bd72-0ca5863f5b22
# ╟─fcf2c92e-bc96-11ee-2232-e9d269857ecd
# ╠═ef205520-e255-468b-ab5e-c8b06480ee00
# ╠═060330c3-fef2-4d99-9f85-7e6a802c5dd0
# ╠═d9e97680-9d32-419d-b08c-78824cdfb886
# ╠═85601ffa-0f4e-4ec4-9e4f-d26303905264
# ╟─4187d88a-a16b-4dd4-ba28-7ca0d52374d4
# ╠═b3c5407b-ed14-40a2-8ae5-eedf034436bd
# ╟─27602c4a-a57b-47ad-846e-7d09530fdf56
# ╠═2be6033c-9e6a-4379-9743-3a6b8c14fbcc
# ╟─9e034970-98a4-4b11-8343-5ab0e87e75ef
# ╟─725c95f4-a30a-4acc-afe7-1f267eea89c2
# ╟─25a88f18-81c2-4299-9fa8-f39ad46ecf2e
# ╠═69051927-03da-4205-b9c0-50fcdba34c92
# ╟─2f8c24a4-1d56-4733-8ab5-50a75185c203
# ╟─2db25eec-7461-4af9-8dfa-d731eb82fdc5
# ╠═fde0a499-b19b-4fb8-9bfc-0e7894c4ae25
# ╠═e5e6e4ec-13db-4bad-be6f-144225e5c0aa
# ╠═83327af2-b0ad-489d-8ab9-433c665de2ec
# ╠═5569c2e2-55d6-44ef-be60-e943dd43bea4
# ╠═3886c3ac-87f7-4604-aac5-b3c1f8191f79
# ╠═d8368947-46e4-48cc-b286-c94451263877
# ╠═7b996bfc-3f79-4994-b39e-fb298835d4bb
# ╟─f8c25c03-8190-4b90-b108-124b937070cb
# ╠═d90f4d01-fb77-445a-97e1-e747db9613c2
# ╟─891eb96b-b563-46f7-9738-16b64caea5bd
# ╠═d8bd602a-b6a5-41ad-82c5-5de48469be19
# ╟─e5f071b8-11cb-4620-b69b-33e0a3cc61f3
# ╠═4b9bac2d-ca94-4789-8007-55a837a4052b
# ╠═9254cc5d-9ea5-4e39-a779-ebed9ad97d47
# ╟─3d4278b4-6b65-405a-9d16-a3154e3ec3b9
# ╠═5456067e-b957-4037-8c09-d96e119624ae
# ╟─b0bcc6e0-82fd-4da4-93fe-37ac82e49bd5
# ╠═980398f5-06c8-43c3-b310-84585271f01a
# ╠═61c51495-e812-4984-b279-39c6efea112c
# ╠═45a20abc-5001-4b57-9ba1-39e2d2cd7f30
# ╟─28e6021b-c936-43d6-a2a5-1a509fa12464
# ╠═41797e54-8956-45d0-a9a2-f9b399e813b8
# ╟─5e881f53-2d12-436b-b1a3-c80466beadf4
# ╠═91d9ab2a-ad4a-44c7-be93-eea84da6d276
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
