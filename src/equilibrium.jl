# ──────────────────────────────────────────────────────────────────────────────
# Equilibrium exact diagonalisation on a periodic ring
#
# Used for the sound velocity of Section 4.2. The ring Hamiltonian is built
# directly from bit manipulation rather than through ITensors, so the
# construction is short and easy to check on its own terms.
# ──────────────────────────────────────────────────────────────────────────────

using SparseArrays

# A basis state is an integer 1..2^N. `bits_of` turns it into a length-N vector
# of 0/1 (bit = 1 means spin down at that site); `index_of` is the inverse.
bits_of(state_index, N) = digits(state_index - 1, base=2, pad=N)
index_of(bits) = 1 + sum(bits[site] * 2^(site - 1) for site in eachindex(bits))

"""
    alcaraz_ring_hamiltonian(N, lambda, p) → SparseMatrixCSC

The model of the thesis on a periodic ring of `N` sites (site N+1 = site 1):
H = -Σ [ σᶻσᶻ₊₁ + λ σˣ + p ( σᶻσᶻ₊₂ + λ σˣσˣ₊₁ ) ].
The σᶻ terms are diagonal in this basis; the σˣ terms flip one or two bits.
"""
function alcaraz_ring_hamiltonian(N, lambda, p)
    dimension = 2^N
    rows, cols, vals = Int[], Int[], Float64[]

    for state_index in 1:dimension
        bits = bits_of(state_index, N)
        spin_at(site) = 1 - 2 * bits[site]        # bit 0 -> +1, bit 1 -> -1

        # diagonal part: the two σᶻσᶻ couplings
        diagonal_energy = 0.0
        for site in 1:N
            diagonal_energy -= spin_at(site) * spin_at(mod1(site + 1, N))
            diagonal_energy -= p * spin_at(site) * spin_at(mod1(site + 2, N))
        end
        push!(rows, state_index); push!(cols, state_index); push!(vals, diagonal_energy)

        # transverse field: flips one bit
        for site in 1:N
            flipped = copy(bits)
            flipped[site] = 1 - flipped[site]
            push!(rows, index_of(flipped)); push!(cols, state_index); push!(vals, -lambda)
        end

        # self-dual σˣσˣ term: flips two neighbouring bits
        for site in 1:N
            flipped = copy(bits)
            flipped[site] = 1 - flipped[site]
            neighbour = mod1(site + 1, N)
            flipped[neighbour] = 1 - flipped[neighbour]
            push!(rows, index_of(flipped)); push!(cols, state_index); push!(vals, -p * lambda)
        end
    end

    return sparse(rows, cols, vals, dimension, dimension)
end

"""
    cyclic_shift_operator(N) → SparseMatrixCSC

Translation by one site, needed to read off the momentum: after the shift, site
k carries whatever site k-1 carried before.
"""
function cyclic_shift_operator(N)
    dimension = 2^N
    rows, cols, vals = Int[], Int[], Float64[]
    for state_index in 1:dimension
        bits = bits_of(state_index, N)
        shifted = [bits[mod1(site - 1, N)] for site in 1:N]
        push!(rows, index_of(shifted)); push!(cols, state_index); push!(vals, 1.0)
    end
    return sparse(rows, cols, vals, dimension, dimension)
end

"""
    global_flip_operator(N) → SparseMatrixCSC

The Z2 parity P = Π σˣ, which flips every spin at once.
"""
function global_flip_operator(N)
    dimension = 2^N
    rows, cols, vals = Int[], Int[], Float64[]
    for state_index in 1:dimension
        push!(rows, index_of(1 .- bits_of(state_index, N)))
        push!(cols, state_index); push!(vals, 1.0)
    end
    return sparse(rows, cols, vals, dimension, dimension)
end

"""
    simultaneous_expectation_values(energies, vectors, operators, N; degeneracy_tol) → Matrix

Expectation value of each operator in each eigenstate. Degenerate eigenstates
come back from the diagonalisation as arbitrary mixtures, which hides their
momentum and parity labels. Within each degenerate cluster we therefore
diagonalise one generic combination of the operators: because they commute they
share an eigenbasis, and the mixture lands on it, restoring the labels.
"""
function simultaneous_expectation_values(sorted_energies, eigenvectors_matrix, operators, N;
                                         degeneracy_tol=1e-8)
    n_states, n_operators = length(sorted_energies), length(operators)
    expectation_values = zeros(ComplexF64, n_states, n_operators)

    cluster_start = 1
    while cluster_start <= n_states
        # energies come out sorted, so a degenerate cluster is a contiguous run
        cluster_end = cluster_start
        while cluster_end < n_states &&
              abs(sorted_energies[cluster_end+1] - sorted_energies[cluster_start]) < degeneracy_tol
            cluster_end += 1
        end

        cluster_vectors = eigenvectors_matrix[:, cluster_start:cluster_end]
        if cluster_end == cluster_start
            joint_basis = cluster_vectors          # non-degenerate: already clean
        else
            weights = [1.0 + 0.37im * j for j in 1:n_operators]   # any generic choice works
            combined = sum(weights[j] * (cluster_vectors' * (operators[j] * cluster_vectors))
                           for j in 1:n_operators)
            _, mixing = eigen(combined)
            joint_basis = cluster_vectors * mixing
        end

        for (local_index, global_index) in enumerate(cluster_start:cluster_end)
            state = joint_basis[:, local_index]
            for j in 1:n_operators
                expectation_values[global_index, j] = state' * (operators[j] * state)
            end
        end
        cluster_start = cluster_end + 1
    end
    return expectation_values
end

"""
    sparse_low_spectrum(N, lambda, p; howmany=10)

The `howmany` lowest eigenpairs by Lanczos, instead of the full 2^N spectrum.
"""
function sparse_low_spectrum(N, lambda, p; howmany=10)
    hamiltonian = alcaraz_ring_hamiltonian(N, lambda, p)
    energies, eigenvectors, info =
        KrylovKit.eigsolve(hamiltonian, randn(2^N), howmany, :SR; ishermitian=true)
    return energies[1:howmany], hcat(eigenvectors[1:howmany]...), info
end

"""
    ground_and_gap(N, lambda, p; howmany=10) → NamedTuple

The two levels the sound velocity needs: the lowest excitation at zero momentum
and the lowest one at the smallest non-zero momentum. Momenta come from the
phase of <T>, which is an N-th root of unity, so m = round[N/(2pi) arg<T>].
"""
function ground_and_gap(N, lambda, p; howmany=10)
    energies, eigenvector_matrix, info = sparse_low_spectrum(N, lambda, p; howmany=howmany)
    ground_state_energy = energies[1]
    excitation_energies = energies[2:end] .- ground_state_energy
    excited_eigenvectors = eigenvector_matrix[:, 2:end]

    momentum_values = simultaneous_expectation_values(excitation_energies, excited_eigenvectors,
                                                      [cyclic_shift_operator(N)], N)
    momenta = round.(Int, angle.(momentum_values[:, 1]) .* N ./ (2 * pi))

    zero_momentum = findall(==(0), momenta)
    nonzero_momentum = findall(!=(0), momenta)
    e1_zero = minimum(excitation_energies[zero_momentum])
    lowest = nonzero_momentum[argmin(abs.(momenta[nonzero_momentum]))]

    return (; N, ground_state_energy, e1_zero,
            e1_nonzero = excitation_energies[lowest],
            smallest_momentum = momenta[lowest],
            converged = info.converged)
end

# ── ground-state helpers (open chains, used by the DMRG calculations) ─────────

"""
    alcaraz_H(sites, lambda, p) → MPO

The Hamiltonian of the thesis as an MPO on an open chain, from `alcaraz_opsum`.
"""
alcaraz_H(sites, lambda, p) = MPO(alcaraz_opsum(length(sites), lambda, p), sites)

"""
    compute_vn_entropy(psi) → Vector{Float64}

Von Neumann entropy across every bond of an MPS, from the Schmidt values at that
bond. Probabilities below 1e-12 are dropped, since they only add numerical noise.
"""
function compute_vn_entropy(psi::MPS)
    N = length(psi)
    S = Float64[]
    for b in 1:(N - 1)
        orthogonalize!(psi, b)
        left_inds = b == 1 ? (siteind(psi, b),) : (linkind(psi, b - 1), siteind(psi, b))
        _, Sv, _ = svd(psi[b], left_inds)
        s = 0.0
        for n in 1:dim(Sv, 1)
            pr = Sv[n, n]^2
            pr > 1e-12 && (s -= pr * log(pr))
        end
        push!(S, s)
    end
    return S
end

# ── open chain, exact in time (used to validate the transfer-matrix column) ───

"""
    sparse_open_hamiltonian(N, lambda, p) → SparseMatrixCSC

The same Hamiltonian on an open chain of `N` sites: the sums stop at the edge
instead of wrapping around. This is the chain the Loschmidt echo is computed on
when we need an answer that owes nothing to the tensor-network machinery.
"""
function sparse_open_hamiltonian(N, lambda, p)
    dimension = 2^N
    rows, cols, vals = Int[], Int[], Float64[]

    for state_index in 1:dimension
        bits = bits_of(state_index, N)
        spin_at(site) = 1 - 2 * bits[site]

        diagonal_energy = 0.0
        for site in 1:(N - 1); diagonal_energy -= spin_at(site) * spin_at(site + 1); end
        for site in 1:(N - 2); diagonal_energy -= p * spin_at(site) * spin_at(site + 2); end
        push!(rows, state_index); push!(cols, state_index); push!(vals, diagonal_energy)

        for site in 1:N
            flipped = copy(bits); flipped[site] = 1 - flipped[site]
            push!(rows, index_of(flipped)); push!(cols, state_index); push!(vals, -lambda)
        end

        for site in 1:(N - 1)
            flipped = copy(bits)
            flipped[site] = 1 - flipped[site]
            flipped[site + 1] = 1 - flipped[site + 1]
            push!(rows, index_of(flipped)); push!(cols, state_index); push!(vals, -p * lambda)
        end
    end

    return sparse(rows, cols, vals, dimension, dimension)
end

"""
    krylov_echo_rate(N, p, T; beta0=0.2, dt=0.5) → Float64

Echo rate -log|<psi|psi(T)>| / N of the open chain, with the initial product
state |X+> cooled by exp(-beta0 H) exactly as the tensor network cools it. The
evolution is done by Krylov exponentiation in steps of `dt`, so the only error
is the Krylov tolerance: no Trotter splitting and no truncation. This is the
reference the two transfer-matrix columns are measured against.
"""
function krylov_echo_rate(N, p, T; lambda=1.0, beta0=0.2, dt=0.5)
    H = sparse_open_hamiltonian(N, lambda, p)
    psi = fill(ComplexF64(1 / sqrt(2.0^N)), 2^N)          # |X+> on every site, in the z basis
    psi, _ = KrylovKit.exponentiate(x -> H * x, -beta0, psi; tol=1e-10, ishermitian=true)

    phi = copy(psi)
    for _ in 1:round(Int, T / dt)
        phi, _ = KrylovKit.exponentiate(x -> H * x, -dt * im, phi; tol=1e-10, ishermitian=true)
    end
    return -log(abs(dot(psi, phi))) / N
end
