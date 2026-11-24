"""Standalone algebra check for a proposed NICAM sedimentation flux factorization.

Run with python3 remap_probe.py. Uses nondimensional synthetic columns and the
standard library. Departure distances are inputs. No NICAM data or GPU needed.
Flux is positive downwards and integrated over one time step.
"""

from bisect import bisect_left
from collections import Counter
from fractions import Fraction
from itertools import accumulate
from math import isclose, nextafter
from random import Random


def reference_flux(walls, face, falls, densities, values, heat):
    """Species-first piecewise-constant integral using the remaining distance."""
    mass = [0] * len(falls)
    moments = [0] * len(values)
    for sp, distance in enumerate(falls):
        target = walls[face] + distance
        # The archived Triton reference retains cell = face on top overshoot.
        last = face if target > walls[-1] else max(
            face, bisect_left(walls, target, lo=face + 1) - 1
        )
        remaining = distance
        for j in range(face, last + 1):
            length = walls[j + 1] - walls[j] if j < last else remaining
            remaining -= length
            carried = length * densities[sp][j]
            mass[sp] += carried
            for p, field in enumerate(values):
                moments[p] += carried * field[j] * (heat[sp] if p == 0 else 1)
    return mass, moments


def joint_flux(walls, face, falls, densities, values, heat, local_layers):
    """Use nested short intervals, dispatch long or exceptional rows separately."""
    end = min(face + local_layers, len(walls) - 1)
    furthest = walls[face] + max(falls)
    if furthest > walls[end]:
        kind = "top overshoot" if furthest > walls[-1] else "long interval"
        return (*reference_flux(walls, face, falls, densities, values, heat), kind)
    mass = [0] * len(falls)
    moments = [0] * len(values)
    for j in range(face, end):
        offset = walls[j] - walls[face]
        width = walls[j + 1] - walls[j]
        joint_mass = joint_heat = 0
        for sp, distance in enumerate(falls):
            length = max(0, min(width, distance - offset))
            carried = length * densities[sp][j]
            mass[sp] += carried
            joint_mass += carried
            joint_heat += carried * heat[sp]
        for p, field in enumerate(values):
            moments[p] += field[j] * (joint_heat if p == 0 else joint_mass)
    return mass, moments, "short interval"


def rational_example():
    walls = list(map(Fraction, [0, 2, 5, 9]))
    falls = [Fraction(1, 3), Fraction(2), Fraction(5), Fraction(7, 2)]
    density = [[Fraction(2 + sp + j, 7) for j in range(3)] for sp in range(4)]
    fields = [[Fraction(3 - j * p, 11) for j in range(3)] for p in range(5)]
    heat = list(map(Fraction, [1, 2, 3, 4]))
    direct = reference_flux(walls, 0, falls, density, fields, heat)
    result = joint_flux(walls, 0, falls, density, fields, heat, 2)
    assert result[:2] == direct and result[2] == "short interval"


def ordered_departures_example():
    """Check the mass update and its boundary budget on an ideal metric-free grid."""
    walls = list(map(Fraction, [0, 1, 3, 6]))
    departures = [Fraction(1, 2), Fraction(2), Fraction(4), Fraction(6)]
    falls = [d - z for d, z in zip(departures, walls)]
    # Unit mass in each possible donor in turn constructs the transfer matrix.
    for donor in range(3):
        density = [Fraction(j == donor, walls[j + 1] - walls[j]) for j in range(3)]
        flux = [
            joint_flux(walls, i, [falls[i]], [density], [], [1], 2)[0][0]
            for i in range(3)
        ] + [0]
        updated = [int(i == donor) + flux[i + 1] - flux[i] for i in range(3)]
        for i, mass in enumerate(updated):
            overlap = max(0, min(departures[i + 1], walls[donor + 1])
                          - max(departures[i], walls[donor]))
            assert mass == overlap * density[donor] and mass >= 0
        assert sum(updated) + flux[0] == 1
    # Crossing endpoints, uniform density: the first updated mass is negative.
    crossed = [Fraction(5, 2), Fraction(2), Fraction(4), Fraction(6)]
    flux = [
        joint_flux(walls, i, [crossed[i] - walls[i]], [[1, 1, 1]], [], [1], 2)[0][0]
        for i in range(3)
    ]
    assert (walls[1] - walls[0]) + flux[1] - flux[0] == Fraction(-1, 2)


def floating_point_examples():
    rng = Random(20261124)
    routes = Counter()
    max_scaled_error = max_budget_residual = 0.0
    for column in range(64):
        n = 78
        walls = [0.0, *accumulate(rng.uniform(0.25, 3.0) for _ in range(n))]
        densities = [
            [0.0 if rng.random() < 0.35 else rng.random() for _ in range(n)]
            for _ in range(4)
        ]
        fields = [[rng.uniform(-2, 2) for _ in range(n)] for _ in range(5)]
        heat = [1.0, 1.2, 1.2, 1.3]
        fluxes = []
        for face in range(n):
            width = walls[face + 1] - walls[face]
            typical = [rng.random() * width * 1.5 for _ in range(4)]
            mode = (column + face) % 8
            if mode == 0:
                falls = [0.0] * 4
            elif mode == 1:
                falls = [width, nextafter(width, 0), 0.0, nextafter(width, float("inf"))]
            elif mode == 2:
                falls = [width * (sp + 1) * 2 for sp in range(4)]
            else:
                falls = typical
            direct = reference_flux(walls, face, falls, densities, fields, heat)
            mass, moments, kind = joint_flux(walls, face, falls, densities, fields, heat, 2)
            routes[kind] += 1
            for actual, ref in zip(mass + moments, direct[0] + direct[1]):
                max_scaled_error = max(max_scaled_error, abs(actual - ref) / max(1, abs(ref)))
                assert isclose(actual, ref, rel_tol=1e-12, abs_tol=1e-12)
            fluxes.append(mass)
        fluxes.append([0.0] * 4)  # closed top, open lower boundary
        for sp in range(4):
            change = sum(fluxes[j + 1][sp] - fluxes[j][sp] for j in range(n))
            residual = change + fluxes[0][sp]
            max_budget_residual = max(max_budget_residual, abs(residual))
            assert abs(residual) < 1e-12
    print(f"64 synthetic columns x 78 faces; route counts: {dict(routes)}")
    print(f"max |difference| / max(1, |reference|) = {max_scaled_error:.3e}")
    print(f"max absolute column mass-budget residual = {max_budget_residual:.3e}")


if __name__ == "__main__":
    rational_example()
    print("Exact rational two-layer example passed (four masses and five moments).")
    ordered_departures_example()
    print("Exact ordered-departure mass budgets and crossing counterexample passed.")
    floating_point_examples()
