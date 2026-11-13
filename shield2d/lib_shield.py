import numpy as np
import math
from scipy.optimize import root_scalar


# -------------------------------------------------------------------------------------------------
def bestrbt(T, h):
    """
    Find the two-stage shield geometry that minimizes the aperture-to-bottom
    radius ratio. Also return the radii sampled during the search.
    """

    # Scan candidate bottom radii to identify the range with valid solutions.
    rbtmin, rbtmax = 0, 0
    rbts = np.arange(20, 160.1, 1)
    for rbt in rbts:
        try:
            rap, r1, h1, Td = solve(T, h, rbt)
            if rbtmin == 0:
                rbtmin = rbt
            rbtmax = rbt
            # print(f"{rbt} {T} OK")
        except RuntimeError:
            pass
            # print(f"{rbt} {T} NG")
    print(f"rbt range {rbtmin} {rbtmax} for {T}")

    rbtbest, raprbt = 0, 1000
    rbts = np.arange(rbtmin, rbtmax, 1)
    raps = []
    for rbt in rbts:
        rap, r1, h1, Td = solve(T, h, rbt)
        print(f"{rap / rbt} for {T}")
        if rap / rbt < raprbt:
            raprbt = rap / rbt
            rbtbest = rbt
            r1best, h1best = r1, h1
            rapbest, rbtbest = rap, rbt
            Tdbest = Td
        raps.append(rap)
    Tdbest = np.degrees(Tdbest)

    return rapbest, rbtbest, r1best, h1best, Tdbest, raps, rbts


# -------------------------------------------------------------------------------------------------
def solve(T, h, rbt, e=0):
    """
    tan(T) = h/(r1+rbt)
    F = T/2,
    tan(F) = (r1-rbt)/h1
    tan(Td) = (h-h1)/(rap + r1)
    tan(G)=h1/(rbt+r1+e)
    G+P = Fd
    Fd+P = Td
    tan(P)=(rap-r1)/(h-h1)

    Solve the geometric constraints above for a two-stage shield.

    Parameters
    ----------
    T : float
        Reference elevation angle, in degrees.
    h : float
        Total shield height.
    rbt : float
        Bottom radius.
    e : float, optional
        Offset included in the geometric constraint for G.

    Returns
    -------
    rap : float
        Aperture radius.
    r1 : float
        Intermediate radius.
    h1 : float
        Intermediate height.
    Td : float
        Derived angle, in radians.
    """
    Tr = math.radians(T)
    # Compute the intermediate radius and height directly from the constraints.
    r1 = h / math.tan(Tr) - rbt
    F = Tr / 2
    h1 = (r1 - rbt) / math.tan(F)

    if h1 <= 0:
        raise RuntimeError("Computed intermediate height h1 is nonpositive")
    if h1 >= h:
        raise RuntimeError("Computed intermediate height h1 is not below the total height h")

    G = math.atan2(h1, rbt + r1 + e)
    dh = h - h1

    def values_from_rap(rap):
        Td = math.atan2(dh, rap + r1)
        P = math.atan2(rap - r1, dh)
        Fd = G + P
        return Td, P, Fd

    def residual(rap):
        Td, P, Fd = values_from_rap(rap)
        # Enforce the angular constraint Fd + P = Td.
        return Fd + P - Td
        # Equivalently: return G + 2 * P - Td.

    # Bracket the solution on the side where rap > r1.
    lo = r1 + 1e-12
    hi = max(2 * r1, r1 + h)

    f_lo = residual(lo)
    f_hi = residual(hi)

    for _ in range(100):
        if f_lo * f_hi < 0:
            break
        hi *= 2
        f_hi = residual(hi)
    else:
        raise RuntimeError("Could not bracket a solution for rap")

    sol = root_scalar(residual, bracket=[lo, hi], method="brentq")
    rap = sol.root

    Td, P, Fd = values_from_rap(rap)

    return rap, r1, h1, Td


# -------------------------------------------------------------------------------------------------
def trace_ray(model, P, U):
    """
    Trace a ray from P along U through the shield until it terminates.

    Returns the ray path as an array of points and the ID of the final
    surface. The aperture and bottom are terminating surfaces.
    """
    max_points = 10        # Maximum number of tracing steps.
    idx = -1
    Ps = [P]
    for i in range(max_points):
        cont, Q, V, idy = trace_step(model, P, U, idx)
        Ps.append(Q)
        if not cont:
            break
        P, U, idx = Q, V, idy

    return np.array(Ps), idy


# -------------------------------------------------------------------------------------------------
def trace_step(model, P, U, idx):
    """
    Find the next surface intersected by the ray from P along U.

    Skip the current surface (idx), then select the nearest intersection.
    Return (cont, Q, V, idy), where Q is the intersection point, V is
    the specularly reflected direction, and idy is the surface ID.
    cont is False at the aperture or bottom, and True otherwise.
    The closed model is expected to yield an intersection.
    """
    # print(f"Trace from {idx} P {P} U {U}")
    # IDs of the starting surface, aperture, and bottom.
    st_id, ap_id, bt_id = model["start_id"], model["aperture_id"], model["bottom_id"]
    # Find the nearest intersection among eligible segments.
    dmin = np.finfo(np.float64).max
    for i in range(model["nlines"]):
        if idx == -1 and (i == st_id or i == ap_id):
            continue
        elif idx == i:
            continue
        else:
            # print(f"from {idx}--{i}")
            hit, d, Qt, Vt = line_cross(model[i]["A"], model[i]["B"], P, U)
            # print(f"hit {hit} d {d}")
            if hit and d < dmin:
                dmin, Q, V, idy = d, Qt, Vt, i
    # print(f"out dmin {dmin}")
    if dmin == np.finfo(np.float64).max:
        raise RuntimeError(f"No surface intersection from {idx} P {P} U {U}")
    if idy == ap_id or idy == bt_id:
        cont = False
    else:
        cont = True
    return cont, Q, V, idy


# -------------------------------------------------------------------------------------------------
def line_cross(A, B, P, U):
    """
    Intersect the ray from P along U with the line segment from A to B.

    Returns (hit, d, Q, V): whether the intersection lies on the
    segment, the ray parameter d, the intersection point Q, and the
    specularly reflected direction V. A parallel ray returns None
    for d, Q, and V.
    """
    # print(f"line cross {A} {B} {P} {U}")
    M = (B - A) / np.linalg.norm(B - A)   # Unit tangent from A to B; assumes a nonzero segment length.
    N = np.array([M[1], -M[0]])         # Unit normal to the segment.
    # Q = P + d U; (Q - A) dot N = 0, so d (U dot N) = (A - P) dot N.
    un = np.dot(U, N)
    wn = np.dot(A - P, N)
    if abs(un) < 1e-8:
        hit = False
        Q, V, d = None, None, None
    else:
        d = wn / un
        Q = P + d * U
        t = np.dot(Q - A, M) / np.linalg.norm(B - A)   # Fractional position along the segment.
        V = U - 2 * np.dot(U, N) * N
        if 0.0 <= t and t <= 1.0:
            hit = True
        else:
            hit = False
        # print(f"t {t} {hit} Q {Q} A {A} B {B}")
    return hit, d, Q, V
