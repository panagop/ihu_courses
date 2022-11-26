import sympy as sym
from sympy import Matrix, Number

def round_expr(expr, num_digits):
    """
    taken from:
    https://stackoverflow.com/questions/48491577/printing-the-output-rounded-to-3-decimals-in-sympy
    """
    return expr.xreplace({n : round(n, num_digits) for n in expr.atoms(Number)})


def K_with_axial(A: float, E: float, I: float, L: float) -> Matrix:
    """Return the stiffness matrix for an element with axial force.

    Args:
        A: Cross-sectional area of the element.
        E: Young's modulus of the element.
        I: Moment of inertia of the element.
        L: Length of the element.

    Returns:
        The stiffness matrix for the element.
    """
    k = Matrix([[A*E/L, 0, 0, -A*E/L, 0, 0],
                [0, 12*E*I/L**3, 6*E*I/L**2, 0, -12*E*I/L**3, 6*E*I/L**2],
                [0, 6*E*I/L**2, 4*E*I/L, 0, -6*E*I/L**2, 2*E*I/L],
                [-A*E/L, 0, 0, A*E/L, 0, 0],
                [0, -12*E*I/L**3, -6*E*I/L**2, 0, 12*E*I/L**3, -6*E*I/L**2],
                [0, 6*E*I/L**2, 2*E*I/L, 0, -6*E*I/L**2, 4*E*I/L]])
    return k


def K_without_axial(E: float, I: float, L: float) -> Matrix:
    """Return the stiffness matrix for an element with axial force.

    Args:
        E: Young's modulus of the element.
        I: Moment of inertia of the element.
        L: Length of the element.

    Returns:
        The stiffness matrix for the element.
    """
    k = Matrix([[12*E*I/L**3, 6*E*I/L**2, -12*E*I/L**3, 6*E*I/L**2],
        [6*E*I/L**2, 4*E*I/L, -6*E*I/L**2, 2*E*I/L],
        [-12*E*I/L**3, -6*E*I/L**2, 12*E*I/L**3, -6*E*I/L**2],
        [6*E*I/L**2, 2*E*I/L, -6*E*I/L**2, 4*E*I/L]])
    return k

def P_uniformly_distributed_load(q: float, L: float) -> Matrix:
    """Return the force vector for a uniformly distributed load.

    Args:
        q: Load per unit length.
        L: Length of the element.

    Returns:
        The force vector for the element.
    """
    p = Matrix([[-q*L/2],
                [-q*L**2/12],
                [-q*L/2],
                [q*L**2/12]])
    return p

def P_point_load_at_distance_a(P: float, L: float, a: float) -> Matrix:
    """Return the force vector for a point load at distance a.

    Args:
        P: Point load.
        L: Length of the element.
        a: Distance from the start of the element.

    Returns:
        The force vector for the element.
    """
    b = L - a
    p = Matrix([[-P*b/L**3*(L**2 - a**2 + a*b)],
                [-P*a*b**2/L**2],
                [-(P*a/L**3)*(L**2 - b**2 + a*b)],
                [P*a**2*b/L**2]])
    return p