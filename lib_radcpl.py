import ctypes
import numpy as np
from pathlib import Path

WORDLENGTH = 64


# C-compatible structure definitions---------------------------------------------------------------------------------
# Three-dimensional vector; matches VEC in lib_vec.h
class VEC(ctypes.Structure):
    _fields_ = [("x", ctypes.c_double), ("y", ctypes.c_double), ("z", ctypes.c_double)]


# Ray path; matches RAY in lib_radcouple.h
class RAY(ctypes.Structure):
    _fields_ = [
        ("n", ctypes.c_int32),              # Number of stored ray path points
        ("f", ctypes.c_int32),              # Ray termination status (HIT_*)
        ("points", ctypes.POINTER(VEC)),    # Ray path points
        ("id", ctypes.POINTER(ctypes.c_int32)),  # Surface ID for each path point
        ("final_dir", ctypes.c_int32),           # Index of the final outgoing direction
    ]


# Surface geometry and reflection properties; matches SURFACE in lib_radcouple.h
class SURFACE(ctypes.Structure):
    _fields_ = [("name", ctypes.c_char * WORDLENGTH), ("type", ctypes.c_int),
                ("reflect", ctypes.c_int), ("r", ctypes.c_double), ("h", ctypes.c_double),
                ("p0", VEC), ("p1", VEC), ("p2", VEC), ("p3", VEC), ("p4", VEC), ]


# Data conversion helpers ----------------------------------------------------------------------------------
def vec3(v):
    v = np.asarray(v, dtype=float)
    return VEC(v[0], v[1], v[2])


def make_surface(d):
    s = SURFACE()

    s.name = d.get("name", "").encode("utf-8")[:WORDLENGTH - 1]

    # s.type = d.get("type", "").encode("utf-8")[:WORDLENGTH - 1]
    s.type = d["type"]
    s.reflect = d.get("reflect", 0)

    s.r = float(d.get("radius", 0.0))
    s.h = float(d.get("height", 0.0))

    s.p0 = vec3(d.get("p0", [0, 0, 0]))
    s.p1 = vec3(d.get("p1", [0, 0, 0]))
    s.p2 = vec3(d.get("p2", [0, 0, 0]))
    s.p3 = vec3(d.get("p3", [0, 0, 0]))
    s.p4 = vec3(d.get("p4", [0, 0, 0]))

    return s


def vec_array(values):
    """
    Convert an (N, 3) NumPy array or sequence into a C VEC array.
    """
    values = np.asarray(values, dtype=float)

    if values.ndim == 1:
        values = values.reshape(1, 3)

    if values.ndim != 2 or values.shape[1] != 3:
        raise ValueError(f"VEC array must have shape (N,3), got {values.shape}")

    n = values.shape[0]
    ArrayType = VEC * n
    arr = ArrayType()

    for i in range(n):
        arr[i] = VEC(values[i, 0], values[i, 1], values[i, 2])

    return arr, n


def _copy_rays(raw_ptr, n_rays):
    """
    Copy a C-allocated RAY array into independent Python/NumPy arrays.
    """
    rays_ptr = ctypes.cast(raw_ptr, ctypes.POINTER(RAY))

    points_list = []                                  # Ray path points for each ray
    ids_list = []                                     # Surface IDs for each ray
    flags = np.empty(n_rays, dtype=np.int32)          # Ray termination status (HIT_NONE or HIT_ABSORB)
    final_dirs = np.empty(n_rays, dtype=np.int32)     # Final outgoing direction ID

    for k in range(n_rays):
        ray = rays_ptr[k]
        n = ray.n

        flags[k] = ray.f
        final_dirs[k] = ray.final_dir

        if n > 0 and bool(ray.points):
            p_arr = ctypes.cast(
                ray.points,
                ctypes.POINTER(VEC * n)
            ).contents

            pts = np.array(
                [[p_arr[i].x, p_arr[i].y, p_arr[i].z] for i in range(n)],
                dtype=float,
            )
        else:
            pts = np.empty((0, 3), dtype=float)

        if n > 0 and bool(ray.id):
            id_arr = ctypes.cast(
                ray.id,
                ctypes.POINTER(ctypes.c_int32 * n)
            ).contents

            ids = np.array(
                [id_arr[i] for i in range(n)],
                dtype=np.int32,
            )
        else:
            ids = np.empty((0,), dtype=np.int32)

        points_list.append(pts)
        ids_list.append(ids)

    return points_list, ids_list, flags, final_dirs


# Load the C shared library and declare function signatures -------------------------------------------------------------------------
_lib_path = Path(__file__).with_name("lib_radcouple.so")
lib = ctypes.CDLL(str(_lib_path))

lib.count_out_dirs.argtypes = [
    ctypes.POINTER(VEC),              # Array of ray origins
    ctypes.c_int32,                   # Number of ray origins
    ctypes.c_int32,                   # Surface ID of ray origins
    ctypes.POINTER(VEC),              # Array of initial ray directions
    ctypes.c_int32,                   # Number of initial directions
    ctypes.POINTER(ctypes.c_double),  # Escaping emission weight by outgoing direction
    ctypes.POINTER(ctypes.c_double),  # Total emitted weight, including absorbed rays
    ctypes.POINTER(ctypes.c_double),  # Weight absorbed by surfaces
]
lib.count_out_dirs.restype = None

lib.gen_rays.argtypes = [
    ctypes.POINTER(VEC),   # StartPoints
    ctypes.c_int32,        # num_start
    ctypes.c_int32,        # start_surf_id
    ctypes.POINTER(VEC),   # Directions
    ctypes.c_int32,        # num_dir
    ctypes.c_int32,        # out_id_flag
    ctypes.c_int32,        # add_shell
    ctypes.c_double,       # shell_radius
]
lib.gen_rays.restype = ctypes.c_void_p

lib.free_rays.argtypes = [ctypes.c_void_p, ctypes.c_int32]
lib.free_rays.restype = None

lib.set_surfaces.argtypes = [ctypes.POINTER(SURFACE), ctypes.c_int, ctypes.c_int,]
lib.set_surfaces.restype = None

lib.fibonacci_sphere.argtypes = [ctypes.c_int]
lib.fibonacci_sphere.restype = ctypes.c_void_p

lib.fibonacci_sphere_half.argtypes = [ctypes.c_int, ctypes.c_int]
lib.fibonacci_sphere_half.restype = ctypes.c_void_p

lib.free_vec_array.argtypes = [ctypes.c_void_p]
lib.free_vec_array.restype = None


# Python interface to the C library ------------------------------------------------------------------
# Accumulate cosine-weighted emission by outgoing direction
def count_out_dirs(Start, start_surface, Direction):
    start_arr, num_start = vec_array(Start)   # Ray origins
    dir_arr, num_dir = vec_array(Direction)   # Initial directions
    dir_emission = np.zeros(num_dir, dtype=float)  # Emission weight by outgoing direction
    total_emission = ctypes.c_double()      # Total emitted weight
    absorbed_emission = ctypes.c_double()   # Absorbed emission weight

    lib.count_out_dirs(
        start_arr,
        num_start,
        int(start_surface),
        dir_arr,
        num_dir,
        dir_emission.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
        ctypes.byref(total_emission),
        ctypes.byref(absorbed_emission),
    )

    return dir_emission, total_emission.value, absorbed_emission.value


# Generate and release ray paths
def free_rays(raw_ptr, n_rays):
    if raw_ptr:
        lib.free_rays(raw_ptr, n_rays)


def gen_rays(Start, start_surface, Direction, id_out=False, shell=True, shell_rad=3000):
    """
    Generate ray paths using the C lib_radcouple.gen_rays() function.

    Parameters
    ----------
    Start : array-like, shape (num_start, 3)
        Ray origin coordinates.
    start_surface : int
        Surface ID associated with the ray origins.
    Direction : array-like, shape (num_dir, 3)
        Initial ray directions.
    id_out : bool, optional
        Assign a Fibonacci direction ID to rays escaping the model.
    shell : bool, optional
        Append intersections with the outer spherical shell.
    shell_rad : float, optional
        Radius of the outer spherical shell.

    Returns
    -------
    points_list, ids_list, flags, final_dirs
        Per-ray path coordinates, surface IDs, termination flags,
        and final outgoing direction IDs. The arrays are copied from
        C-owned memory before that memory is released.
    """
    start_arr, num_start = vec_array(Start)
    dir_arr, num_dir = vec_array(Direction)

    raw_ptr = lib.gen_rays(
        start_arr,
        ctypes.c_int32(num_start),
        ctypes.c_int32(start_surface),
        dir_arr,
        ctypes.c_int32(num_dir),
        ctypes.c_int32(1 if id_out else 0),
        ctypes.c_int32(1 if shell else 0),
        ctypes.c_double(shell_rad),
    )

    n_rays = num_start * num_dir

    try:
        points_list, ids_list, flags, final_dirs = _copy_rays(raw_ptr, n_rays)
    finally:
        free_rays(raw_ptr, n_rays)

    return points_list, ids_list, flags, final_dirs


# Transfer surface geometry to the C library
def set_surface(data, debug=False):
    n = len(data)

    SurfaceArray = SURFACE * n
    arr = SurfaceArray()
    dflag = 0
    if debug:
        dflag = 1

    for i, d in enumerate(data):
        arr[i] = make_surface(d)

    lib.set_surfaces(arr, n, dflag)


# Generate approximately uniform directions on a sphere
def fibonacci_sphere(N):
    raw_ptr = lib.fibonacci_sphere(N)
    try:
        return _vec_ptr_to_array(raw_ptr, N)
    finally:
        lib.free_vec_array(raw_ptr)


def fibonacci_sphere_half(N, keep=False):
    raw_ptr = lib.fibonacci_sphere_half(N, 1 if keep else 0)

    if keep:
        return _vec_ptr_to_array(raw_ptr, N)

    try:
        return _vec_ptr_to_array(raw_ptr, N)
    finally:
        lib.free_vec_array(raw_ptr)


def _vec_ptr_to_array(raw_ptr, n):
    ptr = ctypes.cast(raw_ptr, ctypes.POINTER(VEC))
    return np.array(
        [[ptr[i].x, ptr[i].y, ptr[i].z] for i in range(n)],
        dtype=float,
    )
