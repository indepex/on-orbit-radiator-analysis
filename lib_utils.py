import numpy as np
import json5
import sys
from itertools import product
import warnings
from erfa import ErfaWarning
from astropy.time import Time


# -----------------------------------------------------------------------------------------
def jd2iso(jd):
    # Convert a Julian Date to ISO format while suppressing ERFA warnings.
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", ErfaWarning)
        jdiso = Time(jd, format='jd').iso
    return jdiso


# -----------------------------------------------------------------------------------------
def expand_dict(dic, expand_keys):
    # Expand selected dictionary keys into all combinations of their values.
    values = [np.atleast_1d(dic[k]) for k in expand_keys]

    result = []
    for combo in product(*values):
        d = dic.copy()
        for k, v in zip(expand_keys, combo):
            d[k] = v

        result.append(d)

    return result


# -----------------------------------------------------------------------------------------
def load_scalar(data, key):
    value = data[key]
    if hasattr(value, "shape") and value.shape == ():
        return value.item()
    return value


def parse_list(s, cast):
    if s is None:
        return None
    return [cast(x.strip()) for x in s.split(",") if x.strip()]


# -----------------------------------------------------------------------------------------
def load_json5(path, label):
    try:
        with open(path, 'r') as f:
            ret = json5.load(f)
    except FileNotFoundError as e:
        print(f"{label} file not found: {e}")
        sys.exit(1)
    except Exception as e:
        print(f"Unexpected error: {e}")
        sys.exit(1)
    return ret


# -----------------------------------------------------------------------------------------
def normalize(v):  # Normalize a vector or an array of vectors.
    v = np.asarray(v, dtype=float)
    if v.ndim == 1:
        n = np.linalg.norm(v)
    elif v.ndim == 2:
        n = np.linalg.norm(v, axis=1)[:, None]
    else:
        raise ValueError(f"normalize: unsupported shape {v.shape}")
    if np.any(n < np.finfo(float).eps):
        raise ValueError("normalize: vector length too small")

    return v / n


# -----------------------------------------------------------------------------------------
def list_rewrite(arg, item):
    # Return a list of values parsed from arg, or use item if arg is absent.
    item = np.atleast_1d(item)
    if arg:
        item = parse_list(arg, float)
    return item
