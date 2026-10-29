import ctypes as ct
import numpy as np

from .core import load_library
from os import cpu_count

def space_time_energy_nn(
        origin_time_sec,
        longitude,
        latitude,
        depth_km,
        magnitude,
        t_max_sec,
        fractal_dimension,
        conversion_from_sec,
        conversion_from_km,
        num_threads=1,
        device="cpu"
        ):
    """
    """
    lib = load_library(device)

    # cast to expected data type
    num_events = np.int32(len(longitude))
    t_max_sec= np.float64(t_max_sec)
    fractal_dimension = np.float32(fractal_dimension)
    conversion_from_sec = np.float64(conversion_from_sec)
    conversion_from_km = np.float32(conversion_from_km)
    num_threads = np.int32(num_threads)

    origin_time_sec = origin_time_sec.astype(np.float64)
    longitude = longitude.astype(np.float32)
    latitude = latitude.astype(np.float32)
    depth_km = depth_km.astype(np.float32)
    magnitude = magnitude.astype(np.float32)

    log_nn_r = np.zeros(num_events, dtype=np.float32)
    log_nn_t = np.zeros(num_events, dtype=np.float32)

    lib.space_time_energy_nn(
            origin_time_sec.ctypes.data_as(ct.POINTER(ct.c_double)),
            longitude.ctypes.data_as(ct.POINTER(ct.c_float)),
            latitude.ctypes.data_as(ct.POINTER(ct.c_float)),
            depth_km.ctypes.data_as(ct.POINTER(ct.c_float)),
            magnitude.ctypes.data_as(ct.POINTER(ct.c_float)),
            num_events,
            t_max_sec,
            fractal_dimension,
            conversion_from_sec,
            conversion_from_km,
            num_threads,
            log_nn_r.ctypes.data_as(ct.POINTER(ct.c_float)),
            log_nn_t.ctypes.data_as(ct.POINTER(ct.c_float))
            )

    return log_nn_r, log_nn_t
