# coding: utf-8

import os
import ctypes as ct


DIRPATH_LIBRARIES = os.path.join(os.path.dirname(__file__), "lib")

LIBRARIES = {
    "cpu": {
        "filepath_library": os.path.join(DIRPATH_LIBRARIES, "space_time_energy_cpu.so"),
        "is_loaded": False,
        "lib": None,
        "space_time_energy_nn_argtypes": [
            ct.POINTER(ct.c_double),  # origin_time_sec 
            ct.POINTER(ct.c_float),  # longitude
            ct.POINTER(ct.c_float),  # latitude 
            ct.POINTER(ct.c_float),  # depth_km
            ct.POINTER(ct.c_float),  # magnitude
            ct.c_size_t,  # num_events
            ct.c_double, # t_max_sec
            ct.c_float, # b-value
            ct.c_float, # fractal dimension
            ct.c_double, # conversion_from_sec
            ct.c_float, # conversion_from_km
            ct.c_int, # num_threads
            ct.POINTER(ct.c_float),  # log_nn_r
            ct.POINTER(ct.c_float),  # log_nn_t
        ],
        }
    }

def load_library(device="cpu"):
    """Load library for device.

    This function loads the library only once, that is, if the library is
    successfully loaded a first time, it will be stored as a persistent
    variable in the LIBRARIES dictionary.

    Parameters
    ----------
    device: str, optional
        Device-compilated library, either "cpu" or "gpu".
    """
    # Get device name
    device_name = device.lower()

    # Check device name
    if device.lower() not in LIBRARIES:
        raise NameError(f"Device should be cpu or gpu, not {device.lower()}")

    # Check if shared object exsits
    library_info = LIBRARIES[device_name]
    filepath_library = library_info["filepath_library"]
    if os.path.exists(filepath_library):
        # If library was previously loaded, return it
        if library_info["is_loaded"] is True:
            return library_info["lib"]

        # Otherwise load it
        else:

            # Load
            lib = ct.cdll.LoadLibrary(filepath_library)

            # Declare types
            lib.space_time_energy_nn.argtypes = library_info[
                    "space_time_energy_nn_argtypes"
                    ]

            # Store pre-loaded library
            LIBRARIES[device_name]["is_loaded"] = True
            LIBRARIES[device_name]["lib"] = lib

            return lib

    else:
        raise OSError(
            f"Library for {device_name.upper()} was not compiled. The shared object {filepath_library} does not exist. \nYou may want to check the output of the command: \n\n \t python setup.py build_ext"
        )

