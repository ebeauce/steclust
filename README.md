# steclust
C-based Python library for space-time-energy clustering and declustering analysis of seismic events

## Installation

The simplest way to use `steclust`:

```bash
pip install steclust
```

## For Developers

**Clone and install locally:**

```bash
git clone https://github.com/ebeauce/steclust.git
cd steclust
# option 1:
pip install .
# option 2 (check name of wheel in dist/):
python -m build
pip install dist/steclust-XXXX.whl
```

**If you modify C code:**

```bash
# option 1:
pip install . --no-cache-dir --force-reinstall
# option 2:
# Clean up previous Python and CMake build artifacts
rm -rf build/ dist/
# Re-run the build from scratch
python -m build
# Then reinstall to pick up changes:
pip install --force-reinstall dist/steclust-XXXX.whl
```


## Example


```python
log_rescaled_distance, log_rescaled_time = steclust.steclust.space_time_energy_nn(
    event_timings,
    event_longitudes,
    event_latitudes,
    event_depths,
    event_magnitudes,
    max_interevent_time_sec, # eg, 60 days 
    spatial_fraction_dimension, # eg, 2.6
    conversion_from_sec, # factor to convert inter-event times from seconds to whichever units you want to work with
    conversion_from_km, # factor to convert inter-event distances from kilometers to whichever units you want to work with
)

log_nearest_neighbor_distance = log_rescaled_distance + log_rescaled_time
```
