#include <stddef.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <float.h>
#include <omp.h>
#include <limits.h>
#include "space_time_energy_cpu.h"

#define RADIUS_EARTH_KM 6371.
#define PI 3.1415926535897932384626433

float _pairwise_hypo_distance(
        float lon1,
        float lat1,
        float dep_km1,
        float lon2,
        float lat2,
        float dep_km2
        ){
    float epi_distance = _pairwise_epi_distance(
            lon1, lat1, lon2, lat2
            );
    float hypo_distance = sqrt(pow(epi_distance, 2) + pow(dep_km2 - dep_km1, 2));

    return hypo_distance;
}

float _pairwise_epi_distance(
        float lon1,
        float lat1,
        float lon2,
        float lat2
        ){
    float deg2rad = PI / 180.;
    float av_lat = (lat1 + lat2) / 2. * deg2rad;
    float dx = (lon2 - lon1) * deg2rad * cos(av_lat);
    float dy = (lat2 - lat1) * deg2rad;
    float epi_distance = RADIUS_EARTH_KM * sqrt(dx*dx + dy*dy);

    return epi_distance;
}


void space_time_energy_nn(
        double *origin_time_sec,
        float *longitude,
        float *latitude,
        float *depth_km,
        float *magnitude,
        size_t num_events,
        double t_max_sec,
        float fractal_dimension,
        double conversion_from_sec,
        float conversion_from_km,
        int num_threads,
        float *log_nn_r,
        float *log_nn_t
        ){
    size_t i, j;
    size_t counter_every_n=1000;

    for (i=0; i<num_events-1; i++){
        double ie_time[num_events-i-1];
        float ie_distance[num_events-i-1];
        float rescaling_factor;
        double rescaling_factor_d;
        size_t num_valid_pairs=0;
        float log_nn_eta = FLT_MAX;
        float _log_eta;
        size_t nn_eta_index = 0;

        rescaling_factor = pow(10., -magnitude[i] / 2.);
        rescaling_factor_d = (double)rescaling_factor;

        if (i%counter_every_n == 0){
            printf("---- %zu / %zu ----\n", i, num_events);
        }

        for (j=i+1; j<num_events; j++){
            // loop over all events following event i
            // ---------- inter-event time
            ie_time[j-i-1] = (
                    origin_time_sec[j] - origin_time_sec[i]
                    );
            if (ie_time[j-i-1] < t_max_sec){
                num_valid_pairs += 1;
            }
            else{
                break;
            }
            ie_time[j-i-1] *= conversion_from_sec;
            ie_time[j-i-1] *= rescaling_factor_d;
            // ----------- inter-event distance
            ie_distance[j-i-1] = _pairwise_hypo_distance(
                    longitude[i],
                    latitude[i],
                    depth_km[i],
                    longitude[j],
                    latitude[j],
                    depth_km[j]
                    );
            //ie_distance[j-i-1] = _pairwise_epi_distance(
            //        longitude[i],
            //        latitude[i],
            //        longitude[j],
            //        latitude[j]
            //        );
            ie_distance[j-i-1] *= conversion_from_km;
            ie_distance[j-i-1] = pow(
                    ie_distance[j-i-1], fractal_dimension
                    ) * rescaling_factor;
            // ----------- combine both for space-time-energy distance
            _log_eta = (
                    (float)log(ie_time[j-i-1])
                    + logf(ie_distance[j-i-1])
                    );
            // ----------- update nearest neighbor if necessary
            if (_log_eta < log_nn_eta){
                log_nn_eta = _log_eta;
                nn_eta_index = j;
            }
        }
        log_nn_r[i] = logf(ie_distance[nn_eta_index-i-1]);
        log_nn_t[i] = (float)log(ie_time[nn_eta_index-i-1]);
    }
}

