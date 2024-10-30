#!/bin/bash

# Run ROOT script with different function arguments

#root -l -q -b 'pulser_calibration_v1.C("spectrum_gamma")'

#root -l -q -b 'pulser_calibration_v1.C("spectrum_proton")'

root -l -q -b 'pulser_calibration_v1.C("pulser_calibration")'
