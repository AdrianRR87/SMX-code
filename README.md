A repository that contains a library used for the calibration of the SMX v2.2 ASIC.
It contains a class with multiple functions, most of which are oriented towards calibrating the ADC and FAST discriminator, as well as checking the response of the Analog front-end circuits to injected pulses of variable amplitude.

Calibration procedure steps:
1. **def vrefpn_scan(...)** Scan of the VrefP, N, and Thr2_glb parameters in the selected range of amplitudes. The scan is made for only one channel. Since there is a fundamental difference between odd and even channels, a way to optimize this process can be by making the scan of an odd and an even channel, and averaging between the found settings. The selected scan range is defined by
     - amp_cal_min
     - amp_cal_max

2. The values of VrefP, N, and Thr2_glb found should be written in the corresponding ASIC registers
3. **def get_trim_adc_SA(...)** Calibration of the ADC in the desired range
4. **def get_trim_fast(...)** Calibration of the FAST discriminator at an specific threshold value
   
