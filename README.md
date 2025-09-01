A repository that contains a library used for the calibration of the SMX v2.2 ASIC.
It contains a class with multiple functions, most of which are oriented towards calibrating the ADC and FAST discriminator, as well as checking the response of the Analog front-end circuits to injected pulses of variable amplitude.

Calibration procedure steps:
1. Scan of the VrefP, N, and Thr2_glb parameters in the selected range of amplitudes defined by
     - amp_cal_min
     - amp_cal_max
  The scan is made for only one channel. Since there is a fundamental difference between odd and even channels, a way to optimize this process can be by making the scan of an odd and an even channel, and averaging between the found settings.
2. The values of VrefP, N, and Thr2_glb found should be written in the corresponding ASIC registers
3. Calibration of the ADC in the desired range
   
