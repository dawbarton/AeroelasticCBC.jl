Data collected by Irene Tartaruga on the aero-elastic flutter rig

pts = [1×20 struct]
        forcing_freq: [1×1 double]     # Frequency of oscillation - used to calculate phi_bis
      forcing_coeffs: [1×16 double]    # Additional forcing (Fourier coefficients) - typically zero
           rand_ampl: [1×1 double]     # Additional random forcing - typically zero
       x1_coeffs_ave: [1×16 double]    # Average estimated Fourier coefficients of x1 (over 5 cycles)
       x1_coeffs_var: [1×16 double]    # Variance of estimated Fourier coefficients of x1 (over 5 cycles)
    x1_target_coeffs: [1×16 double]    # Fourier coefficients of the control target 
      out_coeffs_ave: [1×16 double]    # Average estimated Fourier coefficients of out (over 5 cycles)
      out_coeffs_var: [1×16 double]    # Variance of estimated Fourier coefficients of out (over 5 cycles)
         aksim_angle: [1×50000 double] # Time series of pitch angle (radians)
             phi_bis: [1×50000 double] # Time series of estimated instantaneous phase
                  x1: [1×50000 double] # Time series of heave (multiply by 470 to obtain in millimetres)
           x1_target: [1×50000 double] # Time series of control target in heave
             Fshaker: [1×50000 double] # Time series of force at shaker interface
                 out: [1×50000 double] # Time series of control signal
              mean_h: [1×50000 double] # Time series of heave offset (actually a constant)
           timestamp: [1×6 double]     # Time stamp of data capture date

Quantities involving x1 (the heave) should be multiplied by 470 to obtain the
results in millimetres.

Time series are sampled at 5kHz.

Fourier coefficients are stored as [8×sine components, 8×cosine components],
starting with zero frequency and incrementing by one up to a (relative)
frequency of seven.

If array rather than structure
'x1' - heave: position measured by the laser
'x1_target' - heave target
'Fshaker' - Force at the transducer at the shaker
'aksim_angle' - pitch: angle measured by the encoder
'out' - control output defined in the beagle bone
'x1d' - derivative of the heave
'phi_bis' - argument of the rebuilt target
'mean_h' - mean value of the heave
