import numpy as np
import sys
from scipy.signal import iirfilter, lfilter, sosfilt
f01_normalyzed = float(sys.argv[1])
f02_normalyzed = float(sys.argv[2])
display = sys.argv[-3].lower() == "display"
in_arr = sys.argv[-2][2:-1]
out_arr_rs = sys.argv[-1][2:-1]
in_arr = np.array([float(x) for x in in_arr.split(",")])
out_arr_rs = np.array([float(x) for x in out_arr_rs.split(",")])

sos = iirfilter(1, [f01_normalyzed, f02_normalyzed], btype="bandpass", output="sos")
filt_domain = sosfilt(sos, in_arr)
print(str([x for x in filt_domain])[1:-1])

# %%
if display:
    import matplotlib.pyplot as plt
    plt.subplot(311)
    plt.title("Time domain")
    plt.plot(out_arr_rs, label="rsfilt")
    plt.plot(filt_domain, label="pyfilt")
    plt.legend()
    plt.subplot(312)
    plt.title("Error domain")
    plt.plot(out_arr_rs-filt_domain, label="Error")
    plt.legend()
    plt.subplot(313)
    plt.title("freq domain")
    plt.semilogx(np.fft.rfftfreq(len(filt_domain))[2:], np.abs(np.fft.rfft(filt_domain)[2:]), label="fft pyfilt")
    plt.semilogx(np.fft.rfftfreq(len(out_arr_rs))[2:], np.abs(np.fft.rfft(out_arr_rs)[2:]), label="fft rs_filt")
    plt.legend()
    plt.show()
