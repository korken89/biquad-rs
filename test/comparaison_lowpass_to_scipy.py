import numpy as np
import sys
from scipy.signal import iirfilter, lfilter
freq_normalyzed = float(sys.argv[1])
display = sys.argv[-3].lower() == "display"
in_arr = sys.argv[-2][2:-1]
our_arr_rs = sys.argv[-1][2:-1]
in_arr = np.array([float(x) for x in in_arr.split(",")])
our_arr_rs = np.array([float(x) for x in our_arr_rs.split(",")])

b, a = iirfilter(2, freq_normalyzed, btype="lowpass")
filt_domain = lfilter(b, a, in_arr)
# print(in_arr, file=sys.stderr)
print(str([x for x in filt_domain])[1:-1])

# %%
if display:
    import matplotlib.pyplot as plt
    plt.subplot(211)
    plt.title("Time domain")
    plt.plot(our_arr_rs, label="rsfilt")
    plt.plot(filt_domain, label="pyfilt")
    plt.legend()
    plt.subplot(212)
    plt.title("Error")
    plt.plot(our_arr_rs-filt_domain, label="Error")
    plt.legend()
    plt.show()

# %%
f01, f02 = 0.2, 0.4
