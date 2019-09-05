# DecayGenerator
A small project to generate double beta decays (with and without neutrinos). Return energies and angles between electrons in the decay. 

## Installation 
The following projects are required for successfull compilation:
* GSL (GNU Scientific Library)
* Boost with components:
   * python3
   * system
   * numpy-py3
* PythonLibs 
* Python with component NumPy

In order to compile it do: 
~~~~
> cmake 
> make 
~~~~
This sohuld produce "lib" folder with libDecayGenerator.so

## Usage in python 

In order to use the package in Python just add the folder with `libDecayGenerator.so` to your sys path as : 
```python
import sys
sys.path.append("/path/to/the/decayer/lib/")

from libDecayGenerator import DecayGenerator

dec = DecayGenerator("MM")
dec.setSeed(5789) # Random number seed


max_ev_per_try = 50000 # Due to internal limitations not more than 1e5 events should be generated at once
nevents = 1000000
events = {}
for m in ["2vbb", "MM", "RHC"]:
    print("Model : ", m)
    events_left = nevents
    events[m] = np.zeros((0,3))
    gen.setModel(m)
    while events_left > 0:
        if events_left > max_ev_per_try: cur_nevents = max_ev_per_try
        else: cur_nevents = events_left
        events[m] = np.concatenate( (events[m],  
                                     gen.GenerateEvents(cur_nevents) )) 
        ### The generation function returns 3xN array with values of [T1, T2, costheta]
        
        events_left = events_left - cur_nevents
        
#### And making a plot 
bins = np.linspace(-1.0,1.0,101)
fig = plt.figure(figsize=(12,8), facecolor="w")

hist_MM = np.histogram(events['MM'][:,2], bins=bins)[0]
plt.step(bins, np.append(hist_MM, [0]), where = "post", color="r", label="Mass mixing")

hist_RHC = np.histogram(events['RHC'][:,2], bins=bins)[0]
plt.step(bins, np.append(hist_RHC, [0]), where = "post", color="b", label="RHC")

hist_2vbb = np.histogram(events['2vbb'][:,2], bins=bins)[0]
plt.step(bins, np.append(hist_2vbb, [0]), where = "post", color="g", label="2vbb")

plt.xlim(bins[0], bins[-1])
plt.ylim(0.0,plt.ylim()[1])
plt.legend(ncol=3, fontsize = 18, loc=9,framealpha=1)
plt.xlabel(r"$\cos(\theta_{12})$", fontsize=18)
plt.ylabel(r"Frequency [ a.u. ]", fontsize=18)
plt.show()
```

It should produce a figure similar to this one: 
![](https://raw.githubusercontent.com/terliuk/DecayGenerator/master/cos_theta12.png)
## Usage in C++
