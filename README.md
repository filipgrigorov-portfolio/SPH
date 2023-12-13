# SPH simple implementation in 2D

#### *Based on "Particle-Based Fluid Simulation for Interactive Applications" and "SPH Based Shallow Water Simulation"*
---

### (1) Setup

```
sudo apt update
sudo apt install libopengl-dev freeglut3-dev libeigen3-dev
```

### (2) Run

```
mkdir build && cd build
cmake ..
make
```

Run without experimental semi-circular boundary: `./sph`

Run with experimental semi-circular boundary: `./sph on`

Run with experimental semi-circular boundary (with penetration, namely, no restriction of fluid entering the BC, but still a cool visual): `./sph on with_penetration`

Note: **Requires C++14**

---

### Notes and experimental observations:

**-** Somewhat bigger time step results in instability of the simulation

**-** The rest density plays a huge rule in the stability of the simulation

**-** After having added surface tension, there is no such bounciness of the particles

**-** Setting the global fluid parameters is very heuristic for a proper visual

**-** The simulation could be re-implemented in CUDA for accelerating it, as it is kind of slow right now

**-** The simulation parameters could be learned by a neural network based on the initial conditions and boundary conditions to avoid the heurstical setup

**-** If boundaries are used, a good additional handling should be provided for the surface slip/friction
