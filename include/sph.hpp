#pragma once

#include <unordered_map>
#include <vector>

#include <eigen3/Eigen/Dense>

namespace fluid {
struct Particle {
  // Note: (x, y)
  Eigen::Vector2f position;
  Eigen::Vector2f velocityField, forceField;

  float pressure = 0.f;
  float density = 0.f;

  float mass = 2.5f;
};

// Smoothing kernels
} // namespace fluid

namespace fluid {
class SPH {
public:
    SPH(float camWidth, float camHeight, uint32_t cubeSize, float particleSize, 
      std::unordered_map<std::string, float>& fluidProps, std::vector<Particle> &fluidParticles);
    ~SPH() = default;

    float getParticleSize() const { return mParticleSize; }

    void accumulateForces(std::vector<Particle> &fluidParticles);

    void accumulateParticlesDensity(std::vector<Particle> &fluidParticles);

    void timeIntegrate(std::vector<Particle> &fluidParticles, float dt, float dampingScale, float boundaryWidth, float boundaryHeight);

private:
    float mParticleSize = 12.f;

    float mRestDensity;
    float mGasConst;
    float mViscosity;
};
} // namespace fluid
