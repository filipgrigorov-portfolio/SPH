#pragma once

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
  SPH(double camWidth, double camHeight, uint32_t cubeSize,
      float particleSize, std::vector<Particle> &fluidParticles);
  ~SPH() = default;

  float getParticleSize() const { return mParticleSize; }

  void accumulateForces(std::vector<Particle> &fluidParticles, float viscosity);

  void accumulateParticlesDensity(std::vector<Particle> &fluidParticles, float restDensity, float gasConstant);

  void timeIntegrate(std::vector<Particle> &fluidParticles, float dt, float mass, float damping, float containerWidth, float containerHeight);

private:
  float mParticleSize = 12.f;
};
} // namespace fluid
