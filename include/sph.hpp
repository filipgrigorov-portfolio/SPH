#pragma once

#include <vector>

#include <eigen3/Eigen/Dense>

namespace fluid {
struct Particle {
  // Note: (x, y)
  Eigen::Vector2f position;
  Eigen::Vector2f velocityField, forceField;

  float pressure;
  float density;
};

// Smoothing kernels
} // namespace fluid

namespace fluid {
class SPH {
public:
  SPH(double camWidth, double camHeight, uint32_t cubeSize,
      float particleDiameter, std::vector<Particle> &fluidParticles);
  ~SPH() = default;

  void updateParticlesDensity(std::vector<Particle> &fluidParticles,
                              float density, float mass, float h,
                              float gasConstant, float radius);
  void timeIntegrate(std::vector<Particle> &fluidParticles, float dt,
                     float mass, float containerWidth, float containerHeight);

private:
};
} // namespace fluid
