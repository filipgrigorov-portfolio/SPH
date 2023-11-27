#pragma once

#include <vector>

#include <eigen3/Eigen/Dense>

namespace fluid {
struct Particle {
  // Note: (x, y, z)
  Eigen::Vector2d position;
  Eigen::Vector2d velocityField, forceField;

  float pressure;
  float density;
};
} // namespace fluid

namespace fluid {
class SPH {
public:
  SPH(double camWidth, double camHeight, uint32_t cubeSize,
      float particleDiameter, std::vector<Particle> &fluidParticles);
  ~SPH() = default;

private:
};
} // namespace fluid
