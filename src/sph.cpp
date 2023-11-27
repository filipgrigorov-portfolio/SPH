#include "sph.hpp"

#include <iostream>

namespace fluid {
/// @brief This class represents the "average" fluid attributes for each
/// particle based on its neighbours in a given radius R
SPH::SPH(double camWidth, double camHeight, uint32_t cubeSize,
         float particleDiameter, std::vector<Particle> &fluidParticles) {
  fluidParticles.reserve(cubeSize * cubeSize);
  const auto startRowPos = 0.2 * camHeight;
  const auto startColPos = 0.1 * camWidth;
  const float increment = particleDiameter;
  const auto endRowPos = startRowPos + (cubeSize * increment);
  const auto endColPos = startColPos + (cubeSize * increment);
  for (float row = startRowPos; row < endRowPos; row += increment) {
    for (float col = startColPos; col < endColPos; col += increment) {
      const float rowUpdated = row / startRowPos;
      const float colUpdated = col / startColPos;
      const uint32_t idx = uint32_t(rowUpdated * cubeSize + colUpdated);
      std::cout << col << "," << row << std::endl;
      fluidParticles.push_back({Eigen::Vector2d(col, row),
                                Eigen::Vector2d(0.0, 0.0),
                                Eigen::Vector2d(0.0, 0.0), 0.f, 0.f});
    }
  }
  std::cout << "Generated " << fluidParticles.size() << " particles"
            << std::endl;
}
} // namespace fluid
