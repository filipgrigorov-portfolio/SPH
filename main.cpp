#include <iostream>
#include <unordered_map>

#include "sph.hpp"

#include <GL/glut.h>

namespace {
/// @brief Colors of fluid particles
static constexpr float RED = 0.05490196078f;
static constexpr float GREEN = 0.5294117647f;
static constexpr float BLUE = 1.0f;
static constexpr float ALPHA = 1;

/// @brief Properties to setting up openGL primitives
static std::unordered_map<std::string, float> openGLProps = {
    {"window_width", 400},
    {"window_height", 300},
    {"frustrum_width", 800},
    {"frustrum_height", 600}
};

/// @brief Global fluid properties
static std::unordered_map<std::string, float> fluidProps = {
    {"rest_density", 400},
    {"viscosity", 200},
    {"gas_const", 2000}
};

/// @brief Time integration properties (during simulation)
// Note: If the time step is bigger, then simulation is not stable.
static constexpr float TIME_STEP = 0.0006f;
static constexpr float DAMPING_SCALE = 0.2f;

static constexpr uint32_t CUBE_SIZE = 20;

// Note: The particles have to be defined global due to OpenGL
static std::vector<fluid::Particle> fluidParticles;
static fluid::SPH sph(
  openGLProps["frustrum_height"], openGLProps["frustrum_width"], 
  CUBE_SIZE, 18.f, 
  fluidProps,
  fluidParticles);
} // namespace

void initOpenGL() {
  // Note: select background color
  glClearColor(0.f, 0.f, 0.f, 1);

  // Note: Remove aliasing
  glEnable(GL_POINT_SMOOTH);

  glPointSize(sph.getParticleSize() / 2.f);

  // Note: initialize viewing values
  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
  // Note: multiply the current matrix with an orthographic matrix
  glOrtho(0.0, openGLProps["frustrum_width"], 0.0, openGLProps["frustrum_height"], 0.0, 1.0);
}

// Note: Incovenient but it uses global static variables
void renderFunction() {
  // Clear all pixels
  glClear(GL_COLOR_BUFFER_BIT);

  // Note: Color of the fluid particles, (R, G, B, Alpha)
  glColor4f(RED, GREEN, BLUE, ALPHA);

  glBegin(GL_POINTS);
  {
    // Note: Go through the particles
    std::cout << "Rendering " << fluidParticles.size() << " particles" << std::endl;
    for (const auto &fluidParticle : fluidParticles) {
      glVertex2d(fluidParticle.position(0), fluidParticle.position(1));
    }
  }
  glEnd();

  glutSwapBuffers();
}

void updateStep() {
  // Update density
  sph.accumulateParticlesDensity(fluidParticles);
  // Update forces
  sph.accumulateForces(fluidParticles);
  // Time integration
  sph.timeIntegrate(fluidParticles, TIME_STEP, DAMPING_SCALE, openGLProps["frustrum_width"], openGLProps["frustrum_height"]);

  glutPostRedisplay();
}

int main(int argc, char **argv) {
  glutInitWindowSize(openGLProps["window_width"], openGLProps["window_height"]);
  glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB);
  glutInit(&argc, argv);
  glutCreateWindow("SPH 3D");

  // Note: Intialize the openGL parameters
  initOpenGL();

  // Note: Display and refresh of window
  glutDisplayFunc(renderFunction);
  glutIdleFunc(updateStep);

  glutMainLoop();

  return EXIT_SUCCESS;
}
