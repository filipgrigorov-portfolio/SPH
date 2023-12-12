#include <iostream>

#include "sph.hpp"

#include <GL/glut.h>

namespace {
static constexpr int WINDOW_WIDTH = 800;
static constexpr int WINDOW_HEIGHT = 600;
static constexpr double FRUSTRUM_WIDTH = 1.2 * WINDOW_WIDTH;
static constexpr double FRUSTRUM_HEIGHT = 1.2 * WINDOW_HEIGHT;

static constexpr float RED = 0.05490196078f;
static constexpr float GREEN = 0.5294117647f;
static constexpr float BLUE = 1.0f;
static constexpr float ALPHA = 1;

static constexpr float REST_DENSITY = 200; // rest density for water (affects pressure forces)
static constexpr float VISCOSITY = 200.f;
static constexpr float TIME_STEP = 0.0005f; // Note: If the time step is bigger, then simulation is not stable
static constexpr float GAS_CONST = 2000;
static constexpr float DAMPING = 0.6f;

static constexpr uint32_t CUBE_SIZE = 20;
static std::vector<fluid::Particle> fluidParticles;
static fluid::SPH sph(FRUSTRUM_HEIGHT, FRUSTRUM_WIDTH, CUBE_SIZE, 20.f, fluidParticles);
} // namespace

void initOpenGL() {
  // Note: select background color
  glClearColor(0.9f, 0.9f, 0.9f, 1);

  // Note: Remove aliasing
  glEnable(GL_POINT_SMOOTH);

  glPointSize(sph.getParticleSize() / 2.f);

  // Note: initialize viewing values
  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
  // Note: multiply the current matrix with an orthographic matrix
  glOrtho(0.0, FRUSTRUM_WIDTH, 0.0, FRUSTRUM_HEIGHT, 0.0, 1.0);
}

// Note: Incovenient but it uses global static variables
void renderFunction() {
  // Clear all pixels
  glClear(GL_COLOR_BUFFER_BIT);

  // Note: Color of the fluid particles, R,G,B,Alpha
  glColor4f(RED, GREEN, BLUE, ALPHA);

  glBegin(GL_POINTS);
  {
    // Note: Go through the particles
    std::cout << "Rendering " << fluidParticles.size() << " particles"
              << std::endl;
    for (const auto &fluidParticle : fluidParticles) {
      // std::cout << fluidParticle.position(0) << "," << fluidParticle.position(1) << std::endl;
      glVertex2d(fluidParticle.position(0), fluidParticle.position(1));
    }
  }
  glEnd();

  glutSwapBuffers();
}

void updateStep() {
  // Update density
  sph.accumulateParticlesDensity(fluidParticles, REST_DENSITY, GAS_CONST);
  // Update forces
  sph.accumulateForces(fluidParticles, VISCOSITY);
  // Time integration
  sph.timeIntegrate(fluidParticles, TIME_STEP, REST_DENSITY, DAMPING, FRUSTRUM_WIDTH, FRUSTRUM_HEIGHT);

  glutPostRedisplay();
}

int main(int argc, char **argv) {
  glutInitWindowSize(WINDOW_WIDTH, WINDOW_HEIGHT);
  glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB);
  glutInit(&argc, argv);
  glutCreateWindow("SPH 3D");

  // Note: Intialize the openGL parameters
  initOpenGL();

  // Note: Display and refresh of window
  glutDisplayFunc(renderFunction);
  glutIdleFunc(updateStep);
  // glutKeyboardFunc(keyboardEvent); // TODO: move the tank around to slosh the
  // water

  glutMainLoop();

  return EXIT_SUCCESS;
}
