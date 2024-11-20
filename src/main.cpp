#include "geometrycentral/surface/manifold_surface_mesh.h"
#include "geometrycentral/surface/meshio.h"
#include "geometrycentral/surface/vertex_position_geometry.h"

#include "geometrycentral/surface/direction_fields.h"

#include "polyscope/polyscope.h"
#include "polyscope/surface_mesh.h"

#include "args/args.hxx"
#include "imgui.h"

using namespace geometrycentral;
using namespace geometrycentral::surface;

// == Geometry-central data
std::unique_ptr<ManifoldSurfaceMesh> mesh;
std::unique_ptr<VertexPositionGeometry> geometry;

// Polyscope visualization handle, to quickly add data to the surface
polyscope::SurfaceMesh *psMesh;

// Some algorithm parameters
float param1 = 42.0;

// Example computation function -- this one computes and registers a scalar
// quantity
void doWork() {
  polyscope::warning("Computing Gaussian curvature.\nalso, parameter value = " +
                     std::to_string(param1));

  geometry->requireVertexGaussianCurvatures();
  psMesh->addVertexScalarQuantity("curvature",
                                  geometry->vertexGaussianCurvatures,
                                  polyscope::DataType::SYMMETRIC);
}

// A user-defined callback, for creating control panels (etc)
// Use ImGUI commands to build whatever you want here, see
// https://github.com/ocornut/imgui/blob/master/imgui.h
//void myCallback() {
//
//  if (ImGui::Button("do work")) {
//    doWork();
//  }
//
//  ImGui::SliderFloat("param", &param1, 0., 100.);
//}

int main(int argc, char **argv) {

  // Configure the argument parser
  args::ArgumentParser parser("geometry-central & Polyscope example project");
  args::Positional<std::string> inputFilename(parser, "mesh", "A mesh file.");

  // Parse args
  try {
    parser.ParseCLI(argc, argv);
  } catch (args::Help &h) {
    std::cout << parser;
    return 0;
  } catch (args::ParseError &e) {
    std::cerr << e.what() << std::endl;
    std::cerr << parser;
    return 1;
  }

  // Make sure a mesh name was given
  if (!inputFilename) {
    std::cerr << "Please specify a mesh file as argument" << std::endl;
    return EXIT_FAILURE;
  }

  // Initialize polyscope
  polyscope::init();

  // Set the callback function
  //polyscope::state::userCallback = myCallback;

  // Load mesh
  std::tie(mesh, geometry) = readManifoldSurfaceMesh(args::get(inputFilename));

  // Register the mesh with polyscope
  psMesh = polyscope::registerSurfaceMesh(
      polyscope::guessNiceNameFromPath(args::get(inputFilename)),
      geometry->inputVertexPositions, mesh->getFaceVertexList(),
      polyscopePermutations(*mesh));

  // number of verticies, edges and faces
  int nVerts = mesh->nVertices();
  int nEdges = mesh->nEdges();
  int nFaces = mesh->nFaces();

  std::cout << nEdges << "\n";

  // Set vertex tangent spaces
  geometry->requireVertexTangentBasis();
  VertexData<Vector3> vBasisX(*mesh);
  VertexData<Vector3> vBasisY(*mesh);
  for (Vertex v : mesh->vertices()) {
    vBasisX[v] = geometry->vertexTangentBasis[v][0];
    vBasisY[v] = geometry->vertexTangentBasis[v][1];
  }

  auto vField =
      geometrycentral::surface::computeSmoothestVertexDirectionField(*geometry);

  geometry->requireVertexGaussianCurvatures();
  psMesh->addVertexScalarQuantity("Gaussian curvature",
                                  geometry->vertexGaussianCurvatures,
                                  polyscope::DataType::SYMMETRIC);

  // This step implements Discrete Exterior Calculus Operations
  geometry->requireDECOperators();
  VertexData<double> gaussianCurve(*mesh);

  // 1-form of velocity
  Eigen::VectorXd u(nEdges);
  std::vector<char>  orient(nEdges);

  for (int e=0; e<nEdges; e++) {
    u[e]=0.0;
    orient[e]=true;
  }

  u[0]=1.0;
  
  // Computes the linear operator discretizing the exterior derivative of a 0-form, 1-form and Hodge duals
  auto d0 = geometry->d0;
  auto d1 = geometry->d1;
  auto s0 = geometry->hodge0;
  auto s0m1 = geometry->hodge0Inverse;
  auto s1 = geometry->hodge1;
  auto s1m1 = geometry->hodge1Inverse;
  auto s2 = geometry->hodge2;
  auto s2m1 = geometry->hodge2Inverse;

  auto divergence = s0m1*(d0.transpose())*s1;

  Eigen::VectorXd divu(nEdges);
  divu=divergence*u;

  std::cout << divu << "\n";

  gaussianCurve = geometry->vertexGaussianCurvatures;




  psMesh->addVertexTangentVectorQuantity("VF", vField, vBasisX, vBasisY);
  psMesh->addOneFormTangentVectorQuantity("1-form", u, orient);

  // Give control to the polyscope gui
  polyscope::show();

  return EXIT_SUCCESS;
}
