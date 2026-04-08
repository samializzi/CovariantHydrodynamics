#include "geometrycentral/surface/manifold_surface_mesh.h"
#include "geometrycentral/surface/meshio.h"
#include "geometrycentral/surface/vertex_position_geometry.h"


#include "geometrycentral/surface/direction_fields.h"


#include "polyscope/polyscope.h"
#include "polyscope/surface_mesh.h"


#include "args/args.hxx"
#include "imgui.h"


#include <fstream>
#include <iostream>


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


  auto vKg =geometry->vertexGaussianCurvatures;


  // This step implements Discrete Exterior Calculus Operations
  geometry->requireDECOperators();
  VertexData<double> gaussianCurve(*mesh);


  // 1-form of velocity
  //Eigen::VectorXd u(nEdges);
  //std::vector<char>  orient(nEdges);


  //for (int e=0; e<nEdges; e++) {
  //  u[e]=0.0;
  //  orient[e]=true;
  //}


  //u[0]=1.0;
  
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
  //divu=divergence*u;


  //std::cout << divu << "\n";


  // Finds the array of Gaussian curvatures of each vertex 
  gaussianCurve = geometry->vertexGaussianCurvatures;
  
  geometry->requireVertexDualAreas();
  geometry->requireVertexIndices();
  auto area = geometry->vertexDualAreas;
  
  for (Vertex v: mesh->vertices()) {
    int i = geometry->vertexIndices[v];
    gaussianCurve[i] = gaussianCurve[i]/area[i];  
  }


  // Computes gaussian curvature for each edge in the given mesh
  Eigen::VectorXd Kedge(nEdges);
  geometry->requireEdgeIndices();
  geometry->requireVertexIndices();
  geometry->requireHalfedgeIndices();


  for (Edge e: mesh->edges()) {
    int i = geometry->edgeIndices[e];
    Vertex vertex_tail = e.halfedge().tailVertex();
    Vertex vertex_head = e.halfedge().tipVertex();
    int p = geometry->vertexIndices[vertex_tail];
    int q = geometry->vertexIndices[vertex_head];
    Kedge[i] = (0.5)*(gaussianCurve[p] + gaussianCurve[q]);
  }


  std::cout << "Number of Verts: " << nVerts << "\n";
  std::cout << "Number of Edges: " << nEdges << "\n";
  std::cout << "Number of Faces: " << nFaces << "\n";
  std::cout << d0.cols() << "\n";
  std::cout << d0.rows() << "\n";


  // Finding the maximum and minimum radius of the mesh
  geometry->requireVertexIndices();
  geometry->requireVertexPositions();


  double max_radius = 0;
  double min_radius = 1000; // !!!
  int index_max_radius = 0;
  int index_min_radius = 0;


  for (Vertex v: mesh->vertices()) {
    int i = geometry->vertexIndices[v];
    Vector3 vertex_coord = geometry->vertexPositions[v];
    double size = norm(vertex_coord);


    if (size > max_radius) {
      max_radius = size;
      index_max_radius = i;
    }


    if (size < min_radius) {
      min_radius = size;
      index_min_radius = i;
    }
  }


  std::cout << "\n" << "Maximum radius at index " << index_max_radius << " with magnitude " << max_radius << std::endl;
  std::cout << "Minimum radius at index " << index_min_radius << " with magnitude " << min_radius << "\n" << std::endl;


  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // Gaussian curvature L2 norm error (for SPHERE ONLY)
  double average_sphere_radius = (0.5)*(max_radius + min_radius);
  double trueKg = 1.0/(average_sphere_radius*average_sphere_radius);


  Eigen::VectorXd Kvertex_error(nVerts);


  // finds average Kg vertex percentage error (for SPHERE ONLY)
  double Kvertex_error_avg = 0.0;
  for (Vertex v: mesh->vertices()) {
    int i = geometry->vertexIndices[v];
    double calc = (gaussianCurve[i] - trueKg)/trueKg;
    Kvertex_error[i] = sqrt(calc*calc);
    Kvertex_error_avg +=  Kvertex_error[i];
  }


  Kvertex_error_avg = Kvertex_error_avg/nVerts;


  // finds the maximum Kg vertex percentage error (for SPHERE ONLY)
  double Kvertex_error_max = 0;
  for (int i = 0; i < nVerts; i++) {
    if (Kvertex_error[i] > Kvertex_error_max) {
      Kvertex_error_max = Kvertex_error[i];
    }
  }


  std::cout << "L2 Norm of Vertex_Gaussian_curvature is:\n" << Kvertex_error_avg << std::endl;
  std::cout << "Max L2 Error on Gaussian curvature:\n" << Kvertex_error_max << std::endl;


  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // Here we declare our Stokes Equations Linear Operator, A
  printf("\nType 1 to perform matrix calculation and solve for vector field u, 0 otherwise: ");
  int matrix_calc_permission;
  
  Eigen::MatrixXd A = Eigen::MatrixXd::Zero(nEdges+nVerts,nEdges+nVerts);
  Eigen::MatrixXd div1 = Eigen::MatrixXd::Zero(nVerts, nEdges);
  Eigen::MatrixXd l1 = Eigen::MatrixXd::Zero(nEdges, nEdges);
  Eigen::MatrixXd idE = Eigen::MatrixXd::Identity(nEdges, nEdges);
  Eigen::MatrixXd grad0 = Eigen::MatrixXd(d0);

  Eigen::VectorXd u(nEdges);
  Eigen::VectorXd force(nEdges);
  Eigen::VectorXd p(nVerts);
  std::vector<char>  orient(nEdges);
  Eigen::VectorXd x(nEdges+nVerts);
  Eigen::VectorXd b(nEdges+nVerts);

  double relative_error = 42;

  scanf(" %d", &matrix_calc_permission);
  if (matrix_calc_permission == 1) {

    div1 = Eigen::MatrixXd(s0m1*(d0.transpose())*s1);
    l1 = (Eigen::MatrixXd(s1m1)*(Eigen::MatrixXd(d1.transpose()))*(Eigen::MatrixXd(s2))*(Eigen::MatrixXd(d1)));

    for(int i = 0; i < nVerts; i++) {
      for(int j = 0; j < nEdges; j++) {
        A(i, j) = div1(i, j); // top left of A
      }
    }


    for(int i = 0; i < nVerts; i++) {
      for(int j = 0; j < nVerts; j++) { 
        A(i, j + nEdges) = 0.0; // top right of A
      }
    }


    // Fluid parameters
    double eta = 0.0;
    double gamma = 0.0;
    printf("\nGive a value for eta (viscosity): ");
    scanf(" %lf", &eta);
    printf("\nGive a value for gamma (surface tension): ");
    scanf(" %lf", &gamma);

    for (int i = 0; i < nEdges; i++) {
      for( int j = 0; j < nEdges; j++) {
        A(i + nVerts, j) = eta*(-l1(i, j) + 2.0*Kedge(i)*idE(i, j)) - gamma*idE(i, j); // bottom left of A
      }
    }


    for(int i = 0; i < nEdges; i++) {
      for(int j = 0; j < nVerts; j++) {
        A(i + nVerts, j + nEdges) = -grad0(i, j); // bottom right of A
      }
    }


    // velocity and force vectors
    int force_index = -1;
    printf("\nType edge index to apply force (0 to %d inclusive): ", nEdges - 1);
    scanf (" %d", &force_index);  // choosing the edge to apply force
    while ((force_index < 0) || (force_index >= nEdges)) {
      printf("Invalid index, again please: ");
      scanf (" %d", &force_index);  
    }
    
    double force_magnitude = 0.0;
    printf("Force magnitude: ");
    scanf(" %lf", &force_magnitude);  // choosing size of force

    geometry->requireEdgeIndices();
    for (Edge e: mesh->edges()) {
      int i = geometry->edgeIndices[e];
      force(i) = 0.0;
      Vertex vertex_tail = e.halfedge().tailVertex();
      Vertex vertex_head = e.halfedge().tipVertex();
      int p = geometry->vertexIndices[vertex_tail];
      int q = geometry->vertexIndices[vertex_head];
      
      if (p < q) {
        orient[i] = true;
      } else {
        orient[i] = false; 
      }
    
      if (i == force_index) {
        force(i) = force_magnitude;
        
      }
    }


    for (int j = 0; j < nEdges; j++) {  // is this to ignore the edge with the applied force? 
       A(force_index, j) = 0.0;
    }

    A(force_index, nEdges) = 1.0; // from A(nVerts, nEdges) = 1.0;


    for(int i = 0; i < nVerts; i++) {
      b(i) = 0.0;
    }


    for(int i = 0; i < nEdges; i++) {
      b(i + nVerts) = -force(i);
    }


    x = A.colPivHouseholderQr().solve(b); // matrix inversion
  

    for (int i = 0; i < nEdges; i++) {
      u(i) = x(i);
    }

    for(int i = 0; i < nVerts; i++) {
      p(i) = x(i + nEdges);
    }


    relative_error = (A*x - b).norm() / b.norm();
    std::cout << "Flow field solved successfully\n";
  }
  

  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // Append run data into textfile (with errors)
  int print_permission;
  printf("\nType 1 to print data into data.txt, 0 otherwise: ");

  scanf(" %d", &print_permission);
  if (print_permission == 1) {
    std::ofstream outstream;
    outstream.open("data.txt", std::ios::app);
    if (outstream.is_open()) {
      time_t now = time(NULL);

      outstream 
        << "Time run: " << ctime(&now) << "File: " << args::get(inputFilename) << "\n"
        << nVerts << " vertices, " << nEdges << " edges" << "\n"
        << "L2 Norm of Vertex_Gaussian_curvature is: " << Kvertex_error_avg << "\n"
        << "Max L2 Error on Gaussian curvature: " << Kvertex_error_max << "\n";

      if (matrix_calc_permission == 1) {
        outstream
          << "The relative error is: " << relative_error << "\n"
          << "L2 Norm of div1 is: " << div1.norm() << "\n"
          << "L2 Norm of grad1 is: " << grad0.norm() << "\n"
          << "L2 Norm of Hodge star 1 is: " << s1.norm() << "\n"
          << "L2 Norm of Hodge star inverse 1 is: " << s1m1.norm() << "\n"
          << "L2 Norm of Hodge star 2 is: " << s2.norm() << "\n"
          << "L2 Norm of Hodge star inverse 2 is: " << s2m1.norm() << "\n"
          << "L2 Norm of Exterior derivative 1 is: " << d1.norm() << "\n";
      }

      outstream << "\n";
      outstream.close();
      std::cout << "Content appended successfully\n";

    } else {
        std::cerr << "Unable to open file for appending\n";
    }
  }

  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // Checks before feeding into polyscope

  for (int i = 0; i < nEdges; i++) {
   std::cout << u(i) << ".";
  }

  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  psMesh->addVertexTangentVectorQuantity("VF", vField, vBasisX, vBasisY);
  psMesh->addVertexScalarQuantity("Gaussian curvature", gaussianCurve, polyscope::DataType::SYMMETRIC);
  psMesh->addOneFormTangentVectorQuantity("force", force, orient);
  psMesh->addOneFormTangentVectorQuantity("velocity", u, orient);
  psMesh->addVertexScalarQuantity("Pressure", p);
  psMesh->addEdgeScalarQuantity("Edge Kg", Kedge);


  // Give control to the polyscope gui
  polyscope::show();


  return EXIT_SUCCESS;  
}
