#include <iostream>
#include <sstream>

#include <Eigen/Core>
#include <igl/opengl/glfw/Viewer.h>
#include <igl/read_triangle_mesh.h>
#include <igl/per_face_normals.h>
#include <igl/unproject_onto_mesh.h>
#include <igl/edge_topology.h>
#include <igl/cut_mesh.h>
#include <directional/read_raw_field.h>
#include <directional/write_raw_field.h>
#include <directional/curl_matching.h>
#include <directional/combing.h>
#include <directional/setup_integration.h>
#include <directional/integrate.h>
#include <directional/cut_mesh_with_singularities.h>
#include <directional/branched_isolines.h>
#include <directional/mesh_function_isolines.h>
#include <directional/setup_mesh_function_isolines.h>
#include <directional/directional_viewer.h>
#include "polygonal_write_OFF.h"

#define NUM_N 2

int currN = 0;
int N[NUM_N];
int degreeField[NUM_N];
double lengthRatio[NUM_N];
bool roundSeams[NUM_N];
bool integralSeamless[NUM_N];
char *rawFieldName[NUM_N];

Eigen::MatrixXi FMeshWhole, FMeshCut[NUM_N];
Eigen::MatrixXd VMeshWhole, VMeshCut[NUM_N];
Eigen::MatrixXd rawField[NUM_N], combedField[NUM_N];
Eigen::VectorXd effort[NUM_N], combedEffort[NUM_N];
Eigen::VectorXi matching[NUM_N], combedMatching[NUM_N];
Eigen::MatrixXi EV, FE, EF;
Eigen::VectorXi DPolyMesh[NUM_N];
Eigen::MatrixXi FPolyMesh[NUM_N];
Eigen::MatrixXd VPolyMesh[NUM_N];
Eigen::VectorXi singIndices[NUM_N], singVertices[NUM_N];
Eigen::MatrixXd NFunction[NUM_N], NCornerFunction[NUM_N];
directional::DirectionalViewer viewer;

typedef enum {FIELD, INTEGRATION} ViewingModes;
ViewingModes viewingMode=FIELD;


void update_viewer()
{
  for (int i=0; i<NUM_N; i++){
    viewer.toggle_field(false,i);
    viewer.toggle_singularities(false,i);
    viewer.toggle_seams(false,i);
    viewer.toggle_isolines(false,i);
  }
  if (viewingMode==FIELD){
    viewer.toggle_field(true,currN);
    viewer.toggle_singularities(true,currN);
    viewer.toggle_seams(true,currN);
    viewer.toggle_isolines(false,currN);
  } else {
    viewer.toggle_field(false,currN);
    viewer.toggle_singularities(false,currN);
    viewer.toggle_seams(false,currN);
    viewer.toggle_isolines(true,currN);
  }
}


// Handle keyboard input
bool key_down(igl::opengl::glfw::Viewer& viewer, int key, int modifiers)
{
  switch (key)
  {
      // Select vector
    case '1': viewingMode = FIELD; break;
    case '2': viewingMode = INTEGRATION; break;
    case '3': currN=(currN+1)%NUM_N; break;
  }
  update_viewer();
  return true;
}


// Convert string to boolean
bool string_to_bool(char *arg)
{
  if (!strcmp(arg, "0"))
    {
      return false;
    }
    else if (!strcmp(arg, "1"))
    {
      return true;
    }
    else
    {
      return 1;
    }
}


// CLI help menu
int print_help() {
  std::cout << "N-RoSy field to mesh *mandatory* arguments. Type --help for additional info. " << std::endl
            << "<mesh>: The name of the .OFF file that stores the base mesh of interest." << std::endl
            << "<rawfield_1>: The name of the first .rawfield file that stores the N-RoSy field to use to generate a new mesh."<< std::endl
            << "<degree_1>: The degree (N) of the first N-RoSy field."<< std::endl
            << "<length_ratio_1> [0.05]: Controls parametrization and mesh density of the first raw field (Smaller value -> denser mesh)." << std::endl
            << "<round_seams_1> [0]: Boolean for whether to round seams or round singularities." << std::endl
            << "<integral_seamless_1> [1]: Boolean flag for whether do full translational seamless." << std::endl
            << "<rawfield_2>: The name of the second .rawfield file to use to generate another new mesh" << std::endl
            << "<degree_2>: The degree (N) of the first N-RoSy field."<< std::endl
            << "<length_ratio_2> [0.05]: Controls parametrization and mesh density of the second rawfield (Smaller value -> denser mesh)." << std::endl
            << "<round_seams_2> [0]: Boolean for whether to round seams or round singularities." << std::endl
            << "<integral_seamless_2> [1]: Boolean flag for whether do full translational seamless." << std::endl;
  return 0;
}


int main(int argc, char *argv[])
{
  std::cout <<
  "  1  Loaded field" << std::endl <<
  "  2  Show isoline mesh" << std::endl <<
  "  3  change between different N" << std::endl;  // TODO: support multiple rawfield loading

  // CLI
  if (argc == 2 && strcmp(argv[1], "--help")==0)
  {
    print_help();
    return 0;
  }

  // parse CLI arguments
  std::string mesh_name = argv[1];

  // TODO: Currently supports the ingestion of only two hard-coded vector fields
  // field 0
  rawFieldName[0] = argv[2];
  degreeField[0] = atoi(argv[3]);
  lengthRatio[0] = atof(argv[4]);
  roundSeams[0] = string_to_bool(argv[5]);
  integralSeamless[0] = string_to_bool(argv[6]);
  // field 1
  rawFieldName[1] = argv[7];
  degreeField[1] = atoi(argv[8]);
  lengthRatio[1] = atof(argv[9]);
  roundSeams[1] = string_to_bool(argv[10]);
  integralSeamless[1] = string_to_bool(argv[11]);

  // log to console or not
  bool verbose=true;

  // base data folder location
  std::string PATH_DATA = "/Users/arpj/code/libraries/directional_clustering/data/";
  std::string PATH_OFF = PATH_DATA + "off/";
  std::string PATH_RAWFIELD = PATH_DATA + "rawfield/";

  // read input mesh
  std::cout<<"Reading mesh from .off file"<<std::endl;
  igl::readOFF(PATH_OFF + mesh_name + ".off", VMeshWhole, FMeshWhole);

  // build mesh edge topology
  igl::edge_topology(VMeshWhole, FMeshWhole, EV, FE, EF);

  //combing and cutting
  for (int i=0; i<NUM_N; i++){
    std::cout<<"------------"<<std::endl;

    // numbers to string conversions
    std::string str_df = std::to_string(degreeField[i]);
    std::string str_lr = std::to_string(lengthRatio[i]);
    std::string str_rs = std::to_string(roundSeams[i]);
    std::string str_is = std::to_string(integralSeamless[i]);
    // remove trailing zeros from length ratio to string conversion
    str_lr.erase(str_lr.find_last_not_of('0') + 1, std::string::npos);
    std::string mesh_generated_name = mesh_name+"_"+rawFieldName[i]+"_lr"+str_lr+"_sround"+str_rs+"_sintegral"+str_is+"_generated.off";

    // load in raw field
    std::cout<<"Reading N-RoSy field #"<<i<<" from "<<rawFieldName[i]<<"_"<<str_df<<"rosy.rawfield"<<std::endl;
    directional::read_raw_field(PATH_RAWFIELD + mesh_name + "_" + rawFieldName[i] + "_" + str_df + "rosy.rawfield", N[i], rawField[i]);

    // principal matching
    std::cout<<"Principal matching #"<<i<<std::endl;
    directional::principal_matching(VMeshWhole, FMeshWhole,EV, EF, FE, rawField[i], matching[i], effort[i],singVertices[i], singIndices[i]);

    // set up integration
    std::cout<<"Setting up Integration Data #"<<i<<": (lengthRatio="<<str_lr<<", roundSeams="<<roundSeams[i]<<", integralSeamless="<<integralSeamless[i]<<")"<<std::endl;
    directional::IntegrationData intData(N[i]);
    directional::setup_integration(VMeshWhole, FMeshWhole,  EV, EF, FE, rawField[i], matching[i], singVertices[i], intData, VMeshCut[i], FMeshCut[i], combedField[i], combedMatching[i]);

    // intData.lengthRatio Controls parametrization and mesh density
    // Small values mean denser meshes (e.g. intData.lengthRatio=0.02).
    // Large values lead to coarser meshes (e.g. intData.lengthRatio=0.10).
    intData.lengthRatio=lengthRatio[i];  // Defaults to 0.02
    intData.roundSeams=roundSeams[i];  // Whether to round seams or round singularities. Defaults to false
    intData.integralSeamless=integralSeamless[i];  // Whether to do full translational seamless. Defaults to true
    intData.verbose=false;  // Output the integration log. Defaults to fasle
  
    std::cout<<"Solving integration for N="<<N[i]<<std::endl;
    directional::integrate(VMeshWhole, FMeshWhole, FE, combedField[i],  intData, VMeshCut[i], FMeshCut[i], NFunction[i], NCornerFunction[i]);
    std::cout<<"Done!"<<std::endl;
    
    // mesh only if N-RoSy field's degree is greater than 2
    if (degreeField[i] > 2){
      // setting up mesh data from integration data
      std::cout<<"Setting up mesh data from integration data #"<<i<<std::endl;
      directional::MeshFunctionIsolinesData mfiData;
      directional::setup_mesh_function_isolines(VMeshCut[i], FMeshCut[i], intData, mfiData);

      // meshing and saving
      std::cout<<"Meshing"<<std::endl;
      directional::mesh_function_isolines(VMeshWhole, FMeshWhole, EV, EF, FE, mfiData, verbose, VPolyMesh[i], DPolyMesh[i], FPolyMesh[i]);
      hedra::polygonal_write_OFF(PATH_OFF + mesh_generated_name, VPolyMesh[i], DPolyMesh[i], FPolyMesh[i]);
      std::cout<<"Exported .OFF successfully!"<<std::endl;
    }
    
    viewer.set_mesh(VMeshWhole, FMeshWhole, i);
    viewer.set_field(combedField[i], directional::DirectionalViewer::indexed_glyph_colors(combedField[i]), i);
    viewer.set_singularities(singVertices[i], singIndices[i]);
    viewer.set_seams(combedMatching[i], i);
    viewer.set_isolines(VMeshCut[i], FMeshCut[i], NFunction[i], i);
  }

  std::cout<<"------------"<<std::endl;
  std::cout<<"Launching viewer. Enjoy!"<<std::endl;
  update_viewer();
  viewer.callback_key_down = &key_down;
  viewer.launch();
}
