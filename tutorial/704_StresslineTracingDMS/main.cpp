#include <igl/edge_topology.h>
#include <igl/read_triangle_mesh.h>
#include <igl/readOFF.h>
#include <directional/streamlines.h>
#include <directional/power_field.h>
#include <directional/power_to_raw.h>
#include <directional/directional_viewer.h>
#include <directional/read_raw_field.h>


// Mesh
Eigen::MatrixXd V;
Eigen::MatrixXi F;
Eigen::MatrixXcd powerField;
Eigen::MatrixXd rawField;
Eigen::MatrixXd P1,P2;

/*directional::StreamlineData sl_data;
directional::StreamlineState sl_state;*/

int N=3;         // degree of the vector field
int anim_t = 0;
int anim_t_dir = 1;


bool pre_draw(igl::opengl::glfw::Viewer &iglViewer)
{
  igl::opengl::glfw::Viewer* iglViewerPointer=&iglViewer;
  directional::DirectionalViewer* directional_viewer = static_cast<directional::DirectionalViewer*>(iglViewerPointer);
  using namespace Eigen;
  using namespace std;
  
  if (!iglViewer.core().is_animating)
    return false;
  
  directional_viewer->advance_streamlines();
  
  anim_t += anim_t_dir;
  
  return false;
}

bool key_down(igl::opengl::glfw::Viewer &viewer, unsigned char key, int modifier)
{
  if (key == ' ')
  {
    viewer.core().is_animating = !viewer.core().is_animating;
    return true;
  }
  return false;
}


// CLI help menu
int print_help() {
  std::cout << "N-RoSy field to mesh *mandatory* arguments. Type --help for additional info. " << std::endl
            << "<structure>: The name of the structure of interest." << std::endl
            << "<algorithm>: The name of the clustering algorithm used." << std::endl
            << "<n_clusters>: The name of clusters generated." << std::endl
            << "<rawfield>: The name of the .rawfield file that stores the N-RoSy field to use to generate a new mesh."<< std::endl
            << "<degree>: The degree (N) of the first N-RoSy field."<< std::endl;
  return 0;
}


int main(int argc, char *argv[])
{
  using namespace Eigen;
  using namespace std;
  
  directional::DirectionalViewer viewer;

  // CLI
  if (argc == 2 && strcmp(argv[1], "--help")==0)
  {
    print_help();
    return 0;
  }

  // parse CLI arguments
  std::string structure_name = argv[1];
  std::string algorithm_name = argv[2];
  std::string n_clusters = argv[3];
  std::string rawFieldName = argv[4];
  std::string degreeField = argv[5];

  // base data folder location
  std::string MESH_NAME = structure_name + "_" + "k" + n_clusters;
  std::string EXPERIMENT_NAME = structure_name + "/" + algorithm_name + "/" + "k" + n_clusters + "/";
  std::string PATH_DATA = "/Users/arpj/princeton/phd/papers/as_psf/experiments/";
  std::string PATH_OFF = PATH_DATA + EXPERIMENT_NAME + "directional/";
  std::string PATH_RAWFIELD = PATH_OFF;

  // read input mesh
  std::cout<<"Reading mesh from .off file"<<std::endl;
  igl::readOFF(PATH_OFF + structure_name + ".off", V, F);

  // load in raw field
  // numbers to string conversions
  int N = 0;
  std::cout<<"Reading N-RoSy field from "<<rawFieldName<<"_"<<degreeField<<"rosy.rawfield"<<std::endl;
  directional::read_raw_field(PATH_RAWFIELD + rawFieldName + "_" + degreeField + "rosy.rawfield", N, rawField);

  // Load a mesh in OFF format
  // igl::readOFF(TUTORIAL_SHARED_PATH "/lion.off", V, F);
  // Create a Vector Field
  // directional::power_field(V, F, Eigen::VectorXi(), Eigen::MatrixXd(), Eigen::MatrixXd(), 3, powerField);
  
  // Convert it to raw field
  // directional::power_to_raw(V,F,powerField,3,rawField, true);
  
  //triangle mesh
  viewer.set_mesh(V,F);
  viewer.set_field(rawField);
  Eigen::MatrixXd fieldColors=directional::DirectionalViewer::indexed_glyph_colors(rawField);
  viewer.set_field_colors(fieldColors);
  viewer.toggle_field(false);
  viewer.init_streamlines();
  viewer.advance_streamlines();  //to get the initial step

  // Viewer Settings
  viewer.callback_pre_draw = &pre_draw;
  viewer.callback_key_down = &key_down;
  viewer.core().is_animating = false;
  viewer.core().animation_max_fps = 30.;
  
  cout << "Press [space] to toggle animation" << endl;
  viewer.launch();
}
