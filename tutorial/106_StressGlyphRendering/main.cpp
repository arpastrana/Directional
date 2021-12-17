#include <directional/directional_viewer.h>
#include <directional/read_raw_field.h>
// #include <directional/read_singularities.h>

int N;
Eigen::MatrixXi F, EV, FE, EF;
Eigen::MatrixXd V, rawField;
// Eigen::VectorXi singVertices, singIndices;
directional::DirectionalViewer viewer;
bool showField=true;
// showSingularities=true;

bool key_down(igl::opengl::glfw::Viewer& iglViewer, int key, int modifiers)
{
  igl::opengl::glfw::Viewer* iglViewerPointer=&iglViewer;
  directional::DirectionalViewer* directional_viewer = static_cast<directional::DirectionalViewer*>(iglViewerPointer);
  switch (key)
  {
    case '1': showField=!showField; directional_viewer->toggle_field(showField); break;
    // case '2': showSingularities=!showSingularities; directional_viewer->toggle_singularities(showSingularities); break;;
    default: return false;
  }
  return true;
}

// CLI help menu
int print_help() {
  std::cout << "Viewer *mandatory* arguments. Type --help for additional info. " << std::endl
            << "<mesh_name>: The name of the .OFF file that stores the base mesh of interest." << std::endl
            << "<rawfield_name>: The name of the .rawfield file that stores the N-RoSy field to display."<< std::endl;
  return 0;
}

int main(int argc, char *argv[])
{
  std::cout <<"1    Show/hide Field" << std::endl;
  // std::cout <<"2    Show/hide Singularities" << std::endl;

  // CLI
  if (argc == 2 && strcmp(argv[1], "--help")==0)
  {
    print_help();
    return 0;
  }

  std::string mesh_name = argv[1];
  std::string field_name = argv[2];
  // igl::readOFF(TUTORIAL_SHARED_PATH "/bumpy.off", V, F);
  // directional::read_raw_field(TUTORIAL_SHARED_PATH "/bumpy.rawfield", N, rawField);
  // directional::read_singularities(TUTORIAL_SHARED_PATH "/bumpy.sings", N, singVertices, singIndices);

  // base data folder location
  std::string PATH_DATA = "/Users/arpj/code/libraries/directional_clustering/data/";
  std::string PATH_OFF = PATH_DATA + "off/";
  std::string PATH_RAWFIELD = PATH_DATA + "rawfield/";

  std::string PATH_MESH = PATH_OFF + mesh_name + ".off";
  std::string PATH_FIELD = PATH_RAWFIELD + mesh_name + "_" + field_name + ".rawfield";

  std::cout<<"Reading mesh from .off file"<<std::endl;
  igl::readOFF(PATH_MESH, V, F);
  std::cout<<PATH_FIELD<<std::endl;
  directional::read_raw_field(PATH_FIELD, N, rawField);

  std::cout<<"Setting viewer"<<std::endl;
  directional::DirectionalViewer viewer;
  
  viewer.set_mesh(V,F);
  viewer.set_field(rawField);
  // viewer.set_singularities(singVertices, singIndices);
  viewer.toggle_mesh_edges(false);
  
  viewer.callback_key_down = &key_down;
  viewer.launch();
  std::cout<<"Enjoy!"<<std::endl;
}
