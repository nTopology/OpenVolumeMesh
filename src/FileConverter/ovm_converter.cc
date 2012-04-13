#include <fstream>
#include <sstream>

#include <boost/spirit/include/qi.hpp>

#include <OpenVolumeMesh/FileManager/FileManager.hh>

#include "MeshGenerator.hpp"
#include "Grammars/Tetmesh.hpp"

int main(int _argc, char* _argv[]) {

    if(_argc != 3) {
        std::cerr << "You need to specify a source and a target file!" << std::endl;
        return -1;
    }

    std::ifstream iff(_argv[1], std::ios::in);
    //iff.unsetf(std::ios::skipws);

    if(!iff.good()) {
        std::cerr << "Could not open file " << _argv[1] << " for reading!" << std::endl;
        return -1;
    }

    MeshGenerator::PolyhedralMesh mesh;
    MeshGenerator generator(mesh);

    std::ostringstream oss;
    oss << iff.rdbuf();

    std::string fileContent = oss.str();

    // Instantiate grammar object
    tetmesh_grammar<std::string::iterator> grammar(generator);

    // Use iterator to parse file data

    //bool r = boost::spirit::qi::parse(fileContent.begin(), fileContent.end(), grammar);
    bool r = boost::spirit::qi::phrase_parse(fileContent.begin(), fileContent.end(), grammar, qi::space);

    if(r) {
        std::cout << "Parsed all data successfully!" << std::endl;
    } else {
        std::cout << "Parsing failed!" << std::endl;
    }

    std::cerr << "Added " << mesh.n_vertices() << " vertices!" << std::endl;
    std::cerr << "Added " << mesh.n_edges() << " edges!" << std::endl;
    std::cerr << "Added " << mesh.n_faces() << " faces!" << std::endl;
    std::cerr << "Added " << mesh.n_cells() << " cells!" << std::endl;

    OpenVolumeMesh::IO::FileManager fileManager;

    // Write mesh to file
    fileManager.writeFile(_argv[2], mesh);

    return 0;
}
