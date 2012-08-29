#include "GetPot"

#include <string.h>
#include <fstream>
#include <iostream>
#include <ostream>
#include <istream>
#include <sstream>
#include <stdio.h>

int main ( int argc, char* argv[] ) {
    GetPot args( argc, argv );

    std::string partialpath = args.follow( "sphere", 1, "-file" );

    std::string inextension = ".mesh";
    std::string inpath = partialpath + inextension;

    std::string outextension = ".tetra";
    std::string outpath = partialpath + outextension;

    std::ifstream infile( inpath.c_str(), std::ios::in | std::ios::ate );
    if( !infile.is_open() ) return 1;

    std::ofstream outfile( outpath.c_str(), std::ios::out );
    if( !outfile.is_open() ) return 1;

    infile.seekg( 0, std::ios::beg );

    char readbuffer[1024];

    unsigned int verts = 0, tetras = 0;
    unsigned int mode = 0;

    std::string token, delim = " ";

    int value;
    unsigned int counter;

    while( !infile.eof() ) {
        infile.getline( readbuffer, 1024 );

        switch( mode ) {
        case 1:
            if( counter == 0 ) {
                // this line must contain the number of vertices to read
                verts = (unsigned int)atoi( readbuffer );
            } else {
                if( counter <= verts ) {
                    // this must be a vertex definition
                    token = strtok( readbuffer, delim.c_str() );
                    outfile << "v" << delim << token;
                    token = strtok( NULL, delim.c_str() );
                    outfile << delim << token;
                    token = strtok( NULL, delim.c_str() );
                    outfile << delim << token << "\n";
                } else {
                    // otherwise done parsing vertices so change to default mode
                    mode = 0;
                }
            }
            counter++;
            break;
        case 2:
            if( counter == 0 ) {
                // this line must contain the number of tetrahedrons to read
                tetras = (unsigned int)atoi( readbuffer );
            } else {
                if( counter <= tetras ) {
                    // this must be a tetrahedron definition
                    token = strtok( readbuffer, delim.c_str() );
                    value = atoi( token.c_str() ) - 1;
                    outfile << "t" << delim << value;
                    token = strtok( NULL, delim.c_str() );
                    value = atoi( token.c_str() ) - 1;
                    outfile << delim << value;
                    token = strtok( NULL, delim.c_str() );
                    value = atoi( token.c_str() ) - 1;
                    outfile << delim << value;
                    token = strtok( NULL, delim.c_str() );
                    value = atoi( token.c_str() ) - 1;
                    outfile << delim << value << "\n";
                } else {
                    // otherwise done parsing tetrahedrons so change to default mode
                    mode = 0;
                }
            }
            counter++;
            break;
        default:
            break;
        }

        if( strcmp( readbuffer, "Vertices" ) == 0 ) {
            //std::cout << "Vertices Encountered" << "\n";
            counter = 0;
            mode = 1;
        }

        if( strcmp( readbuffer, "Tetrahedra" ) == 0 ) {
            //std::cout << "Tetrahedra Encountered" << "\n";
            counter = 0;
            mode = 2;
        }


        //std::cout << readbuffer << "\n";
    }

    infile.close();
    outfile.close();

    //delete readbuffer;

	return 0;
}
