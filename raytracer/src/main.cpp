#include "parser.hpp"
#include "raytracer.hpp"
#include <cstdio>
#include <cstdlib>
#include <cfloat>
#include <cmath>
#include <algorithm>
#include <string>
#include <fstream>
#include <vector>
#include <iostream>
#include <sstream>
#include <map>
#include <vector>


int main(int argc, char **argv)
{
	{
		GetPot getpot2(argc, argv);
		getpot.absorb(getpot2);
	}

	string scenefile;	//If aren't specified, it will render the default scene with only a cube
	string outfilename("scenes/default_out.bmp");	//NOTE: here we only support .bmp
	string outfilename_depth("scenes/default_depth_out.bmp");

    // Scene filename specified
	if (argc > 1)
	{
		scenefile = string(argv[1]);
		int dot = scenefile.find_last_of('.');
		outfilename = scenefile.substr(0, dot) + "_out.bmp";
		outfilename_depth = scenefile.substr(0, dot) + "_depth_out.bmp";
	}

    std::cout << "Rendering " << (scenefile.empty() ? "default scene" : scenefile) << std::endl;
    std::cout << "Output to " << outfilename << std::endl;
   
	Raytracer raytracer;
	if (scenefile.empty())
	{
		//render default scene
		raytracer.render(
			outfilename.c_str(),
			outfilename_depth.c_str(),
			Scene()
			);
	}
	else
	{
		// Parse the scene file
		std::ifstream  input_stream(scenefile);
		Parser parser(input_stream);
		if (!parser.parse()) {
			puts("Scene file can't be parsed. Use default scene.");
			raytracer.render(
				outfilename.c_str(),
				outfilename_depth.c_str(),
				Scene()
				);
		}
		else
		{
			// Render the input scene with our raytracer.
			raytracer.render(
				outfilename.c_str(),
				outfilename_depth.c_str(),
				parser.scene
				);
		}
	}
	
	// Use if you're running visual studio:
	//system("pause");
	return 0;
}

