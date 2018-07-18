#include "parser.hpp"

#include <algorithm>
#include <cctype>
#include <climits>
#include <deque>
#include <string>
#include <fstream>
#include <sstream>
#include <map>
#include <vector>
#include <iostream>
#include <cmath>
#include <cfloat>


#define LEX_ERROR(error) { \
    std::stringstream ss; \
    ss << error; \
    throw ss.str(); \
}


bool Token::operator==(Token const &other) const {
	if (type != other.type) {
		return false;
	}
	switch (type) {
	case NUMBER:
		return number == other.number;
	case STRING:
	case NAME:
		return string == other.string;
	default:;
	}
	return true;
}


std::ostream& operator<<(std::ostream &out, Token const &token) {
	switch (token.type) {
	case STRING:
		out << "STRING:\"" << token.string << "\"";
		break;
	case NUMBER:
		out << "NUMBER:" << token.number;
		break;
	case NAME:
		out << "NAME:" << token.string;
		break;
	case ARRAY_BEGIN:
		out << "ARRAY_BEGIN";
		break;
	case ARRAY_END:
		out << "ARRAY_END";
		break;
	case END_OF_FILE:
		out << "EOF";
		break;
	case ERROR:
		out << "ERROR";
		break;
	default:
		out << "UNKNOWN";
		break;
	}
	return out;
}


Token Lexer::_processStream(void) {

	// Handle immediate/error conditions.
	if (_input.eof()) {
		return Token(END_OF_FILE);
	}
	if (!_input.good()) {
		return Token(ERROR);
	}

	// The next character in the stream.
	char c = _input.peek();

	// Strip whitespace and comments.
	bool did_strip_something;
	do {
		did_strip_something = false;

		// Strip whitespace.
		while (isspace(c)) {
			_input.ignore(1);
			c = _input.peek();
			did_strip_something = true;
		}

		// Strip comments.
		if (c == '#') {
			do {
				c = _input.get();
			} while (c != '\r' && c != '\n');
			did_strip_something = true;
			c = _input.peek();
		}

	} while (did_strip_something);

	// Arrays.
	switch (c) {
	case '[':
		_input.ignore(1);
		return Token(ARRAY_BEGIN);
	case ']':
		_input.ignore(1);
		return Token(ARRAY_END);
	}

	// Strings.
	if (c == '"') {
		_input.get();
		std::string value;
		bool finished = false;
		while (!finished) {
			c = _input.get();
			switch (c) {
			case '"':
				finished = true;
				break;
				// TODO: handle escapes.
			default:
				value += c;
			}
		}
		return Token(STRING, value);
	}

	// Numbers; try to read one and if it doesn't work reset the error
	// state and carry on.
	double number;
	_input >> number;
	if (!_input.fail()) {
		return Token(number);
	}
	_input.clear();

	// Names.
	std::string name;
	_input >> name;
	if (name.size()) {
		return Token(NAME, name);
	}

	// There is nothing left to read.
	return Token(END_OF_FILE);

}


Token Lexer::peek(unsigned int index) {
	if (_buffer.size() <= index) {
		_buffer.push_back(_processStream());
	}
	return _buffer[index];
}


Token Lexer::next() {
	if (!_buffer.size()) {
		return _processStream();
	}
	Token token = _buffer.front();
	_buffer.pop_front();
	return token;
}


void Lexer::skip(unsigned int count) {
	for (unsigned int i = 0; i < count; i++) {
		next();
	}
}


std::string Lexer::getName() {
	Token token = peek();
	if (token.type != NAME) {
		LEX_ERROR("expected NAME; got " << token)
	}
	skip(1);
	return token.string;
}


std::vector<double> Lexer::getNumbers(unsigned int min, unsigned int max) {

	std::vector<double> values;

	bool is_array = peek().type == ARRAY_BEGIN;
	if (is_array) {
		skip(1);
	}

	// UGLY HACK: Only here so we have something to `continue` to.
	do {

		Token token = next();

		switch (token.type) {
		case NUMBER:
			values.push_back(token.number);
			if (!is_array) {
				break;
			}
			continue;

		case ARRAY_END:

			if (is_array) {
				break;
			}
			// falling through here
		default:
			LEX_ERROR("expected NUMBER; got " << token);
		}

		if (values.size() >= min && values.size() <= max) {
			return values;
		}
		LEX_ERROR("expected " << min << " to " << max << " NUMBERs; got " << values.size())

	} while (true);

}


double Lexer::getNumber() {
	return getNumbers(1, 1)[0];
}


std::string Lexer::getString() {
	Token token = peek();
	if (token.type != STRING) {
		LEX_ERROR("expected STRING; got " << token)
	}
	skip(1);
	return token.string;
}


ParamList Lexer::getParamList(unsigned int min, unsigned int max) {
	ParamList map;
	while (true) {
		std::string key;
		try {
			key = getString();
		}
		catch (std::string) {
			return map;
		}
		map[key] = getNumbers(min, max);
	}
}


bool Parser::parse() {
    
    transformStack.push_back(Matrix::identity());

    while (true) {
        
        Token token = lexer.peek();
        switch (token.type) {
            case END_OF_FILE:
                return true;
            case ERROR:
                std::cerr << "parsing failed due to lexing error" << std::endl;
                return false;
            default:;
        }


        std::string name;
        try {
            name = lexer.getName();
        } catch (std::string e) {
            std::cerr << "parsing failed to get next command due to: " << e << std::endl;
            return false;
        }


        #define HANDLE_NAME(_name) if (name == #_name) { parse##_name(); continue; }

        try {
            HANDLE_NAME(PushMatrix)
            HANDLE_NAME(PopMatrix)
            HANDLE_NAME(Translate)
            HANDLE_NAME(Scale)
            HANDLE_NAME(Rotate)
            HANDLE_NAME(Sphere)
            HANDLE_NAME(Plane)
            HANDLE_NAME(Mesh)
            HANDLE_NAME(Conic)
            HANDLE_NAME(Material)
            HANDLE_NAME(PointLight)
            HANDLE_NAME(Dimensions)
            HANDLE_NAME(Perspective)
            HANDLE_NAME(LookAt)
        } catch (std::string e) {
            std::cerr << "parsing failed on command \"" << name << "\" due to: " << e << std::endl;
            return false;
        }

        #undef HANDLE_NAME

        // Should have `continue`d by this point.
        std::cerr << "parsing failed due to unknown command \"" << name << "\"" << std::endl;
        return false;

    }
}


void Parser::parsePopMatrix() {
    transformStack.pop_back();
    if (!transformStack.size()) {
        throw std::string("popped off bottom of stack");
    }
}


void Parser::parsePushMatrix() {
    transformStack.push_back(transformStack.back());
}


void Parser::parseTranslate() {
    // Need to store these in variables because if we pass them directly to
    // Matrix::whatever(...) we are no guarunteed they are called in order.
    double x = lexer.getNumber();
    double y = lexer.getNumber();
    double z = lexer.getNumber();
    transformStack.back() *= Matrix::translation(x, y, z);
}


void Parser::parseScale() {
    // Need to store these in variables because if we pass them directly to
    // Matrix::whatever(...) we are no guarunteed they are called in order.
    double x = lexer.getNumber();
    double y = lexer.getNumber();
    double z = lexer.getNumber();
    transformStack.back() *= Matrix::scale(x, y, z);
}


void Parser::parseRotate() {
    // Need to store these in variables because if we pass them directly to
    // Matrix::whatever(...) we are no guarunteed they are called in order.
    double a = lexer.getNumber();
    double x = lexer.getNumber();
    double y = lexer.getNumber();
    double z = lexer.getNumber();
    transformStack.back() *= Matrix::rotation(deg2rad(a), Vector(x, y, z));
}


void Parser::parseDimensions() {
    scene.resolution[0] = static_cast<int>(lexer.getNumber());
    scene.resolution[1] = static_cast<int>(lexer.getNumber());
}


void Parser::parsePerspective() {
    scene.camera.fov   = lexer.getNumber();
    scene.camera.aspect = lexer.getNumber();
    scene.camera.zNear  = lexer.getNumber();
    scene.camera.zFar   = lexer.getNumber();
}


void Parser::parseLookAt() {
    scene.camera.position[0] = lexer.getNumber();
    scene.camera.position[1] = lexer.getNumber();
    scene.camera.position[2] = lexer.getNumber();
    scene.camera.center  [0] = lexer.getNumber();
    scene.camera.center  [1] = lexer.getNumber();
    scene.camera.center  [2] = lexer.getNumber();
    scene.camera.up      [0] = lexer.getNumber();
    scene.camera.up      [1] = lexer.getNumber();
    scene.camera.up      [2] = lexer.getNumber();
}


void Parser::parseMaterial() {
    std::string name = lexer.getString();
    ParamList params = lexer.getParamList(1, 4);
    scene.materials[name] = Material(params);
}


void Parser::finishObject(Object *obj) {

    // Get the material name, and make sure that material exists.
    std::string materialName = lexer.getString();
    if (!scene.materials.count(materialName)) {
        std::stringstream ss;
        ss << "no material \"" << materialName << "\"" << std::endl;
        throw(ss.str());
    }
    obj->material = scene.materials[materialName];

    // Set transform, inv transform, and normal transform.
    obj->setup_transform(transformStack.back());

    // Add to the list of objects.
    scene.objects.push_back(obj);
}


void Parser::parseSphere() {
    Sphere *obj = new Sphere;
    obj->radius = lexer.getNumber();
    finishObject(obj);
}


void Parser::parsePlane() {
    Plane *obj = new Plane;
    finishObject(obj);
}


void Parser::parseConic() {
    Conic *obj = new Conic;
    obj->radius1 = lexer.getNumber();
    obj->radius2 = lexer.getNumber();
    obj->zMin = lexer.getNumber();
    obj->zMax = lexer.getNumber();
    finishObject(obj);
}


void Parser::parseMesh() {

    Mesh *obj = new Mesh();
    
    // Normally this comes last, but we want the color on the front.
    finishObject(obj);

    // First try to get a filename, and read an OBJ from it.
    std::string filename;
    try {
        filename = lexer.getString();
        std::cout << "got OBJ filename \"" << filename << "\"" << std::endl;
    } catch (std::string e) {
        // OK.
    }
    
    // Read in the specified OBJ.
    if (filename.size()) {
        if (!obj->readOBJ(filename)) {
            throw std::string("could not open OBJ file");
        }
        std::cout << "Mesh: ";
        std::cout << obj->triangles.size() << " triangles" << std::endl;
        
    // Looks like it is getting embedded directly!
    } else {
       
        // Get a list of vertices counts; 3 or 4 a piece.
        std::vector<double> rawVertCounts = lexer.getNumbers();
        std::vector<int> vertCounts(rawVertCounts.size());
        int totalVertCount = 0;
        for (size_t i = 0; i < rawVertCounts.size(); i++) {
            vertCounts[i] = int(rawVertCounts[i]);
            if (vertCounts[i] < 3 || vertCounts[i] > 4) {
                throw std::string("polygons can only be 3 or 4 vertices");
            }
            totalVertCount += vertCounts[i];
        }
        std::cout << "Mesh: ";
        std::cout << vertCounts.size() << " polygons, ";
        std::cout << totalVertCount << " vertices; ";

        // Get the vertex indices.
        std::vector<double> rawIndices = lexer.getNumbers(totalVertCount);
        std::vector<int> indices(rawIndices.size());
        int maxIndex = 0;
        for (size_t i = 0; i < rawIndices.size(); i++) {
            indices[i] = static_cast<int>(rawIndices[i]);
            maxIndex = std::max(maxIndex, indices[i]);
            if (indices[i] < 0) {
                throw std::string("polygon vertex index out of range; too low");
            }
        }

        std::cout << "(" << (maxIndex + 1) << 
            " independent vertices)" << std::endl;

        // Read the data itself.
        ParamList data = lexer.getParamList(3 * (maxIndex + 1));

        int Pcount = data["P"].size() / 3;

        if (!Pcount) {
            throw std::string("mesh missing P data");
        }

        // Pull in all the position data.
        for (int pi = 0; pi < Pcount; pi++) {
            obj->positions.push_back(Vector(
                data["P"][pi * 3 + 0],
                data["P"][pi * 3 + 1],
                data["P"][pi * 3 + 2]
            ));
        }

        // Pull in all the index data and build the polygons.
        int index_i = 0;
        for (size_t poly_i = 0; poly_i < vertCounts.size(); poly_i++) {
            std::vector<Vertex> polygon;
            for (int vert_i = 0; vert_i < vertCounts[poly_i]; vert_i++) {
                polygon.push_back(Vertex(indices[index_i], -1, -1, -1));
                index_i++;
            }
            if (vertCounts[poly_i] == 3) {
                obj->triangles.push_back(Triangle(polygon[0],
                                                  polygon[1],
                                                  polygon[2]));
            } else if (vertCounts[poly_i] == 4) {
                obj->triangles.push_back(Triangle(polygon[0],
                                                  polygon[1],
                                                  polygon[2]));
                obj->triangles.push_back(Triangle(polygon[0],
                                                  polygon[2],
                                                  polygon[3]));
            }
        }

        obj->updateBBox();
    }
    obj->init();
}


void Parser::parsePointLight() {
    Vector position;
    position[0] = lexer.getNumber();
    position[1] = lexer.getNumber();
    position[2] = lexer.getNumber();

    ParamList params = lexer.getParamList(1, 4);
    
    scene.lights.push_back(PointLight(position, params));
}


